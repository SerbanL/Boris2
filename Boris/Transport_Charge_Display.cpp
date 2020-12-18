#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

//prepare displayVEC_VC ready for calculation of display quantity
bool Transport::PrepareDisplayVEC_VC(DBL3 cellsize)
{
	if (pMesh->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC_VC.assign(cellsize, pMesh->meshRect, DBL3(0.0));

		return true;
	}
	else displayVEC_VC.clear();

	return false;
}

//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
//VERIFIED - CORRECT
VEC_VC<DBL3>& Transport::GetChargeCurrent(void)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC_VC(pMesh->h_e)) return displayVEC_VC;

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		cu_obj<cuVEC_VC<cuReal3>>& cudisplayVEC_VC = GetChargeCurrentCUDA();
		cudisplayVEC_VC()->copy_to_cpuvec(displayVEC_VC);

		return displayVEC_VC;
	}
#endif

	//compute charge current and store result in displayVEC_VC

	if (!pSMesh->disabled_transport_solver) {

		if (stsolve == STSOLVE_NONE) {

			//calculate current density using Jc = -sigma * grad V
#pragma omp parallel for
			for (int idx = 0; idx < displayVEC_VC.linear_size(); idx++) {

				//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
				if (pMesh->V.is_not_empty(idx)) {

					displayVEC_VC[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri(idx);
				}
				else displayVEC_VC[idx] = DBL3(0);
			}
		}
		else {

			bool cppgmr_enabled = IsNZ(pMesh->betaD.get0());
			bool cpump_enabled = IsNZ(pMesh->cpump_eff.get0());
			bool the_enabled = IsNZ(pMesh->the_eff.get0());

			//Current density when contributions from spin accumulation are present
#pragma omp parallel for
			for (int idx = 0; idx < displayVEC_VC.linear_size(); idx++) {

				//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
				if (pMesh->V.is_not_empty(idx)) {

					if (stsolve == STSOLVE_NORMALMETAL) {

						//non-magnetic mesh

						if (IsZ(pMesh->iSHA.get0())) {

							//no iSHE contribution.
							displayVEC_VC[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri(idx);
						}
						else {

							double SHA = pMesh->SHA;
							double iSHA = pMesh->iSHA;
							double De = pMesh->De;
							pMesh->update_parameters_ecoarse(idx, pMesh->SHA, SHA, pMesh->iSHA, iSHA, pMesh->De, De);

							//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
							displayVEC_VC[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri_nneu(idx, (iSHA * De / (MUB_E * pMesh->elC[idx])) * pMesh->S.curl_neu(idx));

							//must also add iSHE contribution -> here we must use non-homogeneous Neumann boundary conditions when calculating S differentials
							displayVEC_VC[idx] += (iSHA * De / MUB_E) * pMesh->S.curl_nneu(idx, epsilon3(pMesh->E[idx]) * (SHA * pMesh->elC[idx] * MUB_E / De));
						}
					}
					else {

						//magnetic mesh

						DBL3 grad_V = pMesh->V.grad_diri(idx);

						//1. principal term : always present
						displayVEC_VC[idx] = -pMesh->elC[idx] * grad_V;

						//additional contributions if enabled
						if (cppgmr_enabled || cpump_enabled || the_enabled) {

							double Ms = pMesh->Ms;
							pMesh->update_parameters_ecoarse(idx, pMesh->Ms, Ms);

							int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(idx));

							DBL3 m = pMesh->M[idx_M] / Ms;
							DBL33 grad_S = pMesh->S.grad_neu(idx);		//homogeneous Neumann since SHA = 0 in magnetic meshes

							//2. CPP-GMR contribution
							if (cppgmr_enabled) {

								double betaD = pMesh->betaD;
								double De = pMesh->De;
								pMesh->update_parameters_ecoarse(idx, pMesh->betaD, betaD, pMesh->De, De);

								displayVEC_VC[idx] += (grad_S * m) * betaD * De / MUB_E;
							}

							//3. topological Hall effect contribution
							//4. charge pumping contribution
							if (cpump_enabled || the_enabled) {

								double P = pMesh->P;
								double n_density = pMesh->n_density;
								pMesh->update_parameters_ecoarse(idx, pMesh->P, P, pMesh->n_density, n_density);

								DBL33 grad_M = pMesh->M.grad_neu(idx_M);
								DBL3 dx_m = grad_M.x / Ms;
								DBL3 dy_m = grad_M.y / Ms;

								//topological Hall effect contribution
								if (the_enabled) {

									double Bz_the = (dx_m ^ dy_m) * m;
									displayVEC_VC[idx] += pMesh->the_eff.get0() * (P * pMesh->elC[idx] * HBAR_E / (ECHARGE * n_density)) * pMesh->elC[idx] * DBL3(grad_V.y * Bz_the, -grad_V.x *Bz_the, 0.0);
								}

								//charge pumping contribution
								if (cpump_enabled) {

									DBL3 dm_dt = dM_dt[idx_M] / Ms;
									displayVEC_VC[idx] += pMesh->cpump_eff.get0() * (P * pMesh->elC[idx] * HBAR_E / 2) * DBL3((dm_dt ^ dx_m) * m, (dm_dt ^ dy_m) * m, 0.0);
								}
							}
						}
					}
				}
				else displayVEC_VC[idx] = DBL3(0);
			}
		}
	}
	else {

		//if transport solver disabled we need to set displayVEC_VC directly from E and elC as Jc = elC * E
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->E.linear_size(); idx++) {

			if (pMesh->elC.is_not_empty(idx)) {

				displayVEC_VC[idx] = pMesh->E[idx] * pMesh->elC[idx];
			}
			else {

				displayVEC_VC[idx] = DBL3(0.0);
			}
		}
	}

	return displayVEC_VC;
}

DBL3 Transport::GetAverageChargeCurrent(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		//update displayVEC_VC with charge current in TransportCUDA
		GetChargeCurrentCUDA();

		//average charge current in displayVEC_VC in TransportCUDA
		return cuReal3(pTransportCUDA->displayVEC_VC()->average_nonempty(pMesh->n_e.dim(), rectangle));
	}
#endif

	//update displayVEC_VC with charge current
	GetChargeCurrent();

	//average charge current in displayVEC_VC
	return displayVEC_VC.average_nonempty_omp(rectangle);
}

#if COMPILECUDA == 1
cu_obj<cuVEC_VC<cuReal3>>& Transport::GetChargeCurrentCUDA(void)
{
	return pTransportCUDA->GetChargeCurrent();
}
#endif

#endif