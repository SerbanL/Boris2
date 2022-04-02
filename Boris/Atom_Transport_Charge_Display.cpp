#include "stdafx.h"
#include "Atom_Transport.h"

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_TransportCUDA.h"
#endif

//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
//VERIFIED - CORRECT
VEC_VC<DBL3>& Atom_Transport::GetChargeCurrent(void)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC_VC(paMesh->h_e)) return displayVEC_VC;

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
				if (paMesh->V.is_not_empty(idx)) {

					displayVEC_VC[idx] = -paMesh->elC[idx] * paMesh->V.grad_diri(idx);
				}
				else displayVEC_VC[idx] = DBL3(0);
			}
		}
		else {

			bool cppgmr_enabled = IsNZ(paMesh->betaD.get0());
			bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0());
			bool the_enabled = IsNZ(paMesh->the_eff.get0());

			//Current density when contributions from spin accumulation are present
#pragma omp parallel for
			for (int idx = 0; idx < displayVEC_VC.linear_size(); idx++) {

				//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
				if (paMesh->V.is_not_empty(idx)) {

					//magnetic mesh

					DBL3 grad_V = paMesh->V.grad_diri(idx);

					//1. principal term : always present
					displayVEC_VC[idx] = -paMesh->elC[idx] * grad_V;

					//additional contributions if enabled
					if (cppgmr_enabled || cpump_enabled || the_enabled) {

						double mu_s = paMesh->mu_s;
						paMesh->update_parameters_ecoarse(idx, paMesh->mu_s, mu_s);

						int idx_M = paMesh->M1.position_to_cellidx(paMesh->S.cellidx_to_position(idx));

						DBL3 m = paMesh->M1[idx_M] / mu_s;
						DBL33 grad_S = paMesh->S.grad_neu(idx);		//homogeneous Neumann since SHA = 0 in magnetic meshes

						//2. CPP-GMR contribution
						if (cppgmr_enabled) {

							double betaD = paMesh->betaD;
							double De = paMesh->De;
							paMesh->update_parameters_ecoarse(idx, paMesh->betaD, betaD, paMesh->De, De);

							displayVEC_VC[idx] += (grad_S * m) * betaD * De / MUB_E;
						}

						//3. topological Hall effect contribution
						//4. charge pumping contribution
						if (cpump_enabled || the_enabled) {

							double P = paMesh->P;
							double n_density = paMesh->n_density;
							paMesh->update_parameters_ecoarse(idx, paMesh->P, P, paMesh->n_density, n_density);

							DBL33 grad_M = paMesh->M1.grad_neu(idx_M);
							DBL3 dx_m = grad_M.x / mu_s;
							DBL3 dy_m = grad_M.y / mu_s;

							//topological Hall effect contribution
							if (the_enabled) {

								double Bz_the = (dx_m ^ dy_m) * m;
								displayVEC_VC[idx] += paMesh->the_eff.get0() * (P * paMesh->elC[idx] * HBAR_E / (ECHARGE * n_density)) * paMesh->elC[idx] * DBL3(grad_V.y * Bz_the, -grad_V.x * Bz_the, 0.0);
							}

							//charge pumping contribution
							if (cpump_enabled) {

								DBL3 dm_dt = dM_dt[idx_M] / mu_s;
								displayVEC_VC[idx] += paMesh->cpump_eff.get0() * (P * paMesh->elC[idx] * HBAR_E / 2) * DBL3((dm_dt ^ dx_m) * m, (dm_dt ^ dy_m) * m, 0.0);
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
		for (int idx = 0; idx < paMesh->E.linear_size(); idx++) {

			if (paMesh->elC.is_not_empty(idx)) {

				displayVEC_VC[idx] = paMesh->E[idx] * paMesh->elC[idx];
			}
			else {

				displayVEC_VC[idx] = DBL3(0.0);
			}
		}
	}

	return displayVEC_VC;
}

#endif