#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

Transport::Transport(Mesh *pMesh_) :
	Modules(),
	TransportBase(pMesh_),
	ProgramStateNames(this, { VINFO(TAMR_conductivity_equation) }, {})
{
	pMesh = pMesh_;

	Set_STSolveType();

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Transport::~Transport()
{
	//free memory used for electrical quantities
	pMesh->V.clear();
	pMesh->elC.clear();
	pMesh->E.clear();
	pMesh->S.clear();
}

//check if dM_dt Calculation should be enabled
bool Transport::Need_dM_dt_Calculation(void)
{
	//enabled for ferromagnetic st solver only
	if (stsolve == STSOLVE_FERROMAGNETIC) {

		//enabled for charge pumping and spin pumping
		if (IsNZ((double)pMesh->cpump_eff.get0()) || IsNZ((double)pMesh->pump_eff.get0())) return true;
	}

	return false;
}

//check if the delsq_V_fixed VEC is needed
bool Transport::Need_delsq_V_fixed_Precalculation(void)
{
	//precalculation only needed in magnetic meshes if CPP-GMR or charge pumping are enabled
	return (stsolve == STSOLVE_FERROMAGNETIC && (IsNZ(pMesh->betaD.get0()) || IsNZ(pMesh->cpump_eff.get0())));
}

//check if the delsq_S_fixed VEC is needed
bool Transport::Need_delsq_S_fixed_Precalculation(void)
{
	//precalculation only needed in magnetic meshes, or in normal metal meshes with SHE enabled
	return (stsolve == STSOLVE_FERROMAGNETIC || (stsolve == STSOLVE_NORMALMETAL && IsNZ(pMesh->SHA.get0())));
}

//-------------------Abstract base class method implementations

BError Transport::Initialize(void)
{
	BError error(CLASS_STR(Transport));

	//is this a "thermoelectric mesh"?
	is_thermoelectric_mesh = (pMesh->Temp.linear_size() && IsNZ(pMesh->Sc.get0()));

	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt.assign(pMesh->h, pMesh->meshRect, DBL3(), pMesh->M)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		dM_dt.clear();
	}

	if (pSMesh->SolveSpinCurrent()) {

		//Poisson equation helper storage allocation
		if (Need_delsq_V_fixed_Precalculation()) {

			if (!delsq_V_fixed.assign(pMesh->h_e, pMesh->meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_V_fixed.clear();

		if (Need_delsq_S_fixed_Precalculation()) {

			if (!delsq_S_fixed.assign(pMesh->h_e, pMesh->meshRect, DBL3())) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_S_fixed.clear();
	}
	else {

		//Poisson equation helper storage not needed
		delsq_V_fixed.clear();
		delsq_S_fixed.clear();
	}

	initialized = true;

	return error;
}

BError Transport::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Transport));

	Uninitialize();

	bool success = true;

	Set_STSolveType();

	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt.assign(pMesh->h, pMesh->meshRect, DBL3(), pMesh->M)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		if (!dM_dt.linear_size()) dM_dt.clear();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		//make sure the cellsize divides the mesh rectangle
		pMesh->n_e = round(pMesh->meshRect / pMesh->h_e);
		if (pMesh->n_e.x < 2) pMesh->n_e.x = 2;
		if (pMesh->n_e.y < 2) pMesh->n_e.y = 2;
		if (pMesh->n_e.z < 2) pMesh->n_e.z = 2;
		pMesh->h_e = pMesh->meshRect / pMesh->n_e;

		//make sure correct memory is assigned for electrical quantities

		//electrical conductivity
		if (pMesh->M.linear_size()) {

			//for magnetic meshes set empty cells using information in M (empty cells have zero electrical conductivity) only on initialization. If already initialized then shape already set.
			if (pMesh->elC.linear_size()) {

				success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
			}
			else {

				success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, pMesh->elecCond, pMesh->M);
			}
		}
		else {

			//for non-ferromagnetic meshes (e.g. simple electrical conductor) then elC contains directly any empty cells information
			if (pMesh->elC.linear_size()) {

				success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
			}
			else {

				success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, pMesh->elecCond);
			}
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMesh->V.linear_size()) {

			success &= pMesh->V.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
		}
		else {

			success &= pMesh->V.assign(pMesh->h_e, pMesh->meshRect, 0.0, pMesh->elC);
		}

		//electric field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMesh->E.linear_size()) {

			success &= pMesh->E.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
		}
		else {

			success &= pMesh->E.assign(pMesh->h_e, pMesh->meshRect, DBL3(0.0), pMesh->elC);
		}
	}

	//here also need to update on ODE change as SolveSpinCurrent state might change
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		//spin accumulation if needed - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (success && stsolve != STSOLVE_NONE) {

			if (pMesh->S.linear_size()) {

				success &= pMesh->S.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
			}
			else {

				success &= pMesh->S.assign(pMesh->h_e, pMesh->meshRect, DBL3(0.0), pMesh->elC);
			}
		}
		else {

			pMesh->S.clear();

			//clear display VEC used for spin currents and torque
			displayVEC.clear();
		}
	}

	if (!success) return error(BERROR_OUTOFMEMORY_CRIT);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Transport::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Transport));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		
		pModuleCUDA = new TransportCUDA(this);
		pTransportCUDA = dynamic_cast<TransportCUDA*>(pModuleCUDA);
		pTransportBaseCUDA = pTransportCUDA;
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Transport::UpdateField(void)
{	
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be interated only at the end of a step or stage
	if (!pMesh->static_transport_solver && !pSMesh->disabled_transport_solver) {

		//update elC (AMR and temperature)
		if (pSMesh->CurrentTimeStepSolved()) CalculateElectricalConductivity();
	}

	//no contribution to total energy density
	return 0.0;
}

//-------------------Calculation Methods

//calculate electric field as the negative gradient of V
//VERIFIED - CORRECT
void Transport::CalculateElectricField(void)
{
	if (stsolve == STSOLVE_NONE || stsolve == STSOLVE_FERROMAGNETIC) {

		//reset counter so we can accumulate total thermoelectric current if needed
		double mesh_thermoelectric_net_current_ = 0.0;

#pragma omp parallel for reduction(+:mesh_thermoelectric_net_current_)
		for (int idx = 0; idx < pMesh->E.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (pMesh->V.is_not_empty(idx)) {

				//include thermoelectric effect?
				if (is_thermoelectric_mesh) {

					double Sc = pMesh->Sc;
					pMesh->update_parameters_ecoarse(idx, pMesh->Sc, Sc);

					//corresponding index in Temp
					int idx_temp = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(idx));

					if (!pSTrans->IsOpenPotential()) {

						//include thermoelectric effect, but no net current generated

						DBL3 shift = pMesh->V.get_shift_to_emptycell(idx);
						if (!shift.IsNull()) {

							//use sided differentials if both neighbors not available
							pMesh->E[idx] = -1 * (pMesh->V.grad_sided(idx) + Sc * pMesh->Temp.grad_sided(idx_temp));
						}
						else pMesh->E[idx] = -1 * (pMesh->V.grad_neu(idx) + Sc * pMesh->Temp.grad_neu(idx_temp));
					}
					else {

						//thermoelectric effect but with open potential, so a net current is generated - must count total current coming out of this mesh

						pMesh->E[idx] = Sc * pMesh->Temp.grad_sided(idx_temp);

						if (pMesh->V.is_cmbnd(idx) || pMesh->V.is_dirichlet(idx)) {

							DBL3 therm_current_density = pMesh->elC[idx] * pMesh->E[idx];
							
							//below divide by 2 since we want net current, not twice the current come into or out of the mesh

							//x
							if (pMesh->V.is_cmbnd_x(idx) || pMesh->V.is_dirichlet_x(idx)) {

								double cmbnd_current_density = therm_current_density.x;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * pMesh->V.h.y * pMesh->V.h.z / 2;
							}

							//y
							else if (pMesh->V.is_cmbnd_y(idx) || pMesh->V.is_dirichlet_y(idx)) {

								double cmbnd_current_density = therm_current_density.y;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * pMesh->V.h.x * pMesh->V.h.z / 2;
							}
							
							//z
							if (pMesh->V.is_cmbnd_z(idx) || pMesh->V.is_dirichlet_z(idx)) {

								double cmbnd_current_density = therm_current_density.z;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * pMesh->V.h.x * pMesh->V.h.y / 2;
							}
						}
					}
				}
				//no thermoelectric effect//include thermoelectric effect?
				if (is_thermoelectric_mesh) {

					double Sc = pMesh->Sc;
					pMesh->update_parameters_ecoarse(idx, pMesh->Sc, Sc);

					//corresponding index in Temp
					int idx_temp = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(idx));

					if (!pSTrans->IsOpenPotential()) {

						//include thermoelectric effect, but no net current generated

						DBL3 shift = pMesh->V.get_shift_to_emptycell(idx);
						if (!shift.IsNull()) {

							//use sided differentials if both neighbors not available
							pMesh->E[idx] = -1 * (pMesh->V.grad_sided(idx) + Sc * pMesh->Temp.grad_sided(idx_temp));
						}
						else pMesh->E[idx] = -1 * (pMesh->V.grad_neu(idx) + Sc * pMesh->Temp.grad_neu(idx_temp));
					}
					else {

						//thermoelectric effect but with open potential, so a net current is generated - must count total current coming out of this mesh

						pMesh->E[idx] = Sc * pMesh->Temp.grad_sided(idx_temp);

						if (pMesh->V.is_cmbnd(idx) || pMesh->V.is_dirichlet(idx)) {

							DBL3 therm_current_density = pMesh->elC[idx] * pMesh->E[idx];
							
							//below divide by 2 since we want net current, not twice the current come into or out of the mesh

							//x
							if (pMesh->V.is_cmbnd_x(idx) || pMesh->V.is_dirichlet_x(idx)) {

								double cmbnd_current_density = therm_current_density.x;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * pMesh->V.h.y * pMesh->V.h.z / 2;
							}

							//y
							else if (pMesh->V.is_cmbnd_y(idx) || pMesh->V.is_dirichlet_y(idx)) {

								double cmbnd_current_density = therm_current_density.y;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * pMesh->V.h.x * pMesh->V.h.z / 2;
							}
							
							//z
							if (pMesh->V.is_cmbnd_z(idx) || pMesh->V.is_dirichlet_z(idx)) {

								double cmbnd_current_density = therm_current_density.z;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * pMesh->V.h.x * pMesh->V.h.y / 2;
							}
						}
					}
				}
				//no thermoelectric effect
				else pMesh->E[idx] = -1.0 * pMesh->V.grad_diri(idx);
			}
			else pMesh->E[idx] = DBL3(0);
		}

		mesh_thermoelectric_net_current = mesh_thermoelectric_net_current_;
	}
	else {

		bool ishe_enabled = IsNZ(pMesh->iSHA.get0()) && (stsolve == STSOLVE_NORMALMETAL);

		//calculate electric field using E = - grad V, but using non-homogeneous Neumann boundary conditions
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->E.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (pMesh->V.is_not_empty(idx)) {

				if (ishe_enabled) {

					//ISHE enabled - use nonhomogeneous Neumann boundary conditions
					double iSHA = pMesh->iSHA;
					double De = pMesh->De;
					pMesh->update_parameters_ecoarse(idx, pMesh->iSHA, iSHA, pMesh->De, De);

					pMesh->E[idx] = -1.0 * pMesh->V.grad_diri_nneu(idx, (iSHA * De / (MUB_E * pMesh->elC[idx])) * pMesh->S.curl_neu(idx));
				}
				else {

					//use homogeneous Neumann boundary conditions - no ISHE
					pMesh->E[idx] = -1.0 * pMesh->V.grad_diri(idx);
				}
			}
			else pMesh->E[idx] = DBL3(0);
		}
	}
}

//-------------------Other Calculation Methods

void Transport::CalculateElectricalConductivity(bool force_recalculate)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pTransportCUDA->CalculateElectricalConductivity(force_recalculate);
		return;
	}
#endif

	//Include AMR or TAMR?
	if (pMesh->M.linear_size() && (IsNZ((double)pMesh->amrPercentage) || IsNZ((double)pMesh->tamrPercentage))) {

		if (IsZ((double)pMesh->tamrPercentage)) {

			//only amr
#pragma omp parallel for
			for (int idx = 0; idx < pMesh->elC.linear_size(); idx++) {

				if (pMesh->elC.is_not_empty(idx)) {

					double elecCond = pMesh->elecCond;
					double amrPercentage = pMesh->amrPercentage;
					pMesh->update_parameters_ecoarse(idx, pMesh->elecCond, elecCond, pMesh->amrPercentage, amrPercentage);

					//get current density value at this conductivity cell
					DBL3 jc_value = normalize(pMesh->elC[idx] * pMesh->E[idx]);

					//get M value (M is on n, h mesh so could be different)
					DBL3 m_value = normalize(pMesh->M[pMesh->elC.cellidx_to_position(idx)]);

					double dotproduct = jc_value * m_value;

					pMesh->elC[idx] = elecCond / (1 + amrPercentage * dotproduct*dotproduct / 100);
				}
			}
		}
		else if (IsZ((double)pMesh->amrPercentage)) {

			bool use_custom_formula = TAMR_conductivity_equation.is_set();

			//only tamr
#pragma omp parallel for
			for (int idx = 0; idx < pMesh->elC.linear_size(); idx++) {

				if (pMesh->elC.is_not_empty(idx)) {

					double elecCond = pMesh->elecCond;
					double tamrPercentage = pMesh->tamrPercentage;
					DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
					DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
					DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
					pMesh->update_parameters_ecoarse(idx, pMesh->elecCond, elecCond, pMesh->tamrPercentage, tamrPercentage, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

					//get M value (M is on n, h mesh so could be different)
					DBL3 m_value = normalize(pMesh->M[pMesh->elC.cellidx_to_position(idx)]);

					double dotprod1 = m_value * mcanis_ea1;
					double dotprod2 = m_value * mcanis_ea2;
					double dotprod3 = m_value * mcanis_ea3;

					if (use_custom_formula) {

						pMesh->elC[idx] = elecCond * TAMR_conductivity_equation.evaluate(tamrPercentage / 100, dotprod1, dotprod2, dotprod3);
					}
					else {

						
						//default formula : ro = ro0 * (1 + TAMR * sin^2(theta)), TAMR is the ratio.
						pMesh->elC[idx] = elecCond / (1 + (tamrPercentage / 100) * (1 - dotprod1*dotprod1));
					}
				}
			}
		}
		else {

			bool use_custom_formula = TAMR_conductivity_equation.is_set();

			//both amr and tamr
#pragma omp parallel for
			for (int idx = 0; idx < pMesh->elC.linear_size(); idx++) {

				if (pMesh->elC.is_not_empty(idx)) {

					double elecCond = pMesh->elecCond;
					double amrPercentage = pMesh->amrPercentage;
					double tamrPercentage = pMesh->tamrPercentage;
					DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
					DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
					DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
					pMesh->update_parameters_ecoarse(idx, pMesh->elecCond, elecCond, pMesh->amrPercentage, amrPercentage, pMesh->tamrPercentage, tamrPercentage, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

					//get current density value at this conductivity cell
					DBL3 jc_value = normalize(pMesh->elC[idx] * pMesh->E[idx]);

					//get M value (M is on n, h mesh so could be different)
					DBL3 m_value = normalize(pMesh->M[pMesh->elC.cellidx_to_position(idx)]);

					double dotproduct = jc_value * m_value;

					double dotprod1 = m_value * mcanis_ea1;
					double dotprod2 = m_value * mcanis_ea2;
					double dotprod3 = m_value * mcanis_ea3;

					if (use_custom_formula) {

						double resistivity_ratio = (1 / TAMR_conductivity_equation.evaluate(tamrPercentage / 100, dotprod1, dotprod2, dotprod3) + (1 + amrPercentage * dotproduct*dotproduct / 100));
						pMesh->elC[idx] = elecCond / resistivity_ratio;
					}
					else {

						double resistivity_ratio = (1 + (tamrPercentage / 100) * (1 - dotprod1 * dotprod1)) + (1 + amrPercentage * dotproduct*dotproduct / 100);
						pMesh->elC[idx] = elecCond / resistivity_ratio;
					}
				}
			}
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
	else if (force_recalculate || (pMesh->elecCond.is_tdep() && pMesh->Temp.linear_size())) {
			
		//no amr but recalculation is being forced, which may include a temperature dependence with non-uniform temperature
		#pragma omp parallel for
		for (int idx = 0; idx < pMesh->elC.linear_size(); idx++) {

			if (pMesh->elC.is_not_empty(idx)) {

				double elecCond = pMesh->elecCond;
				pMesh->update_parameters_ecoarse(idx, pMesh->elecCond, elecCond);

				pMesh->elC[idx] = elecCond;
			}
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
}

#endif