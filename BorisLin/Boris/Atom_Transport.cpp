#include "stdafx.h"
#include "Atom_Transport.h"

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

Atom_Transport::Atom_Transport(Atom_Mesh *paMesh_) :
	Modules(),
	TransportBase(paMesh_),
	ProgramStateNames(this, { VINFO(TAMR_conductivity_equation) }, {})
{
	paMesh = paMesh_;

	Set_STSolveType();

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_Transport::~Atom_Transport()
{
	//free memory used for electrical quantities
	paMesh->V.clear();
	paMesh->elC.clear();
	paMesh->E.clear();
	paMesh->S.clear();
}

//check if dM_dt Calculation should be enabled
bool Atom_Transport::Need_dM_dt_Calculation(void)
{
	//enabled for ferromagnetic st solver only
	if (stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

		//enabled for charge pumping and spin pumping
		if (IsNZ((double)paMesh->cpump_eff.get0()) || IsNZ((double)paMesh->pump_eff.get0())) return true;
	}

	return false;
}

//check if the delsq_V_fixed VEC is needed
bool Atom_Transport::Need_delsq_V_fixed_Precalculation(void)
{
	//precalculation only needed in magnetic meshes if CPP-GMR or charge pumping are enabled
	return (stsolve == STSOLVE_FERROMAGNETIC_ATOM && (IsNZ(paMesh->betaD.get0()) || IsNZ(paMesh->cpump_eff.get0())));
}

//check if the delsq_S_fixed VEC is needed
bool Atom_Transport::Need_delsq_S_fixed_Precalculation(void)
{
	return (stsolve == STSOLVE_FERROMAGNETIC_ATOM);
}

//-------------------Abstract base class method implementations

BError Atom_Transport::Initialize(void)
{
	BError error(CLASS_STR(Atom_Transport));

	//is this a "thermoelectric mesh"?
	is_thermoelectric_mesh = (paMesh->Temp.linear_size() && IsNZ(paMesh->Sc.get0()));

	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt.assign(paMesh->h, paMesh->meshRect, DBL3(), paMesh->M1)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		dM_dt.clear();
	}

	if (pSMesh->SolveSpinCurrent()) {

		//Poisson equation helper storage allocation
		if (Need_delsq_V_fixed_Precalculation()) {

			if (!delsq_V_fixed.assign(paMesh->h_e, paMesh->meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_V_fixed.clear();

		if (Need_delsq_S_fixed_Precalculation()) {

			if (!delsq_S_fixed.assign(paMesh->h_e, paMesh->meshRect, DBL3())) return error(BERROR_OUTOFMEMORY_CRIT);
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

BError Atom_Transport::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Transport));

	Uninitialize();

	bool success = true;

	Set_STSolveType();

	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt.assign(paMesh->h, paMesh->meshRect, DBL3(), paMesh->M1)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		if (!dM_dt.linear_size()) dM_dt.clear();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		//make sure the cellsize divides the mesh rectangle
		paMesh->n_e = round(paMesh->meshRect / paMesh->h_e);
		if (paMesh->n_e.x < 2) paMesh->n_e.x = 2;
		if (paMesh->n_e.y < 2) paMesh->n_e.y = 2;
		if (paMesh->n_e.z < 2) paMesh->n_e.z = 2;
		paMesh->h_e = paMesh->meshRect / paMesh->n_e;

		//make sure correct memory is assigned for electrical quantities

		//electrical conductivity
		if (paMesh->M1.linear_size()) {

			//for magnetic meshes set empty cells using information in M (empty cells have zero electrical conductivity) only on initialization. If already initialized then shape already set.
			if (paMesh->elC.linear_size()) {

				success = paMesh->elC.resize(paMesh->h_e, paMesh->meshRect);
			}
			else {

				success = paMesh->elC.assign(paMesh->h_e, paMesh->meshRect, paMesh->elecCond, paMesh->M1);
			}
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (paMesh->V.linear_size()) {

			success &= paMesh->V.resize(paMesh->h_e, paMesh->meshRect, paMesh->elC);
		}
		else {

			success &= paMesh->V.assign(paMesh->h_e, paMesh->meshRect, 0.0, paMesh->elC);
		}

		//electric field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (paMesh->E.linear_size()) {

			success &= paMesh->E.resize(paMesh->h_e, paMesh->meshRect, paMesh->elC);
		}
		else {

			success &= paMesh->E.assign(paMesh->h_e, paMesh->meshRect, DBL3(0.0), paMesh->elC);
		}
	}

	//here also need to update on ODE change as SolveSpinCurrent state might change
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		//spin accumulation if needed - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (success && stsolve != STSOLVE_NONE) {

			if (paMesh->S.linear_size()) {

				success &= paMesh->S.resize(paMesh->h_e, paMesh->meshRect, paMesh->elC);
			}
			else {

				success &= paMesh->S.assign(paMesh->h_e, paMesh->meshRect, DBL3(0.0), paMesh->elC);
			}
		}
		else {

			paMesh->S.clear();

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

BError Atom_Transport::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Transport));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_TransportCUDA(this);
		pTransportCUDA = dynamic_cast<Atom_TransportCUDA*>(pModuleCUDA);
		pTransportBaseCUDA = pTransportCUDA;
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Transport::UpdateField(void)
{
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be iterated only at the end of a step or stage
	if (!paMesh->static_transport_solver && !pSMesh->disabled_transport_solver) {

		//update elC (AMR and temperature)
		if (pSMesh->CurrentTimeStepSolved()) CalculateElectricalConductivity();
	}

	//no contribution to total energy density
	return 0.0;
}

//-------------------Calculation Methods

//calculate electric field as the negative gradient of V
//VERIFIED - CORRECT
void Atom_Transport::CalculateElectricField(void)
{
	if (stsolve == STSOLVE_NONE || stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

		//reset counter so we can accumulate total thermoelectric current if needed
		double mesh_thermoelectric_net_current_ = 0.0;

#pragma omp parallel for reduction(+:mesh_thermoelectric_net_current_)
		for (int idx = 0; idx < paMesh->E.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (paMesh->V.is_not_empty(idx)) {

				//include thermoelectric effect?
				if (is_thermoelectric_mesh) {

					double Sc = paMesh->Sc;
					paMesh->update_parameters_ecoarse(idx, paMesh->Sc, Sc);

					//corresponding index in Temp
					int idx_temp = paMesh->Temp.position_to_cellidx(paMesh->V.cellidx_to_position(idx));

					if (!pSTrans->IsOpenPotential()) {

						//include thermoelectric effect, but no net current generated

						DBL3 shift = paMesh->V.get_shift_to_emptycell(idx);
						if (!shift.IsNull()) {

							//use sided differentials if both neighbors not available
							paMesh->E[idx] = -1 * (paMesh->V.grad_sided(idx) + Sc * paMesh->Temp.grad_sided(idx_temp));
						}
						else paMesh->E[idx] = -1 * (paMesh->V.grad_neu(idx) + Sc * paMesh->Temp.grad_neu(idx_temp));
					}
					else {

						//thermoelectric effect but with open potential, so a net current is generated - must count total current coming out of this mesh

						paMesh->E[idx] = Sc * paMesh->Temp.grad_sided(idx_temp);

						if (paMesh->V.is_cmbnd(idx) || paMesh->V.is_dirichlet(idx)) {

							DBL3 therm_current_density = paMesh->elC[idx] * paMesh->E[idx];

							//below divide by 2 since we want net current, not twice the current come into or out of the mesh

							//x
							if (paMesh->V.is_cmbnd_x(idx) || paMesh->V.is_dirichlet_x(idx)) {

								double cmbnd_current_density = therm_current_density.x;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * paMesh->V.h.y * paMesh->V.h.z / 2;
							}

							//y
							else if (paMesh->V.is_cmbnd_y(idx) || paMesh->V.is_dirichlet_y(idx)) {

								double cmbnd_current_density = therm_current_density.y;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * paMesh->V.h.x * paMesh->V.h.z / 2;
							}

							//z
							if (paMesh->V.is_cmbnd_z(idx) || paMesh->V.is_dirichlet_z(idx)) {

								double cmbnd_current_density = therm_current_density.z;
								mesh_thermoelectric_net_current_ -= cmbnd_current_density * paMesh->V.h.x * paMesh->V.h.y / 2;
							}
						}
					}
				}
				//no thermoelectric effect
				else paMesh->E[idx] = -1.0 * paMesh->V.grad_diri(idx);
			}
			else paMesh->E[idx] = DBL3(0);
		}

		mesh_thermoelectric_net_current = mesh_thermoelectric_net_current_;
	}
}

//-------------------Other Calculation Methods

void Atom_Transport::CalculateElectricalConductivity(bool force_recalculate)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pTransportCUDA->CalculateElectricalConductivity(force_recalculate);
		return;
	}
#endif

	//Include AMR or TAMR?
	if (IsNZ((double)paMesh->amrPercentage) || IsNZ((double)paMesh->tamrPercentage)) {

		if (IsZ((double)paMesh->tamrPercentage)) {

			//only amr
#pragma omp parallel for
			for (int idx = 0; idx < paMesh->elC.linear_size(); idx++) {

				if (paMesh->elC.is_not_empty(idx)) {

					double elecCond = paMesh->elecCond;
					double amrPercentage = paMesh->amrPercentage;
					paMesh->update_parameters_ecoarse(idx, paMesh->elecCond, elecCond, paMesh->amrPercentage, amrPercentage);

					//get current density value at this conductivity cell
					DBL3 jc_value = normalize(paMesh->elC[idx] * paMesh->E[idx]);

					//get M value (M is on n, h mesh so could be different)
					DBL3 m_value = normalize(paMesh->M1[paMesh->elC.cellidx_to_position(idx)]);

					double dotproduct = jc_value * m_value;

					paMesh->elC[idx] = elecCond / (1 + amrPercentage * dotproduct*dotproduct / 100);
				}
			}
		}
		else if (IsZ((double)paMesh->amrPercentage)) {

			bool use_custom_formula = TAMR_conductivity_equation.is_set();

			//only tamr
#pragma omp parallel for
			for (int idx = 0; idx < paMesh->elC.linear_size(); idx++) {

				if (paMesh->elC.is_not_empty(idx)) {

					double elecCond = paMesh->elecCond;
					double tamrPercentage = paMesh->tamrPercentage;
					DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
					DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
					DBL3 mcanis_ea3 = paMesh->mcanis_ea3;
					paMesh->update_parameters_ecoarse(idx, paMesh->elecCond, elecCond, paMesh->tamrPercentage, tamrPercentage, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2, paMesh->mcanis_ea3, mcanis_ea3);

					//get M value (M is on n, h mesh so could be different)
					DBL3 m_value = normalize(paMesh->M1[paMesh->elC.cellidx_to_position(idx)]);

					double dotprod1 = m_value * mcanis_ea1;
					double dotprod2 = m_value * mcanis_ea2;
					double dotprod3 = m_value * mcanis_ea3;

					if (use_custom_formula) {

						paMesh->elC[idx] = elecCond * TAMR_conductivity_equation.evaluate(tamrPercentage / 100, dotprod1, dotprod2, dotprod3);
					}
					else {

						//default formula : ro = ro0 * (1 + TAMR * sin^2(theta)), TAMR is the ratio.
						paMesh->elC[idx] = elecCond / (1 + (tamrPercentage / 100) * (1 - dotprod1 * dotprod1));
					}
				}
			}
		}
		else {

			bool use_custom_formula = TAMR_conductivity_equation.is_set();

			//both amr and tamr
#pragma omp parallel for
			for (int idx = 0; idx < paMesh->elC.linear_size(); idx++) {

				if (paMesh->elC.is_not_empty(idx)) {

					double elecCond = paMesh->elecCond;
					double amrPercentage = paMesh->amrPercentage;
					double tamrPercentage = paMesh->tamrPercentage;
					DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
					DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
					DBL3 mcanis_ea3 = paMesh->mcanis_ea3;
					paMesh->update_parameters_ecoarse(idx, paMesh->elecCond, elecCond, paMesh->amrPercentage, amrPercentage, paMesh->tamrPercentage, tamrPercentage, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2, paMesh->mcanis_ea3, mcanis_ea3);

					//get current density value at this conductivity cell
					DBL3 jc_value = normalize(paMesh->elC[idx] * paMesh->E[idx]);

					//get M value (M is on n, h mesh so could be different)
					DBL3 m_value = normalize(paMesh->M1[paMesh->elC.cellidx_to_position(idx)]);

					double dotproduct = jc_value * m_value;
					double dotprod1 = m_value * mcanis_ea1;
					double dotprod2 = m_value * mcanis_ea2;
					double dotprod3 = m_value * mcanis_ea3;

					if (use_custom_formula) {

						double resistivity_ratio = (1 / TAMR_conductivity_equation.evaluate(tamrPercentage / 100, dotprod1, dotprod2, dotprod3) + (1 + amrPercentage * dotproduct*dotproduct / 100));
						paMesh->elC[idx] = elecCond / resistivity_ratio;
					}
					else {
						
						double resistivity_ratio = (1 + (tamrPercentage / 100) * (1 - dotprod1 * dotprod1)) + (1 + amrPercentage * dotproduct*dotproduct / 100);
						paMesh->elC[idx] = elecCond / resistivity_ratio;
					}
				}
			}
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
	else if (force_recalculate || (paMesh->elecCond.is_tdep() && paMesh->Temp.linear_size())) {

		//no amr but recalculation is being forced, which may include a temperature dependence with non-uniform temperature
#pragma omp parallel for
		for (int idx = 0; idx < paMesh->elC.linear_size(); idx++) {

			if (paMesh->elC.is_not_empty(idx)) {

				double elecCond = paMesh->elecCond;
				paMesh->update_parameters_ecoarse(idx, paMesh->elecCond, elecCond);

				paMesh->elC[idx] = elecCond;
			}
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
}

#endif