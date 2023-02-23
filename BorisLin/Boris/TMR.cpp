#include "stdafx.h"
#include "TMR.h"

#ifdef MODULE_COMPILATION_TMR

#include "Mesh.h"
#include "MeshParamsControl.h"
#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TMRCUDA.h"
#endif

TMR::TMR(Mesh *pMesh_) :
	RAtmr_p_equation({ "RAp", "RAap", "V" }),
	RAtmr_ap_equation({ "RAp", "RAap", "V" }),
	Modules(),
	TransportBase(pMesh_),
	ProgramStateNames(this, { VINFO(RAtmr_p_equation), VINFO(RAtmr_ap_equation) }, {})
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

TMR::~TMR()
{
	//free memory used for electrical quantities
	pMesh->V.clear();
	pMesh->elC.clear();
	pMesh->E.clear();
	pMesh->S.clear();
}

//check if dM_dt Calculation should be enabled
bool TMR::Need_dM_dt_Calculation(void)
{
	return false;
}

//check if the delsq_V_fixed VEC is needed
bool TMR::Need_delsq_V_fixed_Precalculation(void)
{
	//precalculation only needed in magnetic meshes if CPP-GMR or charge pumping are enabled
	return false;
}

//check if the delsq_S_fixed VEC is needed
bool TMR::Need_delsq_S_fixed_Precalculation(void)
{
	//precalculation only needed in magnetic meshes, or in normal metal meshes with SHE enabled
	return false;
}

//-------------------Abstract base class method implementations

BError TMR::Initialize(void)
{
	BError error(CLASS_STR(TMR));

	//Make list of contacts for TMR calculation
	SuperMesh* pSMesh = pMesh->pSMesh;

	pMeshFM_Bot.clear();
	pMeshFM_Top.clear();
	pMeshAtom_Bot.clear();
	pMeshAtom_Top.clear();

	Rect meshRect = pMesh->GetMeshRect();

	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			Rect candidate_meshRect = (*pSMesh)[idx]->GetMeshRect();

			//candidate mesh at the top
			if (IsE(candidate_meshRect.s.z, meshRect.e.z)) {

				if (candidate_meshRect.intersects(meshRect)) {

					if ((*pSMesh)[idx]->is_atomistic()) pMeshAtom_Top.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
					else if ((*pSMesh)[idx]->GetMeshType() == MESH_FERROMAGNETIC) pMeshFM_Top.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
				}
			}

			//candidate mesh at the botttom
			else if (IsE(candidate_meshRect.e.z, meshRect.s.z)) {

				if (candidate_meshRect.intersects(meshRect)) {

					if ((*pSMesh)[idx]->is_atomistic()) pMeshAtom_Bot.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
					else if ((*pSMesh)[idx]->GetMeshType() == MESH_FERROMAGNETIC) pMeshFM_Bot.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
				}
			}
		}
	}

	//in insulator meshes must mark cells for which magnetic neighbors are not available either side, as empty
	//NOTE, this will require STransport to recalculate CMBND flags, so this is done in STransport initialization which comes after all mesh modules are initialized
#pragma omp parallel for
	for (int j = 0; j < pMesh->n_e.j; j++) {
		for (int i = 0; i < pMesh->n_e.i; i++) {

			int idx = i + j * pMesh->n_e.x;

			DBL3 m_top = DBL3(), m_bot = DBL3();

			//TOP

			//Look at FM meshes first, if any
			for (int mesh_idx = 0; mesh_idx < pMeshFM_Top.size(); mesh_idx++) {

				Rect tmeshRect = pMeshFM_Top[mesh_idx]->GetMeshRect();

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				DBL3 cell_rel_pos = DBL3(
					(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - tmeshRect.s.x,
					(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - tmeshRect.s.y,
					pMeshFM_Top[mesh_idx]->h.z / 2);

				if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

				m_top = normalize(pMeshFM_Top[mesh_idx]->M.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshFM_Top[mesh_idx]->h.z)));
				break;
			}

			if (m_top == DBL3()) {

				//Look at Atomistic meshes, if any
				for (int mesh_idx = 0; mesh_idx < pMeshAtom_Top.size(); mesh_idx++) {

					Rect tmeshRect = pMeshAtom_Top[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - tmeshRect.s.x,
						(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - tmeshRect.s.y,
						pMeshAtom_Top[mesh_idx]->h.z / 2);

					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

					m_top = normalize(pMeshAtom_Top[mesh_idx]->M1.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshAtom_Top[mesh_idx]->h.z)));
					break;
				}
			}

			//BOTTOM

			//Look at FM meshes first, if any
			for (int mesh_idx = 0; mesh_idx < pMeshFM_Bot.size(); mesh_idx++) {

				Rect bmeshRect = pMeshFM_Bot[mesh_idx]->GetMeshRect();

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				DBL3 cell_rel_pos = DBL3(
					(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - bmeshRect.s.x,
					(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - bmeshRect.s.y,
					pMeshFM_Bot[mesh_idx]->meshRect.e.z - pMeshFM_Bot[mesh_idx]->meshRect.s.z - (pMeshFM_Bot[mesh_idx]->h.z / 2));

				//can't couple to an empty cell
				if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

				m_bot = normalize(pMeshFM_Bot[mesh_idx]->M.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshFM_Bot[mesh_idx]->h.z)));
				break;
			}

			if (m_bot == DBL3()) {

				//Look at Atomistic meshes, if any
				for (int mesh_idx = 0; mesh_idx < pMeshAtom_Bot.size(); mesh_idx++) {

					Rect bmeshRect = pMeshAtom_Bot[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - bmeshRect.s.y,
						pMeshAtom_Bot[mesh_idx]->meshRect.e.z - pMeshAtom_Bot[mesh_idx]->meshRect.s.z - (pMeshAtom_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

					m_bot = normalize(pMeshAtom_Bot[mesh_idx]->M1.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshAtom_Bot[mesh_idx]->h.z)));
					break;
				}
			}
			
			if (m_top == DBL3() || m_bot == DBL3()) {

				for (int k = 0; k < pMesh->n_e.k; k++) {

					pMesh->elC.mark_empty(idx + k * pMesh->n_e.x * pMesh->n_e.y);
					pMesh->elC[idx + k * pMesh->n_e.x * pMesh->n_e.y] = 0.0;
					
					pMesh->V.mark_empty(idx + k * pMesh->n_e.x * pMesh->n_e.y);
					pMesh->V[idx + k * pMesh->n_e.x * pMesh->n_e.y] = 0.0;
					
					pMesh->E.mark_empty(idx + k * pMesh->n_e.x * pMesh->n_e.y);
					pMesh->E[idx + k * pMesh->n_e.x * pMesh->n_e.y] = DBL3();
					
					if (pMesh->S.linear_size()) {

						pMesh->S.mark_empty(idx + k * pMesh->n_e.x * pMesh->n_e.y);
						pMesh->S[idx + k * pMesh->n_e.x * pMesh->n_e.y] = DBL3();
					}
				}
			}
		}
	}

	//get a copy of TMR type set and calculate conductivity
	TMR_type = dynamic_cast<InsulatorMesh*>(pMesh)->GetTMRType();

	initialized = true;

	//If CUDA is on, then skip calculating electrical conductivity on CPU below this - we just want to initialize the CPU version here
#if COMPILECUDA == 1
	if (pModuleCUDA) return error;
#endif

	CalculateElectricalConductivity(true);
	return error;
}

BError TMR::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(TMR));

	Uninitialize();

	bool success = true;

	Set_STSolveType();

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		//make sure the cellsize divides the mesh rectangle
		pMesh->n_e = round(pMesh->meshRect / pMesh->h_e);
		if (pMesh->n_e.x < 2) pMesh->n_e.x = 2;
		if (pMesh->n_e.y < 2) pMesh->n_e.y = 2;
		if (pMesh->n_e.z < 2) pMesh->n_e.z = 2;
		pMesh->h_e = pMesh->meshRect / pMesh->n_e;

		//make sure correct memory is assigned for electrical quantities

		//electrical conductivity
		if (pMesh->elC.linear_size()) {

			success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
		}
		else {

			success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, 0.0);
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

			//clear display VEC used for spin currents
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

void TMR::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		UpdateTEquationUserConstants(false);
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (RAtmr_p_equation.is_set()) RAtmr_p_equation.clear();
		if (RAtmr_ap_equation.is_set()) RAtmr_ap_equation.clear();
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError TMR::MakeCUDAModule(void)
{
	BError error(CLASS_STR(TMR));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new TMRCUDA(this);
		pTMRCUDA = dynamic_cast<TMRCUDA*>(pModuleCUDA);
		pTransportBaseCUDA = pTMRCUDA;
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double TMR::UpdateField(void)
{
	if (pSMesh->CurrentTimeStepSolved()) CalculateElectricalConductivity(true);

	//no contribution to total energy density
	return 0.0;
}

//calculate resistance in given rectangle (relative to mesh)
double TMR::GetResistance(Rect rect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		return pTMRCUDA->GetResistance(rect);
	}
#endif

	//To find total resistance, first find total conductance : sum all conductance values from each cell, i.e. elC * hx*hy / thickness
	//finally get resistance as 1 over conductance (i.e. treat all cells as parallel resistors)
	
	//since hx*hy/thickness is a fixed constant, first just find sum of values in elC over given rect (and multiply by constant)
	double conductance = pMesh->elC.sum_nonempty_omp(rect) * pMesh->h_e.x * pMesh->h_e.y / (pMesh->meshRect.height() * pMesh->n_e.k);

	//now return resistance
	if (conductance) return 1.0 / conductance;
	else return 0.0;
}

//-------------------Calculation Methods

//calculate electric field as the negative gradient of V
//VERIFIED - CORRECT
void TMR::CalculateElectricField(void)
{
	//calculate electric field using E = - grad V
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->E.linear_size(); idx++) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (pMesh->V.is_not_empty(idx)) {

			pMesh->E[idx] = -1.0 * pMesh->V.grad_diri(idx);
		}
		else pMesh->E[idx] = DBL3(0);
	}
}

//-------------------Other Calculation Methods

void TMR::CalculateElectricalConductivity(bool force_recalculate)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pTMRCUDA->CalculateElectricalConductivity(force_recalculate);
		return;
	}
#endif

	//Here calculate conductivity based on set TMR formula

	if (force_recalculate) {

#pragma omp parallel for
		for (int j = 0; j < pMesh->n_e.j; j++) {
			for (int i = 0; i < pMesh->n_e.i; i++) {

				int idx = i + j * pMesh->n_e.x;

				if (pMesh->elC.is_not_empty(idx)) {

					//m direction values for top and bottom, used to calculate TMR in this cell
					DBL3 m_top = DBL3(), m_bot = DBL3();

					//TOP

					//Look at FM meshes first, if any
					for (int mesh_idx = 0; mesh_idx < pMeshFM_Top.size(); mesh_idx++) {

						Rect tmeshRect = pMeshFM_Top[mesh_idx]->GetMeshRect();

						//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
						DBL3 cell_rel_pos = DBL3(
							(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - tmeshRect.s.x,
							(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - tmeshRect.s.y,
							pMeshFM_Top[mesh_idx]->h.z / 2);

						if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

						m_top = normalize(pMeshFM_Top[mesh_idx]->M.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshFM_Top[mesh_idx]->h.z)));
						break;
					}

					if (m_top == DBL3()) {

						//Look at Atomistic meshes, if any
						for (int mesh_idx = 0; mesh_idx < pMeshAtom_Top.size(); mesh_idx++) {

							Rect tmeshRect = pMeshAtom_Top[mesh_idx]->GetMeshRect();

							//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
							DBL3 cell_rel_pos = DBL3(
								(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - tmeshRect.s.x,
								(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - tmeshRect.s.y,
								pMeshAtom_Top[mesh_idx]->h.z / 2);

							if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

							m_top = normalize(pMeshAtom_Top[mesh_idx]->M1.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshAtom_Top[mesh_idx]->h.z)));
							break;
						}
					}

					//BOTTOM

					//Look at FM meshes first, if any
					for (int mesh_idx = 0; mesh_idx < pMeshFM_Bot.size(); mesh_idx++) {

						Rect bmeshRect = pMeshFM_Bot[mesh_idx]->GetMeshRect();

						//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
						DBL3 cell_rel_pos = DBL3(
							(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - bmeshRect.s.x,
							(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - bmeshRect.s.y,
							pMeshFM_Bot[mesh_idx]->meshRect.e.z - pMeshFM_Bot[mesh_idx]->meshRect.s.z - (pMeshFM_Bot[mesh_idx]->h.z / 2));

						//can't couple to an empty cell
						if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

						m_bot = normalize(pMeshFM_Bot[mesh_idx]->M.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshFM_Bot[mesh_idx]->h.z)));
						break;
					}

					if (m_bot == DBL3()) {

						//Look at Atomistic meshes, if any
						for (int mesh_idx = 0; mesh_idx < pMeshAtom_Bot.size(); mesh_idx++) {

							Rect bmeshRect = pMeshAtom_Bot[mesh_idx]->GetMeshRect();

							//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
							DBL3 cell_rel_pos = DBL3(
								(i + 0.5) * pMesh->h_e.x + pMesh->meshRect.s.x - bmeshRect.s.x,
								(j + 0.5) * pMesh->h_e.y + pMesh->meshRect.s.y - bmeshRect.s.y,
								pMeshAtom_Bot[mesh_idx]->meshRect.e.z - pMeshAtom_Bot[mesh_idx]->meshRect.s.z - (pMeshAtom_Bot[mesh_idx]->h.z / 2));

							//can't couple to an empty cell
							if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

							m_bot = normalize(pMeshAtom_Bot[mesh_idx]->M1.weighted_average(cell_rel_pos, DBL3(pMesh->h_e.x, pMesh->h_e.y, pMeshAtom_Bot[mesh_idx]->h.z)));
							break;
						}
					}

					//now apply TMR formula to store conductivity value
	
					//cos dependence of RA product (where dRA = RAap - RAp):
					//RA = (RAp + dRA * (1 - m1.m2)/2)
					//so resistivity is (thickness t):
					//ro = (RAp + dRA * (1 - m1.m2)/2) / t
					//so set conductivity as: 1 / ro

					//Slonczewski form : cos dependence of conductivity
					//RA = 2*RAp / [ (1 + RAp/RAap) + (1 - RAp/RAap)cos(theta) ]

					double RAtmr_p = pMesh->RAtmr_p;
					double RAtmr_ap = pMesh->RAtmr_ap;
					double elecCond = pMesh->elecCond;
					pMesh->update_parameters_ecoarse(idx, pMesh->RAtmr_p, RAtmr_p, pMesh->RAtmr_ap, RAtmr_ap, pMesh->elecCond, elecCond);

					//RA bias dependence if set
					if (RAtmr_p_equation.is_set() || RAtmr_ap_equation.is_set()) {

						double bias = 0.0;

						if (pMesh->V.linear_size()) {

							double Vt_1 = pMesh->V[idx + (pMesh->n_e.z - 1) * pMesh->n_e.x * pMesh->n_e.y];
							double Vt_2 = pMesh->V[idx + (pMesh->n_e.z - 2) * pMesh->n_e.x * pMesh->n_e.y];

							double Vb_1 = pMesh->V[idx];
							double Vb_2 = pMesh->V[idx + pMesh->n_e.x * pMesh->n_e.y];

							double Vt = 1.5 * Vt_1 - 0.5 * Vt_2;
							double Vb = 1.5 * Vb_1 - 0.5 * Vb_2;

							bias = Vt - Vb;
						}

						double RAtmr_p0 = RAtmr_p;
						double RAtmr_ap0 = RAtmr_ap;

						if (RAtmr_p_equation.is_set()) RAtmr_p = RAtmr_p_equation.evaluate(RAtmr_p0, RAtmr_ap0, bias);
						if (RAtmr_ap_equation.is_set()) RAtmr_ap = RAtmr_ap_equation.evaluate(RAtmr_p0, RAtmr_ap0, bias);
					}

					for (int k = 0; k < pMesh->n_e.k; k++) {

						if (elecCond > 0.0) {

							//Metallic pinholes
							pMesh->elC[idx + k * pMesh->n_e.x * pMesh->n_e.y] = elecCond;
						}
						else {

							//TMR
	
							//cos form of RA theta dependence
							if (TMR_type == TMR_COS) {

								pMesh->elC[idx + k * pMesh->n_e.x * pMesh->n_e.y] =
									pMesh->meshRect.height() / (RAtmr_p + (RAtmr_ap - RAtmr_p) * (1 - m_top * m_bot) / 2);
							}
							//alternatively Slonczewski form (cos dependence of G)
							else if (TMR_type == TMR_SLONCZEWSKI) {

								pMesh->elC[idx + k * pMesh->n_e.x * pMesh->n_e.y] =
									pMesh->meshRect.height() * ((1 + RAtmr_p / RAtmr_ap) + (1 - RAtmr_p / RAtmr_ap) * m_top * m_bot) / (2 * RAtmr_p);
							}
						}
					}
				}
			}
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
}

//-------------------Other Settings

//Update TEquation object with user constants values
void TMR::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMesh->userConstants.size()) {

		std::vector<std::pair<std::string, double>> constants(pMesh->userConstants.size());
		for (int idx = 0; idx < pMesh->userConstants.size(); idx++) {

			constants[idx] = { pMesh->userConstants.get_key_from_index(idx), pMesh->userConstants[idx] };
		}

		RAtmr_p_equation.set_constants(constants);
		RAtmr_ap_equation.set_constants(constants);

		//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
		if (pModuleCUDA && makeCuda) {

			pTMRCUDA->SetBiasEquation_Parallel(RAtmr_p_equation.get_scalar_fspec());
			pTMRCUDA->SetBiasEquation_AntiParallel(RAtmr_ap_equation.get_scalar_fspec());
		}
#endif
	}
}

BError TMR::SetBiasEquation_Parallel(std::string equation_string)
{
	BError error(CLASS_STR(TMR));

	RAtmr_p_equation.clear();

	if (equation_string.length()) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants(false);

		if (!RAtmr_p_equation.make_from_string(equation_string)) return error(BERROR_INCORRECTSTRING);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		error = pTMRCUDA->SetBiasEquation_Parallel(RAtmr_p_equation.get_scalar_fspec());
	}
#endif

	return error;
}

BError TMR::SetBiasEquation_AntiParallel(std::string equation_string)
{
	BError error(CLASS_STR(TMR));

	RAtmr_ap_equation.clear();

	if (equation_string.length()) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants(false);

		if (!RAtmr_ap_equation.make_from_string(equation_string)) return error(BERROR_INCORRECTSTRING);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		error = pTMRCUDA->SetBiasEquation_AntiParallel(RAtmr_ap_equation.get_scalar_fspec());
	}
#endif

	return error;
}

#endif