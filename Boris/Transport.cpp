#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

Transport::Transport(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	pSMesh = pMesh->pSMesh;

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

//set the stsolve indicator depending on current configuration
void Transport::Set_STSolveType(void)
{
	if (!pSMesh->SolveSpinCurrent()) {

		//no spin transport required
		stsolve = STSOLVE_NONE;
	}
	else {

		switch (pMesh->GetMeshType()) {

		case MESH_FERROMAGNETIC:
		case MESH_DIPOLE:
			stsolve = STSOLVE_FERROMAGNETIC;
			break;

		case MESH_ANTIFERROMAGNETIC:
			//to be implemented
			stsolve = STSOLVE_NONE;
			break;

		case MESH_METAL:
			stsolve = STSOLVE_NORMALMETAL;
			break;

		case MESH_INSULATOR:
			//to be implemented - currently insulator meshes do not even allow a Transport module
			stsolve = STSOLVE_TUNNELING;
			break;
		}
	}

#if COMPILECUDA == 1
	//copy over to TransportCUDA module - need to use it in .cu files
	if (pModuleCUDA) pTransportCUDA->Set_STSolveType();
#endif
}

bool Transport::SetFixedPotentialCells(Rect rectangle, double potential)
{ 
	bool success = true;

#if COMPILECUDA == 1
	if (pModuleCUDA) success &= pTransportCUDA->SetFixedPotentialCells(rectangle, potential);
#endif

	success &= pMesh->V.set_dirichlet_conditions(rectangle, potential);

	return success;
}

void Transport::ClearFixedPotentialCells(void)
{
	pMesh->V.clear_dirichlet_flags();

#if COMPILECUDA == 1
	if (pModuleCUDA) pTransportCUDA->ClearFixedPotentialCells();
#endif
}

//check if dM_dt Calculation should be enabled
bool Transport::Need_dM_dt_Calculation(void)
{
	//enabled for ferromagnetic st solver only
	if (stsolve == STSOLVE_FERROMAGNETIC) {

		//enabled for charge pumping (for spin pumping we can just read it cell by cell, don't need vector calculus in this case).
		if (IsNZ((double)pMesh->cpump_eff.get0())) return true;
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
			if (pMesh->elC.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

				success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
			}
			else {

				success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, pMesh->elecCond, pMesh->M);
			}
		}
		else {

			//for non-ferromagnetic meshes (e.g. simple electrical conductor) then elC contains directly any empty cells information
			if (pMesh->elC.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

				success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
			}
			else {

				success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, pMesh->elecCond);
			}
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMesh->V.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

			success &= pMesh->V.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
		}
		else {

			success &= pMesh->V.assign(pMesh->h_e, pMesh->meshRect, 0.0, pMesh->elC);
		}

		//electric field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMesh->E.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

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

			if (pMesh->S.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

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
		pModuleCUDA = new TransportCUDA(pMesh, pSMesh, this);
		pTransportCUDA = dynamic_cast<TransportCUDA*>(pModuleCUDA);
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
		if (pSMesh->CurrentTimeStepSolved())
			CalculateElectricalConductivity();
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

double Transport::CalculateElectrodeCurrent(Rect &electrode_rect)
{
	if (!pMesh->meshRect.intersects(electrode_rect)) return 0.0;

	Box electrode_box = pMesh->V.box_from_rect_max(electrode_rect);

#if COMPILECUDA == 1
	if (pModuleCUDA) return pTransportCUDA->CalculateElectrodeCurrent(electrode_box);
#endif

	//look at all cells surrounding the electrode

	double current = 0;

	//cells on -x side
	if (electrode_box.s.i > 1) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current += -pMesh->elC[INT3(electrode_box.s.i - 1, j, k)] *
					(pMesh->V[INT3(electrode_box.s.i - 1, j, k)] - pMesh->V[INT3(electrode_box.s.i - 2, j, k)]) * pMesh->h_e.y * pMesh->h_e.z / pMesh->h_e.x;
			}
		}
	}

	//cells on +x side
	if (electrode_box.e.i + 1 < pMesh->n_e.i) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current -= -pMesh->elC[INT3(electrode_box.e.i, j, k)] *
					(pMesh->V[INT3(electrode_box.e.i + 1, j, k)] - pMesh->V[INT3(electrode_box.e.i, j, k)]) * pMesh->h_e.y * pMesh->h_e.z / pMesh->h_e.x;
			}
		}
	}

	//cells on -y side
	if (electrode_box.s.j > 1) {
#pragma omp parallel for reduction(+:current)
		for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current += -pMesh->elC[INT3(i, electrode_box.s.j - 1, k)] *
					(pMesh->V[INT3(i, electrode_box.s.j - 1, k)] - pMesh->V[INT3(i, electrode_box.s.j - 2, k)]) * pMesh->h_e.x * pMesh->h_e.z / pMesh->h_e.y;
			}
		}
	}

	//cells on +y side
	if (electrode_box.e.j + 1 < pMesh->n_e.j) {
#pragma omp parallel for reduction(+:current)
		for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current -= -pMesh->elC[INT3(i, electrode_box.e.j, k)] *
					(pMesh->V[INT3(i, electrode_box.e.j + 1, k)] - pMesh->V[INT3(i, electrode_box.e.j, k)]) * pMesh->h_e.x * pMesh->h_e.z / pMesh->h_e.y;
			}
		}
	}

	//cells on -z side
	if (electrode_box.s.k > 1) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {

				current += -pMesh->elC[INT3(i, j, electrode_box.s.k - 1)] *
					(pMesh->V[INT3(i, j, electrode_box.s.k - 1)] - pMesh->V[INT3(i, j, electrode_box.s.k - 2)]) * pMesh->h_e.x * pMesh->h_e.y / pMesh->h_e.z;
			}
		}
	}

	//cells on +z side
	if (electrode_box.e.k + 1 < pMesh->n_e.k) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {

				current -= -pMesh->elC[INT3(i, j, electrode_box.e.k)] *
					(pMesh->V[INT3(i, j, electrode_box.e.k + 1)] - pMesh->V[INT3(i, j, electrode_box.e.k)]) * pMesh->h_e.x * pMesh->h_e.y / pMesh->h_e.z;
			}
		}
	}

	return current;
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

	//Include AMR?
	if (pMesh->M.linear_size() && (IsNZ((double)pMesh->amrPercentage))) {

		//with amr
		#pragma omp parallel for
		for (int idx = 0; idx < pMesh->elC.linear_size(); idx++) {

			if (pMesh->elC.is_not_empty(idx)) {

				double elecCond = pMesh->elecCond;
				double amrPercentage = pMesh->amrPercentage;
				pMesh->update_parameters_ecoarse(idx, pMesh->elecCond, elecCond, pMesh->amrPercentage, amrPercentage);

				//get current density value at this conductivity cell
				DBL3 Jc_value = pMesh->elC[idx] * pMesh->E[idx];

				//get M value (M is on n, h mesh so could be different)
				DBL3 M_value = pMesh->M[pMesh->elC.cellidx_to_position(idx)];

				double magnitude = Jc_value.norm() * M_value.norm();
				double dotproduct = 0.0;

				if (IsNZ(magnitude)) dotproduct = (Jc_value * M_value) / magnitude;

				pMesh->elC[idx] = elecCond / (1 + amrPercentage*dotproduct*dotproduct / 100);
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

//-------------------Others

//called by MoveMesh method in this mesh - move relevant transport quantities
void Transport::MoveMesh_Transport(double x_shift)
{
	double mesh_end_size = pMesh->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(pMesh->meshRect.s + DBL3(mesh_end_size, 0, 0), pMesh->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->elC()->shift_x(pMesh->elC.linear_size(), x_shift, shift_rect);

		//shift spin accumulation if present
		if (pMesh->S.linear_size()) pMesh->pMeshCUDA->S()->shift_x(pMesh->S.linear_size(), x_shift, shift_rect);

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);

		return;
	}
#endif

	//1. electrical conductivity
	pMesh->elC.shift_x(x_shift, shift_rect);

	//shift spin accumulation if present
	if(pMesh->S.linear_size()) pMesh->S.shift_x(x_shift, shift_rect);

	//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
	pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
}

#endif