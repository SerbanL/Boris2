#include "stdafx.h"
#include "Atom_Transport.h"

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

Atom_Transport::Atom_Transport(Atom_Mesh *paMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	paMesh = paMesh_;
	pMeshBase = paMesh;

	pSMesh = paMesh->pSMesh;

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

//set the stsolve indicator depending on current configuration
void Atom_Transport::Set_STSolveType(void)
{
	if (!pSMesh->SolveSpinCurrent()) {

		//no spin transport required
		stsolve = STSOLVE_NONE;
	}
	else {

		stsolve = STSOLVE_FERROMAGNETIC_ATOM;
	}

#if COMPILECUDA == 1
	//copy over to TransportCUDA module - need to use it in .cu files
	if (pModuleCUDA) pTransportCUDA->Set_STSolveType();
#endif
}

bool Atom_Transport::SetFixedPotentialCells(Rect rectangle, double potential)
{
	bool success = true;

#if COMPILECUDA == 1
	if (pModuleCUDA) success &= pTransportCUDA->SetFixedPotentialCells(rectangle, potential);
#endif

	success &= paMesh->V.set_dirichlet_conditions(rectangle, potential);

	return success;
}

void Atom_Transport::ClearFixedPotentialCells(void)
{
	paMesh->V.clear_dirichlet_flags();

#if COMPILECUDA == 1
	if (pModuleCUDA) pTransportCUDA->ClearFixedPotentialCells();
#endif
}

void Atom_Transport::Set_Linear_PotentialDrop(DBL3 ground_electrode_center, double ground_potential, DBL3 electrode_center, double electrode_potential)
{
	paMesh->V.set_linear(ground_electrode_center, ground_potential, electrode_center, electrode_potential);
}

//check if dM_dt Calculation should be enabled
bool Atom_Transport::Need_dM_dt_Calculation(void)
{
	//TO DO
	//enabled for ferromagnetic st solver only
	//if (stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

		//enabled for charge pumping (for spin pumping we can just read it cell by cell, don't need vector calculus in this case).
		//if (IsNZ((double)paMesh->cpump_eff.get0())) return true;
	//}

	return false;
}

//check if the delsq_V_fixed VEC is needed
bool Atom_Transport::Need_delsq_V_fixed_Precalculation(void)
{
	//TO DO
	//precalculation only needed in magnetic meshes if CPP-GMR or charge pumping are enabled
	//return (stsolve == STSOLVE_FERROMAGNETIC_ATOM && (IsNZ(paMesh->betaD.get0()) || IsNZ(paMesh->cpump_eff.get0())));
	return false;
}

//check if the delsq_S_fixed VEC is needed
bool Atom_Transport::Need_delsq_S_fixed_Precalculation(void)
{
	//TO DO
	//return (stsolve == STSOLVE_FERROMAGNETIC || (stsolve == STSOLVE_NORMALMETAL && IsNZ(paMesh->SHA.get0())));
	return false;
}

//-------------------Abstract base class method implementations

BError Atom_Transport::Initialize(void)
{
	BError error(CLASS_STR(Atom_Transport));

	//TO DO
	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		//if (!dM_dt.assign(paMesh->h, paMesh->meshRect, DBL3(), paMesh->M)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		dM_dt.clear();
	}

	if (pSMesh->SolveSpinCurrent()) {

		//Poisson equation helper storage allocation
		if (Need_delsq_V_fixed_Precalculation()) {

			//if (!delsq_V_fixed.assign(paMesh->h_e, paMesh->meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_V_fixed.clear();

		if (Need_delsq_S_fixed_Precalculation()) {

			//if (!delsq_S_fixed.assign(paMesh->h_e, paMesh->meshRect, DBL3())) return error(BERROR_OUTOFMEMORY_CRIT);
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

	//TO DO
	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		//if (!dM_dt.assign(paMesh->h, paMesh->meshRect, DBL3(), paMesh->M)) return error(BERROR_OUTOFMEMORY_CRIT);
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
			if (paMesh->elC.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

				success = paMesh->elC.resize(paMesh->h_e, paMesh->meshRect);
			}
			else {

				success = paMesh->elC.assign(paMesh->h_e, paMesh->meshRect, paMesh->elecCond, paMesh->M1);
			}
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (paMesh->V.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

			success &= paMesh->V.resize(paMesh->h_e, paMesh->meshRect, paMesh->elC);
		}
		else {

			success &= paMesh->V.assign(paMesh->h_e, paMesh->meshRect, 0.0, paMesh->elC);
		}

		//electric field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (paMesh->E.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

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

			if (paMesh->S.linear_size() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

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
		pModuleCUDA = new Atom_TransportCUDA(paMesh, pSMesh, this);
		pTransportCUDA = dynamic_cast<Atom_TransportCUDA*>(pModuleCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Transport::UpdateField(void)
{
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be interated only at the end of a step or stage
	if (!paMesh->static_transport_solver && !pSMesh->disabled_transport_solver) {

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
void Atom_Transport::CalculateElectricField(void)
{
	if (stsolve == STSOLVE_NONE || stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

		//calculate electric field using E = - grad V
#pragma omp parallel for
		for (int idx = 0; idx < paMesh->E.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (paMesh->V.is_not_empty(idx)) {

				paMesh->E[idx] = -1.0 * paMesh->V.grad_diri(idx);
			}
			else paMesh->E[idx] = DBL3(0);
		}
	}
}

double Atom_Transport::CalculateElectrodeCurrent(Rect &electrode_rect)
{
	if (!paMesh->meshRect.intersects(electrode_rect)) return 0.0;

	Box electrode_box = paMesh->V.box_from_rect_max(electrode_rect);

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

				current += -paMesh->elC[INT3(electrode_box.s.i - 1, j, k)] *
					(paMesh->V[INT3(electrode_box.s.i - 1, j, k)] - paMesh->V[INT3(electrode_box.s.i - 2, j, k)]) * paMesh->h_e.y * paMesh->h_e.z / paMesh->h_e.x;
			}
		}
	}

	//cells on +x side
	if (electrode_box.e.i + 1 < paMesh->n_e.i) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current -= -paMesh->elC[INT3(electrode_box.e.i, j, k)] *
					(paMesh->V[INT3(electrode_box.e.i + 1, j, k)] - paMesh->V[INT3(electrode_box.e.i, j, k)]) * paMesh->h_e.y * paMesh->h_e.z / paMesh->h_e.x;
			}
		}
	}

	//cells on -y side
	if (electrode_box.s.j > 1) {
#pragma omp parallel for reduction(+:current)
		for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current += -paMesh->elC[INT3(i, electrode_box.s.j - 1, k)] *
					(paMesh->V[INT3(i, electrode_box.s.j - 1, k)] - paMesh->V[INT3(i, electrode_box.s.j - 2, k)]) * paMesh->h_e.x * paMesh->h_e.z / paMesh->h_e.y;
			}
		}
	}

	//cells on +y side
	if (electrode_box.e.j + 1 < paMesh->n_e.j) {
#pragma omp parallel for reduction(+:current)
		for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current -= -paMesh->elC[INT3(i, electrode_box.e.j, k)] *
					(paMesh->V[INT3(i, electrode_box.e.j + 1, k)] - paMesh->V[INT3(i, electrode_box.e.j, k)]) * paMesh->h_e.x * paMesh->h_e.z / paMesh->h_e.y;
			}
		}
	}

	//cells on -z side
	if (electrode_box.s.k > 1) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {

				current += -paMesh->elC[INT3(i, j, electrode_box.s.k - 1)] *
					(paMesh->V[INT3(i, j, electrode_box.s.k - 1)] - paMesh->V[INT3(i, j, electrode_box.s.k - 2)]) * paMesh->h_e.x * paMesh->h_e.y / paMesh->h_e.z;
			}
		}
	}

	//cells on +z side
	if (electrode_box.e.k + 1 < paMesh->n_e.k) {
#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {

				current -= -paMesh->elC[INT3(i, j, electrode_box.e.k)] *
					(paMesh->V[INT3(i, j, electrode_box.e.k + 1)] - paMesh->V[INT3(i, j, electrode_box.e.k)]) * paMesh->h_e.x * paMesh->h_e.y / paMesh->h_e.z;
			}
		}
	}

	return current;
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

	//TO DO
	//Include AMR?
	/*
	if (paMesh->M.linear_size() && (IsNZ((double)paMesh->amrPercentage))) {

		//with amr
#pragma omp parallel for
		for (int idx = 0; idx < paMesh->elC.linear_size(); idx++) {

			if (paMesh->elC.is_not_empty(idx)) {

				double elecCond = paMesh->elecCond;
				double amrPercentage = paMesh->amrPercentage;
				paMesh->update_parameters_ecoarse(idx, paMesh->elecCond, elecCond, paMesh->amrPercentage, amrPercentage);

				//get current density value at this conductivity cell
				DBL3 Jc_value = paMesh->elC[idx] * paMesh->E[idx];

				//get M value (M is on n, h mesh so could be different)
				DBL3 M_value = paMesh->M[paMesh->elC.cellidx_to_position(idx)];

				double magnitude = Jc_value.norm() * M_value.norm();
				double dotproduct = 0.0;

				if (IsNZ(magnitude)) dotproduct = (Jc_value * M_value) / magnitude;

				paMesh->elC[idx] = elecCond / (1 + amrPercentage * dotproduct*dotproduct / 100);
			}
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
	else if (force_recalculate || (paMesh->elecCond.is_tdep() && paMesh->Temp.linear_size())) {
	*/
	if (force_recalculate || (paMesh->elecCond.is_tdep() && paMesh->Temp.linear_size())) {

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

//-------------------Properties

bool Atom_Transport::GInterface_Enabled(void)
{
	return paMesh->GInterface_Enabled();
}

//-------------------Others

//called by MoveMesh method in this mesh - move relevant transport quantities
void Atom_Transport::MoveMesh_Transport(double x_shift)
{
	double mesh_end_size = paMesh->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(paMesh->meshRect.s + DBL3(mesh_end_size, 0, 0), paMesh->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		paMesh->paMeshCUDA->elC()->shift_x(paMesh->elC.linear_size(), x_shift, shift_rect);

		//shift spin accumulation if present
		if (paMesh->S.linear_size()) paMesh->paMeshCUDA->S()->shift_x(paMesh->S.linear_size(), x_shift, shift_rect);

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);

		return;
	}
#endif

	//1. electrical conductivity
	paMesh->elC.shift_x(x_shift, shift_rect);

	//shift spin accumulation if present
	if (paMesh->S.linear_size()) paMesh->S.shift_x(x_shift, shift_rect);

	//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
	pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
}

#endif