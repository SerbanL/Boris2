#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

Transport::Transport(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	pSMesh = pMesh->pSMesh;

	error_on_create = UpdateConfiguration();

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

	pMesh->Jc.clear();

	pMesh->S.clear();
}

bool Transport::SetFixedPotentialCells(Rect rectangle, double potential)
{ 
	bool success = true;

#if COMPILECUDA == 1
	if (pModuleCUDA) success = reinterpret_cast<TransportCUDA*>(pModuleCUDA)->SetFixedPotentialCells(rectangle, potential);
#endif

	success &= pMesh->V.set_dirichlet_conditions(rectangle, potential);

	return success;
}

void Transport::ClearFixedPotentialCells(void)
{
	pMesh->V.clear_dirichlet_flags();

#if COMPILECUDA == 1
	if (pModuleCUDA) reinterpret_cast<TransportCUDA*>(pModuleCUDA)->ClearFixedPotentialCells();
#endif
}

//-------------------Abstract base class method implementations

BError Transport::Initialize(void)
{
	BError error(CLASS_STR(Transport));

	if (!initialized) {

		initialized = true;

		CalculateElectricalConductivity(true);
	}

	return error;
}

BError Transport::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Transport));

	Uninitialize();

	bool success = true;

	//make sure the cellsize divides the mesh rectangle
	pMesh->n_e = round(pMesh->meshRect / pMesh->h_e);
	if (pMesh->n_e.x < 2) pMesh->n_e.x = 2;
	if (pMesh->n_e.y < 2) pMesh->n_e.y = 2;
	if (pMesh->n_e.z < 2) pMesh->n_e.z = 2;
	pMesh->h_e = pMesh->meshRect / pMesh->n_e;

	//make sure correct memory is assigned for electrical quantities

	//electrical conductivity
	if (pMesh->M.linear_size()) {

		//for ferromagnetic meshes set empty cells using information in M (empty cells have zero electrical conductivity) only on initialization. If already initialized then shape already set.
		if(pMesh->elC.linear_size()) 
			success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
		else 
			success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, pMesh->elecCond, pMesh->M);
	}
	else {

		//for non-ferromagnetic meshes (e.g. simple electrical conductor) then elC contains directly any empty cells information
		if(pMesh->elC.linear_size()) 
			success = pMesh->elC.resize(pMesh->h_e, pMesh->meshRect);
		else 
			success = pMesh->elC.assign(pMesh->h_e, pMesh->meshRect, pMesh->elecCond);
	}

	//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
	if (pMesh->V.linear_size())
		success = pMesh->V.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
	else
		success = pMesh->V.assign(pMesh->h_e, pMesh->meshRect, 0.0, pMesh->elC);
	
	//electrical current density - set empty cells using information in elC (empty cells have zero electrical conductivity)
	if (pMesh->Jc.linear_size())
		success = pMesh->Jc.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
	else
		success = pMesh->Jc.assign(pMesh->h_e, pMesh->meshRect, DBL3(0.0), pMesh->elC);

	//spin accumulation if needed - set empty cells using information in elC (empty cells have zero electrical conductivity)
	if (success && pSMesh->SolveSpinCurrent()) {

		if (pMesh->S.linear_size())
			success = pMesh->S.resize(pMesh->h_e, pMesh->meshRect, pMesh->elC);
		else
			success = pMesh->S.assign(pMesh->h_e, pMesh->meshRect, DBL3(0.0), pMesh->elC);
	}
	else {

		pMesh->S.clear();

		//clear display VEC used for spin currents and torque
		displayVEC.clear();
	}

	if (!success) return error(BERROR_OUTOFMEMORY_CRIT);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
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
		pModuleCUDA = new TransportCUDA(pMesh, pSMesh);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Transport::UpdateField(void)
{	
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be interated only at the end of a step or stage
	if (!pMesh->static_transport_solver) {

		//update elC (AMR and temperature)
		if (pSMesh->CurrentTimeStepSolved())
			CalculateElectricalConductivity();
	}

	//no contribution to total energy density
	return 0.0;
}

//-------------------Calculation Methods

void Transport::CalculateCurrentDensity(void)
{
	if (!pSMesh->SolveSpinCurrent()) {

		//calculate current density using Jc = -sigma * grad V
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Jc.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (pMesh->V.is_not_empty(idx)) {

				pMesh->Jc[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri(idx);
			}
			else pMesh->Jc[idx] = DBL3(0);
		}
	}
	else {

		//Current density when contributions from spin accumulation are present
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Jc.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (pMesh->V.is_not_empty(idx)) {

				double iSHA = pMesh->iSHA;
				double SHA = pMesh->SHA;
				double De = pMesh->De;
				double betaD = pMesh->betaD;
				pMesh->update_parameters_ecoarse(idx, pMesh->iSHA, iSHA, pMesh->SHA, SHA, pMesh->De, De, pMesh->betaD, betaD);

				//1. Ohm's law contribution
				if (IsZ((double)iSHA) || pMesh->M.linear_size()) {

					//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
					pMesh->Jc[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri(idx);
				}
				else {

					//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
					pMesh->Jc[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri_nneu(idx, (iSHA * De / (MUB_E * pMesh->elC[idx])) * pMesh->S.curl_neu(idx));

					//must also add iSHE contribution -> here we must use non-homogeneous Neumann boundary conditions when calculating S differentials
					pMesh->Jc[idx] += (iSHA * De / MUB_E) * pMesh->S.curl_nneu(idx, epsilon3(pMesh->Jc[idx]) * (SHA * MUB_E / De));
				}

				//2. CPP-GMR contribution
				if (pMesh->M.linear_size() && IsNZ((double)betaD)) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_ecoarse(idx, pMesh->Ms, Ms);

					int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(idx));

					DBL3 M = pMesh->M[idx_M];
					DBL33 grad_S = pMesh->S.grad_neu(idx);		//homogeneous Neumann since SHA = 0 in magnetic meshes
						
					pMesh->Jc[idx] += (grad_S * M) * betaD * De / (MUB_E * Ms);
				}
			}
			else pMesh->Jc[idx] = DBL3(0);
		}
	}
}

double Transport::CalculateElectrodeCurrent(Rect &electrode_rect)
{
	if (!pMesh->meshRect.intersects(electrode_rect)) return 0.0;

	Box electrode_box = pMesh->V.box_from_rect_max(electrode_rect);

#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<TransportCUDA*>(pModuleCUDA)->CalculateElectrodeCurrent(electrode_box);
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

		reinterpret_cast<TransportCUDA*>(pModuleCUDA)->CalculateElectricalConductivity(force_recalculate);
		return;
	}
#endif

	if (pMesh->M.linear_size() && (IsNZ((double)pMesh->amrPercentage))) {

		//with amr
		#pragma omp parallel for
		for (int idx = 0; idx < pMesh->elC.linear_size(); idx++) {

			if (pMesh->elC.is_not_empty(idx)) {

				double elecCond = pMesh->elecCond;
				double amrPercentage = pMesh->amrPercentage;
				pMesh->update_parameters_ecoarse(idx, pMesh->elecCond, elecCond, pMesh->amrPercentage, amrPercentage);

				//get current density value at this conductivity cell
				DBL3 Jc_value = pMesh->Jc[idx];

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