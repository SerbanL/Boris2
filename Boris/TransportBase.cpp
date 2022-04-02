#include "stdafx.h"
#include "TransportBase.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "SuperMesh.h"

TransportBase::TransportBase(MeshBase* pMeshBase_) :
	pMeshBase(pMeshBase_)
{
	pSMesh = pMeshBase->pSMesh;
}

//-------------------Auxiliary

//set the stsolve indicator depending on current configuration
void TransportBase::Set_STSolveType(void)
{
	if (!pSMesh->SolveSpinCurrent()) {

		//no spin transport required
		stsolve = STSOLVE_NONE;
	}
	else {

		switch (pMeshBase->GetMeshType()) {

		case MESH_FERROMAGNETIC:
		case MESH_DIPOLE:
			stsolve = STSOLVE_FERROMAGNETIC;
			break;

		case MESH_ATOM_CUBIC:
			stsolve = STSOLVE_FERROMAGNETIC_ATOM;
			break;

		case MESH_ANTIFERROMAGNETIC:
			//to be implemented. currently treated as NM
			stsolve = STSOLVE_NORMALMETAL;
			break;

		case MESH_METAL:
			stsolve = STSOLVE_NORMALMETAL;
			break;

		case MESH_INSULATOR:
			stsolve = STSOLVE_TUNNELING;
			break;
		}
	}

#if COMPILECUDA == 1
	//copy over to TransportBaseCUDA module - need to use it in .cu files
	if (pTransportBaseCUDA) pTransportBaseCUDA->Set_STSolveType();
#endif
}

//prepare displayVEC ready for calculation of display quantity
bool TransportBase::PrepareDisplayVEC(DBL3 cellsize)
{
	if (pSMesh->SolveSpinCurrent() && pMeshBase->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC.assign(cellsize, pMeshBase->meshRect, DBL3(0.0));

		return true;
	}
	else displayVEC.clear();

	return false;
}

//prepare displayVEC_VC ready for calculation of display quantity
bool TransportBase::PrepareDisplayVEC_VC(DBL3 cellsize)
{
	if (pMeshBase->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC_VC.assign(cellsize, pMeshBase->meshRect, DBL3(0.0));

		return true;
	}
	else displayVEC_VC.clear();

	return false;
}

//-------------------Properties

bool TransportBase::GInterface_Enabled(void)
{
	return pMeshBase->GInterface_Enabled();
}

bool TransportBase::iSHA_nonzero(void)
{
	return pMeshBase->iSHA_nonzero();
}

bool TransportBase::SHA_nonzero(void)
{
	return pMeshBase->SHA_nonzero();
}

//------------------Others

bool TransportBase::SetFixedPotentialCells(Rect rectangle, double potential)
{
	bool success = true;

#if COMPILECUDA == 1
	if (pTransportBaseCUDA) success &= pTransportBaseCUDA->SetFixedPotentialCells(rectangle, potential);
#endif

	success &= pMeshBase->V.set_dirichlet_conditions(rectangle, potential);

	return success;
}

void TransportBase::ClearFixedPotentialCells(void)
{
	pMeshBase->V.clear_dirichlet_flags();

#if COMPILECUDA == 1
	if (pTransportBaseCUDA) pTransportBaseCUDA->ClearFixedPotentialCells();
#endif
}

void TransportBase::Set_Linear_PotentialDrop(Rect electrode1, double potential1, Rect electrode2, double potential2, DBL2 degeneracy)
{
#if COMPILECUDA == 1
	if (pTransportBaseCUDA) {

		pTransportBaseCUDA->Set_Linear_PotentialDrop(electrode1, potential1, electrode2, potential2, degeneracy);
		return;
	}
#endif

	pMeshBase->V.set_linear(electrode1, potential1, electrode2, potential2, degeneracy);
}

//called by MoveMesh method in this mesh - move relevant transport quantities
void TransportBase::MoveMesh_Transport(double x_shift)
{
	double mesh_end_size = pMeshBase->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(pMeshBase->meshRect.s + DBL3(mesh_end_size, 0, 0), pMeshBase->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pMeshBase->pMeshBaseCUDA) {

		pMeshBase->pMeshBaseCUDA->elC()->shift_x(pMeshBase->elC.linear_size(), x_shift, shift_rect);

		//shift spin accumulation if present
		if (pMeshBase->S.linear_size()) pMeshBase->pMeshBaseCUDA->S()->shift_x(pMeshBase->S.linear_size(), x_shift, shift_rect);

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);

		return;
	}
#endif

	//1. electrical conductivity
	pMeshBase->elC.shift_x(x_shift, shift_rect);

	//shift spin accumulation if present
	if (pMeshBase->S.linear_size()) pMeshBase->S.shift_x(x_shift, shift_rect);

	//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
	pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
}

//-------------------Calculation Methods

double TransportBase::CalculateElectrodeCurrent(Rect &electrode_rect)
{
	if (!pMeshBase->meshRect.intersects(electrode_rect)) return 0.0;

	Box electrode_box = pMeshBase->V.box_from_rect_max(electrode_rect);

	//unit vector normal to electrode
	DBL3 ndir = electrode_rect.get_normal();
	//position of electrode_box start cell in absolute coordinates
	DBL3 bpos = ((electrode_box.s + DBL3(0.5)) & pMeshBase->h_e) + pMeshBase->meshRect.s;

	//signs for different directions : 0 of direction not applicable, +1 for cells on -ve side of electrode, -1 for cells on +ve side of electrode
	INT3 sign = INT3(
		(fabs(ndir.x) > 0.0 ? (bpos.x < electrode_rect.s.x ? 1 : -1) : 0),
		(fabs(ndir.y) > 0.0 ? (bpos.y < electrode_rect.s.y ? 1 : -1) : 0),
		(fabs(ndir.z) > 0.0 ? (bpos.z < electrode_rect.s.z ? 1 : -1) : 0));

#if COMPILECUDA == 1
	if (pTransportBaseCUDA) return pTransportBaseCUDA->CalculateElectrodeCurrent(electrode_box, sign);
#endif

	//look at all cells surrounding the electrode

	double current = 0;

	if (sign.x) {

#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current += sign.x * pMeshBase->elC[INT3(electrode_box.s.i, j, k)] * pMeshBase->E[INT3(electrode_box.s.i, j, k)].x * pMeshBase->h_e.y * pMeshBase->h_e.z;
			}
		}
	}

	if (sign.y) {

#pragma omp parallel for reduction(+:current)
		for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {
			for (int k = electrode_box.s.k; k < electrode_box.e.k; k++) {

				current += sign.y * pMeshBase->elC[INT3(i, electrode_box.s.j, k)] * pMeshBase->E[INT3(i, electrode_box.s.j, k)].y * pMeshBase->h_e.x * pMeshBase->h_e.z;
			}
		}
	}

	if (sign.z) {

#pragma omp parallel for reduction(+:current)
		for (int j = electrode_box.s.j; j < electrode_box.e.j; j++) {
			for (int i = electrode_box.s.i; i < electrode_box.e.i; i++) {

				current += sign.z * pMeshBase->elC[INT3(i, j, electrode_box.s.k)] * pMeshBase->E[INT3(i, j, electrode_box.s.k)].z * pMeshBase->h_e.x * pMeshBase->h_e.y;
			}
		}
	}

	return current;
}

//-------------------Display Calculation Methods

DBL3 TransportBase::GetAverageChargeCurrent(Rect rectangle, std::vector<MeshShape> shapes)
{
#if COMPILECUDA == 1
	if (pTransportBaseCUDA) {

		//update displayVEC_VC with charge current in TransportCUDA
		GetChargeCurrentCUDA();

		//average charge current in displayVEC_VC in TransportCUDA
		return cuReal3(pTransportBaseCUDA->displayVEC_VC()->average_nonempty(pMeshBase->n_e.dim(), rectangle));
	}
#endif

	//update displayVEC_VC with charge current
	GetChargeCurrent();

	//average charge current in displayVEC_VC
	if (!shapes.size()) return displayVEC_VC.average_nonempty_omp(rectangle);
	else return displayVEC_VC.shape_getaverage(shapes);
}

DBL3 TransportBase::GetAverageSpinCurrent(int component, Rect rectangle, std::vector<MeshShape> shapes)
{
#if COMPILECUDA == 1
	if (pTransportBaseCUDA) {

		//update displayVEC with spin current in TransportCUDA
		GetSpinCurrentCUDA(component);

		//average spin current in displayVEC in TransportCUDA
		if (pSMesh->SolveSpinCurrent()) return cuReal3(pTransportBaseCUDA->displayVEC()->average_nonempty(pMeshBase->n_e.dim(), rectangle));
		else return cuReal3(0.0);
	}
#endif

	//update displayVEC with spin current
	GetSpinCurrent(component);

	//average spin current in displayVEC
	if (!shapes.size()) return displayVEC.average_nonempty_omp(rectangle);
	else return displayVEC.shape_getaverage(shapes);
}

//Get average bulk spin torque in given rectangle - calculate it first
DBL3 TransportBase::GetAverageSpinTorque(Rect rectangle, std::vector<MeshShape> shapes)
{
	GetSpinTorque();
	if (displayVEC.linear_size()) {

		if (!shapes.size()) return displayVEC.average_nonempty_omp(rectangle);
		else return displayVEC.shape_getaverage(shapes);
	}

	return DBL3();
}

//Get average interfacial spin torque in given rectangle - must have been calculated already in displayVEC through the supermesh transport module
DBL3 TransportBase::GetAverageInterfacialSpinTorque(Rect rectangle, std::vector<MeshShape> shapes)
{
	if (displayVEC.linear_size()) {

		if (!shapes.size()) return displayVEC.average_nonempty_omp(rectangle);
		else return displayVEC.shape_getaverage(shapes);
	}

	return DBL3();
}

#if COMPILECUDA == 1

cu_obj<cuVEC_VC<cuReal3>>& TransportBase::GetChargeCurrentCUDA(void)
{
	return pTransportBaseCUDA->GetChargeCurrent();
}

cu_obj<cuVEC<cuReal3>>& TransportBase::GetSpinCurrentCUDA(int component)
{
	return pTransportBaseCUDA->GetSpinCurrent(component);
}

cu_obj<cuVEC<cuReal3>>& TransportBase::GetSpinTorqueCUDA(void)
{
	return pTransportBaseCUDA->GetSpinTorque();
}

#endif

#endif