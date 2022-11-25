#include "stdafx.h"
#include "TransportBaseCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "SuperMesh.h"
#include "MeshBase.h"
#include "TransportBase.h"

TransportBaseCUDA::TransportBaseCUDA(TransportBase* pTransportBase_) :
	pTransportBase(pTransportBase_)
{
	pMeshBaseCUDA = pTransportBase->pMeshBase->pMeshBaseCUDA;
	pSMesh = pTransportBase->pSMesh;
	pSMeshCUDA = pSMesh->pSMeshCUDA;
}

TransportBaseCUDA::~TransportBaseCUDA()
{
	if (pTransportBase) pTransportBase->pTransportBaseCUDA = nullptr;
}

//-------------------Auxiliary

//set the stsolve indicator depending on current configuration
void TransportBaseCUDA::Set_STSolveType(void)
{
	stsolve = pTransportBase->stsolve;

	poisson_Spin_V()->set_stsolve(stsolve);
	poisson_Spin_S()->set_stsolve(stsolve);
}

//prepare displayVEC ready for calculation of display quantity
bool TransportBaseCUDA::PrepareDisplayVEC(DBL3 cellsize)
{
	if (pSMesh->SolveSpinCurrent() && pMeshBaseCUDA->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC()->assign(cellsize, pMeshBaseCUDA->meshRect, cuReal3(0.0));

		return true;
	}
	else displayVEC()->clear();

	return false;
}

//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
bool TransportBaseCUDA::PrepareDisplayVEC_VC(DBL3 cellsize)
{
	if (pMeshBaseCUDA->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC_VC()->assign(cellsize, pMeshBaseCUDA->meshRect, cuReal3(0.0));

		return true;
	}
	else displayVEC_VC()->clear();

	return false;
}

//-------------------TAMR

BError TransportBaseCUDA::Set_TAMR_Conductivity_Equation(std::vector<std::vector< EqComp::FSPEC >> fspec)
{
	BError error(CLASS_STR(TransportBaseCUDA));

	if (!TAMR_conductivity_equation.make_scalar(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

void TransportBaseCUDA::Clear_TAMR_Conductivity_Equation(void)
{
	TAMR_conductivity_equation.clear();
}

//------------------Others

//set fixed potential cells in this mesh for given rectangle
bool TransportBaseCUDA::SetFixedPotentialCells(cuRect rectangle, cuBReal potential)
{
	return pMeshBaseCUDA->V()->set_dirichlet_conditions(rectangle, potential);
}

void TransportBaseCUDA::ClearFixedPotentialCells(void)
{
	pMeshBaseCUDA->V()->clear_dirichlet_flags();
}

void TransportBaseCUDA::Set_Linear_PotentialDrop(cuRect electrode1, cuBReal potential1, cuRect electrode2, cuBReal potential2, cuReal2 degeneracy)
{
	pMeshBaseCUDA->V()->set_linear(electrode1, potential1, electrode2, potential2, degeneracy);
}

//-------------------Properties

bool TransportBaseCUDA::GInterface_Enabled(void)
{
	return pTransportBase->GInterface_Enabled();
}

bool TransportBaseCUDA::iSHA_nonzero(void)
{
	return pTransportBase->iSHA_nonzero();
}

bool TransportBaseCUDA::SHA_nonzero(void)
{
	return pTransportBase->SHA_nonzero();
}

#endif

#endif