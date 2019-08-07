#include "stdafx.h"
#include "TransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "Transport.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMeshCUDA.h"
#include "SuperMesh.h"

TransportCUDA::TransportCUDA(Mesh* pMesh_, SuperMesh* pSMesh_)
	: ModulesCUDA()
{
	pMesh = pMesh_;
	pMeshCUDA = pMesh->pMeshCUDA;
	pSMesh = pSMesh_;
	pSMeshCUDA = pSMesh->pSMeshCUDA;

	error_on_create = UpdateConfiguration();

	//setup objects with methods used for Poisson equation solvers
	if (!error_on_create) error_on_create = poisson_V()->set_pointers(pMeshCUDA);
	if (!error_on_create) {

		if (pMesh->MComputation_Enabled()) {

			error_on_create = poisson_Spin_S()->set_pointers(pMeshCUDA, reinterpret_cast<FMesh*>(pMesh)->Get_DifferentialEquation().Get_DifferentialEquationCUDA_ptr());
		}
		else error_on_create = poisson_Spin_S()->set_pointers(pMeshCUDA, nullptr);
	}
	if (!error_on_create) error_on_create = poisson_Spin_V()->set_pointers(pMeshCUDA, poisson_Spin_S);
}

TransportCUDA::~TransportCUDA()
{
	//Transport managed quantities must be transferred back to their cpu versions - values only as sizes must match
	if (Holder_Module_Available()) {

		pMeshCUDA->elC()->copy_to_cpuvec(pMesh->elC);
		pMeshCUDA->V()->copy_to_cpuvec(pMesh->V);
		pMeshCUDA->Jc()->copy_to_cpuvec(pMesh->Jc);
		pMeshCUDA->S()->copy_to_cpuvec(pMesh->S);
	}

	//clear mesh quantities as they're not used any more
	pMeshCUDA->elC()->clear();
	pMeshCUDA->V()->clear();
	pMeshCUDA->Jc()->clear();
	pMeshCUDA->S()->clear();
}

//------------------Others

//set fixed potential cells in this mesh for given rectangle
bool TransportCUDA::SetFixedPotentialCells(cuRect rectangle, cuReal potential)
{
	return pMeshCUDA->V()->set_dirichlet_conditions(rectangle, potential);
}

void TransportCUDA::ClearFixedPotentialCells(void)
{
	pMeshCUDA->V()->clear_dirichlet_flags();
}

//-------------------Abstract base class method implementations

BError TransportCUDA::Initialize(void)
{
	BError error(CLASS_STR(TransportCUDA));

	//no energy density contribution here
	ZeroEnergy();

	if (!initialized) {

		initialized = true;

		CalculateElectricalConductivity(true);
	}

	return error;
}

BError TransportCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(TransportCUDA));

	Uninitialize();

	bool success = true;

	//make sure correct memory is assigned for electrical quantities
	
	if (pMeshCUDA->elC()->size_cpu().dim())
		success = pMeshCUDA->elC()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect);
	else
		success = pMeshCUDA->elC()->set_from_cpuvec(pMesh->elC);

	//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
	if (pMeshCUDA->V()->size_cpu().dim())
		success = pMeshCUDA->V()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuReal>&)pMeshCUDA->elC);
	else
		success = pMeshCUDA->V()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuReal)0.0, (cuVEC_VC<cuReal>&)pMeshCUDA->elC);

	//electrical current density - set empty cells using information in elC (empty cells have zero electrical conductivity)
	if (pMeshCUDA->Jc()->size_cpu().dim())
		success = pMeshCUDA->Jc()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuReal>&)pMeshCUDA->elC);
	else
		success = pMeshCUDA->Jc()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, cuReal3(0.0), (cuVEC_VC<cuReal>&)pMeshCUDA->elC);

	//spin accumulation if needed - set empty cells using information in elC (empty cells have zero electrical conductivity)
	if (success && pSMesh->SolveSpinCurrent()) {

		if (pMeshCUDA->S()->size_cpu().dim())
			success = pMeshCUDA->S()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuReal>&)pMeshCUDA->elC);
		else
			success = pMeshCUDA->S()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, cuReal3(0.0), (cuVEC_VC<cuReal>&)pMeshCUDA->elC);
	}
	else {

		pMeshCUDA->S()->clear();
		displayVEC()->clear();
	}

	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	
	return error;
}

void TransportCUDA::UpdateField(void)
{
	//update elC (AMR and temperature)
	if (pSMesh->CurrentTimeStepSolved())
		CalculateElectricalConductivity();
}

//-------------------Public calculation Methods

void TransportCUDA::CalculateElectricalConductivity(bool force_recalculate)
{
	if (pMesh->M.linear_size() && (IsNZ((double)pMesh->amrPercentage))) {

		CalculateElectricalConductivity_AMR();

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
	else if (force_recalculate || (pMesh->elecCond.is_tdep() && pMesh->Temp.linear_size())) {

		//no amr but recalculation is being forced, which may include a temperature dependence with non-uniform temperature
		CalculateElectricalConductivity_NoAMR();

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
}

#endif

#endif