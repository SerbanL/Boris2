#include "stdafx.h"
#include "TransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Transport.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMeshCUDA.h"
#include "SuperMesh.h"

TransportCUDA::TransportCUDA(Transport* pTransport_) : 
	ModulesCUDA(),
	TransportBaseCUDA(pTransport_)
{
	pTransport = pTransport_;
	pMesh = pTransport->pMesh;
	pMeshCUDA = pMesh->pMeshCUDA;
	
	Set_STSolveType();
	
	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
	
	//setup objects with methods used for Poisson equation solvers
	
	if (!error_on_create) error_on_create = poisson_V()->set_pointers_transport(pMeshCUDA);
	if (!error_on_create) error_on_create = poisson_Spin_S()->set_pointers_transport(pMeshCUDA, this);
	if (!error_on_create) error_on_create = poisson_Spin_V()->set_pointers_transport(pMeshCUDA, this);
}

TransportCUDA::~TransportCUDA()
{
	//Transport managed quantities must be transferred back to their cpu versions - values only as sizes must match
	if (Holder_Module_Available()) {

		//If holder module still available, this means the cpu version of this module has not been deleted.
		//The only way this could happen is if CUDA is being switched off. 
		//In this case we want to copy over to cpu vecs, but no need to clear memory explicitly, as this will be done in the cu-obj managed destructor when these cuVECs go out of scope.
		//This is done in the CUDA version of Mesh where these quantities are held.
	}
	else {

		//Holder module not available. This means this module has been deleted entirely, but CUDA must still be switched on.
		//In this case free up GPU memory as these cuVECs will not be going out of scope, but in any case they're not needed anymore.
		pMeshCUDA->elC()->clear();
		pMeshCUDA->V()->clear();
		pMeshCUDA->E()->clear();
		pMeshCUDA->S()->clear();
	}
}

//check if dM_dt Calculation should be enabled
bool TransportCUDA::Need_dM_dt_Calculation(void)
{
	return pTransport->Need_dM_dt_Calculation();
}

//check if the delsq_V_fixed VEC is needed
bool TransportCUDA::Need_delsq_V_fixed_Precalculation(void)
{
	return pTransport->Need_delsq_V_fixed_Precalculation();
}

//check if the delsq_S_fixed VEC is needed
bool TransportCUDA::Need_delsq_S_fixed_Precalculation(void)
{
	return pTransport->Need_delsq_S_fixed_Precalculation();
}

//-------------------Abstract base class method implementations

BError TransportCUDA::Initialize(void)
{
	BError error(CLASS_STR(TransportCUDA));

	//no energy density contribution here
	ZeroEnergy();

	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt()->assign(pMeshCUDA->h, pMeshCUDA->meshRect, cuReal3(), (cuVEC_VC<cuReal3>&)pMeshCUDA->M)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		dM_dt()->clear();
	}

	if (pSMesh->SolveSpinCurrent()) {

		//Poisson equation helper storage allocation
		if (Need_delsq_V_fixed_Precalculation()) {

			if (!delsq_V_fixed()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_V_fixed()->clear();

		if (Need_delsq_S_fixed_Precalculation()) {

			if (!delsq_S_fixed()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, cuReal3())) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_S_fixed()->clear();
	}
	else {

		//Poisson equation helper storage not needed
		delsq_V_fixed()->clear();
		delsq_S_fixed()->clear();
	}

	initialized = true;

	return error;
}

BError TransportCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(TransportCUDA));

	Uninitialize();

	bool success = true;
	
	//make sure correct memory is assigned for electrical quantities
	
	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt()->assign(pMeshCUDA->h, pMeshCUDA->meshRect, cuReal3(), (cuVEC_VC<cuReal3>&)pMeshCUDA->M)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		if (!pTransport->dM_dt.linear_size()) dM_dt()->clear();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		if (pMeshCUDA->elC()->size_cpu().dim() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

			success &= pMeshCUDA->elC()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect);
		}
		else {

			success &= pMeshCUDA->elC()->set_from_cpuvec(pMesh->elC);
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMeshCUDA->V()->size_cpu().dim() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

			success &= pMeshCUDA->V()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
		}
		else {

			success &= pMeshCUDA->V()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuBReal)0.0, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
		}
		
		//electrical field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMeshCUDA->E()->size_cpu().dim() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

			success &= pMeshCUDA->E()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
		}
		else {

			success &= pMeshCUDA->E()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, cuReal3(0.0), (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
		}
	}
	
	//here also need to update on ODE change as SolveSpinCurrent state might change
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		//spin accumulation if needed - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (success && stsolve != STSOLVE_NONE) {

			if (pMeshCUDA->S()->size_cpu().dim() && cfgMessage != UPDATECONFIG_MESHSHAPECHANGE) {

				success &= pMeshCUDA->S()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
			}
			else {

				success &= pMeshCUDA->S()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, cuReal3(0.0), (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
			}
		}
		else {

			pMeshCUDA->S()->clear();
			displayVEC()->clear();
		}
	}
	
	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

void TransportCUDA::UpdateField(void)
{
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be interated only at the end of a step or stage
	if (!pMesh->static_transport_solver && !pSMesh->disabled_transport_solver) {

		//update elC (AMR and temperature)
		if (pSMesh->CurrentTimeStepSolved())
			CalculateElectricalConductivity();
	}
}

//-------------------Public calculation Methods

void TransportCUDA::CalculateElectricalConductivity(bool force_recalculate)
{
	//Include AMR?
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