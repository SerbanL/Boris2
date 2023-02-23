#include "stdafx.h"
#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_Transport.h"
#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "SuperMeshCUDA.h"
#include "SuperMesh.h"

Atom_TransportCUDA::Atom_TransportCUDA(Atom_Transport* pTransport_) : 
	ModulesCUDA(),
	TransportBaseCUDA(pTransport_)
{
	pTransport = pTransport_;
	paMesh = pTransport->paMesh;
	paMeshCUDA = paMesh->paMeshCUDA;

	Set_STSolveType();

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
	
	//setup objects with methods used for Poisson equation solvers
	
	if (!error_on_create) error_on_create = poisson_V()->set_pointers_atomtransport(paMeshCUDA);
	if (!error_on_create) error_on_create = poisson_Spin_S()->set_pointers_atomtransport(paMeshCUDA, this);
	if (!error_on_create) error_on_create = poisson_Spin_V()->set_pointers_atomtransport(paMeshCUDA, this);
}

Atom_TransportCUDA::~Atom_TransportCUDA()
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
		paMeshCUDA->elC()->clear();
		paMeshCUDA->V()->clear();
		paMeshCUDA->E()->clear();
		paMeshCUDA->S()->clear();
	}
}

//------------------Others

//check if dM_dt Calculation should be enabled
bool Atom_TransportCUDA::Need_dM_dt_Calculation(void)
{
	return pTransport->Need_dM_dt_Calculation();
}

//check if the delsq_V_fixed VEC is needed
bool Atom_TransportCUDA::Need_delsq_V_fixed_Precalculation(void)
{
	return pTransport->Need_delsq_V_fixed_Precalculation();
}

//check if the delsq_S_fixed VEC is needed
bool Atom_TransportCUDA::Need_delsq_S_fixed_Precalculation(void)
{
	return pTransport->Need_delsq_S_fixed_Precalculation();
}

//-------------------Abstract base class method implementations

BError Atom_TransportCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_TransportCUDA));

	//no energy density contribution here
	ZeroEnergy();

	//is this a "thermoelectric mesh"?
	is_thermoelectric_mesh = (paMesh->Temp.linear_size() && IsNZ(paMesh->Sc.get0()));
	poisson_V()->set_thermoelectric_mesh_flag(is_thermoelectric_mesh);
	poisson_V()->set_open_potential_flag(pTransport->IsOpenPotential());

	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt()->assign(paMeshCUDA->h, paMeshCUDA->meshRect, cuReal3(), (cuVEC_VC<cuReal3>&)paMeshCUDA->M1)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		dM_dt()->clear();
	}

	if (pSMesh->SolveSpinCurrent()) {

		//Poisson equation helper storage allocation
		if (Need_delsq_V_fixed_Precalculation()) {

			if (!delsq_V_fixed()->assign(paMeshCUDA->h_e, paMeshCUDA->meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else delsq_V_fixed()->clear();

		if (Need_delsq_S_fixed_Precalculation()) {

			if (!delsq_S_fixed()->assign(paMeshCUDA->h_e, paMeshCUDA->meshRect, cuReal3())) return error(BERROR_OUTOFMEMORY_CRIT);
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

BError Atom_TransportCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_TransportCUDA));

	Uninitialize();

	bool success = true;
	
	//make sure correct memory is assigned for electrical quantities
	
	//do we need to enable dM_dt calculation?
	if (Need_dM_dt_Calculation()) {

		if (!dM_dt()->assign(paMeshCUDA->h, paMeshCUDA->meshRect, cuReal3(), (cuVEC_VC<cuReal3>&)paMeshCUDA->M1)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//disabled
		if (!pTransport->dM_dt.linear_size()) dM_dt()->clear();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		if (paMeshCUDA->elC()->size_cpu().dim()) {

			success &= paMeshCUDA->elC()->resize(paMeshCUDA->h_e, paMeshCUDA->meshRect);
		}
		else {

			success &= paMeshCUDA->elC()->set_from_cpuvec(paMesh->elC);
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (paMeshCUDA->V()->size_cpu().dim()) {

			success &= paMeshCUDA->V()->resize(paMeshCUDA->h_e, paMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);
		}
		else {

			success &= paMeshCUDA->V()->assign(paMeshCUDA->h_e, paMeshCUDA->meshRect, (cuBReal)0.0, (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);
		}
		
		//electrical field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (paMeshCUDA->E()->size_cpu().dim()) {

			success &= paMeshCUDA->E()->resize(paMeshCUDA->h_e, paMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);
		}
		else {

			success &= paMeshCUDA->E()->assign(paMeshCUDA->h_e, paMeshCUDA->meshRect, cuReal3(0.0), (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);
		}
	}
	
	//here also need to update on ODE change as SolveSpinCurrent state might change
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		//spin accumulation if needed - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (success && stsolve != STSOLVE_NONE) {

			if (paMeshCUDA->S()->size_cpu().dim()) {

				success &= paMeshCUDA->S()->resize(paMeshCUDA->h_e, paMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);
			}
			else {

				success &= paMeshCUDA->S()->assign(paMeshCUDA->h_e, paMeshCUDA->meshRect, cuReal3(0.0), (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);
			}
		}
		else {

			paMeshCUDA->S()->clear();
			displayVEC()->clear();
		}
	}
	
	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

void Atom_TransportCUDA::UpdateField(void)
{
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be interated only at the end of a step or stage
	if (!paMesh->static_transport_solver && !pSMesh->disabled_transport_solver) {

		//update elC (AMR and temperature)
		if (pSMesh->CurrentTimeStepSolved()) CalculateElectricalConductivity();
	}
}

//-------------------Public calculation Methods

void Atom_TransportCUDA::CalculateElectricalConductivity(bool force_recalculate)
{
	//Include AMR?
	if (paMesh->M1.linear_size() && (IsNZ((double)paMesh->amrPercentage) || IsNZ((double)paMesh->tamrPercentage))) {

		if (IsZ((double)paMesh->tamrPercentage)) {

			//only amr
			CalculateElectricalConductivity_AMR();
		}
		else if (IsZ((double)paMesh->amrPercentage)) {

			//only tamr
			CalculateElectricalConductivity_TAMR();
		}
		else {

			//both amr and tamr
			CalculateElectricalConductivity_TAMR_and_AMR();
		}

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
	else if (force_recalculate || (paMesh->elecCond.is_tdep() && paMesh->Temp.linear_size())) {

		//no amr but recalculation is being forced, which may include a temperature dependence with non-uniform temperature
		CalculateElectricalConductivity_NoAMR();

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
}

#endif

#endif