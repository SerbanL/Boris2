#include "stdafx.h"
#include "TMRCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TMR

#include "TMR.h"
#include "MeshCUDA.h"
#include "Atom_MeshCUDA.h"
#include "Mesh.h"
#include "SuperMeshCUDA.h"
#include "SuperMesh.h"

//--------------- UNIAXIAL

TMRCUDA::TMRCUDA(TMR* pTMR_) : 
	ModulesCUDA(),
	TransportBaseCUDA(pTMR_)
{
	pTMR = pTMR_;
	pMesh = pTMR->pMesh;
	pMeshCUDA = pMesh->pMeshCUDA;
	
	Set_STSolveType();

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//setup objects with methods used for Poisson equation solvers

	if (!error_on_create) error_on_create = poisson_V()->set_pointers_transport(pMeshCUDA);
	if (!error_on_create) error_on_create = poisson_Spin_S()->set_pointers_transport(pMeshCUDA, this);
	if (!error_on_create) error_on_create = poisson_Spin_V()->set_pointers_transport(pMeshCUDA, this);
}

TMRCUDA::~TMRCUDA()
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

//------------------Others

//check if dM_dt Calculation should be enabled
bool TMRCUDA::Need_dM_dt_Calculation(void)
{
	return pTMR->Need_dM_dt_Calculation();
}

//check if the delsq_V_fixed VEC is needed
bool TMRCUDA::Need_delsq_V_fixed_Precalculation(void)
{
	return pTMR->Need_delsq_V_fixed_Precalculation();
}

//check if the delsq_S_fixed VEC is needed
bool TMRCUDA::Need_delsq_S_fixed_Precalculation(void)
{
	return pTMR->Need_delsq_S_fixed_Precalculation();
}

//-------------------Abstract base class method implementations

BError TMRCUDA::Initialize(void)
{
	BError error(CLASS_STR(TMRCUDA));

	//no energy density contribution here
	ZeroEnergy();

	//clear cu_arrs then rebuild them from information in TMR module
	pMeshFM_Bot.clear();
	pMeshFM_Top.clear();
	pMeshAtom_Bot.clear();
	pMeshAtom_Top.clear();

	//dM_dt calculation not needed
	dM_dt()->clear();

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
	
	//make sure information in TMR module is up to date
	error = pTMR->Initialize();

	//now update flags information in gpu versions
	if (!pMeshCUDA->elC()->copyflags_from_cpuvec(pMesh->elC)) error(BERROR_GPUERROR_CRIT);
	if (!pMeshCUDA->V()->copyflags_from_cpuvec(pMesh->V)) error(BERROR_GPUERROR_CRIT);
	if (!pMeshCUDA->E()->copyflags_from_cpuvec(pMesh->E)) error(BERROR_GPUERROR_CRIT);
	if (pMesh->S.linear_size()) {

		if (!pMeshCUDA->S()->copyflags_from_cpuvec(pMesh->S)) error(BERROR_GPUERROR_CRIT);
	}

	if (!error) {

		for (int idx = 0; idx < pTMR->pMeshFM_Bot.size(); idx++) {

			pMeshFM_Bot.push_back(pTMR->pMeshFM_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
		}

		for (int idx = 0; idx < pTMR->pMeshFM_Top.size(); idx++) {

			pMeshFM_Top.push_back(pTMR->pMeshFM_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
		}

		for (int idx = 0; idx < pTMR->pMeshAtom_Bot.size(); idx++) {

			pMeshAtom_Bot.push_back(pTMR->pMeshAtom_Bot[idx]->paMeshCUDA->cuaMesh.get_managed_object());
		}

		for (int idx = 0; idx < pTMR->pMeshAtom_Top.size(); idx++) {

			pMeshAtom_Top.push_back(pTMR->pMeshAtom_Top[idx]->paMeshCUDA->cuaMesh.get_managed_object());
		}

		initialized = true;
	}

	CalculateElectricalConductivity(true);

	return error;
}

BError TMRCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(TMRCUDA));

	Uninitialize();

	bool success = true;

	//make sure correct memory is assigned for electrical quantities

	//dM_dt calculation not needed
	dM_dt()->clear();

	//need this when we switch cuda mode
	if (!RAtmr_p_equation.is_set() && pTMR->RAtmr_p_equation.is_set()) error = SetBiasEquation_Parallel(pTMR->RAtmr_p_equation.get_scalar_fspec());
	if (!RAtmr_ap_equation.is_set() && pTMR->RAtmr_ap_equation.is_set()) error = SetBiasEquation_AntiParallel(pTMR->RAtmr_ap_equation.get_scalar_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		if (pMeshCUDA->elC()->size_cpu().dim()) {

			success &= pMeshCUDA->elC()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect);
		}
		else {

			success &= pMeshCUDA->elC()->set_from_cpuvec(pMesh->elC);
		}

		//electrical potential - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMeshCUDA->V()->size_cpu().dim()) {

			success &= pMeshCUDA->V()->resize(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
		}
		else {

			success &= pMeshCUDA->V()->assign(pMeshCUDA->h_e, pMeshCUDA->meshRect, (cuBReal)0.0, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
		}

		//electrical field - set empty cells using information in elC (empty cells have zero electrical conductivity)
		if (pMeshCUDA->E()->size_cpu().dim()) {

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

			if (pMeshCUDA->S()->size_cpu().dim()) {

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

void TMRCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		if (RAtmr_p_equation.is_set()) {

			SetBiasEquation_Parallel(pTMR->RAtmr_p_equation.get_scalar_fspec());
		}

		if (RAtmr_ap_equation.is_set()) {

			SetBiasEquation_AntiParallel(pTMR->RAtmr_ap_equation.get_scalar_fspec());
		}
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (RAtmr_p_equation.is_set()) RAtmr_p_equation.clear();
		if (RAtmr_ap_equation.is_set()) RAtmr_ap_equation.clear();
	}
}

void TMRCUDA::UpdateField(void)
{
	if (pSMesh->CurrentTimeStepSolved()) CalculateElectricalConductivity(true);
}

//-------------------Public calculation Methods

void TMRCUDA::CalculateElectricalConductivity(bool force_recalculate)
{
	if (force_recalculate) {

		CalculateElectricalConductivity_TMR(pTMR->TMR_type);

		//set flag to force transport solver recalculation (note, the SuperMesh modules always update after the Mesh modules) - elC has changed
		pSMesh->CallModuleMethod(&STransport::Flag_Recalculate_Transport);
	}
}

//calculate resistance in given rectangle (relative to mesh)
cuBReal TMRCUDA::GetResistance(cuRect rect)
{
	cuBReal conductance = pMeshCUDA->elC()->sum_nonempty(pMeshCUDA->n_e.dim(), rect) * pMeshCUDA->h_e.x * pMeshCUDA->h_e.y / (pMeshCUDA->meshRect.height() * pMeshCUDA->n_e.k);

	//now return resistance
	if (conductance) return 1.0 / conductance;
	else return 0.0;
}

#endif

#endif