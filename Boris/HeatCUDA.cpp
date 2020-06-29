#include "stdafx.h"
#include "HeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_HEAT

#include "Heat.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMesh.h"

HeatCUDA::HeatCUDA(Mesh* pMesh_, SuperMesh* pSMesh_, Heat* pHeat_) : 
	ModulesCUDA(),
	HeatBaseCUDA()
{
	pMesh = pMesh_;
	pMeshCUDA = pMesh->pMeshCUDA;
	pSMesh = pSMesh_;
	pHeat = pHeat_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	if (!error_on_create) error_on_create = temp_cmbnd_funcs()->set_pointers(pMeshCUDA);
}

HeatCUDA::~HeatCUDA()
{
	//Transport managed quantities must be transferred back to their cpu versions - values only as sizes must match
	if (Holder_Module_Available()) {

		//If holder module still available, this means the cpu version of this module has not been deleted.
		//The only way this could happen is if CUDA is being switched off. 
		//In this case we want to copy over to cpu vecs, but no need to clear memory explicitly, as this will be done in the cu-obj managed destructor when these cuVECs go out of scope.
		pMeshCUDA->Temp()->copy_to_cpuvec(pMesh->Temp);
		pMeshCUDA->Temp_l()->copy_to_cpuvec(pMesh->Temp_l);
	}
	else {

		//Holder module not available. This means this module has been deleted entirely, but CUDA must still be switched on.
		//In this case free up GPU memory as these cuVECs will not be going out of scope, but in any case they're not needed anymore.
		pMeshCUDA->Temp()->clear();
		pMeshCUDA->Temp_l()->clear();
	}
}


//-------------------Abstract base class method implementations

BError HeatCUDA::Initialize(void)
{
	BError error(CLASS_STR(HeatCUDA));

	//no energy density contribution here
	ZeroEnergy();

	pHeat->SetRobinBoundaryConditions();

	initialized = true;

	return error;
}

BError HeatCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(HeatCUDA));

	Uninitialize();
	
	bool success = true;
	
	//need this when we switch cuda mode
	if (!Q_equation.is_set() && pHeat->Q_equation.is_set()) error = SetQEquation(pHeat->Q_equation.get_scalar_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, UPDATECONFIG_HEAT_MODELTYPE)) {

		//update mesh dimensions in equation constants
		if (Q_equation.is_set()) {

			error = SetQEquation(pHeat->Q_equation.get_scalar_fspec());
		}

		//first clear any VECs which are not needed
		switch (pHeat->tmtype) {

		case TMTYPE_1TM:
			//single temperature only
			pMeshCUDA->Temp_l()->clear();
			break;

		case TMTYPE_2TM:
			//itinerant electrons temperature <-> lattice temperature
			break;
		}

		//if we've just changed the model type, make sure any new temperature VECs start from a reasonable setting - same as the primary Temp VEC.
		bool initialize_Temp_l = !(pMeshCUDA->Temp_l()->linear_size_cpu());

		if (pMeshCUDA->M()->size_cpu().dim()) {

			//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
			if (pMeshCUDA->Temp()->size_cpu().dim()) {

				success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect);

				//lattice temperature for many-temperature models
				if (pHeat->tmtype == TMTYPE_2TM) success &= pMeshCUDA->Temp_l()->resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);

				//lattice temperature for many-temperature models
				if (pHeat->tmtype == TMTYPE_2TM) success &= pMeshCUDA->Temp_l()->set_from_cpuvec(pMesh->Temp_l);
			}
		}
		else if (pMeshCUDA->elC()->size_cpu().dim()) {

			//in a normal metal the electrical conductivity sets the shape

			//before doing this must make sure elC was itself set in the Transport module (it could be this module is being updated before the Transport module)
			//Transport module is guaranteed to be set otherwise elC would have zero size - it does mean Transport has UpdateConfiguration called twice but it doesn't matter.
			if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

				error = (*pMesh)(MOD_TRANSPORT)->UpdateConfiguration(cfgMessage);
				if (error) return error;
			}

			if (success) {

				if (pMeshCUDA->Temp()->size_cpu().dim()) {

					success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);

					//lattice temperature for many-temperature models
					if (pHeat->tmtype == TMTYPE_2TM) success &= pMeshCUDA->Temp_l()->resize(pMesh->h_t, pMesh->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->Temp);
				}
				else {

					success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);

					//lattice temperature for many-temperature models
					if (pHeat->tmtype == TMTYPE_2TM) success &= pMeshCUDA->Temp_l()->set_from_cpuvec(pMesh->Temp_l);
				}
			}
		}
		else {

			//in an insulator (or conductor with Transport module not set) the shape is set directly in Temp
			if (pMeshCUDA->Temp()->size_cpu().dim()) {

				success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect);

				//lattice temperature for many-temperature models
				if (pHeat->tmtype == TMTYPE_2TM) success &= pMeshCUDA->Temp_l()->resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);

				//lattice temperature for many-temperature models
				if (pHeat->tmtype == TMTYPE_2TM) success &= pMeshCUDA->Temp_l()->set_from_cpuvec(pMesh->Temp_l);
			}
		}

		//if we've just changed the model type, make sure any new temperature VECs start from a reasonable setting - same as the primary Temp VEC.
		if (initialize_Temp_l && pMeshCUDA->Temp_l()->linear_size_cpu()) pMeshCUDA->Temp_l()->set_from_cpuvec(pMesh->Temp_l);

		//allocate memory for the heatEq_RHS auxiliary vector
		if (success) {

			if (heatEq_RHS.size() != pMesh->n_t.dim()) {

				success = heatEq_RHS.resize(pMesh->n_t.dim());
				if (success) heatEq_RHS.set(0.0);
			}
		}
	}

	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	
	return error;
}

void HeatCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		if (Q_equation.is_set()) {

			SetQEquation(pHeat->Q_equation.get_scalar_fspec());
		}
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (Q_equation.is_set()) Q_equation.clear();
	}
}

void HeatCUDA::UpdateField(void)
{
}

#endif

#endif