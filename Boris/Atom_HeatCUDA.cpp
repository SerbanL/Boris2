#include "stdafx.h"
#include "Atom_HeatCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_HEAT) && ATOMISTIC == 1

#include "Atom_Heat.h"
#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "SuperMesh.h"

Atom_HeatCUDA::Atom_HeatCUDA(Atom_Mesh* paMesh_, SuperMesh* pSMesh_, Atom_Heat* paHeat_) :
	ModulesCUDA(),
	HeatBaseCUDA()
{
	paMesh = paMesh_;
	paMeshCUDA = paMesh->paMeshCUDA;
	pSMesh = pSMesh_;
	paHeat = paHeat_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	if (!error_on_create) error_on_create = temp_cmbnd_funcs()->set_pointers(paMeshCUDA);
}

Atom_HeatCUDA::~Atom_HeatCUDA()
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
		paMeshCUDA->Temp()->clear();
		paMeshCUDA->Temp_l()->clear();
	}
}


//-------------------Abstract base class method implementations

BError Atom_HeatCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_HeatCUDA));

	//no energy density contribution here
	ZeroEnergy();

	paHeat->SetRobinBoundaryConditions();

	initialized = true;

	return error;
}

BError Atom_HeatCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_HeatCUDA));

	Uninitialize();

	bool success = true;

	//need this when we switch cuda mode
	if (!Q_equation.is_set() && paHeat->Q_equation.is_set()) error = SetQEquation(paHeat->Q_equation.get_scalar_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, UPDATECONFIG_HEAT_MODELTYPE)) {

		//update mesh dimensions in equation constants
		if (Q_equation.is_set()) {

			error = SetQEquation(paHeat->Q_equation.get_scalar_fspec());
		}

		//first clear any VECs which are not needed
		switch (paHeat->tmtype) {

		case TMTYPE_1TM:
			//single temperature only
			paMeshCUDA->Temp_l()->clear();
			break;

		case TMTYPE_2TM:
			//itinerant electrons temperature <-> lattice temperature
			break;
		}

		//if we've just changed the model type, make sure any new temperature VECs start from a reasonable setting - same as the primary Temp VEC.
		bool initialize_Temp_l = !(paMeshCUDA->Temp_l()->linear_size_cpu());

		if (paMeshCUDA->M1()->size_cpu().dim()) {

			//in a magnet the moments set the shape only on initialization. If already initialized then shape already set.
			if (paMeshCUDA->Temp()->size_cpu().dim()) {

				success = paMeshCUDA->Temp()->resize(paMesh->h_t, paMesh->meshRect);

				//lattice temperature for many-temperature models
				if (paHeat->tmtype == TMTYPE_2TM) success &= paMeshCUDA->Temp_l()->resize(paMesh->h_t, paMesh->meshRect);
			}
			else {

				success = paMeshCUDA->Temp()->set_from_cpuvec(paMesh->Temp);

				//lattice temperature for many-temperature models
				if (paHeat->tmtype == TMTYPE_2TM) success &= paMeshCUDA->Temp_l()->set_from_cpuvec(paMesh->Temp_l);
			}
		}
		else if (paMeshCUDA->elC()->size_cpu().dim()) {

			//in a normal metal the electrical conductivity sets the shape

			//before doing this must make sure elC was itself set in the Transport module (it could be this module is being updated before the Transport module)
			//Transport module is guaranteed to be set otherwise elC would have zero size - it does mean Transport has UpdateConfiguration called twice but it doesn't matter.
			//TO DO
			/*
			if (paMesh->IsModuleSet(MOD_TRANSPORT)) {

				error = (*pMesh)(MOD_TRANSPORT)->UpdateConfiguration(cfgMessage);
				if (error) return error;
			}
			*/

			if (success) {

				if (paMeshCUDA->Temp()->size_cpu().dim()) {

					success = paMeshCUDA->Temp()->resize(paMesh->h_t, paMesh->meshRect, (cuVEC_VC<cuBReal>&)paMeshCUDA->elC);

					//lattice temperature for many-temperature models
					if (paHeat->tmtype == TMTYPE_2TM) success &= paMeshCUDA->Temp_l()->resize(paMesh->h_t, paMesh->meshRect, (cuVEC_VC<cuBReal>&)paMeshCUDA->Temp);
				}
				else {

					success = paMeshCUDA->Temp()->set_from_cpuvec(paMesh->Temp);

					//lattice temperature for many-temperature models
					if (paHeat->tmtype == TMTYPE_2TM) success &= paMeshCUDA->Temp_l()->set_from_cpuvec(paMesh->Temp_l);
				}
			}
		}
		else {

			//in an insulator (or conductor with Transport module not set) the shape is set directly in Temp
			if (paMeshCUDA->Temp()->size_cpu().dim()) {

				success = paMeshCUDA->Temp()->resize(paMesh->h_t, paMesh->meshRect);

				//lattice temperature for many-temperature models
				if (paHeat->tmtype == TMTYPE_2TM) success &= paMeshCUDA->Temp_l()->resize(paMesh->h_t, paMesh->meshRect);
			}
			else {

				success = paMeshCUDA->Temp()->set_from_cpuvec(paMesh->Temp);

				//lattice temperature for many-temperature models
				if (paHeat->tmtype == TMTYPE_2TM) success &= paMeshCUDA->Temp_l()->set_from_cpuvec(paMesh->Temp_l);
			}
		}

		//if we've just changed the model type, make sure any new temperature VECs start from a reasonable setting - same as the primary Temp VEC.
		if (initialize_Temp_l && paMeshCUDA->Temp_l()->linear_size_cpu()) paMeshCUDA->Temp_l()->set_from_cpuvec(paMesh->Temp_l);

		//allocate memory for the heatEq_RHS auxiliary vector
		if (success) {

			if (heatEq_RHS.size() != paMesh->n_t.dim()) {

				success = heatEq_RHS.resize(paMesh->n_t.dim());
				if (success) heatEq_RHS.set(0.0);
			}
		}
	}

	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

void Atom_HeatCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		if (Q_equation.is_set()) {

			SetQEquation(paHeat->Q_equation.get_scalar_fspec());
		}
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (Q_equation.is_set()) Q_equation.clear();
	}
}

void Atom_HeatCUDA::UpdateField(void)
{
}

#endif

#endif