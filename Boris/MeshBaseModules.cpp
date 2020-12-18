#include "stdafx.h"
#include "MeshBase.h"
#include "SuperMesh.h"

//----------------------------------- MODULES CONTROL

void MeshBase::DelModule(MOD_ moduleID)
{
	//mesh modules allowed to occur more than once, e.g. SDemag_Demag : calling this method will erase all.
	while (IsModuleSet(moduleID)) {

		//delete memory allocated, then erase from list of modules. 
		if (pMod(moduleID)) delete pMod(moduleID);
		pMod.erase(INT2(moduleID, 0));

		//there may be cases where we want multiple modules of same type
		//we always want the minor ids to start numbering at 0 so regularise them
		pMod.regularise_minor_ids(moduleID);
	}
}

//Delete this particular module, if found.
//This method is intended to be called from within the module itself, asking for it to be deleted - a call to its dtor will be issued.
//When using this method in this way, the calling module should immediately return as any further data access there will result in bad memory access.
//The return addess is still valid as it will be stored on the stack.
void MeshBase::DelModule(Modules* pModule)
{
	for (int idx = 0; idx < pMod.size(); idx++) {

		if (pMod[idx] == pModule) {

			//found it

			//the major ID of this module
			MOD_ moduleID = (MOD_)pMod.get_ID_from_index(idx);

			//delete it - free memory then delete it from the list of active module
			delete pMod[idx];
			pMod.erase(idx);

			//there may be cases where we want multiple modules of same type
			//we always want the minor ids to start numbering at 0 so regularise them
			pMod.regularise_minor_ids(moduleID);
		}
	}
}

BError MeshBase::InitializeAllModules(void)
{
	BError error(std::string(CLASS_STR(MeshBase)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		error = pMod[idx]->Initialize();
		if (error) return error;
	}

	return error;
}

#if COMPILECUDA == 1
BError MeshBase::InitializeAllModulesCUDA(void)
{
	BError error(std::string(CLASS_STR(MeshBase)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		error = pMod[idx]->InitializeCUDA();
		if (error) return error;
	}

	return error;
}
#endif

double MeshBase::UpdateModules(void)
{
	//if base_temperature is set using a text equation then set it now
	if (T_equation.is_set()) {

		SetBaseTemperature(T_equation.evaluate(pSMesh->GetStageTime()), false);
	}

	//total energy density
	double energy = 0;

	//Update effective field by adding in contributions from each set module
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		//if for a module it doesn't make sense to contribute to the total energy density, then it should return zero.
		energy += pMod[idx]->UpdateField();
	}

	return energy;
}

#if COMPILECUDA == 1
void MeshBase::UpdateModulesCUDA(void)
{
	//if base_temperature is set using a text equation then set it now
	if (T_equation.is_set()) {

		SetBaseTemperature(T_equation.evaluate(pSMesh->GetStageTime()), false);
	}

	//Update effective field by adding in contributions from each set module
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		pMod[idx]->UpdateFieldCUDA();
	}
}
#endif