#include "stdafx.h"
#include "SHeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "SHeat.h"
#include "SuperMesh.h"

SHeatCUDA::SHeatCUDA(SuperMesh* pSMesh_, SHeat* pSHeat_) :
	ModulesCUDA()
{
	pSMesh = pSMesh_;
	pSHeat = pSHeat_;

	error_on_create = UpdateConfiguration();
}

SHeatCUDA::~SHeatCUDA()
{

}

//-------------------Abstract base class method implementations

BError SHeatCUDA::Initialize(void)
{
	BError error(CLASS_STR(SHeatCUDA));

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError SHeatCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SHeatCUDA));

	Uninitialize();
	
	//check meshes to set heat boundary flags (NF_CMBND flags for Temp)
	
	//clear everything then rebuild
	pHeat.clear();
	CMBNDcontactsCUDA.clear();
	CMBNDcontacts.clear();
	pTemp.clear();

	//now build pHeat (and pTemp)
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->IsModuleSet(MOD_HEAT)) {

			pHeat.push_back(dynamic_cast<HeatCUDA*>((*pSMesh)[idx]->GetCUDAModule(MOD_HEAT)));
			pTemp.push_back(&(*pSMesh)[idx]->pMeshCUDA->Temp);
		}
	}

	//set cmbnd flags
	for (int idx = 0; idx < (int)pHeat.size(); idx++) {

		//it's easier to just copy the flags entirely from the cpu versions.
		//Notes :
		//1. By design the cpu versions are required to keep size and flags up to date (but not mesh values)
		//2. pHeat in SHeat has exactly the same size and order
		//3. SHeat UpdateConfiguration was called just before, which called this CUDA version at the end.

		if (!(*pTemp[idx])()->copyflags_from_cpuvec(*pSHeat->pTemp[idx])) error(BERROR_GPUERROR_CRIT);
	}
	
	for (int idx = 0; idx < pSHeat->CMBNDcontacts.size(); idx++) {

		vector<cu_obj<CMBNDInfoCUDA>> mesh_contacts;
		vector<CMBNDInfoCUDA> mesh_contacts_cpu;

		for (int idx_contact = 0; idx_contact < pSHeat->CMBNDcontacts[idx].size(); idx_contact++) {

			cu_obj<CMBNDInfoCUDA> contact;

			contact()->copy_from_CMBNDInfo<CMBNDInfo>(pSHeat->CMBNDcontacts[idx][idx_contact]);

			mesh_contacts.push_back(contact);

			mesh_contacts_cpu.push_back(pSHeat->CMBNDcontacts[idx][idx_contact]);
		}

		CMBNDcontactsCUDA.push_back(mesh_contacts);
		CMBNDcontacts.push_back(mesh_contacts_cpu);
	}

	return error;
}

void SHeatCUDA::UpdateField(void)
{
	//only need to update this after an entire magnetisation equation time step is solved
	if (!pSMesh->CurrentTimeStepSolved()) return;
	
	cuReal dT = pSHeat->heat_dT;

	//number of sub_steps to cover magnetic_dT required when advancing in smaller heat_dT steps
	int sub_steps = (int)floor_epsilon(pSHeat->magnetic_dT / pSHeat->heat_dT);

	//any left-over epsilon_dT < heat_dT
	cuReal epsilon_dT = pSHeat->magnetic_dT - pSHeat->heat_dT * sub_steps;

	for (int step_idx = 0; step_idx < sub_steps + 1; step_idx++) {

		//the last step may have a different time step - take this epsilon_dT step (if not zero)
		if (step_idx == sub_steps) {

			if (epsilon_dT) dT = epsilon_dT;
			else continue;
		}

		//1. solve Temp in each mesh separately (1 iteration each) - CMBND cells not set yet
		for (int idx = 0; idx < (int)pHeat.size(); idx++) {

			pHeat[idx]->IterateHeatEquation(dT);
		}

		//2. calculate boundary conditions (i.e. temperature values at CMBND cells)
		set_cmbnd_values();
	}

	//3. update the magnetic dT that will be used next time around to increment the heat solver by
	pSHeat->magnetic_dT = pSMesh->GetTimeStep();
}

#endif

#endif