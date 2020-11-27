#include "stdafx.h"
#include "SHeat.h"

#ifdef MODULE_COMPILATION_HEAT

#include "SuperMesh.h"
#include "HeatBase.h"

SHeat::SHeat(SuperMesh *pSMesh_) :
	Modules(),
	ProgramStateNames(this, {VINFO(heat_dT)}, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

//-------------------Abstract base class method implementations

BError SHeat::Initialize(void)
{
	BError error(CLASS_STR(SHeat));

	//heat_dT must be set correctly using the magnetic time step
	magnetic_dT = pSMesh->GetTimeStep();

	//check meshes to set heat boundary flags (NF_CMBND flags for Temp)

	//clear everything then rebuild
	pHeat.clear();
	pTemp.clear();
	CMBNDcontacts.clear();

	//now build pHeat (and pTemp)
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->IsModuleSet(MOD_HEAT)) {

			pHeat.push_back(dynamic_cast<HeatBase*>((*pSMesh)[idx]->GetModule(MOD_HEAT)));
			pTemp.push_back(&(*pSMesh)[idx]->Temp);
		}
	}

	//set cmbnd flags (also building contacts)
	for (int idx = 0; idx < (int)pHeat.size(); idx++) {

		//build CMBND contacts and set flags for Temp
		CMBNDcontacts.push_back(pTemp[idx]->set_cmbnd_flags(idx, pTemp));
	}

	initialized = true;

	return error;
}

BError SHeat::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SHeat));

	Uninitialize();

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, UPDATECONFIG_SWITCHCUDASTATE)) {

		//clear everything then rebuild
		pHeat.clear();
		pTemp.clear();

		//now build pHeat (and pTemp)
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->IsModuleSet(MOD_HEAT)) {

				pHeat.push_back(dynamic_cast<HeatBase*>((*pSMesh)[idx]->GetModule(MOD_HEAT)));
				pTemp.push_back(&(*pSMesh)[idx]->Temp);
			}
		}
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError SHeat::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SHeat));

#if COMPILECUDA == 1

	pModuleCUDA = new SHeatCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

double SHeat::UpdateField(void)
{
	//only need to update this after an entire magnetization equation time step is solved
	//also if heat_dT is set to zero skip the heat equation solver : this will maintain a fixed temperature
	if (!pSMesh->CurrentTimeStepSolved() || heat_dT < MINTIMESTEP) return 0.0;

	double dT = heat_dT;

	//number of sub_steps to cover magnetic_dT required when advancing in smaller heat_dT steps
	int sub_steps = (int)floor_epsilon(magnetic_dT / heat_dT);

	//any left-over epsilon_dT < heat_dT
	double epsilon_dT = magnetic_dT - heat_dT * sub_steps;

	for (int step_idx = 0; step_idx < sub_steps + 1; step_idx++) {

		//the last step may have a different time step - take this epsilon_dT step (if not zero)
		if (step_idx == sub_steps) {

			if (epsilon_dT) dT = epsilon_dT;
			else continue;
		}

		//1. solve Temp in each mesh separately (1 iteration each) - CMBND cells not set yet
		for (int idx = 0; idx < (int)pHeat.size(); idx++) {

			if (pHeat[idx]->Get_TMType() == TMTYPE_1TM) pHeat[idx]->IterateHeatEquation_1TM(dT);
			else if (pHeat[idx]->Get_TMType() == TMTYPE_2TM) pHeat[idx]->IterateHeatEquation_2TM(dT);
		}
		
		//2. calculate boundary conditions (i.e. temperature values at CMBND cells)
		set_cmbnd_values();
	}

	//3. update the magnetic dT that will be used next time around to increment the heat solver by
	magnetic_dT = pSMesh->GetTimeStep();

	//no contribution to total energy density
	return 0.0;
}

//calculate and set values at composite media boundaries after all other cells have been computed and set
void SHeat::set_cmbnd_values(void)
{
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;

			pTemp[idx_pri]->set_cmbnd_continuous<HeatBase>(
				*pTemp[idx_sec], CMBNDcontacts[idx1][idx2],
				&HeatBase::afunc_sec, &HeatBase::afunc_pri,
				&HeatBase::bfunc_sec, &HeatBase::bfunc_pri,
				&HeatBase::diff2_sec, &HeatBase::diff2_pri,
				*pHeat[idx_sec], *pHeat[idx_pri]);
		}
	}
}

#endif