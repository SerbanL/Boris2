#include "stdafx.h"
#include "MeshBase.h"
#include "SuperMesh.h"

int MeshBase::meshIdCounter = 0;

MeshBase::MeshBase(MESH_ meshType, SuperMesh *pSMesh_) :
	prng(GetSystemTickCount()),
	pSMesh(pSMesh_)
{
	this->meshType = meshType;

	OmpThreads = omp_get_num_procs();
}

MeshBase::~MeshBase()
{
	//cleanup done in (Atom_)Mesh for reasons explained there
}

BError MeshBase::Error_On_Create(void)
{
	for (int idx = 0; idx < pMod.size(); idx++) {

		if (!error_on_create) error_on_create = pMod[idx]->Error_On_Create();
	}

	return error_on_create;
}

//-----------------------------------

void MeshBase::swap_ids(MeshBase* pMesh_swap)
{
	int meshId_temp = meshId;
	meshId = pMesh_swap->meshId;
	pMesh_swap->meshId = meshId_temp;
}

bool MeshBase::IsOutputDataSet_withRect(int datumId)
{ 
	return pSMesh->IsOutputDataSet_withRect(datumId, this); 
}

//return true if data is set (with any rectangle)
bool MeshBase::IsOutputDataSet(int datumId)
{
	return pSMesh->IsOutputDataSet(datumId, this);
}

//check if given stage is set
bool MeshBase::IsStageSet(int stageType)
{
	return pSMesh->IsStageSet(stageType);
}

//set computefields_if_MC flag on SuperMesh
void MeshBase::Set_Force_MonteCarlo_ComputeFields(bool status)
{
	pSMesh->Set_Force_MonteCarlo_ComputeFields(status);
}

//----------------------------------- VALUE GETTERS

//get energy value for given module or one of its exclusive modules (if none active return 0); call it with MOD_ALL to return total energy density in this mesh;
//If avRect is null then get energy density in entire mesh
double MeshBase::GetEnergyDensity(MOD_ moduleType, Rect avRect)
{
	if (moduleType == MOD_ALL) {

		double energy = 0.0;

		//get total energy density in currently set modules
		for (int idx = 0; idx < pMod.size(); idx++) {

			if (avRect.IsNull() || avRect == (meshRect - meshRect.s)) energy += pMod[idx]->GetEnergyDensity();
			else energy += pMod[idx]->GetEnergyDensity(avRect);
		}

		return energy;
	}
	else {

		for (int idx = 0; idx < (int)exclusiveModules[moduleType].size(); idx++) {

			MOD_ moduleID = exclusiveModules[moduleType][idx];
			if (IsModuleSet(moduleID)) {
				
				if (avRect.IsNull() || avRect == (meshRect - meshRect.s)) return pMod(moduleID)->GetEnergyDensity();
				return pMod(moduleID)->GetEnergyDensity(avRect);
			}
		}

		return 0.0;
	}
}

//as above but for the torque
DBL3 MeshBase::GetTorque(MOD_ moduleType, Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect - meshRect.s;

	if (moduleType == MOD_ALL) return DBL3();
	else {

		for (int idx = 0; idx < (int)exclusiveModules[moduleType].size(); idx++) {

			MOD_ moduleID = exclusiveModules[moduleType][idx];
			if (IsModuleSet(moduleID)) {

				return pMod(moduleID)->GetTorque(avRect);
			}
		}
	}

	return DBL3();
}

//get all set modules IDs
std::vector<MOD_> MeshBase::GetModulesIDs(void)
{
	std::vector<MOD_> set_modules;

	for (int idx = 0; idx < pMod.size(); idx++) {

		set_modules.push_back((MOD_)pMod.get_ID_from_index(idx));
	}

	return set_modules;
}