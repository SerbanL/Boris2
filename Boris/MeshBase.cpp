#include "stdafx.h"
#include "MeshBase.h"
#include "SuperMesh.h"

int MeshBase::meshIdCounter = 0;

MeshBase::MeshBase(MESH_ meshType, SuperMesh *pSMesh_) :
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

//----------------------------------- VALUE GETTERS

//get energy value for given module or one of its exclusive modules (if none active return 0); call it with MOD_ALL to return total energy density in this mesh;
//If avRect is null then get energy density in entire mesh
double MeshBase::GetEnergyDensity(MOD_ moduleType, Rect avRect)
{
	if (moduleType == MOD_ALL) {

		double energy = 0.0;

		//get total energy density in currently set modules
		for (int idx = 0; idx < pMod.size(); idx++) {

			if (avRect.IsNull() || avRect == meshRect) energy += pMod[idx]->GetEnergyDensity();
			else energy += pMod[idx]->GetEnergyDensity(avRect);
		}

		return energy;
	}
	else {

		for (int idx = 0; idx < (int)exclusiveModules[moduleType].size(); idx++) {

			MOD_ moduleID = exclusiveModules[moduleType][idx];
			if (IsModuleSet(moduleID)) {
				
				if (avRect.IsNull() || avRect == meshRect) return pMod(moduleID)->GetEnergyDensity();
				return pMod(moduleID)->GetEnergyDensity(avRect);
			}
		}

		return 0.0;
	}
}