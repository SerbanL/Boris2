#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

Mesh::Mesh(MESH_ meshType, SuperMesh *pSMesh_) :
	MeshBase(meshType, pSMesh_),
	//MeshParams constructor after Meshbase, since they both have virtual inheritance from MeshParamsBase : we want MeshParams to control setting values in MeshParamsBase
	MeshParams(params_for_meshtype(meshType))
{
}

Mesh::~Mesh() 
{
	//delete all allocated Modules
	//This has to go here, not in MeshBase destructor, even though pMod is held there
	//The reason for this, some modules can access Mesh data when destructing, and MeshBase destructor is called after Mesh destructor
	//(so Mesh data no longer defined at that point resulting in undefined behaviour)
	clear_vector(pMod);

#if COMPILECUDA == 1
	//free cuda memory by deleting allocated pMeshBaseCUDA
	if (pMeshBaseCUDA) {
		
		//mark implementation of Mesh as destroyed so the CUDA mesh version doesn't attempt to use its data in destructor
		pMeshBaseCUDA->Holder_Mesh_Destroyed();

		delete pMeshBaseCUDA;
		pMeshBaseCUDA = nullptr;
		pMeshCUDA = nullptr;
	}
#endif
}

//calls pSMesh->UpdateConfiguration
BError Mesh::SuperMesh_UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	return pSMesh->UpdateConfiguration(cfgMessage);
}