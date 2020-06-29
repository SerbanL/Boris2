#include "stdafx.h"
#include "Mesh_MetalCUDA.h"

#ifdef MESH_COMPILATION_METAL

#include "Mesh_Metal.h"

#if COMPILECUDA == 1

MetalMeshCUDA::MetalMeshCUDA(MetalMesh* pMesh) :
	MeshCUDA(pMesh)
{
	pMetalMesh = pMesh;
}

MetalMeshCUDA::~MetalMeshCUDA()
{

}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError MetalMeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MetalMeshCUDA));

	return error;
}

#endif
#endif