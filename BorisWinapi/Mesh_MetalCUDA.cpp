#include "stdafx.h"
#include "Mesh_MetalCUDA.h"
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

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
BError MetalMeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MetalMeshCUDA));

	return error;
}

#endif