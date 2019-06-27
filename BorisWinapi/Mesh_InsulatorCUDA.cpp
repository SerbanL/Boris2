#include "stdafx.h"
#include "Mesh_InsulatorCUDA.h"
#include "Mesh_Insulator.h"

#if COMPILECUDA == 1

InsulatorMeshCUDA::InsulatorMeshCUDA(InsulatorMesh* pMesh) :
	MeshCUDA(pMesh)
{
	pInsulatorMesh = pMesh;
}

InsulatorMeshCUDA::~InsulatorMeshCUDA()
{

}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
BError InsulatorMeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(InsulatorMeshCUDA));

	return error;
}

#endif