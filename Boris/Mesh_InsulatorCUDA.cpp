#include "stdafx.h"
#include "Mesh_InsulatorCUDA.h"

#ifdef MESH_COMPILATION_INSULATOR

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

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError InsulatorMeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(InsulatorMeshCUDA));

	return error;
}

#endif
#endif