#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "MeshCUDA.h"

#ifdef MESH_COMPILATION_INSULATOR

#include "BorisCUDALib.h"

class InsulatorMesh;

//Store Mesh quantities as cu_obj managed cuda VECs
class InsulatorMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	InsulatorMesh * pInsulatorMesh;

public:

	//make this object by copying data from the Mesh holding this object
	InsulatorMeshCUDA(InsulatorMesh* pMesh);

	~InsulatorMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS
};

#else

class InsulatorMeshCUDA :
	public MeshCUDA
{

private:

public:

	//make this object by copying data from the Mesh holding this object
	InsulatorMeshCUDA(Mesh* pMesh) :
		MeshCUDA(pMesh)
	{}

	~InsulatorMeshCUDA() {}

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}
};

#endif
#endif

