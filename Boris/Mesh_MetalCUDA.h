#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "MeshCUDA.h"

#ifdef MESH_COMPILATION_METAL

#include "BorisCUDALib.h"

class MetalMesh;

//Store Mesh quantities as cu_obj managed cuda VECs
class MetalMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	MetalMesh * pMetalMesh;

public:

	//make this object by copying data from the Mesh holding this object
	MetalMeshCUDA(MetalMesh* pMesh);

	~MetalMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS
};

#else

class MetalMeshCUDA :
	public MeshCUDA
{

private:

public:

	//make this object by copying data from the Mesh holding this object
	MetalMeshCUDA(Mesh* pMesh) :
		MeshCUDA(pMesh)
	{}

	~MetalMeshCUDA() {}

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}
};

#endif
#endif

