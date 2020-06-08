#pragma once

#include "CompileFlags.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "MeshCUDA.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "BorisCUDALib.h"

class DiaMesh;

//Store Mesh quantities as cu_obj managed cuda VECs
class DiaMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	DiaMesh *pDiaMesh;

public:

	//make this object by copying data from the Mesh holding this object
	DiaMeshCUDA(DiaMesh* pMesh);

	~DiaMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS
};

#else

//Store Mesh quantities as cu_obj managed cuda VECs
class DiaMeshCUDA :
	public MeshCUDA
{

private:

public:

	//make this object by copying data from the Mesh holding this object
	DiaMeshCUDA(Mesh* pMesh) :
		MeshCUDA(pMesh)
	{}

	~DiaMeshCUDA() {}

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS
};

#endif
#endif