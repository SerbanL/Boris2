#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_SURFEXCHANGE

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class FMeshCUDA;
class SurfExchange;
class ManagedMeshCUDA;

class SurfExchangeCUDA :
	public ModulesCUDA
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA* pMeshCUDA;

	//pointer to cpu version of SurfExchange
	SurfExchange* pSurfExch;

	//cu arrays with pointers to other meshes in surface exchange coupling with the mesh holding this module, top and bottom
	cu_arr<ManagedMeshCUDA> pMesh_Bot;
	cu_arr<ManagedMeshCUDA> pMesh_Top;

public:

	SurfExchangeCUDA(FMeshCUDA* pMeshCUDA_, SurfExchange* pSurfExch_);
	~SurfExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------

};

#else

class SurfExchangeCUDA
{
};

#endif

#endif