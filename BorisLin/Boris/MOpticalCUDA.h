#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MOPTICAL

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MOptical;
class MeshCUDA;

class MOpticalCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	MOpticalCUDA(MeshCUDA* pMeshCUDA_);
	~MOpticalCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);
};

#else

class MOpticalCUDA
{
};

#endif

#endif

