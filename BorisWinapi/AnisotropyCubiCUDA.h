#pragma once

#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1

#ifdef MODULE_ANICUBI

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class FMeshCUDA;

//--------------- CUBIC

class Anisotropy_CubicCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA* pMeshCUDA;

public:

	Anisotropy_CubicCUDA(FMeshCUDA* pMeshCUDA_);
	~Anisotropy_CubicCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------

};

#else

class Anisotropy_CubicCUDA 
{
};

#endif

#endif


