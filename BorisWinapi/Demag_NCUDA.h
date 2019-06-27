#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_DEMAG_N

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class FMeshCUDA;

class Demag_NCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA* pMeshCUDA;

public:

	Demag_NCUDA(FMeshCUDA* pMeshCUDA_);
	~Demag_NCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------

};

#else

class Demag_NCUDA
{
};

#endif

#endif

