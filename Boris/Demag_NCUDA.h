#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DEMAG_N

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MeshCUDA;

class Demag_NCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	Demag_NCUDA(MeshCUDA* pMeshCUDA_);
	~Demag_NCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------

};

#else

class Demag_NCUDA
{
};

#endif

#endif

