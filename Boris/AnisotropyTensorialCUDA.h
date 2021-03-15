#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANITENS

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MeshCUDA;

//--------------- UNIAXIAL

class Anisotropy_TensorialCUDA :
	public ModulesCUDA
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	Anisotropy_TensorialCUDA(MeshCUDA* pMeshCUDA_);
	~Anisotropy_TensorialCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class Anisotropy_TensorialCUDA
{
};

#endif

#endif