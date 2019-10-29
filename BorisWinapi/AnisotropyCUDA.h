#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_ANIUNI

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MeshCUDA;

//--------------- UNIAXIAL

class Anisotropy_UniaxialCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	Anisotropy_UniaxialCUDA(MeshCUDA* pMeshCUDA_);
	~Anisotropy_UniaxialCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------

};

#else

class Anisotropy_UniaxialCUDA
{
};

#endif

#endif


