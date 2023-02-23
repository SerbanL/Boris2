#pragma once

#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANICUBI

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MeshCUDA;

//--------------- CUBIC

class Anisotropy_CubicCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	Anisotropy_CubicCUDA(MeshCUDA* pMeshCUDA_);
	~Anisotropy_CubicCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------Torque methods

	cuReal3 GetTorque(cuRect avRect);
};

#else

class Anisotropy_CubicCUDA 
{
};

#endif

#endif


