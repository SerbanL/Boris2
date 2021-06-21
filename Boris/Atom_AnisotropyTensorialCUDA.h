#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ANITENS) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Atom_MeshCUDA;

//--------------- UNIAXIAL

class Atom_Anisotropy_TensorialCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

public:

	Atom_Anisotropy_TensorialCUDA(Atom_MeshCUDA* paMeshCUDA_);
	~Atom_Anisotropy_TensorialCUDA();

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

class Atom_Anisotropy_TensorialCUDA
{
};

#endif

#endif

