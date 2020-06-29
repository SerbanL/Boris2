#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_DEMAG_N) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Atom_MeshCUDA;

class Atom_Demag_NCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

public:

	Atom_Demag_NCUDA(Atom_MeshCUDA* paMeshCUDA_);
	~Atom_Demag_NCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------

};

#else

class Atom_Demag_NCUDA
{
};

#endif

#endif


