#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_DEMAG) && ATOMISTIC == 1

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCUDA.h"

class Atom_MeshCUDA;
class Atom_Demag;

class Atom_DemagCUDA :
	public ModulesCUDA,
	public ConvolutionCUDA<Atom_DemagCUDA, DemagKernelCUDA>
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

	//point to cpu version of this module
	Atom_Demag *paDemag;

	//The demag field and magnetization computed separately at the demag macrocell size.
	//Hdemag has cellsize h_dm (but can be cleared so need to keep this info separate, above).
	cu_obj<cuVEC<cuReal3>> M, Hdemag;

	//when using the evaluation speedup method we must ensure we have a previous Hdemag evaluation available
	bool Hdemag_calculated = false;

public:

	Atom_DemagCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_Demag *paDemag_);
	~Atom_DemagCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class Atom_DemagCUDA
{
};

#endif

#endif