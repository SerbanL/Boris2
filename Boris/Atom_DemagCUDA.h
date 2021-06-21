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
	//Hd has cellsize h_dm (but can be cleared so need to keep this info separate, above).
	cu_obj<cuVEC<cuReal3>> M, Hd;

	//Evaluation speedup mode data

	//1 Hdemag: no extrapolation, just save evaluation and reuse
	//2 Hdemag: linear extrapolation, need 2
	//3 Hdemag: quadratic extrapolation, need 3
	cu_obj<cuVEC<cuReal3>> Hdemag, Hdemag2, Hdemag3;

	//times at which evaluations were done, used for extrapolation
	double time_demag1 = 0.0, time_demag2 = 0.0, time_demag3 = 0.0;

	int num_Hdemag_saved = 0;

	//-Nxx, -Nyy, -Nzz values at r = r0
	cu_obj<cuReal3> selfDemagCoeff;

private:

	//Set pointers in ManagedAtom_MeshCUDA so we can access them in device code. This is used by MonteCarlo algorithm.
	void set_Atom_DemagCUDA_pointers(void);

	void Atom_Demag_EvalSpeedup_SubSelf(cu_obj<cuVEC<cuReal3>>& H);

	void Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3);
	void Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2);
	void Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H);

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