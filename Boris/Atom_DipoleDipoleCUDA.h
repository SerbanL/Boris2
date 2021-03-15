#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DipoleDipoleKernelCUDA.h"

class Atom_MeshCUDA;
class Atom_DipoleDipole;

class Atom_DipoleDipoleCUDA :
	public ModulesCUDA,
	public ConvolutionCUDA<Atom_DipoleDipoleCUDA, DipoleDipoleKernelCUDA>
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

	//point to cpu version of this module
	Atom_DipoleDipole *paDipoleDipole;

	//The dipole-dipole field and moment computed separately at the macrocell size.
	//Hd has cellsize h_dm (but can be cleared so need to keep this info separate, above).
	//These are only used if the macrocell is enabled (i.e. h_dm is greater than h)
	cu_obj<cuVEC<cuReal3>> M, Hd;

	//use macrocell method, or compute dipole-dipole interaction at the atomic unit cell level (i.e. when h_dm == h)?
	bool using_macrocell = true;

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

	//convert value in energy to energy density by dividing by cellsize volume of V
	void Energy_to_EnergyDensity(cu_obj<cuVEC<cuReal3>>& V);

	//Macrocell mode: subtract self contribution from calculated field
	void Atom_DipoleDipole_EvalSpeedup_SubSelf(cu_obj<cuVEC<cuReal3>>& H);

	//Macrocell mode, QUADRATIC: extrapolate field and add self contribution
	void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3);
	//Macrocell mode, LINEAR: extrapolate field and add self contribution
	void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2);
	//Macrocell mode, STEP: extrapolate field and add self contribution
	void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H);

	//QUADRATIC: extrapolate field
	void Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3);
	//LINEAR: extrapolate field
	void Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2);
	//STEP: extrapolate field
	void Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H);

public:

	Atom_DipoleDipoleCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_DipoleDipole *paDipoleDipole_);
	~Atom_DipoleDipoleCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class Atom_DipoleDipoleCUDA
{
};

#endif

#endif
