#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DipoleDipoleKernelCUDA.h"

class Atom_MeshCUDA;
class Atom_DipoleDipole;

class Atom_DipoleDipoleCUDA :
	public ModulesCUDA,
	public ConvolutionCUDA<DipoleDipoleKernelCUDA>
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

	//when using the evaluation speedup method we must ensure we have a previous Hd evaluation available
	bool Hd_calculated = false;

	//use macrocell method, or compute dipole-dipole interaction at the atomic unit cell level (i.e. when h_dm == h)?
	bool using_macrocell = true;

private:

	//convert value in energy to energy density by dividing by cellsize volume of V
	void Energy_to_EnergyDensity(cu_obj<cuVEC<cuReal3>>& V);

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
