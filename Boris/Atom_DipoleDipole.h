#pragma once

#include "BorisLib.h"
#include "Modules.h"
#include "DemagBase.h"



class Atom_Mesh;

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "Convolution.h"
#include "DipoleDipoleKernel.h"

#if COMPILECUDA == 1
#include "Atom_DipoleDipoleCUDA.h"
#endif

class Atom_DipoleDipole :
	public Modules,
	public DemagBase,
	public Convolution<Atom_DipoleDipole, DipoleDipoleKernel>,
	public ProgramState<Atom_DipoleDipole, std::tuple<INT3>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend Atom_DipoleDipoleCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//divide energy by this to obtain energy density : this is the energy density in the entire mesh, which may not be rectangular.
	double non_empty_volume = 0.0;

	//The dipole-dipole field and moment computed separately at the macrocell size.
	//Hd has cellsize h_dm (but can be cleared so need to keep this info separate, above).
	//These are only used if the macrocell is enabled (i.e. h_dm is greater than h), however Hd is also used if eval speedup is enabled so we can store field computations
	VEC<DBL3> M, Hd;

	//when using the evaluation speedup method we must ensure we have a previous Hd evaluation available
	bool Hd_calculated = false;

	//use macrocell method, or compute dipole-dipole interaction at the atomic unit cell level (i.e. when h_dm == h)?
	bool using_macrocell = true;

private:

	//Initialize mesh transfer from atomistic mesh to M if using a macrocell
	BError Initialize_Mesh_Transfer(void);

public:

	Atom_DipoleDipole(Atom_Mesh *paMesh_);
	~Atom_DipoleDipole();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Setters

	//Set PBC
	BError Set_PBC(INT3 demag_pbc_images_);
};

#else

class Atom_DipoleDipole :
	public Modules,
	public DemagBase
{

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh* paMesh;

public:

	Atom_DipoleDipole(Atom_Mesh *paMesh_) {}
	~Atom_DipoleDipole() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Setters

	//Set PBC
	BError Set_PBC(INT3 demag_pbc_images_) { return BError(); }
};

#endif

