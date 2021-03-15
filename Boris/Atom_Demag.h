#pragma once

#include "BorisLib.h"
#include "Modules.h"
#include "DemagBase.h"



class Atom_Mesh;

#if defined(MODULE_COMPILATION_DEMAG) && ATOMISTIC == 1

#include "Convolution.h"
#include "DemagKernel.h"

#if COMPILECUDA == 1
#include "Atom_DemagCUDA.h"
#endif

class Atom_Demag :
	public Modules,
	public DemagBase,
	public Convolution<Atom_Demag, DemagKernel>,
	public ProgramState<Atom_Demag, std::tuple<INT3>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend Atom_DemagCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//need to calculate non-empty cells here so we don't waste time during computations (M is a VEC, not a VEC_VC, which means non-empty cells need to be calculated on every call)
	//obtained at initialization
	int non_empty_cells = 0;

	//The demag field and magnetization computed separately at the demag macrocell size.
	//Hd has cellsize h_dm (but can be cleared so need to keep this info separate, above).
	VEC<DBL3> M, Hd;

private:

	//Initialize mesh transfer from atomistic mesh to micromagnetic mesh for demag field computation
	BError Initialize_Mesh_Transfer(void);

public:

	Atom_Demag(Atom_Mesh *paMesh_);
	~Atom_Demag();

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

class Atom_Demag :
	public Modules,
	public DemagBase
{

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh* paMesh;

public:

	Atom_Demag(Atom_Mesh *paMesh_) {}
	~Atom_Demag() {}

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
