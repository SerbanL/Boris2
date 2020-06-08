#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Atom_Mesh;

#if defined(MODULE_DEMAG) && ATOMISTIC == 1

#include "Convolution.h"
#include "DemagKernel.h"

#if COMPILECUDA == 1
#include "Atom_DemagCUDA.h"
#endif

class Atom_Demag :
	public Modules,
	public Convolution<DemagKernel>,
	public ProgramState<Atom_Demag, tuple<INT3>, tuple<>>
{

#if COMPILECUDA == 1
	friend Atom_DemagCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//number of pbc images in each dimension (set to zero to disable).
	//There is also a copy of this in ConvolutionData inherited from Convolution - we need another copy here to detect changes
	//these pbc images are applicable in individual demag modules only
	INT3 demag_pbc_images = INT3();

	//need to calculate non-empty cells here so we don't waste time during computations (M is a VEC, not a VEC_VC, which means non-empty cells need to be calculated on every call)
	//obtained at initialization
	int non_empty_cells = 0;

	//The demag field and magnetization computed separately at the demag macrocell size.
	//Hdemag has cellsize h_dm (but can be cleared so need to keep this info separate, above).
	VEC<DBL3> M, Hdemag;

	//when using the evaluation speedup method we must ensure we have a previous Hdemag evaluation available
	bool Hdemag_calculated = false;

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

	//-------------------Getters

	//Get PBC images
	INT3 Get_PBC(void) const { return demag_pbc_images; }
};

#else

class Atom_Demag :
	public Modules
{

private:

	//pointer to mesh object holding this effective field module
	Atom_DemagMesh * paMesh;

public:

	Atom_Demag(Atom_DemagMesh *paMesh_) {}
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

	//-------------------Getters

	//Get PBC images
	INT3 Get_PBC(void) { return INT3(); }
};

#endif
