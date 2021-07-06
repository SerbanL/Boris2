#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#include "DemagBase.h"

#ifdef MODULE_COMPILATION_DEMAG

#include "Convolution.h"
#include "DemagKernel.h"

#if COMPILECUDA == 1
#include "DemagCUDA.h"
#endif

class Demag : 
	public Modules,
	public DemagBase,
	public Convolution<Demag, DemagKernel>,
	public ProgramState<Demag, std::tuple<INT3>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend DemagCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

public:

	Demag(Mesh *pMesh_);
	~Demag();

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

	//-------------------Energy methods

	double Get_EnergyChange(int spin_index, DBL3 Mnew);
};

#else

class Demag :
	public Modules,
	public DemagBase
{

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	Demag(Mesh *pMesh_) {}
	~Demag() {}

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