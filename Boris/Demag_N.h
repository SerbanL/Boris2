#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#ifdef MODULE_COMPILATION_DEMAG_N

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag using demagnetizing factors (for uniform magnetization in ellipsoidal shapes)

class Demag_N :
	public Modules,
	public ProgramState<Demag_N, std::tuple<>, std::tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	Demag_N(Mesh *pMesh_);
	~Demag_N();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	//FM mesh
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//AFM mesh
	DBL2 Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B);
};

#else

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag using demagnetizing factors (for uniform magnetization in ellipsoidal shapes)

class Demag_N :
	public Modules
{

private:

public:

	Demag_N(Mesh *pMesh_) {}
	~Demag_N() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif
