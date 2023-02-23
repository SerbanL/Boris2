#pragma once

#include "BorisLib.h"
#include "Modules.h"



class Atom_Mesh;

#if defined(MODULE_COMPILATION_DEMAG_N) && ATOMISTIC == 1

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag using demagnetizing factors (for uniform magnetization in ellipsoidal shapes)

class Atom_Demag_N :
	public Modules,
	public ProgramState<Atom_Demag_N, std::tuple<>, std::tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh * paMesh;

public:

	Atom_Demag_N(Atom_Mesh *paMesh_);
	~Atom_Demag_N();

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

	//For simple cubic mesh spin_index coincides with index in M1
	double Get_EnergyChange(int spin_index, DBL3 Mnew);
};

#else

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag using demagnetizing factors (for uniform magnetization in ellipsoidal shapes)

class Atom_Demag_N :
	public Modules
{

private:

public:

	Atom_Demag_N(Atom_Mesh *paMesh_) {}
	~Atom_Demag_N() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Energy methods
};

#endif

