#pragma once

#include "BorisLib.h"
#include "Modules.h"



class Atom_Mesh;

#if defined(MODULE_COMPILATION_MOPTICAL) && ATOMISTIC == 1

//Magneto-optical module can only be used in a magnetic mesh

class Atom_MOptical :
	public Modules,
	public ProgramState<Atom_MOptical, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend class Atom_MOpticalCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//divide energy by this to obtain energy density : this is the energy density in the entire mesh, which may not be rectangular.
	double non_empty_volume = 0.0;

public:

	Atom_MOptical(Atom_Mesh *paMesh_);
	~Atom_MOptical();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	//For simple cubic mesh spin_index coincides with index in M1
	double Get_EnergyChange(int spin_index, DBL3 Mnew);
};

#else

class Atom_MOptical :
	public Modules
{

private:

public:

	Atom_MOptical(Atom_Mesh *paMesh_) {}
	~Atom_MOptical() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif

