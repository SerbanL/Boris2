#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Atom_Mesh;

#if defined(MODULE_COMPILATION_ANIUNI) && ATOMISTIC == 1

class Atom_Anisotropy_Uniaxial :
	public Modules,
	public ProgramState<Atom_Anisotropy_Uniaxial, tuple<>, tuple<>>
{
private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh * paMesh;

	//divide energy by this to obtain energy density : this is the energy density in the entire mesh, which may not be rectangular.
	double non_empty_volume = 0.0;

public:

	Atom_Anisotropy_Uniaxial(Atom_Mesh *paMesh_);
	~Atom_Anisotropy_Uniaxial() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy density methods

	double GetEnergyDensity(Rect& avRect);
};

#else

class Atom_Anisotropy_Uniaxial :
	public Modules
{
private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh * paMesh;

public:

	Atom_Anisotropy_Uniaxial(Atom_Mesh *paMesh_) {}
	~Atom_Anisotropy_Uniaxial() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Energy density methods

	double GetEnergyDensity(Rect& avRect) { return 0.0; }
};


#endif