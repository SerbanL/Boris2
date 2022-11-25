#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "Atom_STFieldCUDA.h"
#endif

class Atom_Mesh;
class Mesh;

#if defined(MODULE_COMPILATION_STFIELD) && ATOMISTIC == 1

class Atom_STField :
	public Modules,
	public ProgramState<Atom_STField, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend Atom_STFieldCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh * paMesh;

	//If STp is set to zero, these are the magnetic meshes which provide polarization vectors, top and bottom
	std::vector<Atom_Mesh*> paMesh_Bot, paMesh_Top;

	//as above for micromagnetic meshes top and bottom
	std::vector<Mesh*> pMesh_Bot, pMesh_Top;

public:

	Atom_STField(Atom_Mesh *paMesh_);
	~Atom_STField();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);
};

#else

class Atom_STField :
	public Modules
{

private:

public:

	Atom_STField(Atom_Mesh *paMesh_) {}
	~Atom_STField() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif