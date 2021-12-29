#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Atom_Mesh;

#if defined(MODULE_COMPILATION_SOTFIELD) && ATOMISTIC == 1

class Atom_SOTField :
	public Modules,
	public ProgramState<Atom_SOTField, std::tuple<>, std::tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh* paMesh;

public:

	Atom_SOTField(Atom_Mesh* paMesh_);
	~Atom_SOTField();

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

class Atom_SOTField :
	public Modules
{

private:

public:

	Atom_SOTField(Atom_Mesh* paMesh_) {}
	~Atom_SOTField() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif