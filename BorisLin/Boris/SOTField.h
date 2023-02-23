#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#ifdef MODULE_COMPILATION_SOTFIELD

//SOTField module can only be used in a magnetic mesh

class SOTField :
	public Modules,
	public ProgramState<SOTField, std::tuple<>, std::tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	SOTField(Mesh *pMesh_);
	~SOTField();

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

class SOTField :
	public Modules
{

private:

public:

	SOTField(Mesh *pMesh_) {}
	~SOTField() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif