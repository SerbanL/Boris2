#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_DEMAG_N

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag using demagnetizing factors (for uniform magnetization in ellipsoidal shapes)

class Demag_N :
	public Modules,
	public ProgramState<Demag_N, tuple<>, tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	FMesh * pMesh;

public:

	Demag_N(Mesh *pMesh_);
	~Demag_N();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError MakeCUDAModule(void);

	double UpdateField(void);
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

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif
