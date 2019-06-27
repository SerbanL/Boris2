#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_DEMAG

#include "Convolution.h"
#include "DemagKernel.h"

class Demag : 
	public Modules, 
	public Convolution<DemagKernel>,
	public ProgramState<Demag, tuple<>, tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	FMesh *pMesh;

public:

	Demag(Mesh *pMesh_);
	~Demag() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError MakeCUDAModule(void);

	void UpdateField(void);
};

#else

class Demag :
	public Modules
{

private:

	//pointer to mesh object holding this effective field module
	FMesh * pMesh;

public:

	Demag(Mesh *pMesh_) {}
	~Demag() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	void UpdateField(void) {}
};

#endif