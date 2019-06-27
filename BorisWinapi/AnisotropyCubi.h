#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_ANICUBI

class Anisotropy_Cubic :
	public Modules,
	public ProgramState <Anisotropy_Cubic, tuple<>, tuple<>>
{
private:

	//pointer to mesh object holding this effective field module
	FMesh * pMesh;

public:

	Anisotropy_Cubic(Mesh *pMesh_);
	~Anisotropy_Cubic() {}

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

class Anisotropy_Cubic :
	public Modules
{
private:

	//pointer to mesh object holding this effective field module
	FMesh * pMesh;

public:

	Anisotropy_Cubic(Mesh *pMesh_) {}
	~Anisotropy_Cubic() {}

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
