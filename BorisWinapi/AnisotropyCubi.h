#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;

#ifdef MODULE_ANICUBI

class Anisotropy_Cubic :
	public Modules,
	public ProgramState <Anisotropy_Cubic, tuple<>, tuple<>>
{
private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	Anisotropy_Cubic(Mesh *pMesh_);
	~Anisotropy_Cubic() {}

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

class Anisotropy_Cubic :
	public Modules
{
private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	Anisotropy_Cubic(Mesh *pMesh_) {}
	~Anisotropy_Cubic() {}

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
