#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_ANIUNI

//Anisotropy modules can only be used in a ferromagnetic mesh

class Anisotropy_Uniaxial :
	public Modules,
	public ProgramState<Anisotropy_Uniaxial, tuple<>, tuple<>>
{
private:

	//pointer to mesh object holding this effective field module
	FMesh * pMesh;

public:
	
	Anisotropy_Uniaxial(Mesh *pMesh_);
	~Anisotropy_Uniaxial() {}

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

class Anisotropy_Uniaxial :
	public Modules
{
private:

	//pointer to mesh object holding this effective field module
	FMesh * pMesh;

public:

	Anisotropy_Uniaxial(Mesh *pMesh_) {}
	~Anisotropy_Uniaxial() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};


#endif