#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#ifdef MODULE_COMPILATION_ANIUNI

//Anisotropy modules can only be used in a magnetic mesh

class Anisotropy_Uniaxial :
	public Modules,
	public ProgramState<Anisotropy_Uniaxial, std::tuple<>, std::tuple<>>
{
private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:
	
	Anisotropy_Uniaxial(Mesh *pMesh_);
	~Anisotropy_Uniaxial() {}

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

class Anisotropy_Uniaxial :
	public Modules
{
private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	Anisotropy_Uniaxial(Mesh *pMesh_) {}
	~Anisotropy_Uniaxial() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};


#endif