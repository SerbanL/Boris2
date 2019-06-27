#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_EXCHANGE

//Exchange modules can only be used in a ferromagnetic mesh

class Exch_6ngbr_Neu : 
	public Modules,
	public ProgramState<Exch_6ngbr_Neu, tuple<>, tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	FMesh *pMesh;

public:

	Exch_6ngbr_Neu(Mesh *pMesh_);
	~Exch_6ngbr_Neu();

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

class Exch_6ngbr_Neu :
	public Modules
{

public:

	Exch_6ngbr_Neu(Mesh *pMesh_) {}
	~Exch_6ngbr_Neu() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	void UpdateField(void) {}
};

#endif