#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;

#ifdef MODULE_EXCHANGE

#include "ExchangeBase.h"

//Exchange modules can only be used in a magnetic mesh

class Exch_6ngbr_Neu :
	public Modules,
	public ExchangeBase,
	public ProgramState<Exch_6ngbr_Neu, tuple<>, tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

public:

	Exch_6ngbr_Neu(Mesh *pMesh_);
	~Exch_6ngbr_Neu();

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

class Exch_6ngbr_Neu :
	public Modules
{

public:

	Exch_6ngbr_Neu(Mesh *pMesh_) {}
	~Exch_6ngbr_Neu() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif