#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;

#ifdef MODULE_IDMEXCHANGE

#include "ExchangeBase.h"

//Exchange modules can only be used in a ferromagnetic mesh

class iDMExchange :
	public Modules,
	public ExchangeBase,
	public ProgramState<iDMExchange, tuple<>, tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

public:

	iDMExchange(Mesh *pMesh_);
	~iDMExchange();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);
};

#else

class iDMExchange :
	public Modules
{

private:

public:

	iDMExchange(Mesh *pMesh_) {}
	~iDMExchange() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif

