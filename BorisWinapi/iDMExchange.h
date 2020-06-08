#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;

#ifdef MODULE_IDMEXCHANGE

#include "ExchangeBase.h"

//Exchange modules can only be used in a magnetic mesh

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
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy density methods

	double GetEnergyDensity(Rect& avRect);
	double GetEnergy_Max(Rect& rectangle);

	//Compute exchange energy density and store it in displayVEC
	void Compute_Exchange(VEC<double>& displayVEC);
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

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Energy density methods

	double GetEnergyDensity(Rect& avRect) { return 0.0; }
	double GetEnergy_Max(Rect& rectangle) { return 0.0; }
	
	void Compute_Exchange(VEC<double>& displayVEC) {}
};

#endif

