#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#ifdef MODULE_COMPILATION_EXCHANGE

#include "ExchangeBase.h"

//Exchange modules can only be used in a magnetic mesh

class Exch_6ngbr_Neu :
	public Modules,
	public ExchangeBase,
	public ProgramState<Exch_6ngbr_Neu, std::tuple<>, std::tuple<>>
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

	//-------------------Energy methods

	//FM mesh
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//AFM mesh
	DBL2 Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B);

	//-------------------Torque methods

	DBL3 GetTorque(Rect& avRect);
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