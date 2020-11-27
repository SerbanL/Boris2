#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class SuperMesh;
class HeatBase;

#ifdef MODULE_COMPILATION_HEAT

#if COMPILECUDA == 1
#include "SHeatCUDA.h"
#endif

class SHeat :
	public Modules,
	public ProgramState<SHeat, tuple<double>, tuple<>>
{

#if COMPILECUDA == 1
	friend SHeatCUDA;
#endif

private:

	//pointer to supermesh
	SuperMesh* pSMesh;

	//---------------------- CMBND data

	//CMBND contacts for all contacting transport meshes - these are ordered by first vector index; for each mesh there could be multiple contacting meshes and these are ordered by second vector index
	//CMBNDInfo describes the contact between 2 meshes, allowing calculation of values at cmbnd cells based on continuity of a potential and flux (temperature and heat flux)
	vector< vector<CMBNDInfo> > CMBNDcontacts;

	//list of all Heat modules in meshes (same ordering as first vector in CMBNDcontacts)
	vector<HeatBase*> pHeat;

	//vector of pointers to all Temp - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	vector<VEC_VC<double>*> pTemp;

	//----------------------

	//time step for the heat equation - if in a magnetic mesh must always be smaller or equal to dT (the magnetization equation time-step)
	double heat_dT = 0.2e-12;

	//save the last magnetic dT used: when advancing the heat equation this is the time we need to advance by. 
	//Update magnetic_dT after each heat equation advance (in case an adaptive time-step method is used for the magnetic part).
	double magnetic_dT;

private:

	//calculate and set values at composite media boundaries after all other cells have been computed and set
	void set_cmbnd_values(void);

public:

	SHeat(SuperMesh *pSMesh_);
	~SHeat() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Getters

	double get_heat_dT(void) { return heat_dT; }

	//-------------------Setters

	void set_heat_dT(double dT) { heat_dT = dT; }

};

#else

class SHeat :
	public Modules
{

private:

private:

public:

	SHeat(SuperMesh *pSMesh_) {}
	~SHeat() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Getters

	double get_heat_dT(void) { return 0.0; }

	//-------------------Setters

	void set_heat_dT(double dT) {}

};

#endif

