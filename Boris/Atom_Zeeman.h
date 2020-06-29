#pragma once

#include "Boris_Enums_Defs.h"

#include "BorisLib.h"
#include "Modules.h"

#include "ZeemanBase.h"

using namespace std;

class Atom_Mesh;

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

//Zeeman module can only be used in a magnetic mesh

class SuperMesh;

class Atom_Zeeman : 
	public Modules,
	public ZeemanBase,
	public ProgramState<Atom_Zeeman, tuple<DBL3, TEquation<double, double, double, double>>, tuple<>>
{

#if COMPILECUDA == 1
	friend class Atom_ZeemanCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//pointer to supermesh
	SuperMesh* pSMesh;

	//divide energy by this to obtain energy density : this is the energy density in the entire mesh, which may not be rectangular.
	double non_empty_volume = 0.0;

private:

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true);

public:

	Atom_Zeeman(Atom_Mesh *paMesh_);
	~Atom_Zeeman();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy density methods

	double GetEnergyDensity(Rect& avRect);

	//-------------------

	void SetField(DBL3 Hxyz);
	DBL3 GetField(void) { return Ha; }

	BError SetFieldEquation(string equation_string, int step);

	//if base temperature changes we need to adjust Tb in H_equation if it's used.
	void SetBaseTemperature(double Temperature);
};

#else

class Atom_Zeeman :
	public Modules,
	public ZeemanBase
{

private:

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true) {}

public:

	Atom_Zeeman(Atom_Mesh *paMesh_) {}
	~Atom_Zeeman() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Energy density methods

	double GetEnergyDensity(Rect& avRect) { return 0.0; }

	//-------------------

	void SetField(DBL3 Hxyz) {}
	DBL3 GetField(void) { return DBL3(); }

	BError SetFieldEquation(string equation_string, int step = 0) { return BError(); }

	//if base temperature changes we need to adjust Tb in H_equation if it's used.
	void SetBaseTemperature(double Temperature) {}
};

#endif