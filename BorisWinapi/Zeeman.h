#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;

#ifdef MODULE_ZEEMAN

//Zeeman module can only be used in a magnetic mesh

class SuperMesh;

class Zeeman : 
	public Modules,
	public ProgramState<Zeeman, tuple<DBL3, TEquation<double, double, double, double>>, tuple<>>
{

#if COMPILECUDA == 1
	friend class ZeemanCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

	//pointer to supermesh
	SuperMesh* pSMesh;

	//Applied field
	DBL3 Ha;

	//Applied field using user equation, thus allowing simultaneous spatial (x, y, z), stage time (t); base temperature (Tb) and stage step (Ss) are introduced as user constants.
	//A number of constants are always present : mesh dimensions in m (Lx, Ly, Lz)
	TEquation<double, double, double, double> H_equation;

private:

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true);

public:

	Zeeman(Mesh *pMesh_);
	~Zeeman();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------

	void SetField(DBL3 Hxyz);
	DBL3 GetField(void);

	BError SetFieldEquation(string equation_string, int step);

	//if base temperature changes we need to adjust Tb in H_equation if it's used.
	void SetBaseTemperature(double Temperature);
};

#else

class Zeeman :
	public Modules
{

private:

public:

	Zeeman(Mesh *pMesh_) {}
	~Zeeman() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------

	void SetField(DBL3 Hxyz) {}
	DBL3 GetField(void) { return DBL3(); }

	BError SetFieldEquation(string equation_string, int step = 0) {}

	//if base temperature changes we need to adjust Tb in H_equation if it's used.
	void SetBaseTemperature(double Temperature) {}
};

#endif