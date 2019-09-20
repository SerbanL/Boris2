#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class SuperMesh;
class SHeat;

#ifdef MODULE_HEAT

#if COMPILECUDA == 1
#include "HeatCUDA.h"
#endif

class Heat :
	public Modules,
	public ProgramState<Heat, tuple<double, double, bool, bool, bool, bool, bool, bool>, tuple<>>
{
	friend SHeat;

#if COMPILECUDA == 1
	friend HeatCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

	SuperMesh* pSMesh;

	//evaluate heat equation and store result here. After this is done advance time for temperature based on values stored here.
	vector<double> heatEq_RHS;

	//ambient temperature and alpha boundary value used in Robin boundary conditions (Newton's law of cooling):
	//Flux in direction of surface normal = alpha_boundary * (T_boundary - T_ambient)
	//Note : alpha_boundary = 0 results in insulating boundary
	double T_ambient = 297;
	double alpha_boundary = 1e7;

	//alpha_boundary may be non-zero but we can still have insulating boundaries at mesh sides as specified by flags below
	bool insulate_px = false;
	bool insulate_nx = false;
	bool insulate_py = false;
	bool insulate_ny = false;
	bool insulate_pz = false;
	bool insulate_nz = false;

private:

	//-------------------Calculation Methods

	void IterateHeatEquation(double dT);

	//------------------Others

	void SetRobinBoundaryConditions(void);

public:

	Heat(Mesh *pMesh_);
	~Heat();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------CMBND computation methods

	//CMBND values set based on continuity of temperature and heat flux

	//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
	double afunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double afunc_pri(int cell1_idx, int cell2_idx, DBL3 shift) const;

	double bfunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double bfunc_pri(int cell1_idx, int cell2_idx) const;

	//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC
	double diff2_sec(DBL3 relpos_m1, DBL3 stencil) const;
	double diff2_pri(int cell1_idx) const;

	//-------------------Setters

	//set values for Robin boundary conditions.
	void SetAmbientTemperature(double T_ambient_);
	void SetAlphaBoundary(double alpha_boundary_);

	//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
	void SetBaseTemperature(double Temperature);

	//set insulating mesh sides flags. literal can be "x", "-x", "y", "-y", "z", "-z"
	void SetInsulatingSides(string literal, bool status);

	//-------------------Getters

	double GetAmbientTemperature(void) { return T_ambient; }

	double GetAlphaBoundary(void) { return alpha_boundary; }

	//get status of insulating side. literal can be "x", "-x", "y", "-y", "z", "-z"
	bool GetInsulatingSide(string literal);
	//get all at once
	vector<bool> GetInsulatingSides(void) { return { insulate_px, insulate_nx, insulate_py, insulate_ny, insulate_pz, insulate_nz }; }

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant Heat quantities
	void MoveMesh_Heat(double x_shift);
};

#else

class Heat :
	public Modules
{

private:

private:

public:

	Heat(Mesh *pMesh_) {}
	~Heat() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }


	//-------------------Setters

	//set values for Robin boundary conditions.
	void SetAmbientTemperature(double T_ambient_) {}
	void SetAlphaBoundary(double alpha_boundary_) {}

	//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
	void SetBaseTemperature(double Temperature) {}

	//set insulating mesh sides flags. literal can be "x", "-x", "y", "-y", "z", "-z"
	void SetInsulatingSides(string literal, bool status) {}

	//-------------------Getters

	double GetAmbientTemperature(void) { return 0.0; }

	double GetAlphaBoundary(void) { return 0.0; }

	//get status of insulating side. literal can be "x", "-x", "y", "-y", "z", "-z"
	bool GetInsulatingSide(string literal) { return true; }
	//get all at once
	vector<bool> GetInsulatingSides(void) { return {}; }

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant Heat quantities
	void MoveMesh_Heat(double x_shift) {}
};

#endif