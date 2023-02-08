#pragma once

#include "BorisLib.h"
#include "Modules.h"
#include "Boris_Enums_Defs.h"
#include "HeatBase.h"

class Mesh;
class SuperMesh;

#ifdef MODULE_COMPILATION_HEAT

#if COMPILECUDA == 1
#include "HeatCUDA.h"
#endif

class Heat :
	public virtual Modules,
	public HeatBase,
	public ProgramState<Heat, std::tuple<
	int, 
	double, double, 
	bool, bool, bool, bool, bool, bool, 
	TEquation<double, double, double, double>>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend HeatCUDA;
#endif

private:

	//pointer to meshbase object holding this effective field module
	Mesh* pMesh;

	SuperMesh* pSMesh;

private:

	//-------------------Calculation Methods

	//1-temperature model
	void IterateHeatEquation_1TM(double dT);
	
	//2-temperature model : itinerant electrons <-> lattice
	void IterateHeatEquation_2TM(double dT);

public:

	Heat(Mesh *pMesh_);
	~Heat();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------CMBND computation methods

	//CMBND values set based on continuity of temperature and heat flux

	//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
	double afunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double afunc_pri(int cell1_idx, int cell2_idx, DBL3 shift) const;

	double bfunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double bfunc_pri(int cell1_idx, int cell2_idx) const;

	//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC - Q / K - many-temperature model coupling terms / K
	double diff2_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const;
	double diff2_pri(int cell1_idx, DBL3 shift) const;

	//-------------------Setters

	//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
	void SetBaseTemperature(double Temperature);

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

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }


	//-------------------Setters

	//set values for Robin boundary conditions.
	void SetAmbientTemperature(double T_ambient_) {}
	void SetAlphaBoundary(double alpha_boundary_) {}

	//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
	void SetBaseTemperature(double Temperature) {}

	//set insulating mesh sides flags. literal can be "x", "-x", "y", "-y", "z", "-z"
	void SetInsulatingSides(std::string literal, bool status) {}

	//Set Q_equation text equation object
	BError SetQEquation(std::string equation_string, int step) { return BError(); }

	//set temperature solver type
	BError Set_TMType(TMTYPE_ tmtype_ = TMTYPE_DEFAULT) { return BError(); }

	//-------------------Getters

	double GetAmbientTemperature(void) { return 0.0; }

	double GetAlphaBoundary(void) { return 0.0; }

	//get status of insulating side. literal can be "x", "-x", "y", "-y", "z", "-z"
	bool GetInsulatingSide(std::string literal) { return true; }
	//get all at once
	std::vector<bool> GetInsulatingSides(void) { return {}; }

	//get the set temperature model type
	int Get_TMType(void) { return 0; }

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant Heat quantities
	void MoveMesh_Heat(double x_shift) {}
};

#endif