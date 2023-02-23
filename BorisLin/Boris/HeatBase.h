#pragma once

#include "BorisLib.h"



#include "ErrorHandler.h"
#include "Boris_Enums_Defs.h"
#include "Heat_Defs.h"
#include "Modules.h"

class MeshBase;
class SHeat;

#ifdef MODULE_COMPILATION_HEAT

class HeatBase :
	public virtual Modules
{
	friend SHeat;

private:

	MeshBase* pMeshBase;

protected:

	//temperature model type
	int tmtype;

	//evaluate heat equation and store result here. After this is done advance time for temperature based on values stored here.
	std::vector<double> heatEq_RHS;

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

	//Set Q using user equation, thus allowing simultaneous spatial (x, y, z), stage time (t); stage step (Ss) introduced as user constant.
	//A number of constants are always present : mesh dimensions in m (Lx, Ly, Lz)
	TEquation<double, double, double, double> Q_equation;

protected:

	//-------------------Calculation Methods (pure virtual)

	//1-temperature model
	virtual void IterateHeatEquation_1TM(double dT) = 0;

	//2-temperature model : itinerant electrons <-> lattice
	virtual void IterateHeatEquation_2TM(double dT) = 0;

	//------------------Others

	void SetRobinBoundaryConditions(void);

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true);

public:

	HeatBase(MeshBase* pMeshBase_);
	virtual ~HeatBase() {}

	//-------------------CMBND computation methods (pure virtual)

	//CMBND values set based on continuity of temperature and heat flux

	//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
	virtual double afunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double afunc_pri(int cell1_idx, int cell2_idx, DBL3 shift) const = 0;

	virtual double bfunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double bfunc_pri(int cell1_idx, int cell2_idx) const = 0;

	//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC - Q / K - many-temperature model coupling terms / K
	virtual double diff2_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const = 0;
	virtual double diff2_pri(int cell1_idx, DBL3 shift) const = 0;

	//-------------------Setters

	//set values for Robin boundary conditions.
	void SetAmbientTemperature(double T_ambient_);
	void SetAlphaBoundary(double alpha_boundary_);

	//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
	virtual void SetBaseTemperature(double Temperature) = 0;

	//set insulating mesh sides flags. literal can be "x", "-x", "y", "-y", "z", "-z"
	void SetInsulatingSides(std::string literal, bool status);

	//Set Q_equation text equation object
	BError SetQEquation(std::string equation_string, int step);

	//set temperature solver type
	BError Set_TMType(TMTYPE_ tmtype_ = TMTYPE_DEFAULT);

	//-------------------Getters

	double GetAmbientTemperature(void) { return T_ambient; }

	double GetAlphaBoundary(void) { return alpha_boundary; }

	//get status of insulating side. literal can be "x", "-x", "y", "-y", "z", "-z"
	bool GetInsulatingSide(std::string literal)
	{
		if (literal == "x") return insulate_px;
		else if (literal == "-x") return insulate_nx;
		else if (literal == "y") return insulate_py;
		else if (literal == "-y") return insulate_ny;
		else if (literal == "z") return insulate_pz;
		else if (literal == "-z") return insulate_nz;

		return false;
	}

	//get all at once
	std::vector<bool> GetInsulatingSides(void) { return { insulate_px, insulate_nx, insulate_py, insulate_ny, insulate_pz, insulate_nz }; }

	//get the set temperature model type
	int Get_TMType(void) { return tmtype; }

	//-------------------Others (pure virtual)

	//called by MoveMesh method in this mesh - move relevant Heat quantities
	virtual void MoveMesh_Heat(double x_shift) = 0;
};

#else

class HeatBase :
	public virtual Modules
{
	friend SHeat;

private:

	MeshBase* pMeshBase;

protected:

protected:

	//-------------------Calculation Methods (pure virtual)

	//------------------Others


public:

	HeatBase(MeshBase* pMeshBase_) {}
	virtual ~HeatBase() {}

	//-------------------CMBND computation methods (pure virtual)


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
	bool GetInsulatingSide(std::string literal)
	{
		return false;
	}

	//get all at once
	std::vector<bool> GetInsulatingSides(void) { return {}; }

	//get the set temperature model type
	int Get_TMType(void) { return 0; }

	//-------------------Others (pure virtual)

	//called by MoveMesh method in this mesh - move relevant Heat quantities
	void MoveMesh_Heat(double x_shift) {}
};

#endif
