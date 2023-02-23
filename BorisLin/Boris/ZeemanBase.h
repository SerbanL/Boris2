#pragma once

#include "BorisLib.h"
#include "Boris_Enums_Defs.h"
#include "ErrorHandler.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Abstract base class for Zeeman modules (there are more than one: micromagnetic, atomistic)
//
//	All modules are implementations of Modules ABC, so this is an additional ABC (so multiple inheritance).
//	The main reason for this, it allows calling a Zeeman-type method (defined here, implemented in the specific Zeeman type module) using a Modules base pointer.
//	This is very convenient since Meshes hold Modules base pointers, so e.g. if you want to call the SetField method on a mesh, this will only take effect if a Zeeman-type module is active in it.
//  Thus this works irrespective of the particular mesh type (micromagnetic or atomistic), which have different implementations of ZeemanBase, and is accomplished by checking return of 
//	dynamic_cast from Modules* to ZeemanBase*.
//	Also defines common data types and common interface, as these modules have the same structure but do things in a different way.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ZeemanBase {

protected:

	//Applied field
	DBL3 Ha;

	//Applied field using user equation, thus allowing simultaneous spatial (x, y, z), stage time (t); base temperature (Tb) and stage step (Ss) are introduced as user constants.
	//A number of constants are always present : mesh dimensions in m (Lx, Ly, Lz)
	TEquation<double, double, double, double> H_equation;

	//Applied field but as a VEC (e.g. loaded from file), with same resolution as M.
	VEC<DBL3> Havec;

	//global field as obtained in this mesh
	//When a global field is set in supermesh, then globalField VEC is initialized here with mesh transfered values
	VEC<DBL3> globalField;

protected:

	//Update TEquation object with user constants values
	virtual void UpdateTEquationUserConstants(bool makeCuda = true) = 0;

	//setup globalField transfer
	virtual BError InitializeGlobalField(void) = 0;

public:

	ZeemanBase(void) :
		H_equation({ "x", "y", "z", "t" })
	{
		Ha = DBL3(0, 0, 0);
	}

	virtual ~ZeemanBase() {}

	//-------------------

	virtual void SetField(DBL3 Hxyz) = 0;
	virtual DBL3 GetField(void) = 0;

	virtual BError SetFieldEquation(std::string equation_string, int step) = 0;

	virtual BError SetFieldVEC_FromOVF2(std::string fileName) = 0;

	//if base temperature changes we need to adjust Tb in H_equation if it's used.
	virtual void SetBaseTemperature(double Temperature) = 0;
};
