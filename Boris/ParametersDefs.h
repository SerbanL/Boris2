#pragma once

#include <string>

//maximum allowed temperature (K) for scaling array
#define MAX_TEMPERATURE 3000.0

//enum for all mesh material parameters. Add new parameters at the end to keep compatibility with older simulation files.

enum PARAM_ {

	PARAM_ALL = -1,
	PARAM_GREL, PARAM_GDAMPING, PARAM_MS, PARAM_DEMAGXY, PARAM_A, PARAM_D, PARAM_J1, PARAM_J2, PARAM_K1, PARAM_K2, PARAM_EA1, PARAM_EA2, PARAM_SUSREL, PARAM_SUSPREL,
	PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_BETA, PARAM_DE, PARAM_BETAD, PARAM_SHA, PARAM_ISHA, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX, PARAM_TSEFF, PARAM_TSIEFF, PARAM_PUMPEFF,
	PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_FLSOT, PARAM_HA, PARAM_TC, PARAM_MUB, PARAM_T, PARAM_Q,
	PARAM_GREL_AFM, PARAM_GDAMPING_AFM, PARAM_MS_AFM, PARAM_A_AFM, PARAM_A_AFNH, PARAM_D_AFM,
	PARAM_CPUMP_EFF, PARAM_THE_EFF, PARAM_NDENSITY,
	PARAM_MECOEFF, PARAM_YOUNGSMOD, PARAM_POISSONRATIO,
	PARAM_NETADIA,
	PARAM_SHC_E, PARAM_G_E,
	PARAM_A_AFH, PARAM_SUSREL_AFM, PARAM_AFTAU, PARAM_AFTAUCROSS, PARAM_MUB_AFM, PARAM_K1_AFM, PARAM_K2_AFM,
	PARAM_HMO,

	PARAM_ATOM_SC_DAMPING, PARAM_ATOM_SC_MUS, PARAM_ATOM_SC_J, PARAM_ATOM_SC_D, PARAM_ATOM_SC_K, PARAM_ATOM_EA1, PARAM_ATOM_EA2
};

//classification of parameter
enum PARAMTYPE_ { PARAMTYPE_NONE, PARAMTYPE_MAGNETIC, PARAMTYPE_ELECTRIC, PARAMTYPE_THERMAL, PARAMTYPE_MECHANICAL };

//enum for material parameter temperature dependence type
enum MATPTDEP_ { MATPTDEP_NONE = 0, MATPTDEP_ARRAY, MATPTDEP_EQUATION };

//enum for spatial variation types
enum MATPVAR_ { 
	MATPVAR_NONE = 0, MATPVAR_MASK, MATPVAR_RANDOM, MATPVAR_JAGGED, MATPVAR_DEFECTS, MATPVAR_FAULTS, 
	MATPVAR_VORONOI2D, MATPVAR_VORONOI3D, MATPVAR_VORONOIBND2D, MATPVAR_VORONOIBND3D, MATPVAR_VORONOIROT2D, MATPVAR_VORONOIROT3D,
	MATPVAR_EQUATION, MATPVAR_OVF2
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//MeshParamDescriptor used for console display / control of mesh material parameters
//

struct MeshParamDescriptor {

	//the unit used when converting from a string containing units to a numerical value and conversely
	std::string unit;

	PARAMTYPE_ paramType = PARAMTYPE_NONE;

	//display or hide this parameter?
	//e.g. we may enable a parameter for a mesh type, but may want to hide it so it doesn't appear in the usual lists of parameters (params, paramstemp, paramsvar commands).
	bool hidden;

	MeshParamDescriptor(PARAMTYPE_ paramType_, std::string unit_ = "", bool hidden_ = false) :
		unit(unit_), paramType(paramType_), hidden(hidden_)
	{}

	//what parameter is this?
	PARAMTYPE_ get_type(void) { return paramType; }
};
