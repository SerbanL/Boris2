#pragma once

//maximum allowed temperature (K) for scaling array
#define MAX_TEMPERATURE 3000.0

//enum for all mesh material parameters. Add new parameters at the end to keep compatibility with older simulation files.
//
//Some parameters are "special", i.e. they do not appear in the list of parameters in meshes (thus cannot be assigned temperature or spatial dependences).
//They are flagged as hidden in the Mesh constructor.
//We still want to save/load them from materials database however.
//These are designated as SPARAM_...
//
enum PARAM_ {

	PARAM_ALL = -1,
	PARAM_GREL, PARAM_GDAMPING, PARAM_MS, PARAM_DEMAGXY, PARAM_A, PARAM_D, PARAM_J1, PARAM_J2, PARAM_K1, PARAM_K2, PARAM_EA1, PARAM_EA2, PARAM_SUSREL, PARAM_SUSPREL,
	PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_BETA, PARAM_DE, PARAM_BETAD, PARAM_SHA, PARAM_ISHA, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX, PARAM_TSEFF, PARAM_TSIEFF, PARAM_PUMPEFF,
	PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_FLSOT, PARAM_HA, PARAM_TC, PARAM_MUB, PARAM_T, PARAM_Q,
	PARAM_GREL_AFM, PARAM_GDAMPING_AFM, PARAM_MS_AFM, PARAM_A_AFM, PARAM_A12, PARAM_D_AFM,
	PARAM_CPUMP_EFF, PARAM_THE_EFF, PARAM_NDENSITY
};

//classification of parameter
enum PARAMTYPE_ { PARAMTYPE_NONE, PARAMTYPE_MAGNETIC, PARAMTYPE_ELECTRIC, PARAMTYPE_THERMAL, PARAMTYPE_MECHANICAL };

//enum for all analytical formulas for parameter temperature scaling
enum MATPFORM_ { MATPFORM_NONE = 0, MATPFORM_LINEAR, MATPFORM_PARABOLIC, MATPFORM_INVERSELINEAR };

//enum for spatial variation types
enum MATPVAR_ { 
	MATPVAR_NONE = 0, MATPVAR_CUSTOM, MATPVAR_RANDOM, MATPVAR_JAGGED, MATPVAR_DEFECTS, MATPVAR_FAULTS, 
	MATPVAR_VORONOI2D, MATPVAR_VORONOI3D, MATPVAR_VORONOIBND2D, MATPVAR_VORONOIBND3D, MATPVAR_VORONOIROT2D, MATPVAR_VORONOIROT3D	
};

//maximum number of coefficients used for MATPFORM
#define MAXFORMULACOEFFICIENTS	2
