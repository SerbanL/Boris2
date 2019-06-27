#pragma once

#include "BorisLib.h"

#include "MaterialParameter.h"
#include "ErrorHandler.h"
#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1
#include "MeshParamsCUDA.h"
#endif

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////
//MeshParamDescriptor used for console display / control of mesh material parameters
//

struct MeshParamDescriptor {

	//the unit used when converting from a string containing units to a numerical value and conversely
	string unit;

	PARAMTYPE_ paramType = PARAMTYPE_NONE;

	//display or hide this parameter?
	//e.g. we may enable a parameter for a mesh type, but may want to hide it so it doesn't appear in the usual lists of parameters (params, paramstemp, paramsvar commands).
	bool hidden;

	MeshParamDescriptor(PARAMTYPE_ paramType_, string unit_ = "", bool hidden_ = false) :
		unit(unit_), paramType(paramType_), hidden(hidden_)
	{}

	//what parameter is this?
	PARAMTYPE_ get_type(void) { return paramType; }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//The Mesh material parameters
//

class MeshParams
{

#if COMPILECUDA == 1
	friend MeshParamsCUDA;
#endif

private:

	//list of all mesh parameters storing the units, handles and parameter type : used for console info display
	vector_key_lut<MeshParamDescriptor> meshParams;

public:

	//-------------------------------- LIST ALL MESH PARAMETERS HERE

	//Relative electron gyromagnetic ratio
	MatP<double, double> grel = 1.0;

	//Gilbert damping
	MatP<double, double> alpha = 0.02;

	//Saturation magnetisation (A/m)
	MatP<double, double> Ms = 8e5;

	//in-plane demagnetizing factors (used for Demag_N module)
	MatP<DBL2, double> Nxy = DBL2(0);

	//Exchange stifness (J/m)
	MatP<double, double> A = 1.3e-11;

	//Dzyaloshinskii-Moriya exchange constant (J/m^2)
	MatP<double, double> D = 3e-3;

	//bilinear surface exchange coupling (J/m^2) : J1
	//biquadratic surface exchange coupling (J/m^2) : J2
	//For coupled meshes it is the top mesh that sets the J values.
	MatP<double, double> J1 = -1e-3;
	MatP<double, double> J2 = 0;

	//Magneto-crystalline anisotropy K1 and K2 constants (J/m^3) and easy axes directions. For uniaxial anisotropy only ea1 is needed, for cubic ea1 and ea2 should be orthogonal.
	MatP<double, double> K1 = 1e4;
	MatP<double, double> K2 = 0;
	MatP<DBL3, DBL3> mcanis_ea1 = DBL3(1, 0, 0);
	MatP<DBL3, DBL3> mcanis_ea2 = DBL3(0, 1, 0);

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	MatP<double, double> susrel = 1.0;

	//perpendicular (transverse) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	MatP<double, double> susprel = 1.0;

	//applied field spatial variation coefficient (unitless)
	MatP<double, double> cHA = 1.0;

	//electrical conductivity (units S/m).
	//this is the value at RT for Ni80Fe20.
	MatP<double, double> elecCond = 7e6;

	//anisotropic magnetoresistance as a percentage (of base resistance)
	MatP<double, double> amrPercentage = 0.0;

	//spin current polarization (also the charge current spin polarization for the spin current solver) and non-adiabaticity (for Zhang-Li STT). (unitless)
	MatP<double, double> P = 0.4;
	MatP<double, double> beta = 0.04;

	//parameters for spin current solver

	//electron diffusion constant (m^2/s)
	MatP<double, double> De = 1e-2;

	//diffusion spin polarisation (unitless)
	MatP<double, double> betaD = 0.5;

	//spin Hall angle (unitless)
	MatP<double, double> SHA = 0.1;

	//"inverse" spin Hall angle (unitless) -> should normally be the same as SHA but e.g. can be set to zero to turn off the inverse SHE in the spin transport equation
	MatP<double, double> iSHA = 0.1;

	//field-like spin torque coefficient (unitless)
	MatP<double, double> flSOT = 0.0;

	//spin-flip length (m)
	MatP<double, double> l_sf = 10e-9;

	//spin exchange rotation length (m)
	MatP<double, double> l_ex = 2e-9;

	//spin dephasing length (m)
	MatP<double, double> l_ph = 4e-9;

	//interface spin-dependent conductivity (spin-up and spin-down) (S/m^2)
	MatP<DBL2, double> Gi = DBL2(1e15, 1e14);

	//interface spin-mixing conductivity (real and imaginary parts) (S/m^2)
	MatP<DBL2, double> Gmix = DBL2(1e15, 1e14);

	//spin accumulation torque efficiency in the bulk (unitless, varies from 0 : no torque, up to 1 : full torque)
	MatP<double, double> ts_eff = 1;

	//spin accumulation torque efficiency at interfaces (unitless, varies from 0 : no torque, up to 1 : full torque)
	MatP<double, double> tsi_eff = 1;

	//spin pumping efficiency (unitless, varies from 0 : no spin pumping, up to 1 : full strength)
	MatP<double, double> pump_eff = 1;

	//the mesh base temperature (K)
	double base_temperature = 0.0;

	//Curie temperature - 870K for permalloy but turn it off by default. If LLG is the default equation we don't want temperature dependencies to be updated every time the applied field changes.
	//This is the actually set value
	double T_Curie = 0.0;

	//This is the indicative Curie temperature of the material, but not used in any calculations.
	//If you want to turn default temperature dependeces on, this is the value you should set in T_Curie.
	MatP<double, double> T_Curie_material = 870;

	//The atomic magnetic moment as a multiple of the Bohr magneton - default 1 ub for permalloy.
	MatP<double, double> atomic_moment = 1.0;

	//thermal conductivity (W/mK) - default for permalloy
	MatP<double, double> thermCond = 46.4;

	//mass density (kg/m^3) - default for permalloy
	MatP<double, double> density = 8740;

	//specific heat capacity (J/kgK) - default for permalloy
	MatP<double, double> shc = 430;

private:

	//-------------------------Parameter control

	//run set code (run_this) using set parameters (run_this_args) on a MatP object selected through the paramID selector (majorID in meshParams identifying the required material parameter).
	template <typename RType, typename Lambda, typename ... PType>
	RType run_on_param(PARAM_ paramID, Lambda& run_this, PType& ... run_this_args);

protected:

	//-------------------------Setters

	//call this to update given parameter output value to current base_temperature
	void update_parameters(PARAM_ paramID = PARAM_ALL);

	//-------------------------Info (Getters)

	//get reference to mesh parameter spatial scaling VEC
	//use void* since the scaling is templated. Caller must recast it correctly - see is_s_scaling_scalar method
	void* get_meshparam_s_scaling(PARAM_ paramID);

	//check if scaling array is scalar (or else vectorial)
	bool is_s_scaling_scalar(PARAM_ paramID);

public:

	//------------------------CTOR/DTOR

	//inherited by Mesh implementations
	MeshParams(vector<PARAM_>& enabledParams);
	virtual ~MeshParams() {}

	//-------------------------Getters

	//return number of mesh parameters
	int get_num_meshparams(void) { return meshParams.size(); }

	//get a reference to the mesh param
	template <typename RType>
	RType* get_meshparam_pointer(PARAM_ paramID) 
	{ 
		auto code = [](auto& MatP_object) -> RType* { return reinterpret_cast<RType*>(&MatP_object); };
		return run_on_param<RType*>(paramID, code);
	}

	//get id of indexed mesh parameter (this is the value from PARAM_ enum)
	int get_meshparam_id(int index) { return meshParams.get_ID_from_index(index); }
	int get_meshparam_id(string paramHandle) { return meshParams.get_ID_from_key(paramHandle); }

	//get handle of indexed mesh parameter
	string get_meshparam_handle(int index) { return meshParams.get_key_from_index(index); }
	string get_meshparam_handle(PARAM_ paramID) { return meshParams.get_key_from_ID(paramID); }

	//get unit of indexed mesh parameter
	string get_meshparam_unit(int index) { return meshParams[index].unit; }
	string get_meshparam_unit(PARAM_ paramID) { return meshParams(paramID).unit; }

	bool contains_param(PARAM_ paramID) { return meshParams.is_ID_set(paramID); }
	bool contains_param(string paramHandle) { return meshParams.has_key(paramHandle); }

	//get value of indexed mesh parameter as a string (with unit)
	string get_meshparam_value(int index);
	string get_meshparam_value(PARAM_ paramID);

	//get value of indexed mesh parameter as a string (without unit)
	string get_meshparam_value_sci(int index);
	string get_meshparam_value_sci(PARAM_ paramID);

	PARAMTYPE_ get_meshparam_type(PARAM_ paramID) { return meshParams(paramID).get_type(); }

	//returns a string describing the set temperature dependence ("none", "array" or set formula : "name parameters...") 
	string get_paraminfo_string(PARAM_ paramID);

	//returns a string describing the set spatial dependence with any parameters
	string get_paramvarinfo_string(PARAM_ paramID);

	//check if the given parameter has a temperature dependence set
	bool is_paramtemp_set(PARAM_ paramID);
	//check if the given parameter has a  spatial variation set
	bool is_paramvar_set(PARAM_ paramID);
	//check if the given parameter has a  temperature dependence or a spatial variation set
	bool is_param_nonconst(PARAM_ paramID);

	//get mesh parameter temperature scaling up to max_temperature : return a vector from 0K up to and including max_temperature with scaling coefficients
	vector<double> get_meshparam_tempscaling(PARAM_ paramID, double max_temperature);

	//is this param hidden or can we display it?
	bool is_param_hidden(PARAM_ paramID) { return meshParams(paramID).hidden; }

	//-------------------------Setters

	//copy all parameters from another Mesh
	void copy_parameters(MeshParams& copy_this);

	//-------------------------Setters : value and temperature dependence

	//set value from string for named parameter (units allowed in string)
	void set_meshparam_value(PARAM_ paramID, string value_text);

	//set the mesh parameter formula with given coefficients
	void set_meshparam_formula(PARAM_ paramID, MATPFORM_ formulaID, vector<double> coefficients);

	//set mesh parameter array scaling
	bool set_meshparam_tscaling_array(PARAM_ paramID, vector<double>& temp, vector<double>& scaling);

	//-------------------------Setters : spatial variation

	//clear mesh parameter spatial variation (all if paramID == PARAM_ALL)
	void clear_meshparam_variation(PARAM_ paramID);

	//update mesh parameter spatial variation (e.g. cellsize or rectangle could have changed)
	bool update_meshparam_var(PARAM_ paramID, DBL3 h, Rect rect);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
	BError set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, DBL3 h, Rect rect, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader);
};

//-------------------------Parameter control

template <typename RType, typename Lambda, typename ... PType>
RType MeshParams::run_on_param(PARAM_ paramID, Lambda& run_this, PType& ... run_this_args)
{
	switch (paramID) {

	case PARAM_GREL:
		return run_this(grel, run_this_args...);
		break;

	case PARAM_GDAMPING:
		return run_this(alpha, run_this_args...);
		break;

	case PARAM_MS:
		return run_this(Ms, run_this_args...);
		break;

	case PARAM_DEMAGXY:
		return run_this(Nxy, run_this_args...);
		break;

	case PARAM_A:
		return run_this(A, run_this_args...);
		break;

	case PARAM_D:
		return run_this(D, run_this_args...);
		break;

	case PARAM_J1:
		return run_this(J1, run_this_args...);
		break;

	case PARAM_J2:
		return run_this(J2, run_this_args...);
		break;

	case PARAM_K1:
		return run_this(K1, run_this_args...);
		break;

	case PARAM_K2:
		return run_this(K2, run_this_args...);
		break;

	case PARAM_EA1:
		return run_this(mcanis_ea1, run_this_args...);
		break;

	case PARAM_EA2:
		return run_this(mcanis_ea2, run_this_args...);
		break;

	case SPARAM_TC:
		return run_this(T_Curie_material, run_this_args...);
		break;

	case SPARAM_MUB:
		return run_this(atomic_moment, run_this_args...);
		break;

	case PARAM_SUSREL:
		return run_this(susrel, run_this_args...);
		break;

	case PARAM_SUSPREL:
		return run_this(susprel, run_this_args...);
		break;

	case PARAM_HA:
		return run_this(cHA, run_this_args...);
		break;

	case PARAM_ELC:
		return run_this(elecCond, run_this_args...);
		break;

	case PARAM_AMR:
		return run_this(amrPercentage, run_this_args...);
		break;

	case PARAM_P:
		return run_this(P, run_this_args...);
		break;

	case PARAM_BETA:
		return run_this(beta, run_this_args...);
		break;

	case PARAM_DE:
		return run_this(De, run_this_args...);
		break;

	case PARAM_BETAD:
		return run_this(betaD, run_this_args...);
		break;

	case PARAM_SHA:
		return run_this(SHA, run_this_args...);
		break;

	case PARAM_FLSOT:
		return run_this(flSOT, run_this_args...);
		break;

	case PARAM_ISHA:
		return run_this(iSHA, run_this_args...);
		break;

	case PARAM_LSF:
		return run_this(l_sf, run_this_args...);
		break;

	case PARAM_LEX:
		return run_this(l_ex, run_this_args...);
		break;

	case PARAM_LPH:
		return run_this(l_ph, run_this_args...);
		break;

	case PARAM_GI:
		return run_this(Gi, run_this_args...);
		break;

	case PARAM_GMIX:
		return run_this(Gmix, run_this_args...);
		break;

	case PARAM_TSEFF:
		return run_this(ts_eff, run_this_args...);
		break;

	case PARAM_TSIEFF:
		return run_this(tsi_eff, run_this_args...);
		break;

	case PARAM_PUMPEFF:
		return run_this(pump_eff, run_this_args...);
		break;

	case PARAM_THERMCOND:
		return run_this(thermCond, run_this_args...);
		break;

	case PARAM_DENSITY:
		return run_this(density, run_this_args...);
		break;

	case PARAM_SHC:
		return run_this(shc, run_this_args...);
		break;

	default:
		//this is needed to stop the "not all control paths return a value" error, but should never get here
		return run_this(grel, run_this_args...);
		break;
	}
}