#include "stdafx.h"
#include "MeshParams.h"

MeshParams::MeshParams(vector<PARAM_>& enabledParams)
{
	//store all simulation parameters in meshParams : allows easy handling in the console, as well as saving and loading. In the computational routines use the parmeters directly, not through simParams
	for (int idx = 0; idx < (int)enabledParams.size(); idx++) {

		//List all parameters here, but only add them if configured to
		switch (enabledParams[idx]) {

		case PARAM_GREL:
			meshParams.push_back("grel", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GREL);
			break;

		case PARAM_GDAMPING:
			meshParams.push_back("damping", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GDAMPING);
			break;

		case PARAM_MS:
			meshParams.push_back("Ms", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "A/m"), PARAM_MS);
			break;

		case PARAM_DEMAGXY:
			meshParams.push_back("Nxy", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_DEMAGXY);
			break;

		case PARAM_A:
			meshParams.push_back("A", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m"), PARAM_A);
			break;

		case PARAM_D:
			meshParams.push_back("D", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_D);
			break;

		case PARAM_J1:
			meshParams.push_back("J1", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_J1);
			break;

		case PARAM_J2:
			meshParams.push_back("J2", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_J2);
			break;

		case PARAM_K1:
			meshParams.push_back("K1", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_K1);
			break;

		case PARAM_K2:
			meshParams.push_back("K2", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_K2);
			break;

		case PARAM_TC:
			//add it to mesh but mark it as hidden
			meshParams.push_back("Tc", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "K", true), PARAM_TC);
			break;

		case PARAM_MUB:
			//add it to mesh but mark it as hidden
			meshParams.push_back("muB", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "K", true), PARAM_MUB);
			break;

		case PARAM_EA1:
			meshParams.push_back("ea1", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_EA1);
			break;

		case PARAM_EA2:
			meshParams.push_back("ea2", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_EA2);
			break;

		case PARAM_SUSREL:
			meshParams.push_back("susrel", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "As2/kg"), PARAM_SUSREL);
			break;

		case PARAM_SUSPREL:
			meshParams.push_back("susprel", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "As2/kg"), PARAM_SUSPREL);
			break;

		case PARAM_HA:
			meshParams.push_back("cHa", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_HA);
			break;

		case PARAM_T:
			meshParams.push_back("cT", MeshParamDescriptor(PARAMTYPE_THERMAL), PARAM_T);
			break;

		case PARAM_Q:
			meshParams.push_back("Q", MeshParamDescriptor(PARAMTYPE_THERMAL, "W/m3"), PARAM_Q);
			break;

		case PARAM_ELC:
			meshParams.push_back("elC", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "S/m"), PARAM_ELC);
			break;

		case PARAM_AMR:
			meshParams.push_back("amr", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "%"), PARAM_AMR);
			break;

		case PARAM_P:
			meshParams.push_back("P", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_P);
			break;

		case PARAM_BETA:
			meshParams.push_back("beta", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_BETA);
			break;

		case PARAM_DE:
			meshParams.push_back("De", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m2/s"), PARAM_DE);
			break;

		case PARAM_BETAD:
			meshParams.push_back("betaD", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_BETAD);
			break;

		case PARAM_SHA:
			meshParams.push_back("SHA", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_SHA);
			break;

		case PARAM_FLSOT:
			meshParams.push_back("flST", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_FLSOT);
			break;

		case PARAM_ISHA:
			meshParams.push_back("iSHA", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_ISHA);
			break;

		case PARAM_LSF:
			meshParams.push_back("l_sf", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m"), PARAM_LSF);
			break;

		case PARAM_LEX:
			meshParams.push_back("l_J", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m"), PARAM_LEX);
			break;

		case PARAM_LPH:
			meshParams.push_back("l_phi", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m"), PARAM_LPH);
			break;

		case PARAM_GI:
			meshParams.push_back("Gi", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "S/m2"), PARAM_GI);
			break;

		case PARAM_GMIX:
			meshParams.push_back("Gmix", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "S/m2"), PARAM_GMIX);
			break;

		case PARAM_TSEFF:
			meshParams.push_back("ts_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_TSEFF);
			break;

		case PARAM_TSIEFF:
			meshParams.push_back("tsi_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_TSIEFF);
			break;

		case PARAM_PUMPEFF:
			meshParams.push_back("pump_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_PUMPEFF);
			break;

		case PARAM_THERMCOND:
			meshParams.push_back("thermK", MeshParamDescriptor(PARAMTYPE_THERMAL, "W/mK"), PARAM_THERMCOND);
			break;

		case PARAM_DENSITY:
			meshParams.push_back("density", MeshParamDescriptor(PARAMTYPE_MECHANICAL, "kg/m3"), PARAM_DENSITY);
			break;

		case PARAM_SHC:
			meshParams.push_back("shc", MeshParamDescriptor(PARAMTYPE_THERMAL, "J/kgK"), PARAM_SHC);
			break;
		}
	}
}

//-------------------------Parameter control

void MeshParams::update_parameters(PARAM_ paramID)
{
	auto update_param = [&](PARAM_ update_paramID) {

		auto code = [&](auto& MatP_object) -> void {

			MatP_object.update(base_temperature);
		};

		run_on_param<void>(update_paramID, code);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			update_param((PARAM_)meshParams.get_ID_from_index(index));
		}
	}
	else update_param(paramID);
}

//copy all parameters from another Mesh
void MeshParams::copy_parameters(MeshParams& copy_this)
{
	auto copy_param = [&](PARAM_ paramID) {
		
		auto copy_from = [&](auto& MatP_object, MeshParams& copy_this) {

			auto get_meshparam = [&](decltype(MatP_object)& MatP_object2) -> decltype(MatP_object)& {

				return MatP_object2;
			};

			MatP_object = copy_this.run_on_param<decltype(MatP_object)>(paramID, get_meshparam);
		};

		run_on_param<void>(paramID, copy_from, copy_this);
	};

	//copy all MatP parameters
	for (int index = 0; index < meshParams.size(); index++) {

		PARAM_ paramID = (PARAM_)meshParams.get_ID_from_index(index);

		copy_param(paramID);
	}

	//now copy special values which are not MatP
	base_temperature = copy_this.base_temperature;
	T_Curie = copy_this.T_Curie;
}

//-------------------------Getters

//get value of indexed mesh parameter as a string (with unit)
string MeshParams::get_meshparam_value(int index)
{
	PARAM_ paramID = (PARAM_)get_meshparam_id(index);

	auto code = [&](auto& MatP_object, string unit) -> string {

		return ToString(MatP_object.get0(), unit);
	};

	return run_on_param<string>(paramID, code, meshParams(paramID).unit);
}

string MeshParams::get_meshparam_value(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object, string unit) -> string {

		return ToString(MatP_object.get0(), unit);
	};

	return run_on_param<string>(paramID, code, meshParams(paramID).unit);
}

//get value of indexed mesh parameter as a string (without unit)
string MeshParams::get_meshparam_value_sci(int index)
{
	PARAM_ paramID = (PARAM_)get_meshparam_id(index);

	auto code = [&](auto& MatP_object) -> string {

		return ToString(MatP_object.get0());
	};

	return run_on_param<string>(paramID, code);
}

string MeshParams::get_meshparam_value_sci(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object) -> string {

		return ToString(MatP_object.get0());
	};

	return run_on_param<string>(paramID, code);
}

string MeshParams::get_paraminfo_string(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object, PARAM_ paramID) -> string {

		return get_meshparam_handle(paramID) + ": " + MatP_object.get_info_string();
	};

	return run_on_param<string>(paramID, code, paramID);
}

//returns a string describing the set spatial dependence with any parameters
string MeshParams::get_paramvarinfo_string(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object, PARAM_ paramID) -> string {

		return get_meshparam_handle(paramID) + ": " + MatP_object.get_varinfo_string();
	};

	return run_on_param<string>(paramID, code, paramID);
}

bool MeshParams::is_paramtemp_set(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return MatP_object.is_tdep();
	};

	return run_on_param<bool>(paramID, code);
}

//check if the given parameter has a  spatial variation set
bool MeshParams::is_paramvar_set(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return MatP_object.is_sdep();
	};

	return run_on_param<bool>(paramID, code);
}

//check if the given parameter has a  temperature dependence or a spatial variation set
bool MeshParams::is_param_nonconst(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return (MatP_object.is_tdep() || MatP_object.is_sdep());
	};

	return run_on_param<bool>(paramID, code);
}

vector<double> MeshParams::get_meshparam_tempscaling(PARAM_ paramID, double max_temperature)
{
	auto code = [](auto& MatP_object, double max_temperature) -> vector<double> {

		return MatP_object.get_temperature_scaling(max_temperature);
	};

	return run_on_param<vector<double>>(paramID, code, max_temperature);
}

//get reference to mesh parameter spatial scaling VEC
void* MeshParams::get_meshparam_s_scaling(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> void* {

		return &(MatP_object.s_scaling_ref());
	};

	return run_on_param<void*>(paramID, code);
}

//check if scaling array is scalar (or else vectorial)
bool MeshParams::is_s_scaling_scalar(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		typedef contained_type<decltype(MatP_object.s_scaling_ref())>::type SType;

		return std::is_same<SType, double>::value;
	};

	return run_on_param<bool>(paramID, code);
}

//-------------------------Setters

//set value from string for named parameter (units allowed in string)
void MeshParams::set_meshparam_value(PARAM_ paramID, string value_text)
{
	auto code = [](auto& MatP_object, string value_text, string unit) -> void {

		MatP_object = (decltype(MatP_object.get0()))ToNum(value_text, unit);
	};

	return run_on_param<void>(paramID, code, value_text, meshParams(paramID).unit);

	update_parameters(paramID);
}

//set the mesh parameter named formula (see handles in MaterialsParameterFormulas.h) with given coefficients
void MeshParams::set_meshparam_formula(PARAM_ paramID, MATPFORM_ formulaID, vector<double> coefficients)
{
	auto code = [](auto& MatP_object, MATPFORM_ formulaID, vector<double>& coefficients) -> void {

		MatP_object.set_scaling_formula(formulaID, coefficients);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param<void>((PARAM_)meshParams.get_ID_from_index(index), code, formulaID, coefficients);
		}
	}
	else run_on_param<void>(paramID, code, formulaID, coefficients);
}

//set mesh parameter array scaling
bool MeshParams::set_meshparam_tscaling_array(PARAM_ paramID, vector<double>& temp, vector<double>& scaling)
{
	auto code = [](auto& MatP_object, vector<double>& temp, vector<double>& scaling) -> bool {

		return MatP_object.set_scaling_array(temp, scaling);
	};

	return run_on_param<bool>(paramID, code, temp, scaling);
}

//-------------------------Setters : spatial variation

//clear mesh parameter spatial variation (all if paramID == PARAM_ALL)
void MeshParams::clear_meshparam_variation(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> void {

		MatP_object.clear_s_scaling();
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param<void>((PARAM_)meshParams.get_ID_from_index(index), code);
		}
	}
	else run_on_param<void>(paramID, code);
}

//update mesh parameter spatial variation (e.g. cellsize or rectangle could have changed)
bool MeshParams::update_meshparam_var(PARAM_ paramID, DBL3 h, Rect rect)
{
	auto code = [](auto& MatP_object, DBL3 h, Rect rect) -> bool {

		return MatP_object.update_s_scaling(h, rect);
	};

	return run_on_param<bool>(paramID, code, h, rect);
}

//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
BError MeshParams::set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, DBL3 h, Rect rect, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	auto code = [](auto& MatP_object, DBL3 h, Rect rect, MATPVAR_ generatorID, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader) -> BError {

		return MatP_object.set_s_scaling(h, rect, generatorID, generatorArgs, bitmap_loader);
	};

	error = run_on_param<BError>(paramID, code, h, rect, generatorID, generatorArgs, bitmap_loader);

	return error;
}