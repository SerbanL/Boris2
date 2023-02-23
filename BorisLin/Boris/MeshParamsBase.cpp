#include "stdafx.h"
#include "MeshParamsBase.h"
#include "MeshParams.h"
#include "Atom_MeshParams.h"

//-------------------------Parameter control

template <typename RType, typename Lambda, typename ... PType>
RType MeshParamsBase::run_on_param_switch(PARAM_ paramID, Lambda& run_this, PType& ... run_this_args)
{
	switch (meshBaseType)
	{
	default:
	case MESHTYPE_MICROMAGNETIC:
		return reinterpret_cast<MeshParams*>(pImplementation)->run_on_param<RType>(paramID, run_this, run_this_args...);
		break;

#if ATOMISTIC == 1
	case MESHTYPE_ATOMISTIC:
		return reinterpret_cast<Atom_MeshParams*>(pImplementation)->run_on_param<RType>(paramID, run_this, run_this_args...);
		break;
#endif
	}
}

//-------------------------Getters

//get value of indexed mesh parameter as a std::string (with unit)
std::string MeshParamsBase::get_meshparam_value(int index)
{
	PARAM_ paramID = (PARAM_)get_meshparam_id(index);

	auto code = [&](auto& MatP_object, std::string unit) -> std::string {

		return ToString(MatP_object.get0(), unit);
	};

	return run_on_param_switch<std::string>(paramID, code, meshParams(paramID).unit);
}

std::string MeshParamsBase::get_meshparam_value(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object, std::string unit) -> std::string {

		return ToString(MatP_object.get0(), unit);
	};

	return run_on_param_switch<std::string>(paramID, code, meshParams(paramID).unit);
}

//get value of indexed mesh parameter as a std::string (without unit)
std::string MeshParamsBase::get_meshparam_value_sci(int index)
{
	PARAM_ paramID = (PARAM_)get_meshparam_id(index);

	auto code = [&](auto& MatP_object) -> std::string {

		return ToString(MatP_object.get0());
	};

	return run_on_param_switch<std::string>(paramID, code);
}

std::string MeshParamsBase::get_meshparam_value_sci(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object) -> std::string {

		return ToString(MatP_object.get0());
	};

	return run_on_param_switch<std::string>(paramID, code);
}

std::string MeshParamsBase::get_paraminfo_string(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object, PARAM_ paramID) -> std::string {

		return get_meshparam_handle(paramID) + ": " + MatP_object.get_info_string();
	};

	return run_on_param_switch<std::string>(paramID, code, paramID);
}

//returns a std::string describing the set spatial dependence with any parameters
std::string MeshParamsBase::get_paramvarinfo_string(PARAM_ paramID)
{
	auto code = [&](auto& MatP_object, PARAM_ paramID) -> std::string {

		return get_meshparam_handle(paramID) + ": " + MatP_object.get_varinfo_string();
	};

	return run_on_param_switch<std::string>(paramID, code, paramID);
}

bool MeshParamsBase::is_paramtemp_set(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return MatP_object.is_tdep();
	};

	return run_on_param_switch<bool>(paramID, code);
}

//check if the given parameters has a temperature dependence specified using a text equation (equation may still not be set due to missing constants)
bool MeshParamsBase::is_paramtempequation_set(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return MatP_object.is_t_equation_set();
	};

	return run_on_param_switch<bool>(paramID, code);
}

//check if the given parameter has a  spatial variation set
bool MeshParamsBase::is_paramvar_set(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return MatP_object.is_sdep();
	};

	return run_on_param_switch<bool>(paramID, code);
}

//check if the given parameter has a  temperature dependence or a spatial variation set
bool MeshParamsBase::is_param_nonconst(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return (MatP_object.is_tdep() || MatP_object.is_sdep());
	};

	return run_on_param_switch<bool>(paramID, code);
}

bool MeshParamsBase::get_meshparam_tempscaling(PARAM_ paramID, double max_temperature, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{
	auto code = [](auto& MatP_object, double max_temperature, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) -> bool {

		return MatP_object.get_temperature_scaling(max_temperature, x, y, z);
	};

	return run_on_param_switch<bool>(paramID, code, max_temperature, x, y, z);
}

//check if scaling array is scalar (or else vectorial)
bool MeshParamsBase::is_paramvar_scalar(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		typedef typename contained_type<decltype(MatP_object.s_scaling_ref())>::type SType;

		return std::is_same<SType, double>::value;
	};

	return run_on_param_switch<bool>(paramID, code);
}

//check if the spatial dependence is set using a text equation
bool MeshParamsBase::is_paramvarequation_set(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> bool {

		return MatP_object.is_s_equation_set();
	};

	return run_on_param_switch<bool>(paramID, code);
}

//-------------------------Spatial scaling VEC get / calculate methods

//get reference to mesh parameter spatial scaling VEC
void* MeshParamsBase::get_meshparam_s_scaling(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> void* {

		return &(MatP_object.s_scaling_ref());
	};

	return run_on_param_switch<void*>(paramID, code);
}

//get value of mesh parameter spatial scaling coefficient at given position (and time)
Any MeshParamsBase::get_meshparam_s_scaling_value(PARAM_ paramID, DBL3 rel_pos, double stime)
{
	auto code = [](auto& MatP_object, DBL3 rel_pos, double stime) -> Any {

		return (Any)MatP_object.get_s_scaling_value(rel_pos, stime);
	};

	return run_on_param_switch<Any>(paramID, code, rel_pos, stime);
}

//calculate spatial variation into the provided VECs - intended to be used when the spatial variation is set using a text equation
void MeshParamsBase::calculate_meshparam_s_scaling(PARAM_ paramID, VEC<double>& displayVEC_SCA, double stime)
{
	auto code = [](auto& MatP_object, VEC<double>& displayVEC_SCA, double stime) -> void {

		MatP_object.calculate_s_scaling(displayVEC_SCA, stime);
	};

	return run_on_param_switch<void>(paramID, code, displayVEC_SCA, stime);
}

void MeshParamsBase::calculate_meshparam_s_scaling(PARAM_ paramID, VEC<DBL3>& displayVEC_VEC, double stime)
{
	auto code = [](auto& MatP_object, VEC<DBL3>& displayVEC_VEC, double stime) -> void {

		MatP_object.calculate_s_scaling(displayVEC_VEC, stime);
	};

	return run_on_param_switch<void>(paramID, code, displayVEC_VEC, stime);
}

//-------------------------Setters : value and temperature dependence

//set value from std::string for named parameter (units allowed in std::string)
void MeshParamsBase::set_meshparam_value(PARAM_ paramID, std::string value_text)
{
	auto code = [](auto& MatP_object, std::string value_text, std::string unit) -> void {

		decltype(MatP_object.get0()) value = ToNum(value_text, unit);
		MatP_object = value;
		//C++17:
		//MatP_object = (decltype(MatP_object.get0()))ToNum(value_text, unit);
	};

	run_on_param_switch<void>(paramID, code, value_text, meshParams(paramID).unit);

	update_parameters(paramID);
}

//clear mesh parameter temperature dependence
void MeshParamsBase::clear_meshparam_temp(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> void {

		MatP_object.clear_t_scaling();
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param_switch<void>((PARAM_)meshParams.get_ID_from_index(index), code);
		}
	}
	else run_on_param_switch<void>(paramID, code);
}

//set mesh parameter array scaling
bool MeshParamsBase::set_meshparam_tscaling_array(PARAM_ paramID, std::vector<double>& temp, std::vector<double>& scaling_x, std::vector<double>& scaling_y, std::vector<double>& scaling_z)
{
	auto code = [](auto& MatP_object, std::vector<double>& temp, std::vector<double>& scaling_x, std::vector<double>& scaling_y, std::vector<double>& scaling_z) -> bool {

		return MatP_object.set_t_scaling_array(temp, scaling_x, scaling_y, scaling_z);
	};

	return run_on_param_switch<bool>(paramID, code, temp, scaling_x, scaling_y, scaling_z);
}

//set temperature dependence info std::string for console display purposes
void MeshParamsBase::set_meshparam_tscaling_info(PARAM_ paramID, std::string info_text)
{
	auto code = [](auto& MatP_object, std::string info_text) -> void {

		MatP_object.set_t_scaling_info(info_text);
	};

	run_on_param_switch<void>(paramID, code, info_text);
}

//-------------------------Setters : spatial variation

//set the mesh parameter spatial variation equation with given user constants
void MeshParamsBase::set_meshparam_s_equation(PARAM_ paramID, std::string& equationText, vector_key<double>& userConstants, DBL3 meshDimensions)
{
	auto code = [](auto& MatP_object, std::string& equationText, vector_key<double>& userConstants, DBL3 meshDimensions) -> void {

		MatP_object.set_s_scaling_equation(equationText, userConstants, meshDimensions);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param_switch<void>((PARAM_)meshParams.get_ID_from_index(index), code, equationText, userConstants, meshDimensions);
		}
	}
	else run_on_param_switch<void>(paramID, code, equationText, userConstants, meshDimensions);
}

//clear mesh parameter spatial variation (all if paramID == PARAM_ALL)
void MeshParamsBase::clear_meshparam_variation(PARAM_ paramID)
{
	auto code = [](auto& MatP_object) -> void {

		MatP_object.clear_s_scaling();
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param_switch<void>((PARAM_)meshParams.get_ID_from_index(index), code);
		}
	}
	else run_on_param_switch<void>(paramID, code);
}

//update mesh parameter spatial variation (e.g. cellsize or rectangle could have changed)
bool MeshParamsBase::update_meshparam_var(PARAM_ paramID, DBL3 h, Rect rect)
{
	auto code = [](auto& MatP_object, DBL3 h, Rect rect) -> bool {

		return MatP_object.update_s_scaling(h, rect);
	};

	return run_on_param_switch<bool>(paramID, code, h, rect);
}

//set parameter spatial variation using a given generator and arguments (arguments passed as a std::string to be interpreted and converted using ToNum)
BError MeshParamsBase::set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, DBL3 h, Rect rect, std::string generatorArgs, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	auto code = [](auto& MatP_object, DBL3 h, Rect rect, MATPVAR_ generatorID, std::string generatorArgs, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader) -> BError {

		return MatP_object.set_s_scaling(h, rect, generatorID, generatorArgs, bitmap_loader);
	};

	error = run_on_param_switch<BError>(paramID, code, h, rect, generatorID, generatorArgs, bitmap_loader);

	return error;
}

//set parameter spatial variation using a shape : set value in given shape only
BError MeshParamsBase::set_meshparam_shape(PARAM_ paramID, DBL3 h, Rect rect, std::vector<MeshShape> shapes, std::string value_text)
{
	BError error(__FUNCTION__);

	auto code = [](auto& MatP_object, DBL3 h, Rect rect, std::vector<MeshShape> shapes, std::string value_text) -> BError {

		return MatP_object.set_s_scaling_shape(h, rect, shapes, ToNum(value_text));
	};

	error = run_on_param_switch<BError>(paramID, code, h, rect, shapes, value_text);

	return error;
}

//-------------------------General Updaters

void MeshParamsBase::update_parameters(PARAM_ paramID)
{
	auto update_param = [&](PARAM_ update_paramID) {

		auto code = [&](auto& MatP_object) -> void {

			MatP_object.update(base_temperature);
		};

		return run_on_param_switch<void>(update_paramID, code);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			update_param((PARAM_)meshParams.get_ID_from_index(index));
		}
	}
	else update_param(paramID);
}