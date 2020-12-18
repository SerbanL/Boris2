#include "stdafx.h"
#include "MeshBase.h"

#include "SuperMesh.h"

//set text equation for base temperature : when iterating, the base temperature will be evaluated and set using the text equation
BError MeshBase::SetBaseTemperatureEquation(std::string equation_string, int step)
{
	BError error(CLASS_STR(Atom_Mesh));

	//set equation if not already set, or this is the first step in a stage
	if (!T_equation.is_set() || step == 0) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants();

		if (!T_equation.make_from_string(equation_string, { {"Ss", 0} })) return error(BERROR_INCORRECTSTRING);
	}
	else {

		//equation set and not the first step : adjust Ss constant
		T_equation.set_constant("Ss", step);
		UpdateTEquationUserConstants();
	}

	return error;
}

void MeshBase::UpdateTEquationUserConstants(void)
{
	if (userConstants.size()) {

		std::vector<std::pair<std::string, double>> constants(userConstants.size());
		for (int idx = 0; idx < userConstants.size(); idx++) {

			constants[idx] = { userConstants.get_key_from_index(idx), userConstants[idx] };
		}

		T_equation.set_constants(constants);
	}
}

//----------------------------------- PARAMETERS SPATIAL VARIATION

//update all mesh parameters spatial variation (if needed)
bool MeshBase::update_meshparam_var(void)
{
	bool success = true;

	//update all mesh parameters in case cellsize or mesh rectangle has changed
	for (int index = 0; index < get_num_meshparams(); index++) {

		PARAM_ paramID = (PARAM_)get_meshparam_id(index);

		DBL3 cellsize = get_paramtype_cellsize(paramID);

		success &= MeshParamsBase::update_meshparam_var(paramID, cellsize, meshRect);
	}

	return success;
}

//update text equations for mesh parameters with user constants, mesh dimensions, Curie temperature, base temperature
bool MeshBase::update_all_meshparam_equations(void)
{
	bool success = true;

	for (int index = 0; index < get_num_meshparams(); index++) {

		PARAM_ paramID = (PARAM_)get_meshparam_id(index);

		success &= update_meshparam_equations(paramID, userConstants, meshRect.size());
	}

	return success;
}

//set parameter spatial variation using a given generator and arguments (arguments passed as a std::string to be interpreted and converted using ToNum)
BError MeshBase::set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, std::string generatorArgs, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	DBL3 cellsize = get_paramtype_cellsize(paramID);

	error = MeshParamsBase::set_meshparam_var(paramID, generatorID, cellsize, meshRect, generatorArgs, bitmap_loader);

	return error;
}

//set parameter spatial variation using a shape : set value in given shape only
BError MeshBase::set_meshparam_shape(PARAM_ paramID, std::vector<MeshShape> shapes, std::string value_text)
{
	BError error(__FUNCTION__);

	DBL3 cellsize = get_paramtype_cellsize(paramID);

	error = MeshParamsBase::set_meshparam_shape(paramID, cellsize, meshRect, shapes, value_text);

	return error;
}

//----------------------------------- OTHERS

//each parameter has a certain type (PARAMTYPE_) - return cellsize in this mesh associated with this parameter type (e.g. ferromagnetic, electric, thermal, mechanical)
DBL3 MeshBase::get_paramtype_cellsize(PARAM_ paramID)
{
	PARAMTYPE_ pType = get_meshparam_type(paramID);

	switch (pType) {

	case PARAMTYPE_MAGNETIC:
		return h;

	case PARAMTYPE_ELECTRIC:
		return h_e;

	case PARAMTYPE_THERMAL:
		return h_t;

	case PARAMTYPE_MECHANICAL:
		return h_m;

	default:
		return h;
	}
}
