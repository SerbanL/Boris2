#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

//----------------------------------- PARAMETERS TEMPERATURE

//set mesh base temperature
void Mesh::SetBaseTemperature(double Temperature, bool clear_equation)
{
	if (clear_equation) T_equation.clear();

	//new base temperature and adjust parameter output values for new base temperature
	base_temperature = Temperature;

#if COMPILECUDA == 1
	if (pMeshCUDA) pMeshCUDA->base_temperature.from_cpu(base_temperature);
#endif

	//update any text equations used in mesh parameters (base temperature used as a constant - Tb)
	//need to update equations before parameters : if a parameter has a temperature dependence set using an equation by Temp.size() is zero, then updating parameters sets values based on the mesh base temperature
	update_meshparam_equations();

	//update parameter current values : if they have a temperature dependence set the base temperature will change their values
	update_parameters();

	//1a. reset Temp VEC to base temperature
	CallModuleMethod(&Heat::SetBaseTemperature, Temperature);

	//1b. Zeeman module might set field using a custom user equation where the base temperature is a parameter
	CallModuleMethod(&Zeeman::SetBaseTemperature, Temperature);

	//NOTE : do not call UpdateConfiguration here - this is to allow Temperature sequences in simulation stages
	//Instead deal with any adjustments required on a module by module basis (e.g. StrayField)

	//2. If this Mesh is the base for a DipoleMesh then need to force recalculateStrayField to be set to true, in case StrayField super-mesh module is set : dipole magnetisation value could have changed now.
	if (dynamic_cast<DipoleMesh*>(this)) dynamic_cast<DipoleMesh*>(this)->Reset_Mdipole();

	//3. electrical conductivity might also need updating so force it here - if Transport module not set then nothing happens (note elC will have zero size in this case)
	CallModuleMethod(&Transport::CalculateElectricalConductivity, true);
}

//set text equation for base temperature : when iterating, the base temperature will be evaluated and set using the text equation
BError Mesh::SetBaseTemperatureEquation(string equation_string, int step)
{
	BError error(CLASS_STR(Mesh));

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

void Mesh::UpdateTEquationUserConstants(void)
{
	if (userConstants.size()) {

		vector<pair<string, double>> constants(userConstants.size());
		for (int idx = 0; idx < userConstants.size(); idx++) {

			constants[idx] = { userConstants.get_key_from_index(idx), userConstants[idx] };
		}

		T_equation.set_constants(constants);
	}
}

//----------------------------------- PARAMETERS SPATIAL VARIATION

//update all mesh parameters spatial variation (if needed)
bool Mesh::update_meshparam_var(void)
{
	bool success = true;

	//update all mesh parameters in case cellsize or mesh rectangle has changed
	for (int index = 0; index < get_num_meshparams(); index++) {

		PARAM_ paramID = (PARAM_)get_meshparam_id(index);

		DBL3 cellsize = get_paramtype_cellsize(paramID);

		success &= MeshParams::update_meshparam_var(paramID, cellsize, meshRect);
	}

	return success;
}

//update text equations for mesh parameters with user constants, mesh dimensions, Curie temperature, base temperature
bool Mesh::update_meshparam_equations(void)
{
	bool success = true;

	for (int index = 0; index < get_num_meshparams(); index++) {

		PARAM_ paramID = (PARAM_)get_meshparam_id(index);

		success &= MeshParams::update_meshparam_equations(paramID, userConstants, meshRect.size());
	}

	return success;
}

//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
BError Mesh::set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	DBL3 cellsize = get_paramtype_cellsize(paramID);

	error = MeshParams::set_meshparam_var(paramID, generatorID, cellsize, meshRect, generatorArgs, bitmap_loader);

	return error;
}

//----------------------------------- OTHERS

//copy all parameters from another Mesh
BError Mesh::copy_mesh_parameters(Mesh& copy_this)
{
	BError error(__FUNCTION__);

	if (meshType != copy_this.GetMeshType()) return error(BERROR_INCORRECTVALUE);
	else copy_parameters(copy_this);

	//make sure to update the spatial variation as any copied spatial variation may not have the correct cellsize now
	update_meshparam_var();

	return error;
}

//each parameter has a certain type (PARAMTYPE_) - return cellsize in this mesh associated with this parameter type (e.g. ferromagnetic, electric, thermal, mechanical)
DBL3 Mesh::get_paramtype_cellsize(PARAM_ paramID)
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