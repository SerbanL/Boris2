#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- MESH PARAMETERS

//these set parameter values and temperature dependence in the indicated mesh - call through these since it's important to call UpdateConfiguration also
BError SuperMesh::set_meshparam_value(string meshName, string paramHandle, string value_text)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	pMesh[meshName]->set_meshparam_value(paramID, value_text);

	error = UpdateConfiguration();

	return error;
}

//get named parameter value from given mesh. Set value as a string in value_text, without units
BError SuperMesh::get_meshparam_value(string meshName, string paramHandle, string& value_text)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	value_text = pMesh[meshName]->get_meshparam_value_sci(paramID);

	return error;
}

//--------------------------------------------------------- temperature dependence

BError SuperMesh::set_meshparam_formula(string meshName, string paramHandle, string formulaName, vector<double> coefficients)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle) || !formula_descriptor.has_key(formulaName)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);
	MATPFORM_ formulaID = (MATPFORM_)formula_descriptor.get_ID_from_key(formulaName);

	pMesh[meshName]->set_meshparam_formula(paramID, formulaID, coefficients);

	error = UpdateConfiguration();

	return error;
}

BError SuperMesh::set_meshparam_tscaling_array(string meshName, string paramHandle, vector<double>& temp, vector<double>& scaling)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	if (pMesh[meshName]->set_meshparam_tscaling_array(paramID, temp, scaling)) {

		error = UpdateConfiguration();

		return error;
	}
	else return error(BERROR_INCORRECTARRAYS);
}

//clear parameters temperature dependence in given mesh (all meshes if empty string)
BError SuperMesh::clear_meshparam_temp(string meshName)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->set_meshparam_formula(PARAM_ALL, MATPFORM_NONE, {});
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->set_meshparam_formula(PARAM_ALL, MATPFORM_NONE, {});
		}
	}

	error = UpdateConfiguration();

	return error;
}

//--------------------------------------------------------- spatial dependence

//clear parameters spatial dependence (variation) in given mesh (all meshes if empty string)
BError SuperMesh::clear_meshparam_variation(string meshName)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->clear_meshparam_variation(PARAM_ALL);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->clear_meshparam_variation(PARAM_ALL);
		}
	}

	error = UpdateConfiguration();

	return error;
}

//clear parameter spatial dependence (variation) in given mesh for named parameter only
BError SuperMesh::clear_meshparam_variation(string meshName, string paramHandle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	pMesh[meshName]->clear_meshparam_variation(paramID);

	error = UpdateConfiguration();

	return error;
}

//set parameter to display in given mesh when ParamVar spatial variation display is enabled
BError SuperMesh::set_meshparamvar_display(string meshName, string paramHandle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetDisplayedParamVar(pMesh[meshName]->get_meshparam_id(paramHandle));

	return error;
}

//set parameter spatial variation using a given generator
BError SuperMesh::set_meshparam_var(string meshName, string paramHandle, string generatorHandle, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle) || !vargenerator_descriptor.has_key(generatorHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);
	MATPVAR_ generatorID = (MATPVAR_)vargenerator_descriptor.get_ID_from_key(generatorHandle);

	if (!generatorArgs.length()) generatorArgs = vargenerator_descriptor[generatorHandle];

	error = pMesh[meshName]->set_meshparam_var(paramID, generatorID, generatorArgs, bitmap_loader);

	if (!error) error = UpdateConfiguration();

	return error;
}

//--------------------------------------------------------- others

BError SuperMesh::SetBaseTemperature(string meshName, double Temperature)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetBaseTemperature(Temperature);

			//also set ambient temperature (for heat equation Robin boundary conditions) if heat module set
			pMesh[idx]->CallModuleMethod(&Heat::SetAmbientTemperature, Temperature);
		}
	}
	else {

		//named mesh only
		pMesh[meshName]->SetBaseTemperature(Temperature);

		//also set ambient temperature (for heat equation Robin boundary conditions) if heat module set
		pMesh[meshName]->CallModuleMethod(&Heat::SetAmbientTemperature, Temperature);
	}

	return error;
}

//ambient and alpha boundary coefficient for Robin boundary conditions - set in Heat module if active
BError SuperMesh::SetAmbientTemperature(string meshName, double T_ambient)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->CallModuleMethod(&Heat::SetAmbientTemperature, T_ambient);
		}
	}
	else {

		//named mesh only
		pMesh[meshName]->CallModuleMethod(&Heat::SetAmbientTemperature, T_ambient);
	}

	return error;
}

BError SuperMesh::SetAlphaHeatBoundary(string meshName, double alpha_boundary)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->CallModuleMethod(&Heat::SetAlphaBoundary, alpha_boundary);
		}
	}
	else {

		//named mesh only
		pMesh[meshName]->CallModuleMethod(&Heat::SetAlphaBoundary, alpha_boundary);
	}

	return error;
}

BError SuperMesh::SetInsulatingSides(string meshName, string literal, bool status)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->CallModuleMethod(&Heat::SetInsulatingSides, literal, status);

	return error;
}

BError SuperMesh::SetCurieTemperature(string meshName, double T_Curie)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetCurieTemperature(T_Curie);
		}
	}
	else {

		//named mesh only
		pMesh[meshName]->SetCurieTemperature(T_Curie);
	}

	return error;
}

BError SuperMesh::SetCurieTemperatureMaterial(string meshName, double T_Curie_material)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetCurieTemperatureMaterial(T_Curie_material);

	return error;
}

BError SuperMesh::SetAtomicMagneticMoment(string meshName, double atomic_moment)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetAtomicMoment(atomic_moment);
		}
	}
	else {

		//named mesh only
		pMesh[meshName]->SetAtomicMoment(atomic_moment);
	}

	return error;
}

//copy all parameters from another Mesh
BError SuperMesh::copy_mesh_parameters(string meshName_from, string meshName_to)
{
	BError error(__FUNCTION__);

	if (!contains(meshName_from) || !contains(meshName_to)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName_to]->copy_mesh_parameters(*pMesh[meshName_from]);

	return error;
}