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

	error = UpdateConfiguration(UPDATECONFIG_PARAMVALUECHANGED);

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

BError SuperMesh::set_meshparam_t_equation(string meshName, string paramHandle, string equationText)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	pMesh[meshName]->set_meshparam_t_equation(paramID, equationText, userConstants);

	error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

	return error;
}

BError SuperMesh::set_meshparam_tscaling_array(string meshName, string paramHandle, vector<double>& temp, vector<double>& scaling_x, vector<double>& scaling_y, vector<double>& scaling_z, string fileName_info)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	if (pMesh[meshName]->set_meshparam_tscaling_array(paramID, temp, scaling_x, scaling_y, scaling_z)) {

		//set fileName as the temperature scaling info for console displayuu purposes
		ExtractFilenameDirectory(fileName_info);
		ExtractFilenameTermination(fileName_info);
		pMesh[meshName]->set_meshparam_tscaling_info(paramID, fileName_info);

		error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

		return error;
	}
	else return error(BERROR_INCORRECTARRAYS);
}

//clear parameters temperature dependence in given mesh for given parameter (all meshes and parameters if empty string, all parameters in given mesh if empty string)
BError SuperMesh::clear_meshparam_temp(string meshName, string paramHandle)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = PARAM_ALL;

	if (meshName.length() && paramHandle.length()) {

		if (!pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

		paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);
	}

	if (meshName.length()) {

		pMesh[meshName]->clear_meshparam_temp(paramID);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->clear_meshparam_temp(paramID);
		}
	}

	error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

	return error;
}

//--------------------------------------------------------- spatial dependence

BError SuperMesh::set_meshparam_s_equation(string meshName, string paramHandle, string equationText)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	pMesh[meshName]->set_meshparam_s_equation(paramID, equationText, userConstants, pMesh[meshName]->meshRect.size());

	error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

	return error;
}

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

	error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

	return error;
}

//clear parameter spatial dependence (variation) in given mesh for named parameter only
BError SuperMesh::clear_meshparam_variation(string meshName, string paramHandle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);

	pMesh[meshName]->clear_meshparam_variation(paramID);

	error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

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
BError SuperMesh::set_meshparam_var(string meshName, string paramHandle, string generatorHandle, string generatorArgs, function<vector<unsigned char>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->contains_param(paramHandle) || !vargenerator_descriptor.has_key(generatorHandle)) return error(BERROR_INCORRECTNAME);

	PARAM_ paramID = (PARAM_)pMesh[meshName]->get_meshparam_id(paramHandle);
	MATPVAR_ generatorID = (MATPVAR_)vargenerator_descriptor.get_ID_from_key(generatorHandle);

	if (!generatorArgs.length()) generatorArgs = vargenerator_descriptor[generatorHandle];

	error = pMesh[meshName]->set_meshparam_var(paramID, generatorID, generatorArgs, bitmap_loader);

	if (!error) error = UpdateConfiguration(UPDATECONFIG_PARAMCHANGED);

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