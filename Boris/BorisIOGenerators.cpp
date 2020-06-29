#include "stdafx.h"
#include "Simulation.h"
#include "BorisIO_Make.h"
#include "Boris.h"

//-------------------------------------Program update checker (in BorisIOGenerators.cpp)

//Check with "www.boris-spintronics.uk" if program version is up to date
void Simulation::CheckUpdate(void)
{
	string response_message;

	string console_update_message = "[tc0,0.5,0,1/tc]Program version update status : "  + MakeIO(IOI_PROGRAMUPDATESTATUS);

	console_update_message += "</c>[tc0,0.5,0,1/tc] Check for updates on startup : " + MakeIO(IOI_UPDATESTATUSCHECKSTARTUP, start_check_updates);

	BD.DisplayFormattedConsoleMessage(console_update_message);

	//trying to connect
	version_checking = -1;

	HTTP httpComms(domain_name);

	if (httpComms.is_available() && httpComms.http_post(version_checker, string("version=") + ToString(Program_Version), response_message)) {

		int latest_version = ToNum(response_message);

		if (latest_version == Program_Version) {

			//program up to date
			version_checking = 1;
		}
		else {

			//update available
			version_checking = 2;
		}

		if (httpComms.is_available() && httpComms.http_get(mdb_lastupdate, response_message)) {

			BD.DisplayConsoleMessage("Online materials database last updated on (PT) : " + response_message);
		}
	}
	else {

		//connection failure
		version_checking = 0;
	}

	RefreshScreen();
}

void Simulation::OpenDownloadPage(void)
{
	open_web_page(domain_name + download_page);
}

//---------------------------------------------------- MESH LIST

void Simulation::Print_Mesh_List(void) 
{
	string mesh_list = "[tc1,1,1,1/tc]Available meshes ([tc0,0.5,0,1/tc]green [tc1,1,1,1/tc]in focus) :\n";

	//super-mesh
	//1. magnetic
	mesh_list += "</c>Magnetic <b>" + SMesh.superMeshHandle + "</b> " + MakeIO(IOI_FMSMESHRECTANGLE) + "</c> with cell " + MakeIO(IOI_FMSMESHCELLSIZE);
	//2. electric
	mesh_list += "</c> Electric <b>" + SMesh.superMeshHandle + "</b> " + MakeIO(IOI_ESMESHRECTANGLE) + "</c> with cell " + MakeIO(IOI_ESMESHCELLSIZE);

	mesh_list += "\n";

	//normal meshes
	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		mesh_list += Build_Mesh_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(mesh_list);
}

string Simulation::Build_Mesh_ListLine(int meshIndex) 
{
	string mesh_line =
		"Mesh " + MakeIO(IOI_MESH_FORMESHLIST, meshIndex) +
		"</c> [sa0/sa]with rectangle " + MakeIO(IOI_MESHRECTANGLE, meshIndex) +
		"</c> [sa1/sa]Magnetic cell " + MakeIO(IOI_MESHCELLSIZE, meshIndex) +
		"</c> [sa2/sa]Electric cell " + MakeIO(IOI_MESHECELLSIZE, meshIndex) +
		"</c> [sa3/sa]Thermal cell " + MakeIO(IOI_MESHTCELLSIZE, meshIndex) +
		"</c> [sa4/sa]Mechanical cell " + MakeIO(IOI_MESHMCELLSIZE, meshIndex);

	return mesh_line;
}

//---------------------------------------------------- MESH DISPLAY

void Simulation::Print_MeshDisplay_List(void)
{
	string mesh_display_list = "[tc1,1,1,1/tc]Mesh display thresholds (minimum, maximum) : " + MakeIO(IOI_MESHDISPLAYTHRESHOLDS, ToString(displayThresholds)) + "</c>[tc1,1,1,1/tc] For vectors trigger on: " + MakeIO(IOI_MESHDISPLAYTHRESHOLDTRIG) + "</c>\n";
	mesh_display_list += "[tc1,1,1,1/tc]Dual mesh display transparency (foreground, background) : " + MakeIO(IOI_MESHDISPLAYTRANSPARENCY, ToString(displayTransparency)) + "</c>\n";
	
	mesh_display_list += "[tc1,1,1,1/tc]<b>" + SMesh.superMeshHandle + "</b> [a0/a]display options ";

	//super-mesh
	for (int idx = 0; idx < (int)meshAllowedDisplay(MESH_SUPERMESH).size(); idx++) {

		mesh_display_list += "[a" + ToString(idx + 1) + "/a]" + MakeIO(IOI_SMESHDISPLAY, meshAllowedDisplay(MESH_SUPERMESH)[idx]) + "</c> ";
	}

	mesh_display_list += " (Vector display: " + MakeIO(IOI_SMESHVECREP) + "</c>)\n";

	//normal meshes
	for (int meshIndex = 0; meshIndex < (int)SMesh().size(); meshIndex++) {

		mesh_display_list += Build_MeshDisplay_ListLine(meshIndex) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(mesh_display_list);
}

string Simulation::Build_MeshDisplay_ListLine(int meshIndex)
{
	MESH_ meshType = SMesh[meshIndex]->GetMeshType();

	string mesh_display_line = MakeIO(IOI_MESH_FORDISPLAYOPTIONS, meshIndex) + "</c> [sa0/sa]display options ";

	//iterate through all the possible display options for this mesh type
	for (int idx = 0; idx < (int)meshAllowedDisplay(meshType).size(); idx++) {

		mesh_display_line += "[sa" + ToString(idx + 1) + "/sa]" + MakeIO(IOI_MESHDISPLAY, meshIndex, meshAllowedDisplay(meshType)[idx]) + "</c> ";
	}

	mesh_display_line += " (Vector display: " + MakeIO(IOI_MESHVECREP, meshIndex) + "</c>)";

	return mesh_display_line;
}

//---------------------------------------------------- MODULES LIST

void Simulation::Print_Modules_List(void) 
{
	string modules_list = "[tc1,1,1,1/tc]Available modules ([tc0,0.5,0,1/tc]green [tc1,1,1,1/tc]added, [tc1,0,0,1/tc]red [tc1,1,1,1/tc]not added) :\n";

	modules_list += "[tc1,1,1,1/tc]<b>" + SMesh.superMeshHandle + "</b> [a0/a]modules ";

	//first show modules for super-mesh
	for (int idx = 0; idx < (int)modules_for_meshtype(MESH_SUPERMESH).size(); idx++) {

		//super-mesh module
		MOD_ module = modules_for_meshtype(MESH_SUPERMESH)[idx];

		//only display this super-mesh module if it's not a super-mesh companion module
		if (!superMeshCompanionModules[module].size()) {

			modules_list += MakeIO(IOI_SMODULE, modules_for_meshtype(MESH_SUPERMESH)[idx]) + "</c> ";
		}
	}

	modules_list += "</c>\n";

	//next show for normal modules
	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		modules_list += Build_Modules_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(modules_list);
}

string Simulation::Build_Modules_ListLine(int meshIndex)
{
	string modulesListLine = "[tc1,1,1,1/tc]" + MakeIO(IOI_MESH_FORMODULES, meshIndex) + "</c> [sa0/sa]modules ";

	//get the mesh type so we can index modules_for_meshtype with it (modules_for_meshtype is a vector_lut with MESH_ values as the major id - it contains the available modules for this mesh type)
	MESH_ meshType = SMesh[meshIndex]->GetMeshType();

	for (int idx = 0; idx < (int)modules_for_meshtype(meshType).size(); idx++) {

		//get next available module for this mesh type
		MOD_ moduleID = modules_for_meshtype(meshType)[idx];

		//some modules are not set in moduleHandles, as they are not intended to be displayed in the console ("silent" background modules, e.g. SDemag_Demag managed by through SDemag)
		if (moduleHandles.is_ID_set(moduleID)) {

			modulesListLine += "[sa" + ToString(idx + 1) + "/sa]" + MakeIO(IOI_MODULE, meshIndex, modules_for_meshtype(meshType)[idx]) + "</c> ";
		}
	}

	return modulesListLine;
}

//---------------------------------------------------- ODES LIST

void Simulation::Print_ODEs(void) 
{
	string odes_eval_list = "[tc1,1,1,1/tc]Available ODEs and associated evaluation methods ([tc0,0.5,0,1/tc]green [tc1,1,1,1/tc]set, [tc1,0,0,1/tc]red [tc1,1,1,1/tc]not set, [tc0.5,0.5,0.5,1/tc]gray [tc1,1,1,1/tc]not available) : \n";

	odes_eval_list += "[tc1,1,1,1/tc]Set Micromagnetic Equation: ";

	for (int odeIdx = 0; odeIdx < (int)odeHandles.size(); odeIdx++) {

		ODE_ odeId = (ODE_)odeHandles.get_ID_from_index(odeIdx);

		odes_eval_list += "[tc1,1,1,1/tc][sa" + ToString(odeIdx) + "/sa]" + MakeIO(IOI_ODE, odeId) + "</c> ";
	}

	odes_eval_list += "\n[tc1,1,1,1/tc]Set Atomistic Equation: ";

	for (int odeIdx = 0; odeIdx < (int)atom_odeHandles.size(); odeIdx++) {

		ODE_ odeId = (ODE_)atom_odeHandles.get_ID_from_index(odeIdx);

		odes_eval_list += "[tc1,1,1,1/tc][sa" + ToString(odeIdx) + "/sa]" + MakeIO(IOI_ATOMODE, odeId) + "</c> ";
	}

	odes_eval_list += "\n";
	odes_eval_list += "[tc1,1,1,1/tc]Evaluation:   ";

	for (int evalIdx = 0; evalIdx < (int)odeEvalHandles.size(); evalIdx++) {

		EVAL_ evalId = (EVAL_)odeEvalHandles.get_ID_from_index(evalIdx);

		odes_eval_list += "[tc1,1,1,1/tc]" + MakeIO(IOI_ODE_EVAL, evalId) + "</c> ";
	}

	odes_eval_list += "\n";
	odes_eval_list += "[tc1,1,1,1/tc]Time Step: " + MakeIO(IOI_ODEDT) + "</c>";

	BD.DisplayFormattedConsoleMessage(odes_eval_list);
}

//---------------------------------------------------- SHOW DATA LIST (for console display)

void Simulation::Print_ShowData(void) 
{
	string showData = "[tc1,1,1,1/tc]Choose data to show value for (double-click to show, right click to pin it to data box):\n";

	for (int idx = 0; idx < (int)dataDescriptor.size(); idx++) {

		DATA_ dataID = (DATA_)dataDescriptor.get_ID_from_index(idx);

		showData += MakeIO(IOI_SHOWDATA, dataID) + "</c> ";
	}

	BD.DisplayFormattedConsoleMessage(showData);
}

//---------------------------------------------------- OUTPUT DATA LIST (for saving during a simulation)

void Simulation::Print_AvailableOutputData(void) 
{
	string showData = "[tc1,1,1,1/tc]Add data to output list from the following :\n";

	for (int idx = 0; idx < (int)dataDescriptor.size(); idx++) {

		DATA_ dataID = (DATA_)dataDescriptor.get_ID_from_index(idx);

		showData += MakeIO(IOI_DATA, dataID) + "</c> ";
	}

	BD.DisplayFormattedConsoleMessage(showData);
}

void Simulation::Print_SetOutputData_List(void) 
{
	string showData =
		string("[tc1,1,1,1/tc]Current file for output data : ") + 
		"[a0/a]" + MakeIO(IOI_DIRECTORY, directory) + "</c>" +
		MakeIO(IOI_SAVEDATAFILE, savedataFile) +
		"</c>[tc1,1,1,1/tc]. [a1/a]Data save " +
		"[a2/a]" + MakeIO(IOI_SAVEDATAFLAG) + "</c>\n";

	showData +=
		string("</c>[tc1,1,1,1/tc]Current file base for images : ") +
		"[sa0/sa]" + MakeIO(IOI_SAVEIMAGEFILEBASE, imageSaveFileBase) + "</c>" +
		"</c>[tc1,1,1,1/tc]. [sa1/sa]Image save " +
		"[sa2/sa]" + MakeIO(IOI_SAVEIMAGEFLAG) + "</c>\n";

	showData += "</c>[tc1,1,1,1/tc]List of output data as: dataname (<meshname>, ((box)) :\n";

	for (int idx = 0; idx < (int)saveDataList.size(); idx++) {

		showData += Build_SetOutputData_ListLine(idx) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(showData);
}

string Simulation::Build_SetOutputData_ListLine(int index_in_list) 
{
	return MakeIO(IOI_OUTDATA, index_in_list);
}

string Simulation::Build_SetOutputData_Text(int index_in_list)
{
	string meshName, dataBoxString;

	if (!dataDescriptor(saveDataList[index_in_list].datumId).meshless) {

		meshName = " <" + saveDataList[index_in_list].meshName + ">";
	}

	if (!dataDescriptor(saveDataList[index_in_list].datumId).boxless) {

		dataBoxString = " (" + ToString(saveDataList[index_in_list].rectangle, "m") + ")";
	}

	return (ToString(index_in_list) + ". " + dataDescriptor.get_key_from_ID(saveDataList[index_in_list].datumId) + meshName + dataBoxString);
}

//---------------------------------------------------- SIMULATION STAGES

void Simulation::Print_SetStages_List(void) 
{
	//Current set stages
	string simStagesList = "[tc1,1,1,1/tc]Simulation schedule stages :\n";

	for (int idx = 0; idx < (int)simStages.size(); idx++) {

		simStagesList += Build_SetStages_ListLine(idx) + "\n";
	}
	
	//Available configuration options
	string showStageTypes = "</c>\n[tc1,1,1,1/tc]Stage/step stopping conditions, apply to all stages : ";
	
	for (int stopIdx = 0; stopIdx < stageStopDescriptors.size(); stopIdx++) {
		
		showStageTypes += "[a" + ToString(stopIdx) + "/a]" + MakeIO(IOI_STAGESTOPCONDITIONALL, stageStopDescriptors.get_ID_from_index(stopIdx)) + "</c> ";
	}

	showStageTypes += "\n[tc1,1,1,1/tc]Data saving conditions, apply to all stages : ";

	for (int dsaveIdx = 0; dsaveIdx < dataSaveDescriptors.size(); dsaveIdx++) {

		showStageTypes += "[sa" + ToString(dsaveIdx) + "/sa]" + MakeIO(IOI_DSAVETYPEALL, dataSaveDescriptors.get_ID_from_index(dsaveIdx)) + "</c> ";
	}

	showStageTypes += "\n\n[tc1,1,1,1/tc]Add a stage to simulation schedule from the following :\n";

	for (int idx = 0; idx < stageDescriptors.size(); idx++) {

		if (!stageDescriptors[idx].invisible) showStageTypes += MakeIO(IOI_STAGE, stageDescriptors.get_ID_from_index(idx)) + "</c> ";
	}
	
	BD.DisplayFormattedConsoleMessage(simStagesList + showStageTypes);
}

string Simulation::Build_SetStages_ListLine(int index_in_list) 
{
	string stageLineText;

	//the stage tpye
	stageLineText = MakeIO(IOI_SETSTAGE, index_in_list);

	//the stage value (if set)
	if (simStages[index_in_list].IsValueSet()) stageLineText += "</c> [sa0/sa]" + MakeIO(IOI_SETSTAGEVALUE, index_in_list);
	else stageLineText += "</c>[sa0/sa]";

	//next is the stopping condition
	stageLineText += "</c> [sa1/sa]Stop: " + MakeIO(IOI_STAGESTOPCONDITION, index_in_list);

	//the data save types - list all, only the active one will be set with the active color
	stageLineText += "</c> [sa2/sa]Save: ";

	for (int dsaveIdx = 0; dsaveIdx < dataSaveDescriptors.size(); dsaveIdx++) {

		stageLineText += MakeIO(IOI_DSAVETYPE, index_in_list, dataSaveDescriptors.get_ID_from_index(dsaveIdx)) + "</c> ";
	}

	return stageLineText;
}

string Simulation::Build_SetStages_Text(int index_in_list) 
{
	string meshName;
	if (simStages[index_in_list].meshname().length()) meshName = " <" + simStages[index_in_list].meshname() + ">";

	string text = ToString(index_in_list) + ". " + stageDescriptors.get_key_from_ID(simStages[index_in_list].stage_type()) + meshName;

	return text;
}

string Simulation::Build_SetStages_ValueText(int index_in_list) 
{
	string text;

	if (simStages[index_in_list].IsValueSet()) text = simStages[index_in_list].get_value_string();

	return text;
}

string Simulation::Build_SetStages_StopConditionText(int index_in_list) 
{
	//the stop condition and value
	string stopConditionString = stageStopDescriptors.get_key_from_ID(simStages[index_in_list].stop_condition());

	if (simStages[index_in_list].IsStopValueSet())
		stopConditionString += ": " + simStages[index_in_list].get_stopvalue_string();

	return stopConditionString;
}

string Simulation::Build_SetStages_SaveConditionText(int index_in_list, int dsaveIdx) 
{
	string dsaveConditionString = dataSaveDescriptors.get_key_from_index(dsaveIdx);

	if (simStages[index_in_list].IsdSaveValueSet())
		dsaveConditionString += ": " + simStages[index_in_list].get_dsavevalue_string();

	return dsaveConditionString;
}

//---------------------------------------------------- MESH PARAMETERS

void Simulation::Print_MeshParams(string meshName)
{
	if (!SMesh.contains(meshName)) return;

	BD.DisplayFormattedConsoleMessage(Build_MeshParams_Line(SMesh().index_from_key(meshName)));
}

string Simulation::Build_MeshParams_Line(int meshIndex)
{
	string meshParamsString = "[tc1,1,1,1/tc]Parameters for " + MakeIO(IOI_MESH_FORPARAMS, meshIndex) + "</c>\n";

	string ferromagneticParameters, econductionParameters, tconductionParameters, mechanicalParameters;

	int f_index = 0;
	int e_index = 0;
	int t_index = 0;
	int m_index = 0;

	for (int paramIdx = 0; paramIdx < SMesh[meshIndex]->get_num_meshparams(); paramIdx++) {

		PARAM_ paramId = (PARAM_)SMesh[meshIndex]->get_meshparam_id(paramIdx);
		PARAMTYPE_ paramType = (PARAMTYPE_)SMesh[meshIndex]->get_meshparam_type(paramId);

		if (!SMesh[meshIndex]->is_param_hidden(paramId)) {

			switch (paramType) {

			case PARAMTYPE_MAGNETIC:
				//ferromagneticParameters += "[sa" + ToString(f_index) + "/sa]" + MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				ferromagneticParameters += MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				f_index++;
				break;

			case PARAMTYPE_ELECTRIC:
				//econductionParameters += "[sa" + ToString(e_index) + "/sa]" + MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				econductionParameters += MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				e_index++;
				break;

			case PARAMTYPE_THERMAL:
				//tconductionParameters += "[sa" + ToString(t_index) + "/sa]" + MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				tconductionParameters += MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				t_index++;
				break;

			case PARAMTYPE_MECHANICAL:
				//mechanicalParameters += "[sa" + ToString(m_index) + "/sa]" + MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				mechanicalParameters += MakeIO(IOI_MESHPARAM, meshIndex, paramId) + "</c> ";
				m_index++;
				break;
			}
		}
	}

	meshParamsString +=
		"Magnetic " + ferromagneticParameters + "\n" +
		"Electric " + econductionParameters + "\n" +
		"Thermal " + tconductionParameters + "\n" +
		"Mechanical " + mechanicalParameters;

	return meshParamsString;
}

string Simulation::Build_MeshParams_Text(int meshIdx, PARAM_ paramId)
{
	return SMesh[meshIdx]->get_meshparam_handle(paramId) + ": " + SMesh[meshIdx]->get_meshparam_value(paramId);
}

//---------------------------------------------------- MESH PARAMETERS TEMPERATURE DEPENDENCE

//print all mesh parameters temperature dependence for the given mesh name
void Simulation::Print_MeshParamsTemperature(string meshName)
{
	if (!SMesh.contains(meshName)) return;

	string meshParamsTemp = Build_MeshParamsTemp_Text(SMesh().index_from_key(meshName));

	meshParamsTemp += "[tc1,1,1,1/tc]\nAvailable temperature dependence type (drag to param, right click for info): ";

	for (int idx = 0; idx < temperature_dependence_type.size(); idx++) {

		MATPTDEP_ type = (MATPTDEP_)temperature_dependence_type.get_ID_from_index(idx);

		meshParamsTemp += "</c> " + MakeIO(IOI_MESHPARAMTEMPTYPE, type);
	}

	BD.DisplayFormattedConsoleMessage(meshParamsTemp);
}

string Simulation::Build_MeshParamsTemp_Text(int meshIndex)
{
	string meshParamsString = "[tc1,1,1,1/tc]Parameters temperature dependence for " + MakeIO(IOI_MESH_FORPARAMSTEMP, meshIndex) + "</c>\n";

	string ferromagneticParameters, econductionParameters, tconductionParameters, mechanicalParameters;

	int f_index = 0;
	int e_index = 0;
	int t_index = 0;
	int m_index = 0;

	for (int paramIdx = 0; paramIdx < SMesh[meshIndex]->get_num_meshparams(); paramIdx++) {

		PARAM_ paramId = (PARAM_)SMesh[meshIndex]->get_meshparam_id(paramIdx);
		PARAMTYPE_ paramType = (PARAMTYPE_)SMesh[meshIndex]->get_meshparam_type(paramId);

		bool temperature_dependence_enabled = params_enabled_props(paramId).first;

		if (!SMesh[meshIndex]->is_param_hidden(paramId) && temperature_dependence_enabled) {

			switch (paramType) {

			case PARAMTYPE_MAGNETIC:
				//ferromagneticParameters += "[sa" + ToString(f_index) + "/sa]" + MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				ferromagneticParameters += MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				f_index++;
				break;

			case PARAMTYPE_ELECTRIC:
				//econductionParameters += "[sa" + ToString(e_index) + "/sa]" + MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				econductionParameters += MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				e_index++;
				break;

			case PARAMTYPE_THERMAL:
				//tconductionParameters += "[sa" + ToString(t_index) + "/sa]" + MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				tconductionParameters += MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				t_index++;
				break;

			case PARAMTYPE_MECHANICAL:
				//mechanicalParameters += "[sa" + ToString(m_index) + "/sa]" + MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				mechanicalParameters += MakeIO(IOI_MESHPARAMTEMP, meshIndex, paramId) + "</c> ";
				m_index++;
				break;
			}
		}
	}

	meshParamsString +=
		"Magnetic " + ferromagneticParameters + "\n" +
		"Electric " + econductionParameters + "\n" +
		"Thermal " + tconductionParameters + "\n" +
		"Mechanical " + mechanicalParameters;

	return meshParamsString;
}

//---------------------------------------------------- MESH PARAMETERS SPATIAL DEPENDENCE

//print all mesh parameters temperature dependence for the given mesh name
void Simulation::Print_MeshParamsVariation(string meshName)
{
	if (!SMesh.contains(meshName)) return;

	string meshParamsVar = Build_MeshParamsVariation_Text(SMesh().index_from_key(meshName));

	meshParamsVar += "[tc1,1,1,1/tc]\nAvailable spatial variation generators (drag to param, right click for info): ";

	for (int idx = 0; idx < vargenerator_descriptor.size(); idx++) {

		MATPVAR_ generatorID = (MATPVAR_)vargenerator_descriptor.get_ID_from_index(idx);

		meshParamsVar += "</c> " + MakeIO(IOI_MESHPARAMVARGENERATOR, generatorID);
	}

	BD.DisplayFormattedConsoleMessage(meshParamsVar);
}

string Simulation::Build_MeshParamsVariation_Text(int meshIndex)
{
	string meshParamsString = "[tc1,1,1,1/tc]Parameters spatial dependence for " + MakeIO(IOI_MESH_FORPARAMSVAR, meshIndex) + "</c>\n";

	string ferromagneticParameters, econductionParameters, tconductionParameters, mechanicalParameters;

	int f_index = 0;
	int e_index = 0;
	int t_index = 0;
	int m_index = 0;

	for (int paramIdx = 0; paramIdx < SMesh[meshIndex]->get_num_meshparams(); paramIdx++) {

		PARAM_ paramId = (PARAM_)SMesh[meshIndex]->get_meshparam_id(paramIdx);
		PARAMTYPE_ paramType = (PARAMTYPE_)SMesh[meshIndex]->get_meshparam_type(paramId);

		bool spatial_variation_enabled = params_enabled_props(paramId).second;

		if (!SMesh[meshIndex]->is_param_hidden(paramId) && spatial_variation_enabled) {

			switch (paramType) {

			case PARAMTYPE_MAGNETIC:
				//ferromagneticParameters += "[sa" + ToString(f_index) + "/sa]" + MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				ferromagneticParameters += MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				f_index++;
				break;

			case PARAMTYPE_ELECTRIC:
				//econductionParameters += "[sa" + ToString(e_index) + "/sa]" + MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				econductionParameters += MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				e_index++;
				break;

			case PARAMTYPE_THERMAL:
				//tconductionParameters += "[sa" + ToString(t_index) + "/sa]" + MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				tconductionParameters += MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				t_index++;
				break;

			case PARAMTYPE_MECHANICAL:
				//mechanicalParameters += "[sa" + ToString(m_index) + "/sa]" + MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				mechanicalParameters += MakeIO(IOI_MESHPARAMVAR, meshIndex, paramId) + "</c> ";
				m_index++;
				break;
			}
		}
	}

	meshParamsString +=
		"Magnetic " + ferromagneticParameters + "\n" +
		"Electric " + econductionParameters + "\n" +
		"Thermal " + tconductionParameters + "\n" +
		"Mechanical " + mechanicalParameters;

	return meshParamsString;
}

//---------------------------------------------------- MOVING MESH SETTINGS

void Simulation::PrintMovingMeshSettings(void)
{
	string movingmesh_text = "[tc1,1,1,1/tc]Moving mesh trigger : " + MakeIO(IOI_MOVINGMESH);

	movingmesh_text += "</c>[tc1,1,1,1/tc] Moving mesh type : " + MakeIO(IOI_MOVINGMESHASYM);
	movingmesh_text += "</c>[tc1,1,1,1/tc] Mesh shift threshold : " + MakeIO(IOI_MOVINGMESHTHRESH);

	BD.DisplayFormattedConsoleMessage(movingmesh_text);
}

//---------------------------------------------------- ELECTRODES and TRANSPORT SETTINGS

void Simulation::Print_Electrodes_List(void)
{

	string electrode_list = "[tc1,1,1,1/tc]Using " + MakeIO(IOI_CONSTANTCURRENTSOURCE) + "</c>\n";

	if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

		int num_electrodes = SMesh.CallModuleMethod(&STransport::GetNumberofElectrodes);

		//go through all electrodes in this mesh
		for (int el_index = 0; el_index < num_electrodes; el_index++) {

			electrode_list += Build_Electrodes_ListLine(el_index) + "\n";
		}
	}

	BD.DisplayFormattedConsoleMessage(electrode_list);
}

string Simulation::Build_Electrodes_ListLine(int el_index)
{
	//electrode rectangle
	string electrode_list_line = "[tc1,1,1,1/tc]Electrode rectangle " + MakeIO(IOI_ELECTRODERECT, el_index) + "</c>";

	//electrode potential
	electrode_list_line += "[tc1,1,1,1/tc] [sa0/sa]Potential " + MakeIO(IOI_ELECTRODEPOTENTIAL, el_index) + "</c>";

	//electrode ground setting
	electrode_list_line += "[tc1,1,1,1/tc] [sa1/sa]" + MakeIO(IOI_ELECTRODEGROUND, el_index) + "</c>";

	return electrode_list_line;
}

void Simulation::PrintTransportSolverConfig(void)
{
	string tsolver_text;

	tsolver_text = "[tc1,1,1,1/tc]Charge-solver convergence error : " + MakeIO(IOI_TSOLVERCONVERROR) + "</c>";
	tsolver_text += " with iterations timeout : " + MakeIO(IOI_TSOLVERTIMEOUT) + "</c>";
	tsolver_text += " Spin-solver convergence error : " + MakeIO(IOI_SSOLVERCONVERROR) + "</c>";
	tsolver_text += " with iterations timeout : " + MakeIO(IOI_SSOLVERTIMEOUT) + "</c>\n";
	tsolver_text += " SOR damping values (V, S) : " + MakeIO(IOI_SORDAMPING) + "</c>\n";
	tsolver_text += "Static transport solver : " + MakeIO(IOI_STATICTRANSPORT) + "</c>\n";

	BD.DisplayFormattedConsoleMessage(tsolver_text);
}

//---------------------------------------------------- TEMPERATURE

void Simulation::Print_MeshTemperature_List(void)
{
	string mesh_temperature_list;

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		mesh_temperature_list += Build_MeshTemperature_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(mesh_temperature_list);
}

string Simulation::Build_MeshTemperature_ListLine(int meshIndex)
{
	string mesh_temperature_line = MakeIO(IOI_MESH_FORTEMPERATURE, meshIndex) + "</c>[tc1,1,1,1/tc] [sa0/sa]base temperature : " + MakeIO(IOI_BASETEMPERATURE, meshIndex);

	return mesh_temperature_line;
}

void Simulation::Print_HeatBoundaries_List(void)
{
	string mesh_heatboundaries_list;

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		mesh_heatboundaries_list += Build_HeatBoundaries_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(mesh_heatboundaries_list);
}

string Simulation::Build_HeatBoundaries_ListLine(int meshIndex)
{

	string mesh_heatboundaries_line = 
		MakeIO(IOI_MESH_FORHEATBOUNDARIES, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa0/sa]Ambient temperature : " + MakeIO(IOI_AMBIENT_TEMPERATURE, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa1/sa]Robin coefficient : " + MakeIO(IOI_ROBIN_ALPHA, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa2/sa]Insulating sides : " +
		MakeIO(IOI_INSULATINGSIDE, meshIndex, "x") +
		"</c>[sa3/sa]" + MakeIO(IOI_INSULATINGSIDE, meshIndex, "-x") +
		"</c>[sa4/sa]" + MakeIO(IOI_INSULATINGSIDE, meshIndex, "y") +
		"</c>[sa5/sa]" + MakeIO(IOI_INSULATINGSIDE, meshIndex, "-y") +
		"</c>[sa6/sa]" + MakeIO(IOI_INSULATINGSIDE, meshIndex, "z") +
		"</c>[sa7/sa]" + MakeIO(IOI_INSULATINGSIDE, meshIndex, "-z");

	return mesh_heatboundaries_line;
}

//---------------------------------------------------- CURIE TEMPERATURE and ATOMIC MOMENT

void Simulation::Print_CurieandMoment_List(void)
{
	string curie_and_moment_list;

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		curie_and_moment_list += Build_CurieandMoment_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(curie_and_moment_list);
}

string Simulation::Build_CurieandMoment_ListLine(int meshIndex)
{
	if (meshIndex < SMesh.size()) {

		if (SMesh[meshIndex]->GetMeshType() != MESH_ANTIFERROMAGNETIC) {

			string curie_and_moment_line = MakeIO(IOI_MESH_FORCURIEANDMOMENT, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa0/sa]Set Curie temperature : [sa1/sa]" + MakeIO(IOI_CURIETEMP, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa2/sa]Indicative material Curie temperature : [sa3/sa]" + MakeIO(IOI_CURIETEMPMATERIAL, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa4/sa]Atomic moment (multiple of Bohr magneton) : [sa5/sa]" + MakeIO(IOI_ATOMICMOMENT, meshIndex);

			return curie_and_moment_line;
		}
		else {

			string curie_and_moment_line = MakeIO(IOI_MESH_FORCURIEANDMOMENT, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa0/sa]Set Neel temperature : [sa1/sa]" + MakeIO(IOI_CURIETEMP, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa2/sa]Indicative material Neel temperature : [sa3/sa]" + MakeIO(IOI_CURIETEMPMATERIAL, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa4/sa]Atomic moments (multiples of Bohr magneton) : [sa5/sa]" + MakeIO(IOI_ATOMICMOMENT_AFM, meshIndex) +
				"</c>[tc1,1,1,1/tc] [sa6/sa]Tc coupling tau11, tau22, tau12, tau21 : [sa7/sa]" + MakeIO(IOI_TAU, meshIndex);

			return curie_and_moment_line;
		}
	}

	return "";
}

//---------------------------------------------------- TEMPERATURE MODEL TYPE

void Simulation::Print_TemperatureModel_List(void)
{
	string tmodel_list;

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		tmodel_list += Build_TemperatureModel_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(tmodel_list);
}

string Simulation::Build_TemperatureModel_ListLine(int meshIndex)
{
	string tmodel_line = MakeIO(IOI_MESH_FORTMODEL, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa0/sa]Temperature model type : " + MakeIO(IOI_TMODEL, meshIndex);

	return tmodel_line;
}

//---------------------------------------------------- STOCHASTICITY SETIINGS

void Simulation::Print_Stochasticity_List(void)
{
	string stochasticity_list = "[tc1,1,1,1/tc]Stochastic time-step : " + MakeIO(IOI_STOCHDT) + "</c> Linked to ODE dT : " + MakeIO(IOI_LINKSTOCHDT) + "</c>\n";

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		stochasticity_list += Build_Stochasticity_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(stochasticity_list);
}

string Simulation::Build_Stochasticity_ListLine(int meshIndex)
{
	string stochasticity_line = MakeIO(IOI_MESH_FORSTOCHASTICITY, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa0/sa]Stochastic cellsize : " + MakeIO(IOI_MESHSCELLSIZE, meshIndex) + "</c>[tc1,1,1,1/tc] [sa1/sa] Linked to magnetic cellize : " + MakeIO(IOI_LINKSTOCHASTIC, meshIndex);

	return stochasticity_line;
}

//---------------------------------------------------- EVALUATION SPEEDUP SETIINGS

void Simulation::Print_Speedup_List(void)
{
	string speedup_list = "[tc1,1,1,1/tc]Evaluation speedup mode : " + MakeIO(IOI_SPEEDUPMODE) + "</c> Speedup time-step (for extreme mode only): " + MakeIO(IOI_SPEEDUPDT) + "</c> Linked to ODE dT : " + MakeIO(IOI_LINKSPEEDUPDT) + "</c>\n";

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		speedup_list += Build_Speedup_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(speedup_list);
}

string Simulation::Build_Speedup_ListLine(int meshIndex)
{
	string speedup_line = MakeIO(IOI_MESH_FORSPEEDUP, meshIndex) + "</c>[tc1,1,1,1/tc] Macrocell size : " + MakeIO(IOI_MESHDMCELLSIZE, meshIndex);

	return speedup_line;
}

//---------------------------------------------------- CUDA and MEMORY INFO

void Simulation::Print_CUDAStatus(void)
{
	string cuda_info = "[tc1,1,1,1/tc]CUDA status : " + MakeIO(IOI_CUDASTATE, cudaEnabled);

	BD.DisplayFormattedConsoleMessage(cuda_info);
}

void Simulation::Print_MemoryInfo(void)
{
	string memory_info;

	memory_info += "[tc1,1,1,1/tc]CPU free memory (MB) : " + MakeIO(IOI_CPUMEMFREE, cpuMemFree_MB);
	memory_info += "</c>[tc1,1,1,1/tc] [a0/a]out of (MB) : " + MakeIO(IOI_CPUMEMTOTAL, cpuMemTotal_MB) + "\n";

	memory_info += "</c>[tc1,1,1,1/tc]GPU free memory (MB) : " + MakeIO(IOI_GPUMEMFREE, gpuMemFree_MB);
	memory_info += "</c>[tc1,1,1,1/tc] [sa0/sa]out of (MB) : " + MakeIO(IOI_GPUMEMTOTAL, gpuMemTotal_MB) + "\n";

	BD.DisplayFormattedConsoleMessage(memory_info);
}

//---------------------------------------------------- SCALE RECTS STATUS

void Simulation::Print_Scale_Rects_Status(void)
{
	string scale_rects_info = "[tc1,1,1,1/tc]Scale mesh rectangles : " + MakeIO(IOI_SCALERECTSSTATUS, SMesh.Get_Scale_Rects());

	BD.DisplayFormattedConsoleMessage(scale_rects_info);
}

//---------------------------------------------------- COUPLED-To-DIPOLES STATUS

void Simulation::Print_CoupledToDipoles_Settings(void)
{
	string coupled_to_dipoles_info = "[tc1,1,1,1/tc]Exchange coupling to dipoles : " + MakeIO(IOI_COUPLEDTODIPOLESSTATUS, SMesh.Get_Coupled_To_Dipoles());

	BD.DisplayFormattedConsoleMessage(coupled_to_dipoles_info);
}

//---------------------------------------------------- ERROR LOG STATUS and other STARTUP OPTIONS

void Simulation::Print_ErrorLogStatus(void)
{
	string error_log_info = "[tc1,1,1,1/tc]Error log status : " + MakeIO(IOI_ERRORLOGSTATUS, log_errors);

	BD.DisplayFormattedConsoleMessage(error_log_info);
}

void Simulation::Print_StartupUpdateCheckStatus(void)
{
	string startup_info = "[tc1,1,1,1/tc]Check for updates on startup : " + MakeIO(IOI_UPDATESTATUSCHECKSTARTUP, start_check_updates);

	BD.DisplayFormattedConsoleMessage(startup_info);
}

void Simulation::Print_StartupScriptServerStatus(void)
{
	string startup_info = "[tc1,1,1,1/tc]Start script server on startup : " + MakeIO(IOI_SCRIPTSERVERSTARTUP, start_scriptserver);

	BD.DisplayFormattedConsoleMessage(startup_info);
}

//print all startup options
void Simulation::Print_StartupOptions(void)
{
	string startup_info = "[tc1,1,1,1/tc]Check for updates on startup : [sa0/sa]" + MakeIO(IOI_UPDATESTATUSCHECKSTARTUP, start_check_updates);
	startup_info += "</c>\n[tc1,1,1,1/tc]Start script server on startup : [sa0/sa]" + MakeIO(IOI_SCRIPTSERVERSTARTUP, start_scriptserver);
	startup_info += "</c>\n[tc1,1,1,1/tc]Error log status : [sa0/sa]" + MakeIO(IOI_ERRORLOGSTATUS, log_errors);

	BD.DisplayFormattedConsoleMessage(startup_info);
}

//---------------------------------------------------- NEIGHBORING MESHES EXCHANGE COUPLING

void Simulation::Print_ExchangeCoupledMeshes_List(void)
{
	string ec_info = "[tc1,1,1,1/tc]Neighboring meshes exchange coupling.\n";

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		ec_info += Build_ExchangeCoupledMeshes_ListLine(idxMesh);
	}

	BD.DisplayFormattedConsoleMessage(ec_info);
}

string Simulation::Build_ExchangeCoupledMeshes_ListLine(int meshIndex)
{
	string ec_line = MakeIO(IOI_MESH_FOREXCHCOUPLING, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa0/sa]: " + MakeIO(IOI_MESHEXCHCOUPLING, meshIndex) + "\n";

	return ec_line;
}

//---------------------------------------------------- MESH ROUGHNESS REFINEMENT

void Simulation::Print_MeshRoughnessRefinement(string meshName)
{
	string roughness_refinement_info = "[tc1,1,1,1/tc]Roughness refinement cells multiplier : " + MakeIO(IOI_REFINEROUGHNESS, meshName);

	BD.DisplayFormattedConsoleMessage(roughness_refinement_info);
}
//---------------------------------------------------- MULTILAYERED CONVOLUTION CONFIGURATION

void Simulation::Print_MultiConvolution_Config(void)
{
	string multiconv_info = "[tc1,1,1,1/tc]Multi-layered convolution : " + MakeIO(IOI_MULTICONV, SMesh.Get_Multilayered_Convolution_Status());

	multiconv_info += "</c> [tc1,1,1,1/tc]Force 2D : " + MakeIO(IOI_2DMULTICONV, SMesh.Get_2D_Multilayered_Convolution_Status());
	multiconv_info += "</c> [tc1,1,1,1/tc]Use default n : " + MakeIO(IOI_NCOMMONSTATUS, SMesh.Use_Default_n_Status());
	multiconv_info += "</c> [tc1,1,1,1/tc]Common n : " + MakeIO(IOI_NCOMMON, SMesh.Get_n_common());

	BD.DisplayFormattedConsoleMessage(multiconv_info);
}

//---------------------------------------------------- MATERIALS DATABASE

void Simulation::Print_MaterialsDatabase(void)
{
	string mdb_info = "[tc1,1,1,1/tc]Local materials database in use : " + MakeIO(IOI_LOCALMDB, mdb.GetDataBaseName());

	BD.DisplayFormattedConsoleMessage(mdb_info);
}

//---------------------------------------------------- ADAPTIVE TIME STEP CONTROL

void Simulation::Print_AStepCtrl(void)
{
	string astep_ctrl_info = "[tc1,1,1,1/tc]Adaptive time step control. err_fail : " + MakeIO(IOI_ODERELERRFAIL, SMesh.Get_AStepRelErrCtrl().i);
	astep_ctrl_info += "</c> [tc1,1,1,1/tc]err_high : " + MakeIO(IOI_ODERELERRHIGH, SMesh.Get_AStepRelErrCtrl().j);
	astep_ctrl_info += "</c> [tc1,1,1,1/tc]err_low : " + MakeIO(IOI_ODERELERRLOW, SMesh.Get_AStepRelErrCtrl().k);
	astep_ctrl_info += "</c> [tc1,1,1,1/tc]dT_incr : " + MakeIO(IOI_ODEDTINCR, SMesh.Get_AStepdTCtrl().i);
	astep_ctrl_info += "</c> [tc1,1,1,1/tc]dT_min : " + MakeIO(IOI_ODEDTMIN, SMesh.Get_AStepdTCtrl().j);
	astep_ctrl_info += "</c> [tc1,1,1,1/tc]dT_max : " + MakeIO(IOI_ODEDTMAX, SMesh.Get_AStepdTCtrl().k);

	BD.DisplayFormattedConsoleMessage(astep_ctrl_info);
}

//---------------------------------------------------- PERIODIC BOUNDARY CONDITIONS

void Simulation::Print_PBC(void)
{
	string pbc_info = "[tc1,1,1,1/tc]Periodic boundary conditions for magnetization.\n";

	pbc_info += string("[tc1,1,1,1/tc]Supermesh PBC: ") +
		"</c>[tc1,1,1,1/tc] [sa0/sa] x: " + MakeIO(IOI_SPBC_X) +
		"</c>[tc1,1,1,1/tc] [sa1/sa] y: " + MakeIO(IOI_SPBC_Y) +
		"</c>[tc1,1,1,1/tc] [sa2/sa] z: " + MakeIO(IOI_SPBC_Z) + "\n";

	for (int idxMesh = 0; idxMesh < (int)SMesh().size(); idxMesh++) {

		pbc_info += Build_PBC_ListLine(idxMesh) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(pbc_info);
}

string Simulation::Build_PBC_ListLine(int meshIndex)
{
	string pbc_line = MakeIO(IOI_MESH_FORPBC, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa0/sa] x: " + MakeIO(IOI_PBC_X, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa1/sa] y: " + MakeIO(IOI_PBC_Y, meshIndex) +
		"</c>[tc1,1,1,1/tc] [sa2/sa] z: " + MakeIO(IOI_PBC_Z, meshIndex);

	return pbc_line;
}

//---------------------------------------------------- INDIVIDUAL SHAPE CONTROL

void Simulation::Print_IndividualShapeStatus(void)
{
	string individualshape_info = "[tc1,1,1,1/tc]Individual shape status flag : " + MakeIO(IOI_INDIVIDUALSHAPE, shape_change_individual);

	BD.DisplayFormattedConsoleMessage(individualshape_info);
}

//---------------------------------------------------- USER EQUATION CONSTANTS

//Print currently set equation constants
void Simulation::Print_EquationConstants(void)
{
	string showUserConstants = "[tc1,1,1,1/tc]Defined user constants :\n";

	for (int idx = 0; idx < (int)userConstants.size(); idx++) {

		showUserConstants += Build_EquationConstants_ListLine(idx) + "\n";
	}

	BD.DisplayFormattedConsoleMessage(showUserConstants);
}

//build formatted string for interactive objects describing user constants at index_in_list from userConstants (also helper method to build the actual unformatted display text in the object)
string Simulation::Build_EquationConstants_ListLine(int index_in_list)
{
	string name = userConstants.get_key_from_index(index_in_list);
	double value = userConstants[index_in_list];

	return MakeIO(IOI_USERCONSTANT, name, value, index_in_list);
}

string Simulation::Build_EquationConstants_Text(int index_in_list)
{
	string name = userConstants.get_key_from_index(index_in_list);
	double value = userConstants[index_in_list];

	return name + ": " + ToString(value);
}
