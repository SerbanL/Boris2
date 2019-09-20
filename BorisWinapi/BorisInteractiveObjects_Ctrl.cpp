#include "stdafx.h"

//Interactive Objects : 
//
//handle user interactions using ConsoleActionHandler
//update their state depending on current program state using ConsoleInteractiveObjectState

//These are not needed in non-graphical mode
#include "CompileFlags.h"
#if GRAPHICS == 1

#include "Simulation.h"

InteractiveObjectActionOutcome Simulation::ConsoleActionHandler(int actionCode, InteractiveObjectProperties iop, TextObject *pTO) {

	//!!!IMPORTANT!!! Do not access BorisDisplay through thread-safe entry points from here. This method was called from within BorisDisplay (through a function pointer), which was thread-safe accessed so the mutex is now locked.
	//In fact it's better not to make any calls to BorisDisplay here at all (e.g. DisplayConsoleMessage). There should always be a better way. e.g. an interactive object will typically require a console
	//command to be issued - in this case just launch the command on the THREAD_HANDLEMESSAGE thread. If that command needs any console text displayed, it will wait for the display mutex to become available.

	//Makes changes in the program depending on the interactive object properties and the action type with which it was interacted with. No changes to the object are made here, that is handled by ConsoleInteractiveObjectState during Refresh cycles

	InteractiveObjectActionOutcome actionOutcome = AO_NOTHING;

	//action codes which can be handled in the same way for all objects
	if (actionCode == AC_HOVERCHECK) {

		//mouse is hovering over interactive object - display hover info?
		
		//the index in ioInfo to obtain info string
		INT2 id = INT2(iop.majorId, iop.minorId);
		
		//1. special cases first (require modification of minor id)
		switch (iop.majorId) {

		case IOI_SETSTAGEVALUE:
		{
			id.minor = simStages[INT2(0, iop.minorId)].stage_type();
		}
		break;

		case IOI_STAGESTOPCONDITION:
		{
			id.minor = simStages[INT2(0, iop.minorId)].stop_condition();
		}
		break;

		case IOI_DSAVETYPE:
		{
			id.minor = iop.auxId;
		}
		break;
		}

		//2. general method applicable to all interactive objects
		if (ioInfo.is_id_set(id)) {

			//full id entry set - get it
			actionOutcome.text = ioInfo[id];
		}
		else {

			//at most one entry for this majorId, so this could be a generic info string applicable for any minorId - if no entry at all then nothing to display
			if (ioInfo.is_ID_set(iop.majorId)) {

				actionOutcome.text = ioInfo(iop.majorId);
			}
		}

		//only ask to display if we managed to retrieve a non-empty string
		if (actionOutcome.text.length()) {

			actionOutcome = AO_SHOWHOVERINFO;
		}

		return actionOutcome;
	}

	//------------------------------------------------------ DEFINE LAMBDA CLOSURES FOR REUSABLE CODE USED ONLY IN THIS METHOD

	//INTERACTING OBJECTS : make change in their rerpresenting data structures (generic lambda closure needs C++14 minimum)
	auto interaction = [](auto &dataStructure, int io_id, int this_id) {

		//make sure this is not just the object interacting with itself - that would be pointless!
		if (io_id != this_id) {

			// move the data in list referred to by the interacting object ...
			int srcIdx = dataStructure.get_index_from_id(INT2(0, io_id));
			//... to this object's data list entry
			int dstIdx = dataStructure.get_index_from_id(INT2(0, this_id));
			//make the actual move
			dataStructure.move(srcIdx, dstIdx);
		}
	};

	//SEND COMMAND TO COMMAND HANDLER (not verbose) : makes code neater using this quick lambda
	auto sendCommand = [&](CMD_ commandCode, auto ...params) {

		string command = "~" + commands.get_key_from_index(commandCode);

		for (string param : { ToString(params)... })
			command += " " + param;

		//launch command
		single_call_launch<string>(&Simulation::HandleCommand, command, THREAD_HANDLEMESSAGE);
	};

	//SEND COMMAND TO COMMAND HANDLER (verbose, but erase console line entry) : makes code neater using this quick lambda
	auto sendCommand_verbose = [&](CMD_ commandCode, auto ...params) {

		string command = commands.get_key_from_index(commandCode);

		for (string param : { ToString(params)... })
			command += " " + param;

		//launch command
		single_call_launch<string>(&Simulation::HandleCommand, command, THREAD_HANDLEMESSAGE);
		single_call_launch<string>(&Simulation::SetConsoleEntryLineText, "", THREAD_HANDLEMESSAGE2);
	};

	//SET CONSOLE ENTRY LINE : form a command syntax and set it as the console entry line for further editing by user
	auto setConsoleEntry = [&](CMD_ commandCode, auto ...params) {

		string text = commands.get_key_from_index(commandCode);

		for (string param : { ToString(params)... })
			text += " " + param;

		single_call_launch<string>(&Simulation::SetConsoleEntryLineText, text, THREAD_HANDLEMESSAGE2);
	};

	//------------------------------------------------------ SWITCH FOR HANDLING THE DIFFERENT INTERACTIVE OBJECTS

	//take different action depending on the major interactive object identifier (this is a value from IOI_ enum)
	switch (iop.majorId) {

		//Shows program version update status : auxIdis the status as -1: attempting to connect, 0: connection failure, 1: program up to date, 2: update available
	case IOI_PROGRAMUPDATESTATUS:
	{
		if (actionCode == AC_MOUSELEFTDOWN) {

			single_call_launch(&Simulation::OpenDownloadPage, THREAD_HANDLEMESSAGE);
		}
	}
	break;

		//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). textId is the mesh name (if associated with this data type)
		//Note this entry must always represent the entry in Simulation::dataBoxList with the index in auxId.
	case IOI_DATABOXFIELDLABEL:
	{
		//parameters from iop
		int dataBoxList_idminor = iop.minorId;
		int DataBox_index = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN) {

			sendCommand(CMD_DELPINNEDDATA, DataBox_index);

			//need to update values as well as labels (values must correspond to labels)
			single_call_launch(&Simulation::UpdateDataBox_Refresh, THREAD_HANDLEMESSAGE);
		}

		//user might be trying to re-arrange data order : start interaction by creating a moveable pop-up window which holds this data entry
		else if (actionCode == AC_MOUSELEFTDOWN) { actionOutcome = AO_STARTINTERACTION; }

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_DATABOXFIELDLABEL
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == IOI_DATABOXFIELDLABEL) {

			interaction(dataBoxList, iop.interactingObjectId.minor, dataBoxList_idminor);

			//need to update values as well as labels (values must correspond to labels)
			single_call_launch(&Simulation::UpdateDataBox_Refresh, THREAD_HANDLEMESSAGE);
		}
	}
	break;

	//A set or available module for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to
	case IOI_MODULE:
	{
		//parameters from iop
		MOD_ module = (MOD_)iop.minorId;
		int meshId = iop.auxId;

		//try to add module
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_ADDMODULE, SMesh.key_from_meshId(meshId), moduleHandles(module));

		//try to remove module
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELMODULE, SMesh.key_from_meshId(meshId), moduleHandles(module));
	}
	break;

	//super-mesh module : minorId is an entry from MOD_ enum
	case IOI_SMODULE:
	{
		//parameters from iop
		MOD_ module = (MOD_)iop.minorId;

		//try to add module
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_ADDMODULE, SMesh.superMeshHandle, moduleHandles(module));

		//try to remove module
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELMODULE, SMesh.superMeshHandle, moduleHandles(module));
	}
	break;

	//Available/set ode : minorId is an entry from ODE_ (the equation)
	case IOI_ODE:
	{
		//parameters from iop
		ODE_ odeID = (ODE_)iop.minorId;

		//try to set ODE and its default evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_SETODE, odeHandles(odeID), odeEvalHandles(odeDefaultEval(odeID)));
	}
	break;

	//Set ODE time step: textId is the value
	case IOI_ODEDT:
	{
		//parameters from iop
		string dT_string = ToNum(iop.textId);

		//try to set ODE time step
		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SETDT, dT_string);
	}
	break;

	//Set heat equation time step: textId is the value
	case IOI_HEATDT:
	{
		//parameters from iop
		string dT_string = ToNum(iop.textId);

		//try to set ODE time step
		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SETHEATDT, dT_string);
	}
	break;

	//Available/set evaluation method for ode : minorId is an entry from ODE_ (the equation), auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	case IOI_ODE_EVAL:
	{
		//parameters from iop
		ODE_ odeID = (ODE_)iop.minorId;
		string evalHandle = iop.textId;

		//try to set ODE and its evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_SETODE, odeHandles(odeID), evalHandle);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMS:
	case IOI_MESH_FORPARAMSTEMP:
	case IOI_MESH_FORPARAMSVAR:
	case IOI_MESH_FORMODULES:
	case IOI_MESH_FORMESHLIST:
	case IOI_MESH_FORDISPLAYOPTIONS:
	case IOI_MESH_FORTEMPERATURE:
	case IOI_MESH_FORHEATBOUNDARIES:
	case IOI_MESH_FORCURIEANDMOMENT:
	case IOI_MESH_FORPBC:
	case IOI_MESH_FOREXCHCOUPLING:
	{
		//parameters from iop
		string meshName = iop.textId;

		//try to change mesh focus
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_MESHFOCUS, meshName);

		//try to delete mesh
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELMESH, meshName);

		//rename mesh : bring up console command
		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_RENAMEMESH, meshName);
	}
	break;

	//Shows mesh rectangle (units m) : minorId is the unique mesh id number, textId is the mesh rectangle
	case IOI_MESHRECTANGLE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		string meshRect_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) { 
			
			sendCommand(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			setConsoleEntry(CMD_MESHRECT, combine(split(meshRect_string, ", ", "; "), " "));
		}
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		string cellsize_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			setConsoleEntry(CMD_CELLSIZE, combine(split(cellsize_string, ", ", "; "), " "));
		}
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHECELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		string cellsize_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			setConsoleEntry(CMD_ECELLSIZE, combine(split(cellsize_string, ", ", "; "), " "));
		}
	}
	break;

	//Shows thermal mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHTCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		string cellsize_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			setConsoleEntry(CMD_TCELLSIZE, combine(split(cellsize_string, ", ", "; "), " "));
		}
	}
	break;

	//Shows ferromagnetic super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_FMSMESHCELLSIZE:
	{
		string cellsizeValue = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_FMSCELLSIZE, combine(split(cellsizeValue, ", ", "; "), " "));
	}
	break;

	//Shows electric super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_ESMESHCELLSIZE:
	{
		string cellsizeValue = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_ESCELLSIZE, combine(split(cellsizeValue, ", ", "; "), " "));
	}
	break;

	//Simulation output data, specifically used for showing values in console : minorId is the DATA_ id, textId is the data handle
	case IOI_SHOWDATA:
	{
		//parameters from iop
		string dataHandle = iop.textId;

		if (actionCode == AC_DOUBLECLICK) sendCommand_verbose(CMD_SHOWDATA, dataHandle);

		else if (actionCode == AC_MOUSELEFTDOWN) {

			setConsoleEntry(CMD_SHOWDATA, dataHandle);

			//this object can be dragged into the data box to display there
			actionOutcome = AO_STARTINTERACTION;
		}

		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_ADDPINNEDDATA, dataHandle);

		else if (actionCode == AC_INTERACTOBJECTWITHWINDOW) {

			if (iop.interactingObjectId == INT2(WIN_DATABOX, 0)) {

				sendCommand(CMD_ADDPINNEDDATA, dataHandle);

				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}
	}
	break;

	//Simulation output data, specifically used to construct output data list : minorId is the DATA_ id, textId is the data handle
	case IOI_DATA:
	{
		//parameters from iop
		string dataHandle = iop.textId;

		if (actionCode == AC_DOUBLECLICK) sendCommand(CMD_ADDDATA, dataHandle);

		else if (actionCode == AC_MOUSELEFTDOWN) actionOutcome = AO_STARTINTERACTION;

		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_ADDPINNEDDATA, dataHandle);

		else if (actionCode == AC_INTERACTOBJECTWITHWINDOW) {

			if (iop.interactingObjectId == INT2(WIN_DATABOX, 0)) {

				sendCommand(CMD_ADDPINNEDDATA, dataHandle);

				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}
	}
	break;

	//Show currently set directory : textId is the directory
	case IOI_DIRECTORY:
	{
		//parameters from iop
		string directory = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_CHDIR, directory);
	}
	break;

	//Show currently set save data file : textId is the file name
	case IOI_SAVEDATAFILE:
	{
		//parameters from iop
		string savedataFile = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SAVEDATAFILE, savedataFile);
	}
	break;

	//Show currently set image file base : textId is the file name
	case IOI_SAVEIMAGEFILEBASE:
	{
		//parameters from iop
		string imageSaveFileBase = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SAVEIMAGEFILE, imageSaveFileBase);
	}
	break;

	//Show flag status for data/image saving during a simulation : minorId is the flag value (boolean)
	case IOI_SAVEDATAFLAG:
	{
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_DATASAVEFLAG, !saveDataFlag);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DATASAVEFLAG, false);
	}
	break;

	case IOI_SAVEIMAGEFLAG:
	{
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_IMAGESAVEFLAG, !saveImageFlag);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_IMAGESAVEFLAG, false);
	}
	break;

	//Show set output data : minorId is the minor id of elements in Simulation::saveDataList (major id there is always 0), auxId is the number of the interactive object in the list as it appears in the console, textId is the configured output data.
	//Note this entry must always represent the entry in Simulation::saveDataList with the index in auxId.
	case IOI_OUTDATA:
	{
		//parameters from iop
		int outDataId = iop.minorId;
		int io_index = iop.auxId;

		string dataHandle = dataDescriptor.get_key_from_ID(saveDataList[INT2(0, outDataId)].datumId);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELDATA, io_index);

		else if (actionCode == AC_DOUBLECLICK) {

			string extraOptions;
			if (saveDataList[io_index].meshName.length()) extraOptions += saveDataList[io_index].meshName;
			if (!saveDataList[io_index].rectangle.IsNull()) extraOptions += " " + trim(ToString(saveDataList[io_index].rectangle, "m"), ",", ";");

			setConsoleEntry(CMD_EDITDATA, io_index, dataHandle, extraOptions);
		}

		else if (actionCode == AC_MOUSELEFTDOWN) actionOutcome = AO_STARTINTERACTION;

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_OUTDATA
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == IOI_OUTDATA) interaction(saveDataList, iop.interactingObjectId.minor, outDataId);

		else if (actionCode == AC_INTERACTOBJECTWITHWINDOW) {

			if (iop.interactingObjectId == INT2(WIN_DATABOX, 0)) {

				sendCommand(CMD_ADDPINNEDDATA, dataHandle, saveDataList[io_index].meshName);

				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}
	}
	break;

	//Shows a possible stage type, used for adding generic stages to the simulation schedule : minorId is the stage type (SS_ enum value, which is the majorId from stageDescriptors), textId is the stage setting handle
	case IOI_STAGE:
	{
		//parameters from iop
		string stageHandle = iop.textId;

		if (actionCode == AC_DOUBLECLICK) sendCommand(CMD_ADDSTAGE, stageHandle);
	}
	break;

	//Shows a stage added to the simulation schedule : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the configured stage text
	//Note this entry must always represent the entry in Simulation::simStages with the index in auxId.
	case IOI_SETSTAGE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		int io_index = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELSTAGE, io_index);

		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITSTAGE, io_index, stageDescriptors.get_key_from_ID(simStages[io_index].stage_type()), simStages[io_index].meshname());

		else if (actionCode == AC_MOUSELEFTDOWN) actionOutcome = AO_STARTINTERACTION;

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_SETSTAGE
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == IOI_SETSTAGE) interaction(simStages, iop.interactingObjectId.minor, stageId_minor);
	}
	break;

	//Shows the value to set for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the value as a string
	case IOI_SETSTAGEVALUE:
	{
		//parameters from iop
		int io_index = iop.auxId;
		string stageValueText = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITSTAGEVALUE, io_index, trim(stageValueText, ",", ";"));
	}
	break;

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a string
	case IOI_STAGESTOPCONDITION:
	{
		//parameters from iop
		int io_index = iop.auxId;
		string stopConditionText = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITSTAGESTOP, io_index, trim(stopConditionText, ":", ",", ";"));
	}
	break;

	//Shows the saving condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the DSAVE_ value for this data save type, textId is the save type and value as a string
	case IOI_DSAVETYPE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		string saveConditionText = iop.textId;

		int io_index = simStages.get_index_from_id(INT2(0, stageId_minor));

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_EDITDATASAVE, io_index, trim(saveConditionText, ":", ",", ";"));

		else if (actionCode == AC_DOUBLECLICK) {

			if (simStages[INT2(0, stageId_minor)].IsdSaveValueSet()) setConsoleEntry(CMD_EDITDATASAVE, io_index, trim(saveConditionText, ":", ",", ";"));
		}
	}
	break;

	//Shows a stop condition, used to apply the same condition to all simulation stages : minorId is the STOP_ value, textId is the stop type handle
	case IOI_STAGESTOPCONDITIONALL:
	{
		//parameters from iop
		string stopTypeHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_EDITSTAGESTOP, -1, stopTypeHandle);

		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITSTAGESTOP, -1, stopTypeHandle);
	}
	break;

	//Shows a data save condition, used to apply the same condition to all simulation stages : minorId is the DSAVE_ value, textId is the save type handle
	case IOI_DSAVETYPEALL:
	{
		//parameters from iop
		string saveTypeHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_EDITDATASAVE, -1, saveTypeHandle);

		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITDATASAVE, -1, saveTypeHandle);
	}
	break;

	//Shows parameter and value for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter handle and value
	case IOI_MESHPARAM:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		string paramText = iop.textId;

		string meshName = SMesh.key_from_meshId(meshId);

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SETPARAM, meshName, trim(paramText, ":", ",", ";"));
	}
	break;

	//Shows parameter temperature dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter temperature dependence setting
	case IOI_MESHPARAMTEMP:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		string paramText = iop.textId;

		string meshName = SMesh.key_from_meshId(meshId);

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SETPARAMTEMP, meshName, trim(paramText, ":", ",", ";", "(", ")"));
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_SETPARAMTEMP, meshName, SMesh[meshName]->get_meshparam_handle(paramId), "none");
		else if (actionCode == AC_DROPINTERACTOBJECTS) {

			if (iop.interactingObjectId.major == IOI_MESHPARAMTEMPFORMULA) {
				
				//when IOI_MESHPARAMTEMPFORMULA called for AO_STARTINTERACTION then iop.interactingObjectId.minor became the minorId of IOI_MESHPARAMTEMPFORMULA, i.e. the MATPFORM_ enum value
				sendCommand(CMD_SETPARAMTEMP, meshName, SMesh[meshName]->get_meshparam_handle(paramId), formula_descriptor.get_key_from_ID(iop.interactingObjectId.minor));
			}
		}
	}
	break;

	//Shows parameter spatial dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter spatial dependence setting
	case IOI_MESHPARAMVAR:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		string paramText = iop.textId;

		string meshName = SMesh.key_from_meshId(meshId);
		string paramName = SMesh[meshName]->get_meshparam_handle(paramId);

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_SETDISPLAYEDPARAMSVAR, meshName, paramName);
		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SETPARAMVAR, meshName, trim(paramText, ":", ";", ","));
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_CLEARPARAMVAR, meshName, paramName);
		else if (actionCode == AC_DROPINTERACTOBJECTS) {

			if (iop.interactingObjectId.major == IOI_MESHPARAMVARGENERATOR) {

				//when IOI_MESHPARAMVARGENERATOR called for AO_STARTINTERACTION then iop.interactingObjectId.minor became the minorId of IOI_MESHPARAMVARGENERATOR, i.e. the MATPVAR_ enum value
				sendCommand(CMD_SETPARAMVAR, meshName, paramName, trim(vargenerator_descriptor.get_key_from_ID(iop.interactingObjectId.minor), ";", ","));
			}
		}
	}
	break;

	//Shows a possible formula name for mesh parameter temperature dependence : minorId is the MATPFORM_ enum value, textId is the formula name
	case IOI_MESHPARAMTEMPFORMULA:
	{
		//parameters from iop
		MATPFORM_ formulaID = (MATPFORM_)iop.minorId;
		string formulaHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) { actionOutcome = AO_STARTINTERACTION; }
	}
	break;

	//Shows a possible generator name for mesh parameter spatial dependence : minorId is the MATPVAR_ enum value, textId is the generator name
	case IOI_MESHPARAMVARGENERATOR:
	{
		if (actionCode == AC_MOUSELEFTDOWN) { actionOutcome = AO_STARTINTERACTION; }
	}
	break;

	//Shows mesh display option for a given mesh : minorId is the MESHDISPLAY_ value, auxId is the unique mesh id number, textId is the MESHDISPLAY_ handle
	case IOI_MESHDISPLAY:
	{
		//parameters from iop
		MESHDISPLAY_ displayOption = (MESHDISPLAY_)iop.minorId;
		int meshId = iop.auxId;
		string displayHandle = iop.textId;

		string meshName = SMesh.key_from_meshId(meshId);

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_DISPLAY, displayHandle, meshName);
	}
	break;

	//Shows super-mesh display option : minorId is the MESHDISPLAY_ value, textId is the MESHDISPLAY_ handle
	case IOI_SMESHDISPLAY:
	{
		//parameters from iop
		MESHDISPLAY_ displayOption = (MESHDISPLAY_)iop.minorId;
		string displayHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_DISPLAY, displayHandle, SMesh.superMeshHandle);
	}
	break;

	//Shows mesh vectorial quantity display option : minorId is the unique mesh id number, auxId is the display option
	case IOI_MESHVECREP:
	{
		//parameters from iop
		int meshId = iop.minorId;
		VEC3REP_ option = (VEC3REP_)iop.auxId;

		string meshName = SMesh.key_from_meshId(meshId);

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_VECREP, meshName, (option + 1) % VEC3REP_NUMOPTIONS);
	}
	break;

	//Shows supermesh vectorial quantity display option : auxId is the display option
	case IOI_SMESHVECREP:
	{
		//parameters from iop
		VEC3REP_ option = (VEC3REP_)iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_VECREP, SMesh.superMeshHandle, (option + 1) % VEC3REP_NUMOPTIONS);
	}
	break;

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	case IOI_MOVINGMESH:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool moving_mesh = iop.auxId;
		string meshName = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_MOVINGMESH, meshName);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_MOVINGMESH, false);
	}
	break;

	//Shows movingmesh symmetry : auxId is the asymmetry status (1: asymmetric, 0: symmetric)
	case IOI_MOVINGMESHASYM:
	{
		//parameters from iop
		bool asymmetric = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN || actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_MOVINGMESHASYM, !asymmetric);
	}
	break;

	//Shows movingmesh threshold : textId is the threshold value as a string
	case IOI_MOVINGMESHTHRESH:
	{
		//parameters from iop
		string threshold_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_MOVINGMESHTHRESH, threshold_string);
	}
	break;

	//Shows electrode box. minorId is the minor Id in STransport::electrode_boxes, auxId is the number of the interactive object in the list (electrode index), textId is the electrode rect as a string
	case IOI_ELECTRODERECT:
	{
		//parameters from iop
		int electrodeId_minor = iop.minorId;
		int io_index = iop.auxId;

		if (actionCode == AC_DOUBLECLICK) {

			Rect electrode_rect = SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, io_index).first;

			setConsoleEntry(CMD_SETELECTRODERECT, io_index, combine(split(ToString(electrode_rect, "m"), ", ", "; "), " "));
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELELECTRODE, io_index);
	}
	break;

	//Shows electrode potential. minorId is the electrode index, textId is potential value as a string
	case IOI_ELECTRODEPOTENTIAL:
	{
		//parameters from iop
		int el_index = iop.minorId;
		string potential_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_SETELECTRODEPOTENTIAL, el_index, potential_string);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_DELELECTRODE, el_index);
	}
	break;

	//Shows electrode ground setting. minorId is the electrode index, auxId is the setting (0 : not ground, 1 : ground)
	case IOI_ELECTRODEGROUND:
	{
		//parameters from iop
		int el_index = iop.minorId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_DESIGNATEGROUND, el_index);
	}
	break;

	//Shows constant current source setting. auxId is the setting.
	case IOI_CONSTANTCURRENTSOURCE:
	{
		//parameters from iop
		bool is_constant_current = (bool)iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) {

			if (is_constant_current) sendCommand(CMD_SETPOTENTIAL, 0);
			else sendCommand(CMD_SETCURRENT, 0);
		}
	}
	break;

	//Shows transport solver convergence error. textId is the convergence error value.
	case IOI_TSOLVERCONVERROR:
	{
		//parameters from iop
		double conv_error = ToNum(iop.textId);

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_TSOLVERCONFIG, conv_error, SMesh.CallModuleMethod(&STransport::GetConvergenceTimeout));
	}
	break;

	//Shows transport solver timeout iterations. auxId is the timeout value.
	case IOI_TSOLVERTIMEOUT:
	{
		//parameters from iop
		int timeout = iop.auxId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_TSOLVERCONFIG, SMesh.CallModuleMethod(&STransport::GetConvergenceError), timeout);
	}
	break;

	//Shows spin transport solver convergence error. textId is the convergence error value.
	case IOI_SSOLVERCONVERROR:
	{
		//parameters from iop
		double conv_error = ToNum(iop.textId);

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SSOLVERCONFIG, conv_error, SMesh.CallModuleMethod(&STransport::GetSConvergenceTimeout));
	}
	break;

	//Shows spin transport solver timeout iterations. auxId is the timeout value.
	case IOI_SSOLVERTIMEOUT:
	{
		//parameters from iop
		int timeout = iop.auxId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SSOLVERCONFIG, SMesh.CallModuleMethod(&STransport::GetSConvergenceError), timeout);
	}
	break;

	//Shows Poisson solver SOR damping type : true for adaptive, false for fixed. auxId is enabled (1)/disabled(0) status.
	case IOI_SORFIXEDDAMPING:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_SETFIXEDSOR, !status);
	}
	break;

	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a string. (DBL2 since we need different damping values for V and S solvers)
	case IOI_SORDAMPING:
	{
		//parameters from iop
		string SOR_damping = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_SETSORDAMPING, combine(split(SOR_damping, ", ", "; "), " "));
	}
	break;

	//Static transport solver state. auxId is the value (0/1)
	case IOI_STATICTRANSPORT:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_STATICTRANSPORTSOLVER, !status);
	}
	break;

	//Shows mesh temperature. minorId is the unique mesh id number, textId is the temperature value
	case IOI_BASETEMPERATURE:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_TEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_TEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows ambient temperature for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the temperature value
	case IOI_AMBIENT_TEMPERATURE:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;
		bool enabled = iop.auxId;

		if (enabled) {

			if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_AMBIENTTEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_AMBIENTTEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows alpha value (W/m^2K) for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the value
	case IOI_ROBIN_ALPHA:
	{
		int meshId = iop.minorId;
		string alpha_string = iop.textId;
		bool enabled = iop.auxId;

		if (enabled) {

			if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_ROBINALPHA, alpha_string, SMesh.key_from_meshId(meshId));
			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_ROBINALPHA, alpha_string, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows temperature insulating side setting for heat equation. minorId is the unique mesh id number, auxId is the status (Heat module must be active) : -1 disabled (gray), 0 not insulating (green), 1 insulating (red), textId represents the side : "x", "-x", "y", "-y", "z", "-z"
	case IOI_INSULATINGSIDE:
	{
		int meshId = iop.minorId;
		string literal = iop.textId;
		int status = iop.auxId;

		if (status >= 0) {

			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_INSULATINGSIDES, literal, !(bool)status, SMesh.key_from_meshId(meshId));
			if (actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_INSULATINGSIDES, literal, !(bool)status, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows mesh Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	case IOI_CURIETEMP:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_CURIETEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_CURIETEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows indicative material Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	case IOI_CURIETEMPMATERIAL:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_CURIETEMPERATUREMATERIAL, temp_string, SMesh.key_from_meshId(meshId));
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_CURIETEMPERATUREMATERIAL, temp_string, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows atomic moment multiple of Bohr magneton. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the value
	case IOI_ATOMICMOMENT:
	{
		int meshId = iop.minorId;
		string amoment_string = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_ATOMICMOMENT, amoment_string, SMesh.key_from_meshId(meshId));
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_ATOMICMOMENT, amoment_string, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows cuda enabled/disabled or n/a state. auxId is enabled (1)/disabled(0)/not available(-1) status.
	case IOI_CUDASTATE:
	{
		int status = iop.auxId;
		if (status >= 0) {

			if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_CUDA, !status);
		}
	}
	break;
	
	//Shows scale_rects enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCALERECTSSTATUS:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_SCALEMESHRECTS, !status);
	}
	break;

	//Shows coupled_to_dipoles enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_COUPLEDTODIPOLESSTATUS:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_COUPLETODIPOLES, !status);
	}
	break;

	//Shows neighboring meshes exchange coupling setting for this mesh. minorId is the unique mesh id number, auxId is the status (1/0 : on/off, -1 : not available: must be ferromagnetic mesh)
	case IOI_MESHEXCHCOUPLING:
	{
		bool status = iop.auxId;
		int meshId = iop.minorId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_EXCHANGECOUPLEDMESHES, !status, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, auxId is enabled (1)/disabled(0) status. textId is the value
	case IOI_REFINEROUGHNESS:
	{
		int meshId = iop.minorId;
		string refine = iop.textId;
		bool status = iop.auxId;

		if (actionCode == AC_DOUBLECLICK && status) {

			setConsoleEntry(CMD_REFINEROUGHNESS, trim(refine, ","), SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows status of multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_MULTICONV:
	{
		int status = iop.auxId;

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand(CMD_MULTICONV, !status);
	}
	break;

	//Shows status of force 2D multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_2DMULTICONV:
	{
		int status = iop.auxId;

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand(CMD_2DMULTICONV, (status + 1) % 3);
	}
	break;

	//Shows status of use default n for multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_NCOMMONSTATUS:
	{
		int status = iop.auxId;

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand(CMD_NCOMMONSTATUS, !status);
	}
	break;

	//Shows n_common for multi-layered convolution. auxId is the status (-1 : N/A, otherwise available). textId is the value as a SZ3.
	case IOI_NCOMMON:
	{
		int status = iop.auxId;
		string n_common = iop.textId;

		if (status >= 0 && actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_NCOMMON, trim(n_common, ","));
		}
	}
	break;

	//Shows materials database in use. textId is the name of the database, including the path.
	case IOI_LOCALMDB:
	{
		//parameters from iop
		string mdbFile = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_MATERIALSDATABASE, mdbFile);
	}
	break;

	//Shows relative error fail threshold for ode eval. textId is the value.
	case IOI_ODERELERRFAIL:
	{
		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_ASTEPCTRL, trim(ToString(SMesh.Get_AStepRelErrCtrl()), ","), trim(ToString(SMesh.Get_AStepdTCtrl(), "s"), ","));
		}
	}
	break;

	//Shows relative error high threshold for decreasing dT. textId is the value.
	case IOI_ODERELERRHIGH:
	{
		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_ASTEPCTRL, trim(ToString(SMesh.Get_AStepRelErrCtrl()), ","), trim(ToString(SMesh.Get_AStepdTCtrl(), "s"), ","));
		}
	}
	break;

	//Shows relative error low threshold for increasing dT. textId is the value.
	case IOI_ODERELERRLOW:
	{
		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_ASTEPCTRL, trim(ToString(SMesh.Get_AStepRelErrCtrl()), ","), trim(ToString(SMesh.Get_AStepdTCtrl(), "s"), ","));
		}
	}
	break;

	//Shows dT increase factor. textId is the value.
	case IOI_ODEDTINCR:
	{
		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_ASTEPCTRL, trim(ToString(SMesh.Get_AStepRelErrCtrl()), ","), trim(ToString(SMesh.Get_AStepdTCtrl(), "s"), ","));
		}
	}
	break;

	//Shows minimum dT value. textId is the value.
	case IOI_ODEDTMIN:
	{
		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_ASTEPCTRL, trim(ToString(SMesh.Get_AStepRelErrCtrl()), ","), trim(ToString(SMesh.Get_AStepdTCtrl(), "s"), ","));
		}
	}
	break;

	//Shows maximum dT value. textId is the value.
	case IOI_ODEDTMAX:
	{
		if (actionCode == AC_DOUBLECLICK) {

			setConsoleEntry(CMD_ASTEPCTRL, trim(ToString(SMesh.Get_AStepRelErrCtrl()), ","), trim(ToString(SMesh.Get_AStepdTCtrl(), "s"), ","));
		}
	}
	break;

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc; -1 means setting is not available) (must be ferromagnetic mesh)
	case IOI_PBC_X:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand(CMD_PBC, SMesh.key_from_meshId(meshId), "x 10");
		else if (actionCode == AC_DOUBLECLICK && images > 0) setConsoleEntry(CMD_PBC, SMesh.key_from_meshId(meshId), "x", images);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_PBC, SMesh.key_from_meshId(meshId), "x 0");
	}
	break;

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc; -1 means setting is not available) (must be ferromagnetic mesh)
	case IOI_PBC_Y:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand(CMD_PBC, SMesh.key_from_meshId(meshId), "y 10");
		else if (actionCode == AC_DOUBLECLICK && images > 0) setConsoleEntry(CMD_PBC, SMesh.key_from_meshId(meshId), "y", images);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_PBC, SMesh.key_from_meshId(meshId), "y 0");
	}
	break;

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc; -1 means setting is not available) (must be ferromagnetic mesh)
	case IOI_PBC_Z:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand(CMD_PBC, SMesh.key_from_meshId(meshId), "z 10");
		else if (actionCode == AC_DOUBLECLICK && images > 0) setConsoleEntry(CMD_PBC, SMesh.key_from_meshId(meshId), "z", images);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_PBC, SMesh.key_from_meshId(meshId), "z 0");
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc; -1 means setting is not available)
	case IOI_SPBC_X:
	{
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand(CMD_PBC, SMesh.superMeshHandle, "x 10");
		else if (actionCode == AC_DOUBLECLICK && images > 0) setConsoleEntry(CMD_PBC, SMesh.superMeshHandle, "x", images);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_PBC, SMesh.superMeshHandle, "x 0");
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc; -1 means setting is not available)
	case IOI_SPBC_Y:
	{
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand(CMD_PBC, SMesh.superMeshHandle, "y 10");
		else if (actionCode == AC_DOUBLECLICK && images > 0) setConsoleEntry(CMD_PBC, SMesh.superMeshHandle, "y", images);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_PBC, SMesh.superMeshHandle, "y 0");
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc; -1 means setting is not available)
	case IOI_SPBC_Z:
	{
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand(CMD_PBC, SMesh.superMeshHandle, "z 10");
		else if (actionCode == AC_DOUBLECLICK && images > 0) setConsoleEntry(CMD_PBC, SMesh.superMeshHandle, "z", images);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand(CMD_PBC, SMesh.superMeshHandle, "z 0");
	}
	break;

	//Shows individual shape control flag. auxId is the value (0/1)
	case IOI_INDIVIDUALSHAPE:
	{
		bool status = (bool)iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand(CMD_INDIVIDUALMASKSHAPE, !status);
	}
	break;

	//Shows image cropping settings : textId has the DBL4 value as text
	case IOI_IMAGECROPPING:
	{
		string value_text = iop.textId;

		if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_IMAGECROPPING, trim(value_text, ","));
	}
	break;

	default:
		break;
	}

	return actionOutcome;
}

#endif