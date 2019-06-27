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

	InteractiveObjectActionOutcome actionOutcome = AO_NONE;

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

	//Available/set ode and evaluation method for magnetization : minorId is an entry from ODE_ (the equation), auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	case IOI_ODE:
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

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand(CMD_2DMULTICONV, !status);
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

	default:
		break;
	}

	return actionOutcome;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

InteractiveObjectStateChange Simulation::ConsoleInteractiveObjectState(InteractiveObjectProperties &iop, TextObject *pTO) {

	//!!!IMPORTANT!!!: Do not call for a Refresh in this method, as it is called during a Refresh() : causes infinite loop! 
	//Also, this method was called from within BorisDisplay (through a function pointer), which was thread-safe accessed so the mutex is now locked.

	//return true if TextObject was changed in any way (including its state). Set iop.state = IOS_DELETING if this object needs to be deleted.

	InteractiveObjectStateChange stateChanged;

	//------------------------------------------------------ DEFINE LAMBDA CLOSURES FOR IMPLEMENTING REUSABLE CODE USED ONLY IN THIS METHOD

	//used for interactive object lists : index_in_list is the index of the interactive object, lastIndex is the last entry index in the data structure represented by the interactive objects, 
	//simMethod_BuildListEntry is a Simulation method which builds a formatted text line to represent a given data structure entry at a given index
	auto updateList = [&](int index_in_list, int lastIndex, auto simMethod_BuildListEntry) {

		//if this is the last in list, make sure it is marked by setting its state IOS_ISLASTINLIST (e.g. could happen last element which did have IOS_ISLASTINLIST state set, was deleted)
		if (index_in_list == lastIndex) {

			if (iop.state != IOS_ISLASTINLIST)
				iop.state = IOS_ISLASTINLIST;
		}
		else if (index_in_list > lastIndex) {

			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
		else {

			//if not last in list, but marked as being last in list then further elements must be inserted
			if (iop.state == IOS_ISLASTINLIST) {

				stateChanged = true;
				//set IOS_WASLASTINLIST so the caller knows to insert the object below
				iop.state = IOS_WASLASTINLIST;
				//insert a new output data interactive object after this : simMethod_BuildListEntry is a Simulation method which takes an integer argument (the index for which to build the formatted text string) and returns a string
				stateChanged.textMessage = "";
				for (int idx = index_in_list + 1; idx <= lastIndex; idx++) {

					//Note, we need to allow for the possibility of inserting more than one list line at a time (can happen if two or more elements are added before calling for a screen refresh)
					//Use new-line separators, and the caller checking for IOS_WASLASTINLIST will then split the text message using the newline separators two add 2 or more paragraphs at a time
					stateChanged.textMessage += CALLFP(this, simMethod_BuildListEntry)(idx) + "\n";
				}
			}
		}
	};

	//display a mesh interactive object which is tagged to a list line (e.g. list of meshes, modules for meshes, etc.) - does book-keeping like delete the line when the mesh is deleted, update list etc.
	//there are many objects of this type and they only differ in the method used to build the list line
	auto display_meshIO = [&](auto simMethod_Build_Mesh_ListLine) {

		//parameters from iop
		int meshId = iop.minorId;
		bool update = (bool)iop.auxId;
		string meshName = iop.textId;

		//from unique mesh id number get index in pMesh (held in SMesh)
		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0 && meshName == SMesh().get_key_from_index(meshIdx)) {

			//mark currently active mesh
			if (meshName == SMesh.GetMeshFocus()) pTO->SetBackgroundColor(ONCOLOR);
			else pTO->SetBackgroundColor(OFFCOLOR);
		}
		else {

			//mismatch found : either the mesh name has changed or entry has been deleted.
			if (meshIdx >= 0) {

				//mesh still exists, it's just the name that has changed
				meshName = SMesh().get_key_from_index(meshIdx);
				iop.textId = meshName;

				pTO->set(" " + meshName + " ");
				stateChanged = true;
			}
			else {

				//mesh no longer exists : delete the entire paragraph containing this object
				stateChanged = true;
				iop.state = IOS_DELETINGPARAGRAPH;
			}
		}

		//this object is part of a list : make sure this list is updated
		if (update && !stateChanged) updateList(SMesh().index_from_key(meshName), SMesh.size() - 1, simMethod_Build_Mesh_ListLine);
	};

	//------------------------------------------------------ SWITCH FOR HANDLING THE DIFFERENT INTERACTIVE OBJECTS

	//take different action depending on the major interactive object identifier (this is a value from IOI_ enum)
	switch (iop.majorId) {

	//Shows program version update status : auxIdis the status as -1: attempting to connect, 0: connection failure, 1: program up to date, 2: update available
	case IOI_PROGRAMUPDATESTATUS:
	{
		//parameters from iop
		int status = iop.auxId;

		if (status != version_checking) {

			iop.auxId = version_checking;

			switch (version_checking) {
			case -1:
				pTO->set(" checking for updates... ");
				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				break;
			case 0:
				pTO->set(" couldn't connect ");
				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				break;
			case 1:
				pTO->set(" updated ");
				pTO->SetBackgroundColor(ONCOLOR);
				break;
			case 2:
				pTO->set(" new version available - click here ");
				pTO->SetBackgroundColor(OFFCOLOR);
				break;
			}

			stateChanged = true;
		}
	}
	break;

		//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). textId is the mesh name (if associated with this data type)
		//Note this entry must always represent the entry in Simulation::dataBoxList with the index in auxId.
	case IOI_DATABOXFIELDLABEL:
	{
		//parameters from iop
		int dataBoxList_idminor = iop.minorId;
		int DataBox_index = iop.auxId;			//this data box entry should represent the element with this index in dataBoxList
		string meshName = iop.textId;

		//this is the index corresponding to the dataBoxList_idminor - on any mismatch just reconstruct the data box entry to correspond to the element with DataBox_index index in dataBoxList
		int index_in_list = dataBoxList.get_index_from_id(INT2(0, dataBoxList_idminor));

		if (DataBox_index > dataBoxList.last()) {

			//too many fields : get rid of excess fields.
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
			break;
		}

		string actualmeshName = dataBoxList[DataBox_index].meshName;

		//if displayed meshname doesn't match the actual mesh name, or if indexes don't match, update Label (the n-th entry in the data box should represent the n-th entry in dataBoxList).
		if ((actualmeshName.length() && meshName != actualmeshName) || index_in_list != DataBox_index) {

			//meshname is set but doesn't match displayes name: update it.
			string newObjectText;
			if (actualmeshName.length()) {

				iop.textId = actualmeshName;
				newObjectText = "<" + actualmeshName + "> " + dataDescriptor(dataBoxList[DataBox_index].datumId).Label;
			}
			else newObjectText = dataDescriptor(dataBoxList[DataBox_index].datumId).Label;

			iop.minorId = dataBoxList.get_id_from_index(DataBox_index).minor;
			pTO->set(newObjectText);
			stateChanged = true;
		}
	}
	break;

	//A set or available module for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to
	case IOI_MODULE:
	{
		//parameters from iop
		MOD_ module = (MOD_)iop.minorId;
		int meshId = iop.auxId;

		//from unique mesh id number get index in pMesh (held in SMesh)
		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0 && SMesh[meshIdx]->IsModuleSet(module)) {

			//the mesh is contained and the module is set : ON color
			pTO->SetBackgroundColor(ONCOLOR);
		}
		else if (meshIdx < 0) {

			//the mesh is not contained : must have been deleted - delete this object
			stateChanged = true;
			iop.state = IOS_DELETING;
		}
		else pTO->SetBackgroundColor(OFFCOLOR);		//the mesh is contained but module not active : OFF color
	}
	break;

	//super-mesh module : minorId is an entry from MOD_ enum
	case IOI_SMODULE:
	{
		//parameters from iop
		MOD_ module = (MOD_)iop.minorId;

		if (SMesh.IsSuperMeshModuleSet(module))
			pTO->SetBackgroundColor(ONCOLOR);
		else pTO->SetBackgroundColor(OFFCOLOR);
	}
	break;

	//Available/set ode and evaluation method for magnetization : minorIdis an entry from ODE_ (the equation), auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	case IOI_ODE:
	{
		ODE_ actual_odeID;
		EVAL_ actual_evalID;
		SMesh.QueryODE(actual_odeID, actual_evalID);

		//parameters from iop
		ODE_ odeID = (ODE_)iop.minorId;
		EVAL_ evalID = (EVAL_)iop.auxId;

		if (actual_odeID == odeID && actual_evalID == evalID) {

			pTO->SetBackgroundColor(ONCOLOR);
		}
		else pTO->SetBackgroundColor(OFFCOLOR);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMS:
	{
		display_meshIO(&Simulation::Build_MeshParams_Line);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMSTEMP:
	{
		display_meshIO(&Simulation::Build_MeshParamsTemp_Text);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMSVAR:
	{
		display_meshIO(&Simulation::Build_MeshParamsVariation_Text);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORMODULES:
	{
		display_meshIO(&Simulation::Build_Modules_ListLine);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORMESHLIST:
	{
		display_meshIO(&Simulation::Build_Mesh_ListLine);
	}
	break;

	case IOI_MESH_FORDISPLAYOPTIONS:
	{
		display_meshIO(&Simulation::Build_MeshDisplay_ListLine);
	}
	break;

	case IOI_MESH_FORTEMPERATURE:
	{
		display_meshIO(&Simulation::Build_MeshTemperature_ListLine);
	}
	break;

	case IOI_MESH_FORHEATBOUNDARIES:
	{
		display_meshIO(&Simulation::Build_HeatBoundaries_ListLine);
	}
	break;

	case IOI_MESH_FORCURIEANDMOMENT:
	{
		display_meshIO(&Simulation::Build_CurieandMoment_ListLine);
	}
	break;

	//Shows mesh rectangle (units m) : minorId is the unique mesh id number, textId is the mesh rectangle
	case IOI_MESHRECTANGLE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		string rectValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			//update mesh rectangle if not matching
			Rect meshRect = SMesh[meshIdx]->GetMeshRect();
			if (ToString(meshRect, "m") != rectValue) {

				iop.textId = ToString(meshRect, "m");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}
		}
		else {

			//mesh no longer exists : delete the entire paragraph containing this object
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
	}
	break;

	//Shows ferromagnetic super-mesh rectangle (unit m) : textId is the mesh rectangle for the ferromagnetic super-mesh
	case IOI_FMSMESHRECTANGLE:
	{
		//parameters from iop
		string rectValue = iop.textId;

		//update mesh rectangle if not matching
		Rect meshRect = SMesh.GetFMSMeshRect();
		if (ToString(meshRect, "m") != rectValue) {

			iop.textId = ToString(meshRect, "m");
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows electric super-mesh rectangle (unit m) : textId is the mesh rectangle for the ferromagnetic super-mesh
	case IOI_ESMESHRECTANGLE:
	{
		//parameters from iop
		string rectValue = iop.textId;

		//update mesh rectangle if not matching
		Rect meshRect = SMesh.GetESMeshRect();
		
		if (meshRect.IsNull()) {

			if (rectValue != "N/A") {
				
				iop.textId = "N/A";
				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(OFFCOLOR);

				stateChanged = true;
			}
		}
		else if (ToString(meshRect, "m") != rectValue) {

			iop.textId = ToString(meshRect, "m");
			pTO->set(" " + iop.textId + " ");
			pTO->SetBackgroundColor(ONCOLOR);

			stateChanged = true;
		}
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		string cellsizeValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (enabled) {

				if (!SMesh[meshIdx]->MComputation_Enabled()) {

					iop.textId = "N/A";
					iop.auxId = 0;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(OFFCOLOR);
					stateChanged = true;
				}
				else {
					//update mesh cellsize if not matching
					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshCellsize();
					if (ToString(meshCellsize, "m") != cellsizeValue) {

						iop.textId = ToString(meshCellsize, "m");
						pTO->set(" " + iop.textId + " ");
						stateChanged = true;
					}
				}
			}
			else {

				if (SMesh[meshIdx]->MComputation_Enabled()) {

					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshCellsize();
					iop.textId = ToString(meshCellsize, "m");
					iop.auxId = 1;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(ONCOLOR);
					stateChanged = true;
				}
			}
		}
		else {

			//mesh no longer exists : delete the entire paragraph containing this object
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHECELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		string cellsizeValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (enabled) {

				if (!SMesh[meshIdx]->EComputation_Enabled()) {

					iop.textId = "N/A";
					iop.auxId = 0;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(OFFCOLOR);
					stateChanged = true;
				}
				else {
					//update mesh cellsize if not matching
					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshECellsize();
					if (ToString(meshCellsize, "m") != cellsizeValue) {

						iop.textId = ToString(meshCellsize, "m");
						pTO->set(" " + iop.textId + " ");
						stateChanged = true;
					}
				}
			}
			else {

				if (SMesh[meshIdx]->EComputation_Enabled()) {

					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshECellsize();
					iop.textId = ToString(meshCellsize, "m");
					iop.auxId = 1;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(ONCOLOR);
					stateChanged = true;
				}
			}
		}
		else {

			//mesh no longer exists : delete the entire paragraph containing this object
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHTCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		string cellsizeValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (enabled) {

				if (!SMesh[meshIdx]->TComputation_Enabled()) {

					iop.textId = "N/A";
					iop.auxId = 0;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(OFFCOLOR);
					stateChanged = true;
				}
				else {
					//update mesh cellsize if not matching
					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshTCellsize();
					if (ToString(meshCellsize, "m") != cellsizeValue) {

						iop.textId = ToString(meshCellsize, "m");
						pTO->set(" " + iop.textId + " ");
						stateChanged = true;
					}
				}
			}
			else {

				if (SMesh[meshIdx]->TComputation_Enabled()) {

					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshTCellsize();
					iop.textId = ToString(meshCellsize, "m");
					iop.auxId = 1;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(ONCOLOR);
					stateChanged = true;
				}
			}
		}
		else {

			//mesh no longer exists : delete the entire paragraph containing this object
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
	}
	break;

	//Shows ferromagnetic super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_FMSMESHCELLSIZE:
	{
		//parameters from iop
		string cellsizeValue = iop.textId;

		//update mesh cellsize if not matching
		DBL3 meshCellsize = SMesh.GetFMSMeshCellsize();
		if (ToString(meshCellsize, "m") != cellsizeValue) {

			iop.textId = ToString(meshCellsize, "m");
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows electric super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_ESMESHCELLSIZE:
	{
		//parameters from iop
		string cellsizeValue = iop.textId;

		//update mesh cellsize if not matching
		DBL3 meshCellsize = SMesh.GetESMeshCellsize();
		Rect meshRect = SMesh.GetESMeshRect();

		if (meshRect.IsNull()) {

			if (cellsizeValue != "N/A") {

				iop.textId = "N/A";
				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(OFFCOLOR);

				stateChanged = true;
			}
		}
		else if (ToString(meshCellsize, "m") != cellsizeValue) {

			iop.textId = ToString(meshCellsize, "m");
			pTO->set(" " + iop.textId + " ");
			pTO->SetBackgroundColor(ONCOLOR);

			stateChanged = true;
		}
	}
	break;

	//Show currently set directory : textId is the directory
	case IOI_DIRECTORY:
	{
		//parameters from iop
		string directory_fromio = iop.textId;

		//update name if not matching
		if (directory != directory_fromio) {

			iop.textId = directory;
			pTO->set(" " + directory + " ");

			stateChanged = true;
		}
	}
	break;

	//Show currently set save data file : textId is the file name
	case IOI_SAVEDATAFILE:
	{
		//parameters from iop
		string savedataFile_fromiop = iop.textId;

		//update name if not matching
		if (savedataFile != savedataFile_fromiop) {

			iop.textId = savedataFile;
			pTO->set(" " + savedataFile + " ");

			stateChanged = true;
		}
	}
	break;

	//Show currently set image save file base : textId is the file name
	case IOI_SAVEIMAGEFILEBASE:
	{
		//parameters from iop
		string savedataFile_fromiop = iop.textId;

		//update name if not matching
		if (imageSaveFileBase != savedataFile_fromiop) {

			iop.textId = imageSaveFileBase;
			pTO->set(" " + imageSaveFileBase + " ");

			stateChanged = true;
		}
	}
	break;

	//Show flag status for data/image saving during a simulation : minorId is the flag value (boolean)
	case IOI_SAVEDATAFLAG:
	{
		//parameters from iop
		int status = iop.minorId;

		if (status != (int)saveDataFlag) {

			iop.minorId = (int)saveDataFlag;

			if (saveDataFlag) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}

			stateChanged = true;
		}
	}
	break;

	case IOI_SAVEIMAGEFLAG:
	{
		//parameters from iop
		int status = iop.minorId;

		if (status != (int)saveImageFlag) {

			iop.minorId = (int)saveImageFlag;

			if (saveImageFlag) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}

			stateChanged = true;
		}
	}
	break;

	//Show set output data : minorId is the minor id of elements in Simulation::saveDataList (major id there is always 0), auxId is the number of the interactive object in the list as it appears in the console, textId is the configured output data. 
	//Note this entry must always represent the entry in Simulation::saveDataList with the index in auxId.
	case IOI_OUTDATA:
	{
		//parameters from iop
		int outDataId = iop.minorId;
		int io_index = iop.auxId;
		string configuredOutData = iop.textId;

		int index_in_list = saveDataList.get_index_from_id(INT2(0, outDataId));

		if (io_index <= saveDataList.last() && (index_in_list != io_index || configuredOutData != Build_SetOutputData_Text(io_index))) {

			iop.minorId = saveDataList.get_id_from_index(io_index).minor;
			iop.textId = Build_SetOutputData_Text(io_index);

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}

		//this object is part of a list : make sure this list is updated
		updateList(io_index, saveDataList.last(), &Simulation::Build_SetOutputData_ListLine);
	}
	break;

	//Shows a stage added to the simulation schedule : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the configured stage text
	//Note this entry must always represent the entry in Simulation::simStages with the index in auxId.
	case IOI_SETSTAGE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		int io_index = iop.auxId;
		string configuredSetStage = iop.textId;

		int index_in_list = simStages.get_index_from_id(INT2(0, stageId_minor));

		//if there's a mismatch between the object number and the actual index in saveDataList then updating is needed (also needed if meshname or box are mismatched) - update the entire object so that it corresponds to the entry in saveDataList at io_index.
		if (io_index <= simStages.last() && (index_in_list != io_index || configuredSetStage != Build_SetStages_Text(io_index))) {

			//because there are multiple objects on this line, all of them must be replaced. The caller must do this.
			iop.state = IOS_REPLACINGPARAGRAPH;
			stateChanged.textMessage = Build_SetStages_ListLine(io_index);
			stateChanged = true;
			break;
		}

		//this object is part of a list : make sure this list is updated
		updateList(io_index, simStages.last(), &Simulation::Build_SetStages_ListLine);
	}
	break;


	//Shows the value to set for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the value as a string
	case IOI_SETSTAGEVALUE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		string stageValueText = iop.textId;

		//this is the value as a string
		string actualValuestring = simStages[INT2(0, stageId_minor)].get_value_string();

		if (stageValueText != actualValuestring) {

			iop.textId = actualValuestring;
			pTO->set(" " + actualValuestring + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a string
	case IOI_STAGESTOPCONDITION:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		int io_index = iop.auxId;
		string stopConditionText = iop.textId;

		if (stopConditionText != Build_SetStages_StopConditionText(io_index)) {

			iop.textId = Build_SetStages_StopConditionText(io_index);
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows the saving condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the DSAVE_ value for this data save type, textId is the save type and value as a string
	case IOI_DSAVETYPE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		DSAVE_ dSaveType = (DSAVE_)iop.auxId;
		string saveConditionText = iop.textId;

		//this is the actual save type set
		DSAVE_ dsaveTypeSet = simStages[INT2(0, stageId_minor)].dsave_type();

		//set on or off color
		if (dsaveTypeSet != dSaveType) {

			if (iop.state == IOS_ON) {

				//this data save type not enabled anymore - reset background color and text
				pTO->SetBackgroundColor(OFFCOLOR);
				iop.textId = dataSaveDescriptors.get_key_from_ID(dSaveType);
				pTO->set(" " + iop.textId + " ");

				iop.state = IOS_OFF;
				stateChanged = true;
			}
		}
		else {

			//this save type is active
			if (iop.state == IOS_OFF) {

				//show it as enabled now
				pTO->SetBackgroundColor(ONCOLOR);
				iop.state = IOS_ON;

				stateChanged = true;
			}

			//check if object text matches actual data save condition including value
			int io_index = simStages.get_index_from_id(INT2(0, stageId_minor));
			int saveType_index = dataSaveDescriptors.get_index_from_ID(dSaveType);

			if (saveConditionText != Build_SetStages_SaveConditionText(io_index, saveType_index)) {

				iop.textId = Build_SetStages_SaveConditionText(io_index, saveType_index);
				pTO->set(" " + iop.textId + " ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows parameter and value for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter handle and value
	case IOI_MESHPARAM:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		string paramText = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (paramText != Build_MeshParams_Text(meshIdx, paramId)) {

				iop.textId = Build_MeshParams_Text(meshIdx, paramId);
				pTO->set(" " + iop.textId + " ");
				stateChanged = true;
			}
		}
		else {

			//this mesh no longer exists, so delete all associated interactive object parameters
			iop.state = IOS_DELETINGPARAGRAPH;
			stateChanged = true;
		}
	}
	break;

	//Shows parameter temperature dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter temperature dependence setting
	case IOI_MESHPARAMTEMP:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		string paramText = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (paramText != SMesh[meshIdx]->get_paraminfo_string(paramId)) {

				iop.textId = SMesh[meshIdx]->get_paraminfo_string(paramId);
				pTO->set(" " + iop.textId + " ");

				if (SMesh[meshIdx]->is_paramtemp_set(paramId)) pTO->SetBackgroundColor(ONCOLOR);
				else pTO->SetBackgroundColor(OFFCOLOR);

				stateChanged = true;
			}
		}
		else {

			//this mesh no longer exists, so delete all associated interactive object parameters
			iop.state = IOS_DELETINGPARAGRAPH;
			stateChanged = true;
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

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->GetDisplayedParamVar() == paramId) pTO->SetBackgroundColor(ONCOLOR);
			else pTO->SetBackgroundColor(OFFCOLOR);

			if (paramText != SMesh[meshIdx]->get_paramvarinfo_string(paramId)) {

				iop.textId = SMesh[meshIdx]->get_paramvarinfo_string(paramId);
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}
		}
		else {

			//this mesh no longer exists, so delete all associated interactive object parameters
			iop.state = IOS_DELETINGPARAGRAPH;
			stateChanged = true;
		}
	}
	break;

	//Shows mesh display option for a given mesh : minorId is the MESHDISPLAY_ value, auxId is the unique mesh id number, textId is the MESHDISPLAY_ handle
	case IOI_MESHDISPLAY:
	{
		//parameters from iop
		MESHDISPLAY_ displayOption = (MESHDISPLAY_)iop.minorId;
		int meshId = iop.auxId;
		string displayHandle = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->GetDisplayedPhysicalQuantity() == displayOption) {

				//this display option enabled
				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				//this display option disabled
				pTO->SetBackgroundColor(OFFCOLOR);
			}
		}
	}
	break;

	//Shows super-mesh display option : minorId is the MESHDISPLAY_ value, textId is the MESHDISPLAY_ handle
	case IOI_SMESHDISPLAY:
	{
		//parameters from iop
		MESHDISPLAY_ displayOption = (MESHDISPLAY_)iop.minorId;
		string displayHandle = iop.textId;

		if (SMesh.GetDisplayedPhysicalQuantity() == displayOption) {

			//this display option enabled
			pTO->SetBackgroundColor(ONCOLOR);
		}
		else {

			//this display option disabled
			pTO->SetBackgroundColor(OFFCOLOR);
		}
	}
	break;

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	case IOI_MOVINGMESH:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool moving_mesh = iop.auxId;
		string meshName = iop.textId;

		//is there a state mismatch?
		if (moving_mesh != SMesh.IsMovingMeshSet()) {

			moving_mesh = SMesh.IsMovingMeshSet();
			iop.auxId = moving_mesh;

			if (!moving_mesh) {

				iop.minorId = -1;
				iop.textId = "";
				pTO->set(" None ");
				pTO->SetBackgroundColor(OFFCOLOR);
				stateChanged = true;
			}
			else {

				iop.minorId = SMesh.GetId_of_MoveMeshTrigger();
				iop.textId = SMesh.key_from_meshId(iop.minorId);

				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(ONCOLOR);
				stateChanged = true;
			}
		}
		else if (moving_mesh) {

			if (meshName != SMesh.key_from_meshId(SMesh.GetId_of_MoveMeshTrigger())) {

				iop.minorId = SMesh.GetId_of_MoveMeshTrigger();
				iop.textId = SMesh.key_from_meshId(iop.minorId);

				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(ONCOLOR);
				stateChanged = true;
			}
		}
	}
	break;

	//Shows movingmesh symmetry : auxId is the asymmetry status (1: asymmetric, 0: symmetric)
	case IOI_MOVINGMESHASYM:
	{
		//parameters from iop
		bool asymmetric = iop.auxId;

		if (asymmetric != SMesh.MoveMeshAntisymmetric()) {

			iop.auxId = !asymmetric;

			if (iop.auxId) {

				pTO->set(" Antisymmetric ");
				stateChanged = true;
			}
			else {

				pTO->set(" Symmetric ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows movingmesh threshold : textId is the threshold value as a string
	case IOI_MOVINGMESHTHRESH:
	{
		//parameters from iop
		string threshold_string = iop.textId;

		if (threshold_string != ToString(SMesh.MoveMeshThreshold())) {

			iop.textId = ToString(SMesh.MoveMeshThreshold());

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows electrode box. minorId is the minor Id in STransport::electrode_boxes, auxId is the number of the interactive object in the list (electrode index), textId is the electrode rect as a string
	case IOI_ELECTRODERECT:
	{
		//parameters from iop
		int electrodeId_minor = iop.minorId;
		int io_index = iop.auxId;
		string rect_string = iop.textId;

		//actual index in electrodes list for the electrode identifier (should normally be the same as io_index)
		int index_in_list = SMesh.CallModuleMethod(&STransport::GetElectrodeIndex, electrodeId_minor);
		int el_last_index = SMesh.CallModuleMethod(&STransport::GetNumberofElectrodes) - 1;

		//if there's a mismatch between the object number and the actual index then updating is needed - update the entire object so that it corresponds to the entry at io_index.
		if (io_index <= el_last_index && index_in_list != io_index) {

			//because there are multiple objects on this line, all of them must be replaced. The caller must do this.
			iop.state = IOS_REPLACINGPARAGRAPH;
			stateChanged.textMessage = Build_Electrodes_ListLine(io_index);
			stateChanged = true;
			break;
		}

		//this object is part of a list : make sure this list is updated
		updateList(io_index, el_last_index, &Simulation::Build_Electrodes_ListLine);

		if (rect_string != ToString(SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, io_index).first, "m")) {

			iop.textId = ToString(SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, io_index).first, "m");
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows electrode potential. minorId is the electrode index, textId is potential value as a string
	case IOI_ELECTRODEPOTENTIAL:
	{
		//parameters from iop
		int el_index = iop.minorId;
		string potential_string = iop.textId;

		if (ToString(SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, el_index).second, "V") != potential_string) {

			iop.textId = ToString(SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, el_index).second, "V");
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows electrode ground setting. minorId is the electrode index, auxId is the setting (0 : not ground, 1 : ground)
	case IOI_ELECTRODEGROUND:
	{
		//parameters from iop
		int el_index = iop.minorId;
		bool is_ground = (bool)iop.auxId;

		if (SMesh.CallModuleMethod(&STransport::IsGroundElectrode, el_index) != is_ground) {

			iop.auxId = SMesh.CallModuleMethod(&STransport::IsGroundElectrode, el_index);

			if (iop.auxId) pTO->SetBackgroundColor(ONCOLOR);
			else pTO->SetBackgroundColor(OFFCOLOR);

			stateChanged = true;
		}
	}
	break;

	//Shows constant current source setting. auxId is the setting.
	case IOI_CONSTANTCURRENTSOURCE:
	{
		//parameters from iop
		bool is_constant_current = (bool)iop.auxId;

		if (SMesh.CallModuleMethod(&STransport::UsingConstantCurrentSource) != is_constant_current) {

			iop.auxId = (int)SMesh.CallModuleMethod(&STransport::UsingConstantCurrentSource);

			if (iop.auxId) pTO->set(" constant current ");
			else pTO->set(" constant voltage ");

			stateChanged = true;
		}
	}
	break;

	//Shows transport solver convergence error. textId is the convergence error value.
	case IOI_TSOLVERCONVERROR:
	{
		//parameters from iop
		double conv_error = ToNum(iop.textId);

		if (conv_error != SMesh.CallModuleMethod(&STransport::GetConvergenceError)) {

			iop.textId = ToString(SMesh.CallModuleMethod(&STransport::GetConvergenceError));
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows transport solver timeout iterations. auxId is the timeout value.
	case IOI_TSOLVERTIMEOUT:
	{
		//parameters from iop
		int timeout = iop.auxId;

		if (timeout != SMesh.CallModuleMethod(&STransport::GetConvergenceTimeout)) {

			iop.auxId = SMesh.CallModuleMethod(&STransport::GetConvergenceTimeout);
			pTO->set(" " + ToString(iop.auxId) + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows spin-transport solver convergence error. textId is the convergence error value.
	case IOI_SSOLVERCONVERROR:
	{
		//parameters from iop
		double conv_error = ToNum(iop.textId);

		if (conv_error != SMesh.CallModuleMethod(&STransport::GetSConvergenceError)) {

			iop.textId = ToString(SMesh.CallModuleMethod(&STransport::GetSConvergenceError));
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows spin-transport solver timeout iterations. auxId is the timeout value.
	case IOI_SSOLVERTIMEOUT:
	{
		//parameters from iop
		int timeout = iop.auxId;

		if (timeout != SMesh.CallModuleMethod(&STransport::GetSConvergenceTimeout)) {

			iop.auxId = SMesh.CallModuleMethod(&STransport::GetSConvergenceTimeout);
			pTO->set(" " + ToString(iop.auxId) + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows Poisson solver SOR damping type : true for adaptive, false for fixed. auxId is enabled (1)/disabled(0) status.
	case IOI_SORFIXEDDAMPING:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (status != SMesh.CallModuleMethod(&STransport::IsFixedSORdamping)) {

			iop.auxId = !status;

			if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" Fixed ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Adaptive ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a string. (DBL2 since we need different damping values for V and S solvers)
	case IOI_SORDAMPING:
	{
		//parameters from iop
		string SOR_damping = iop.textId;

		if (SOR_damping != ToString(SMesh.CallModuleMethod(&STransport::GetSORDamping))) {

			iop.textId = ToString(SMesh.CallModuleMethod(&STransport::GetSORDamping));
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows mesh base temperature. minorId is the unique mesh id number, textId is the temperature value
	case IOI_BASETEMPERATURE:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0 && ToString(SMesh[meshIdx]->GetBaseTemperature(), "K") != temp_string) {

			iop.textId = ToString(SMesh[meshIdx]->GetBaseTemperature(), "K");
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows ambient temperature for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the temperature value
	case IOI_AMBIENT_TEMPERATURE:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->IsModuleSet(MOD_HEAT) != status) {

				iop.auxId = SMesh[meshIdx]->IsModuleSet(MOD_HEAT);

				if (iop.auxId) {

					pTO->SetBackgroundColor(ONCOLOR);
					iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&Heat::GetAmbientTemperature), "K");
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
			else if (status && ToString(SMesh[meshIdx]->CallModuleMethod(&Heat::GetAmbientTemperature), "K") != temp_string) {

				iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&Heat::GetAmbientTemperature), "K");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}
		}
	}
	break;

	//Shows alpha value (W/m^2K) for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the value
	case IOI_ROBIN_ALPHA:
	{
		int meshId = iop.minorId;
		string alpha_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->IsModuleSet(MOD_HEAT) != status) {

				iop.auxId = SMesh[meshIdx]->IsModuleSet(MOD_HEAT);

				if (iop.auxId) {

					pTO->SetBackgroundColor(ONCOLOR);
					iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&Heat::GetAlphaBoundary), "W/m2K");
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
			else if (status && ToString(SMesh[meshIdx]->CallModuleMethod(&Heat::GetAlphaBoundary), "W/m2K") != alpha_string) {

				iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&Heat::GetAlphaBoundary), "W/m2K");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}
		}
	}
	break;

	//Shows temperature insulating side setting for heat equation. minorId is the unique mesh id number, auxId is the status (Heat module must be active) : -1 disabled (gray), 0 not insulating (green), 1 insulating (red), textId represents the side : "x", "-x", "y", "-y", "z", "-z"
	case IOI_INSULATINGSIDE:
	{
		int meshId = iop.minorId;
		string literal = iop.textId;
		int status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->IsModuleSet(MOD_HEAT) != bool(status + 1)) {

				bool heat_set = SMesh[meshIdx]->IsModuleSet(MOD_HEAT);
				if (!heat_set) iop.auxId = -1;
				else iop.auxId = SMesh[meshIdx]->CallModuleMethod(&Heat::GetInsulatingSide, literal);

				if (heat_set) {

					if (iop.auxId == 0) {

						pTO->SetBackgroundColor(ONCOLOR);
						pTO->set(" " + literal + ": No ");
					}
					else {

						pTO->SetBackgroundColor(OFFCOLOR);
						pTO->set(" " + literal + ": Yes ");
					}	
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
			else if (status >= 0 && SMesh[meshIdx]->CallModuleMethod(&Heat::GetInsulatingSide, literal) != (bool)status) {

				iop.auxId = SMesh[meshIdx]->CallModuleMethod(&Heat::GetInsulatingSide, literal);
				pTO->set(" " + iop.textId + " ");

				if (iop.auxId == 0) {

					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" " + literal + ": No ");
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
					pTO->set(" " + literal + ": Yes ");
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows mesh Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	case IOI_CURIETEMP:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {
			
			if (ToString(SMesh[meshIdx]->GetCurieTemperature(), "K") != temp_string) {

				iop.textId = ToString(SMesh[meshIdx]->GetCurieTemperature(), "K");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->Magnetisation_Enabled()) != status) {

				iop.auxId = (SMesh[meshIdx]->Magnetisation_Enabled());

				if (iop.auxId == 1) {

					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows indicative material Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	case IOI_CURIETEMPMATERIAL:
	{
		int meshId = iop.minorId;
		string temp_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (ToString(SMesh[meshIdx]->GetCurieTemperatureMaterial(), "K") != temp_string) {

				iop.textId = ToString(SMesh[meshIdx]->GetCurieTemperatureMaterial(), "K");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->Magnetisation_Enabled()) != status) {

				iop.auxId = (SMesh[meshIdx]->Magnetisation_Enabled());

				if (iop.auxId == 1) {

					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows atomic moment multiple of Bohr magneton. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the value
	case IOI_ATOMICMOMENT:
	{
		int meshId = iop.minorId;
		string amoment_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (ToString(SMesh[meshIdx]->GetAtomicMoment(), "uB") != amoment_string) {

				iop.textId = ToString(SMesh[meshIdx]->GetAtomicMoment(), "uB");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->GetMeshType() == MESH_FERROMAGNETIC) != status) {

				iop.auxId = (SMesh[meshIdx]->GetMeshType() == MESH_FERROMAGNETIC);

				if (iop.auxId == 1) {

					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows cuda enabled/disabled or n/a state. auxId is enabled (1)/disabled(0)/not available(-1) status.
	case IOI_CUDASTATE:
	{
		int status = iop.auxId;

		//if status was set to -1 then cuda is not available and will not be for the duration of this program execution, so nothing to do
		if (status >= 0) {

			if (status != (int)cudaEnabled) {

				iop.auxId = cudaEnabled;

				if (iop.auxId == 1) {

					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" On ");
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
					pTO->set(" Off ");
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows scale_rects enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCALERECTSSTATUS:
	{
		bool status = iop.auxId;

		if (status != SMesh.Get_Scale_Rects()) {

			iop.auxId = SMesh.Get_Scale_Rects();

			if (iop.auxId) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows coupled_to_dipoles enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_COUPLEDTODIPOLESSTATUS:
	{
		bool status = iop.auxId;

		if (status != SMesh.Get_Coupled_To_Dipoles()) {

			iop.auxId = SMesh.Get_Coupled_To_Dipoles();

			if (iop.auxId) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, auxId is enabled (1)/disabled(0) status. textId is the value
	case IOI_REFINEROUGHNESS:
	{
		int meshId = iop.minorId;
		string refine_string = iop.textId;
		bool status = iop.auxId;
		
		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (status != SMesh[meshIdx]->IsModuleSet(MOD_ROUGHNESS)) {

				status = SMesh[meshIdx]->IsModuleSet(MOD_ROUGHNESS);
				iop.auxId = status;

				if (iop.auxId) {

					pTO->SetBackgroundColor(ONCOLOR);
					iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&Roughness::get_refine));
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}

			if (status && refine_string != ToString(SMesh[meshIdx]->CallModuleMethod(&Roughness::get_refine))) {

				iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&Roughness::get_refine));
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}
		}
	}
	break;

	//Shows status of multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_MULTICONV:
	{
		int status = iop.auxId;

		if (status != SMesh.Get_Multilayered_Convolution_Status()) {

			iop.auxId = SMesh.Get_Multilayered_Convolution_Status();

			if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else if (iop.auxId == 0) {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}
			else {

				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				pTO->set(" N/A ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows status of force 2D multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_2DMULTICONV:
	{
		int status = iop.auxId;

		if (status != SMesh.Get_2D_Multilayered_Convolution_Status()) {

			iop.auxId = SMesh.Get_2D_Multilayered_Convolution_Status();

			if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else if (iop.auxId == 0) {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}
			else {

				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				pTO->set(" N/A ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows status of use default n for multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_NCOMMONSTATUS:
	{
		int status = iop.auxId;

		if (status != SMesh.Use_Default_n_Status()) {

			iop.auxId = SMesh.Use_Default_n_Status();

			if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else if (iop.auxId == 0) {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
			}
			else {

				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				pTO->set(" N/A ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows n_common for multi-layered convolution. auxId is the status (-1 : N/A, otherwise available). textId is the value as a SZ3.
	case IOI_NCOMMON:
	{
		string common_n = iop.textId;

		if (common_n != ToString(SMesh.Get_n_common())) {

			iop.textId = ToString(SMesh.Get_n_common());

			if (SMesh.Get_n_common() == SZ3()) {

				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				pTO->set(" N/A ");

				iop.auxId = -1;
			}
			else {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" " + iop.textId + " ");

				iop.auxId = 0;
			}

			stateChanged = true;
		}
	}
	break;

	//Shows materials database in use. textId is the name of the database, including the path.
	case IOI_LOCALMDB:
	{
		//parameters from iop
		string mdbFile = iop.textId;

		//update name if not matching
		if (mdbFile != mdb.GetDataBaseName()) {

			iop.textId = mdb.GetDataBaseName();
			pTO->set(" " + mdb.GetDataBaseName() + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows gpu free memory. auxId is the value
	case IOI_GPUMEMFREE:
	{
		size_t mem_size = iop.auxId;

		if (mem_size != gpuMemFree_MB) {

			iop.auxId = gpuMemFree_MB;
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows gpu total memory. auxId is the value
	case IOI_GPUMEMTOTAL:
	{
		size_t mem_size = iop.auxId;

		if (mem_size != gpuMemTotal_MB) {

			iop.auxId = gpuMemTotal_MB;
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows cpu free memory. auxId is the value
	case IOI_CPUMEMFREE:
	{
		size_t mem_size = iop.auxId;

		if (mem_size != cpuMemFree_MB) {

			iop.auxId = cpuMemFree_MB;
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows cpu total memory. auxId is the value
	case IOI_CPUMEMTOTAL:
	{
		size_t mem_size = iop.auxId;

		if (mem_size != cpuMemTotal_MB) {

			iop.auxId = cpuMemTotal_MB;
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	default:
		break;
	}

	return stateChanged;
}

#endif