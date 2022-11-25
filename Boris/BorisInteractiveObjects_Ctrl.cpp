#include "stdafx.h"

//Interactive Objects : 
//
//handle user interactions using ConsoleActionHandler
//update their state depending on current program state using ConsoleInteractiveObjectState

//These are not needed in non-graphical mode
#include "CompileFlags.h"
#if GRAPHICS == 1

#include "Simulation.h"

InteractiveObjectActionOutcome Simulation::ConsoleActionHandler(int actionCode, InteractiveObjectProperties& iop, TextObject *pTO)
{
	//!!!IMPORTANT!!! Do not access BorisDisplay through thread-safe entry points from here. This method was called from within BorisDisplay (through a function pointer), which was thread-safe accessed so the std::mutex is now locked.
	//In fact it's better not to make any calls to BorisDisplay here at all (e.g. DisplayConsoleMessage). There should always be a better way. e.g. an interactive object will typically require a console
	//command to be issued - in this case just launch the command on the THREAD_HANDLEMESSAGE thread. If that command needs any console text displayed, it will wait for the display std::mutex to become available.

	//Makes changes in the program depending on the interactive object properties and the action type with which it was interacted with. No changes to the object are made here, that is handled by ConsoleInteractiveObjectState during Refresh cycles

	InteractiveObjectActionOutcome actionOutcome = AO_NOTHING;

	//action codes which can be handled in the same way for all objects
	if (actionCode == AC_HOVERCHECK) {

		//mouse is hovering over interactive object - display hover info?
		
		//the index in ioInfo to obtain info std::string
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

			//at most one entry for this majorId, so this could be a generic info std::string applicable for any minorId - if no entry at all then nothing to display
			if (ioInfo.is_ID_set(iop.majorId)) {

				actionOutcome.text = ioInfo(iop.majorId);
			}
		}

		//only ask to display if we managed to retrieve a non-empty std::string
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

		std::string command = "~" + commands.get_key_from_index(commandCode);

		for (std::string param : { ToString(params)... })
			command += " " + param;

		//launch command
		single_call_launch<std::string>(&Simulation::HandleCommand, command, THREAD_HANDLEMESSAGE);
	};

	//SEND COMMAND TO COMMAND HANDLER (verbose) : makes code neater using this quick lambda
	auto sendCommand_verbose = [&](CMD_ commandCode, auto ...params) {

		std::string command = commands.get_key_from_index(commandCode);

		for (std::string param : { ToString(params)... })
			command += " " + param;

		//launch command
		single_call_launch<std::string>(&Simulation::HandleCommand, command, THREAD_HANDLEMESSAGE);
		single_call_launch<std::string>(&Simulation::SetConsoleEntryLineText, "", THREAD_HANDLEMESSAGE2);
	};

	//SET CONSOLE ENTRY LINE : form a command syntax and set it as the console entry line for further editing by user
	auto setConsoleEntry = [&](CMD_ commandCode, auto ...params) {

		std::string text = commands.get_key_from_index(commandCode);

		for (std::string param : { ToString(params)... })
			text += " " + param;

		single_call_launch<std::string>(&Simulation::SetConsoleEntryLineText, text, THREAD_HANDLEMESSAGE2);
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

			sendCommand_verbose(CMD_DELPINNEDDATA, DataBox_index);

			//need to update values as well as labels (values must correspond to labels)
			single_call_launch(&Simulation::UpdateDataBox_Refresh, THREAD_HANDLEMESSAGE);
		}

		//user might be trying to re-arrange data order : start interaction by creating a moveable pop-up window which holds this data entry
		else if (actionCode == AC_MOUSELEFTDOWN) { actionOutcome = AO_STARTINTERACTION; }

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_DATABOXFIELDLABEL
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == IOI_DATABOXFIELDLABEL) {

			if (iop.interactingObjectId.minor != dataBoxList_idminor) {

				interaction(dataBoxList, iop.interactingObjectId.minor, dataBoxList_idminor);

				//need to update values as well as labels (values must correspond to labels)
				single_call_launch(&Simulation::UpdateDataBox_Refresh, THREAD_HANDLEMESSAGE);

				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}
	}
	break;

	//A set or available module for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to
	case IOI_MODULE:
	{
		//parameters from iop
		MOD_ moduleID = (MOD_)iop.minorId;
		int meshId = iop.auxId;

		//try to add module
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_ADDMODULE, SMesh.key_from_meshId(meshId), moduleHandles(moduleID));

		//try to remove module
		else if (actionCode == AC_MOUSERIGHTDOWN) {

			//special treatment of MOD_DEMAG to control if excluded from multilayered demag convolution
			if (moduleID == MOD_DEMAG && SMesh.IsSuperMeshModuleSet(MODS_SDEMAG)) {

				int meshIdx = SMesh.contains_id(meshId);
				if (meshIdx >= 0) {

					sendCommand_verbose(CMD_EXCLUDEMULTICONVDEMAG, !SMesh[meshIdx]->Get_Demag_Exclusion(), SMesh.key_from_meshId(meshId));
				}
			}
			else sendCommand_verbose(CMD_DELMODULE, SMesh.key_from_meshId(meshId), moduleHandles(moduleID));
		}
	}
	break;

	//super-mesh module : minorId is an entry from MOD_ enum
	case IOI_SMODULE:
	{
		//parameters from iop
		MOD_ moduleID = (MOD_)iop.minorId;

		//try to add module
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_ADDMODULE, SMesh.superMeshHandle, moduleHandles(moduleID));

		//try to remove module
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELMODULE, SMesh.superMeshHandle, moduleHandles(moduleID));
	}
	break;

	//Module used for effective field display for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to, textId is the MOD_ value used for display
	case IOI_DISPLAYMODULE:
	{
		//parameters from iop
		MOD_ moduleID = (MOD_)iop.minorId;
		int meshId = iop.auxId;

		//try to add module to display
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DISPLAYMODULE, moduleHandles(moduleID), SMesh.key_from_meshId(meshId));

		//remove module from display
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DISPLAYMODULE, "none", SMesh.key_from_meshId(meshId));
	}
	break;

	//Available/set ode : minorId is an entry from ODE_ (the equation)
	case IOI_ODE:
	{
		//parameters from iop
		ODE_ odeID = (ODE_)iop.minorId;

		//try to set ODE and its default evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) {

			//choose a default eval applicable for both atomistic and micromagnetic solvers
			EVAL_ defaultEval = odeDefaultEval(odeID);
			ODE_ atom_odeID;
			SMesh.QueryAtomODE(atom_odeID);
			if (!vector_contains(odeAllowedEvals(atom_odeID), defaultEval)) defaultEval = EVAL_AHEUN;

			sendCommand_verbose(CMD_SETODE, odeHandles(odeID), odeEvalHandles(defaultEval));
		}
	}
	break;

	//Available/set ode for atomistic meshes: minorId is an entry from ODE_ (the equation)
	case IOI_ATOMODE:
	{
		//parameters from iop
		ODE_ atom_odeID = (ODE_)iop.minorId;

		//try to set ODE and its default evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) {

			//choose a default eval applicable for both atomistic and micromagnetic solvers
			EVAL_ defaultEval = odeDefaultEval(atom_odeID);
			ODE_ odeID;
			SMesh.QueryODE(odeID);
			if (!vector_contains(odeAllowedEvals(odeID), defaultEval)) defaultEval = EVAL_AHEUN;

			sendCommand_verbose(CMD_SETATOMODE, odeHandles(atom_odeID), odeEvalHandles(defaultEval));
		}
	}
	break;

	//Set ODE time step: textId is the value
	case IOI_ODEDT:
	{
		//parameters from iop
		std::string dT_string = ToNum(iop.textId);

		//try to set ODE time step
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SETDT, trimspaces(to_text));
		}
	}
	break;

	//Set stochastic time-step: textId is the value
	case IOI_STOCHDT:
	{
		//parameters from iop
		std::string dT_string = ToNum(iop.textId);

		//try to set ODE time step
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SETDTSTOCH, trimspaces(to_text));
		}
	}
	break;

	//Link stochastic time-step to ODE dT flag : auxId is the value
	case IOI_LINKSTOCHDT:
	{
		//parameters from iop
		bool state = (bool)iop.auxId;

		//try to set ODE and its default evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_LINKDTSTOCHASTIC, !state);
	}
	break;

	//Set evaluation speedup time-step: textId is the value
	case IOI_SPEEDUPDT:
	{
		//parameters from iop
		std::string dT_string = ToNum(iop.textId);

		//try to set ODE time step
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SETDTSPEEDUP, trimspaces(to_text));
		}
	}
	break;

	//Link evaluation speedup time-step to ODE dT flag : auxId is the value
	case IOI_LINKSPEEDUPDT:
	{
		//parameters from iop
		bool state = (bool)iop.auxId;

		//try to set ODE and its default evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_LINKDTSPEEDUP, !state);
	}
	break;

	//Set heat equation time step: textId is the value
	case IOI_HEATDT:
	{
		//try to set ODE time step
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SETHEATDT, trimspaces(to_text));
		}
	}
	break;

	//Set elastodynamics equation time step: textId is the value
	case IOI_ELDT:
	{
		//try to set ODE time step
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SETELDT, trimspaces(to_text));
		}
	}
	break;

	//Link elastodynamics time-step to ODE dT flag : auxId is the value
	case IOI_LINKELDT:
	{
		//parameters from iop
		bool state = (bool)iop.auxId;

		//try to set ODE and its default evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_LINKDTELASTIC, !state);
	}
	break;

	//Available/set evaluation method for ode : minorId is an entry from ODE_ as : micromagnetic equation value + 100 * atomistic equation value, auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	case IOI_ODE_EVAL:
	{
		//parameters from iop
		std::string evalHandle = iop.textId;

		//try to set ODE and its evaluation method
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_SETODEEVAL, evalHandle);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMS:
	case IOI_MESH_FORPARAMSTEMP:
	case IOI_MESH_FORPARAMSVAR:
	{
		//parameters from iop
		std::string meshName = iop.textId;

		//try to change mesh focus
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_MESHFOCUS, meshName);

		//try to delete mesh
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELMESH, meshName);

		//rename mesh : bring up console command
		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_RENAMEMESH, trimspaces(to_text));
		};
	}
	break;

	case IOI_MESH_FORCURIEANDMOMENT:
	case IOI_MESH_FORPBC:
	case IOI_MESH_FOREXCHCOUPLING:
	case IOI_MESH_FORTMODEL:
	case IOI_MESH_FORHEATBOUNDARIES:
	case IOI_MESH_FORTEMPERATURE:
	case IOI_MESH_FORDISPLAYOPTIONS:
	case IOI_MESH_FORMODULES:
	case IOI_MESH_FORDISPLAYMODULES:
	case IOI_MESH_FORMESHLIST:
	case IOI_MESH_FORSTOCHASTICITY:
	case IOI_MESH_FORELASTICITY:
	case IOI_MESH_FORSPEEDUP:
	case IOI_MESH_FORSKYPOSDMUL:
	case IOI_MESH_FORMC:
	case IOI_MESH_FORDIPOLESHIFT:
	case IOI_MESH_FORTMR:
	{
		//parameters from iop
		std::string meshName = iop.textId;
		int meshId = iop.minorId;

		//try to change mesh focus
		if (actionCode == AC_MOUSELEFTDOWN) {

			sendCommand_verbose(CMD_MESHFOCUS, meshName);
	
			actionOutcome = AO_STARTINTERACTION;
		}

		//try to delete mesh
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELMESH, meshName);

		//rename mesh : bring up console command
		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_RENAMEMESH, trimspaces(to_text));
		}

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_MESH_FORMESHLIST
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == iop.majorId) {

			//swapping mesh positions in mesh list
			if (iop.interactingObjectId.minor != meshId) {

				int idxSrc = SMesh.contains_id(iop.interactingObjectId.minor);
				int idxDst = SMesh.contains_id(meshId);

				if (idxSrc >= 0 && idxDst >= 0) {

					//swap mesh positions in SMesh.pMesh list, and also swap their ids.
					//swapping the ids is easier, otherwise we have to do some book keeping. 
					SMesh[idxSrc]->swap_ids(SMesh[idxDst]);
					SMesh().move(idxSrc, idxDst);
				}

				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}
	}
	break;

	//Shows mesh rectangle (units m) : minorId is the unique mesh id number, textId is the mesh rectangle
	case IOI_MESHRECTANGLE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		std::string meshRect_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { 
			
			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX; 
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_MESHRECT, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsize_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_CELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHECELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsize_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_ECELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows thermal mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHTCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsize_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_TCELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHMCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsize_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_MCELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows stochastic cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHSCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsize_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SCELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows link stochastic flag : minorId is the unique mesh id number, auxId is the value off (0), on (1), N/A (-1)
	case IOI_LINKSTOCHASTIC:
	{
		int meshId = iop.minorId;
		int status = iop.auxId;

		if (status >= 0) {

			if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_LINKSTOCHASTIC, (status + 1) % 2, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows macrocell size (units m) for atomistic meshes: minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHDMCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsize_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && enabled) {

			sendCommand_verbose(CMD_MESHFOCUS, SMesh.key_from_meshId(meshId));
			actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_ATOMDMCELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows evaluation speedup type: auxId is the type value.
	case IOI_SPEEDUPMODE:
	{
		//parameters from iop
		int option = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_EVALSPEEDUP, (option + 1) % EVALSPEEDUP_NUMENTRIES);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_EVALSPEEDUP, (option + EVALSPEEDUP_NUMENTRIES - 1) % EVALSPEEDUP_NUMENTRIES);
	}
	break;

	//Shows ferromagnetic super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_FMSMESHCELLSIZE:
	{
		std::string cellsizeValue = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_FMSCELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Shows electric super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_ESMESHCELLSIZE:
	{
		std::string cellsizeValue = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_ESCELLSIZE, combine(split(trimspaces(to_text), ",", ";"), " "));
		};
	}
	break;

	//Simulation output data, specifically used for showing values in console : minorId is the DATA_ id, textId is the data handle
	case IOI_SHOWDATA:
	{
		//parameters from iop
		std::string dataHandle = iop.textId;

		if (actionCode == AC_DOUBLECLICK) sendCommand_verbose(CMD_SHOWDATA, dataHandle);

		else if (actionCode == AC_MOUSELEFTDOWN) {

			setConsoleEntry(CMD_SHOWDATA, dataHandle);

			//this object can be dragged into the data box to display there
			actionOutcome = AO_STARTINTERACTION;
		}

		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_ADDPINNEDDATA, dataHandle);

		else if (actionCode == AC_INTERACTOBJECTWITHWINDOW) {

			if (iop.interactingObjectId == INT2(WIN_DATABOX, 0)) {

				sendCommand_verbose(CMD_ADDPINNEDDATA, dataHandle);

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
		std::string dataHandle = iop.textId;

		if (actionCode == AC_DOUBLECLICK) sendCommand_verbose(CMD_ADDDATA, dataHandle);

		else if (actionCode == AC_MOUSELEFTDOWN) actionOutcome = AO_STARTINTERACTION;

		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_ADDPINNEDDATA, dataHandle);

		else if (actionCode == AC_INTERACTOBJECTWITHWINDOW) {

			if (iop.interactingObjectId == INT2(WIN_DATABOX, 0)) {

				sendCommand_verbose(CMD_ADDPINNEDDATA, dataHandle);

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
		std::string directory = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_CHDIR, trimendspaces(to_text));
		}
	}
	break;

	//Show currently set save data file : textId is the file name
	case IOI_SAVEDATAFILE:
	{
		//parameters from iop
		std::string savedataFile = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SAVEDATAFILE, trimendspaces(to_text));
		}
	}
	break;

	//Show currently set image file base : textId is the file name
	case IOI_SAVEIMAGEFILEBASE:
	{
		//parameters from iop
		std::string imageSaveFileBase = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();
			sendCommand_verbose(CMD_SAVEIMAGEFILE, trimendspaces(to_text));
		}
	}
	break;

	//Show flag status for data/image saving during a simulation : minorId is the flag value (boolean)
	case IOI_SAVEDATAFLAG:
	{
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DATASAVEFLAG, !saveDataFlag);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DATASAVEFLAG, false);
	}
	break;

	case IOI_SAVEIMAGEFLAG:
	{
		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_IMAGESAVEFLAG, !saveImageFlag);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_IMAGESAVEFLAG, false);
	}
	break;

	//Show set output data : minorId is the minor id of elements in Simulation::saveDataList (major id there is always 0), auxId is the number of the interactive object in the list as it appears in the console, textId is the configured output data.
	//Note this entry must always represent the entry in Simulation::saveDataList with the index in auxId.
	case IOI_OUTDATA:
	{
		//parameters from iop
		int outDataId = iop.minorId;
		int io_index = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELDATA, io_index);

		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimspaces(pTO->GetText());

			//After trimming all spaces has the form (last 2 are optional):
			//index.dataName<meshname>(rectangle)

			std::string index, dataName;

			std::string meshName = remove_last_contained_substring(to_text, "<", ">");
			std::string rectangle = remove_last_contained_substring(to_text, "(", ")");
			replaceall(rectangle, ",", " ");
			replaceall(rectangle, ";", " ");

			auto fields = split(to_text, ".");
			if (fields.size() == 2) {

				index = fields[0];
				dataName = fields[1];

				if (!meshName.length() && !rectangle.length()) sendCommand_verbose(CMD_EDITDATA, index, dataName);
				else if (!rectangle.length()) sendCommand_verbose(CMD_EDITDATA, index, dataName, meshName);
				else sendCommand_verbose(CMD_EDITDATA, index, dataName, meshName, rectangle);
			}
		}

		else if (actionCode == AC_MOUSELEFTDOWN) actionOutcome = AO_STARTINTERACTION;

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_OUTDATA
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == IOI_OUTDATA) {

			if (iop.interactingObjectId.minor != outDataId) {

				interaction(saveDataList, iop.interactingObjectId.minor, outDataId);
				
				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}

		else if (actionCode == AC_INTERACTOBJECTWITHWINDOW) {

			if (iop.interactingObjectId == INT2(WIN_DATABOX, 0)) {

				if (saveDataList.is_id_set(INT2(0, outDataId))) {

					std::string dataHandle = dataDescriptor.get_key_from_ID(saveDataList[INT2(0, outDataId)].datumId);
					sendCommand_verbose(CMD_ADDPINNEDDATA, dataHandle, saveDataList[io_index].meshName, saveDataList[io_index].rectangle);
				}

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
		std::string stageHandle = iop.textId;

		if (actionCode == AC_DOUBLECLICK) sendCommand_verbose(CMD_ADDSTAGE, stageHandle);
	}
	break;

	//Shows a stage added to the simulation schedule : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the configured stage text
	//Note this entry must always represent the entry in Simulation::simStages with the index in auxId.
	case IOI_SETSTAGE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		int io_index = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELSTAGE, io_index);

		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimspaces(pTO->GetText());

			//After trimming all spaces has the form (the last one is optional):
			//index.stageName<meshname>

			std::string index, stageName;
			std::string meshName = remove_last_contained_substring(to_text, "<", ">");

			auto fields = split(to_text, ".");
			if (fields.size() == 2) {

				index = fields[0];
				stageName = fields[1];

				sendCommand_verbose(CMD_EDITSTAGE, index, stageName, meshName);
			}
		}

		else if (actionCode == AC_MOUSELEFTDOWN) actionOutcome = AO_STARTINTERACTION;

		//moveable pop-up window is trying to interact with this object : interact them if the interacting object is also a IOI_SETSTAGE
		else if (actionCode == AC_INTERACTOBJECTS && iop.interactingObjectId.major == IOI_SETSTAGE) {

			if (iop.interactingObjectId.minor != stageId_minor) {

				interaction(simStages, iop.interactingObjectId.minor, stageId_minor);

				//call for interaction to end as purpose achieved
				actionOutcome = AO_ENDINTERACTION;
			}
		}
	}
	break;

	//Shows the value to set for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the value as a std::string
	case IOI_SETSTAGEVALUE:
	{
		//parameters from iop
		int io_index = iop.auxId;
		std::string stageValueText = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_EDITSTAGEVALUE, io_index, to_text);
		}
	}
	break;

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a std::string
	case IOI_STAGESTOPCONDITION:
	{
		//parameters from iop
		int io_index = iop.auxId;
		std::string stopConditionText = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ":", ",", ";");

			sendCommand_verbose(CMD_EDITSTAGESTOP, io_index, to_text);
		}
	}
	break;

	//Shows the saving condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the DSAVE_ value for this data save type, textId is the save type and value as a std::string
	case IOI_DSAVETYPE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		std::string saveConditionText = iop.textId;

		int io_index = simStages.get_index_from_id(INT2(0, stageId_minor));

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_EDITDATASAVE, io_index, trim(saveConditionText, ":", ",", ";"));

		else if (actionCode == AC_DOUBLECLICK) {

			//on double-click make popup edit box to edit the currently displayed value
			if (simStages[INT2(0, stageId_minor)].IsdSaveValueSet()) { actionOutcome = AO_STARTPOPUPEDITBOX; }
		}

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ":", ",", ";");

			sendCommand_verbose(CMD_EDITDATASAVE, io_index, to_text);
		}
	}
	break;

	//Shows a stop condition, used to apply the same condition to all simulation stages : minorId is the STOP_ value, textId is the stop type handle
	case IOI_STAGESTOPCONDITIONALL:
	{
		//parameters from iop
		std::string stopTypeHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_EDITSTAGESTOP, -1, stopTypeHandle);

		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITSTAGESTOP, -1, stopTypeHandle);
	}
	break;

	//Shows a data save condition, used to apply the same condition to all simulation stages : minorId is the DSAVE_ value, textId is the save type handle
	case IOI_DSAVETYPEALL:
	{
		//parameters from iop
		std::string saveTypeHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_EDITDATASAVE, -1, saveTypeHandle);

		else if (actionCode == AC_DOUBLECLICK) setConsoleEntry(CMD_EDITDATASAVE, -1, saveTypeHandle);
	}
	break;

	//Shows parameter and value for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter handle and value
	case IOI_MESHPARAM:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		std::string paramText = iop.textId;

		std::string meshName = SMesh.key_from_meshId(meshId);

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();

			//meshName of holding mesh
			std::string meshName = SMesh.key_from_meshId(meshId);

			//the text should be of the form "parameter_handle: value_with_units"; extract value_with_units if possible
			std::vector<std::string> fields = split(to_text, { ": " });
			if (fields.size() > 1) {
				
				sendCommand_verbose(CMD_SETPARAM, meshName, SMesh[meshName]->get_meshparam_handle(paramId), fields[1]);
			}
		}
	}
	break;

	//Shows parameter temperature dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter temperature dependence setting
	case IOI_MESHPARAMTEMP:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		std::string paramText = iop.textId;

		std::string meshName = SMesh.key_from_meshId(meshId);

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { 
			
			//only allow editing of value if this is a text equation or an  array (i.e. if it has a temperature dependence set)
			if (SMesh[meshName]->is_paramtemp_set(paramId) || SMesh[meshName]->is_paramtempequation_set(paramId)) actionOutcome = AO_STARTPOPUPEDITBOX;
		}

		//popup edit box has returned some text - try to set value from it (only applicable if the parameter has a equation temperature dependence
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			auto fields = split(trimendspaces(pTO->GetText()), ":");

			if (fields.size() == 2) {

				//set equation
				if (SMesh[meshName]->is_paramtempequation_set(paramId)) sendCommand_verbose(CMD_SETPARAMTEMPEQUATION, meshName, SMesh[meshName]->get_meshparam_handle(paramId), fields[1]);
				//set array from file
				else sendCommand_verbose(CMD_SETPARAMTEMPARRAY, meshName, SMesh[meshName]->get_meshparam_handle(paramId), trimendspaces(fields[1]));
			}
		}
		
		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_CLEARPARAMSTEMP, meshName, SMesh[meshName]->get_meshparam_handle(paramId));
		else if (actionCode == AC_DROPINTERACTOBJECTS) {

			if (iop.interactingObjectId.major == IOI_MESHPARAMTEMPTYPE) {
				
				//when IOI_MESHPARAMTEMPTYPE called for AO_STARTINTERACTION then iop.interactingObjectId.minor became the minorId of IOI_MESHPARAMTEMPTYPE, i.e. the MATPTDEP_ enum value
				switch ((MATPTDEP_)iop.interactingObjectId.minor) {

				case MATPTDEP_NONE:
					sendCommand_verbose(CMD_CLEARPARAMSTEMP, meshName, SMesh[meshName]->get_meshparam_handle(paramId));
					break;

				case MATPTDEP_ARRAY:
					sendCommand_verbose(CMD_SETPARAMTEMPARRAY, meshName, SMesh[meshName]->get_meshparam_handle(paramId));
					break;

				case MATPTDEP_EQUATION:
					sendCommand_verbose(CMD_SETPARAMTEMPEQUATION, meshName, SMesh[meshName]->get_meshparam_handle(paramId), "1");
					break;
				}
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
		std::string paramText = iop.textId;

		std::string meshName = SMesh.key_from_meshId(meshId);
		std::string paramName = SMesh[meshName]->get_meshparam_handle(paramId);

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_SETDISPLAYEDPARAMSVAR, meshName, paramName);
		
		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {
			
			//the actual text returned by the popup edit box
			std::string to_text;

			if (SMesh[meshName]->is_paramvarequation_set(paramId)) {

				//don't trim "," as these are used to specify vector equations
				to_text = trim(trimendspaces(pTO->GetText()), ":", ";");
			}
			else {

				//the actual text returned by the popup edit box
				to_text = trim(trimendspaces(pTO->GetText()), ":", ",", ";");
			}

			sendCommand_verbose(CMD_SETPARAMVAR, meshName, to_text);
		}

		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_CLEARPARAMSVAR, meshName, paramName);
		else if (actionCode == AC_DROPINTERACTOBJECTS) {

			if (iop.interactingObjectId.major == IOI_MESHPARAMVARGENERATOR) {

				//when IOI_MESHPARAMVARGENERATOR called for AO_STARTINTERACTION then iop.interactingObjectId.minor became the minorId of IOI_MESHPARAMVARGENERATOR, i.e. the MATPVAR_ enum value
				sendCommand_verbose(CMD_SETPARAMVAR, meshName, paramName, trim(vargenerator_descriptor.get_key_from_ID(iop.interactingObjectId.minor), ";", ","));
			}
		}
	}
	break;

	//Shows a possible temperature dependence : minorId is the type (entry from MATPTDEP_ enum)
	case IOI_MESHPARAMTEMPTYPE:
	{
		//parameters from iop
		MATPTDEP_ type = (MATPTDEP_)iop.minorId;

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
		std::string displayHandle = iop.textId;

		std::string meshName = SMesh.key_from_meshId(meshId);

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DISPLAY, displayHandle, meshName);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DISPLAYBACKGROUND, displayHandle, meshName);
	}
	break;

	//Shows dual mesh display transparency values : textId is the DBL2 value as a std::string
	case IOI_MESHDISPLAYTRANSPARENCY:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::vector<std::string> fields = split(trimendspaces(pTO->GetText()), ", ");

			if (fields.size() == 2) {

				sendCommand_verbose(CMD_DISPLAYTRANSPARENCY, fields[0], fields[1]);
			}
		}
	}
	break;

	//Shows mesh display threshold values : textId is the DBL2 value as a std::string
	case IOI_MESHDISPLAYTHRESHOLDS:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		if (actionCode == AC_MOUSERIGHTDOWN) { sendCommand_verbose(CMD_DISPLAYTHRESHOLDS, 0, 0); }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::vector<std::string> fields = split(trimendspaces(pTO->GetText()), ", ");

			if (fields.size() == 2) {

				sendCommand_verbose(CMD_DISPLAYTHRESHOLDS, fields[0], fields[1]);
			}
		}
	}
	break;

	//Shows mesh display threshold trigger type : auxId is the trigger option
	case IOI_MESHDISPLAYTHRESHOLDTRIG:
	{
		//parameters from iop
		int option = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN) {

			do { option = (option + 1) % VEC3REP_NUMOPTIONS;
			} while (option != (int)VEC3REP_X && option != (int)VEC3REP_Y && option != (int)VEC3REP_Z && option != (int)VEC3REP_MAGNITUDE);

			sendCommand_verbose(CMD_DISPLAYTHRESHOLDTRIGGER, option);
		}
		if (actionCode == AC_MOUSERIGHTDOWN) {

			do {
				option = (option + VEC3REP_NUMOPTIONS - 1) % VEC3REP_NUMOPTIONS;
			} while (option != (int)VEC3REP_X && option != (int)VEC3REP_Y && option != (int)VEC3REP_Z && option != (int)VEC3REP_MAGNITUDE);

			sendCommand_verbose(CMD_DISPLAYTHRESHOLDTRIGGER, option);
		}
	}
	break;

	//Shows super-mesh display option : minorId is the MESHDISPLAY_ value, textId is the MESHDISPLAY_ handle
	case IOI_SMESHDISPLAY:
	{
		//parameters from iop
		MESHDISPLAY_ displayOption = (MESHDISPLAY_)iop.minorId;
		std::string displayHandle = iop.textId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DISPLAY, displayHandle, SMesh.superMeshHandle);
	}
	break;

	//Shows mesh vectorial quantity display option : minorId is the unique mesh id number, auxId is the display option
	case IOI_MESHVECREP:
	{
		//parameters from iop
		int meshId = iop.minorId;
		int option = iop.auxId;

		std::string meshName = SMesh.key_from_meshId(meshId);

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_VECREP, meshName, (option + 1) % VEC3REP_NUMOPTIONS);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_VECREP, meshName, (option + VEC3REP_NUMOPTIONS - 1) % VEC3REP_NUMOPTIONS);
	}
	break;

	//Shows supermesh vectorial quantity display option : auxId is the display option
	case IOI_SMESHVECREP:
	{
		//parameters from iop
		int option = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_VECREP, SMesh.superMeshHandle, (option + 1) % VEC3REP_NUMOPTIONS);
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_VECREP, SMesh.superMeshHandle, (option + VEC3REP_NUMOPTIONS - 1) % VEC3REP_NUMOPTIONS);
	}
	break;

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	case IOI_MOVINGMESH:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool moving_mesh = iop.auxId;
		std::string meshName = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_MOVINGMESH, to_text);
		}

		else if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_MOVINGMESH, false);
	}
	break;

	//Shows movingmesh symmetry : auxId is the asymmetry status (1: asymmetric, 0: symmetric)
	case IOI_MOVINGMESHASYM:
	{
		//parameters from iop
		bool asymmetric = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN || actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_MOVINGMESHASYM, !asymmetric);
	}
	break;

	//Shows movingmesh threshold : textId is the threshold value as a std::string
	case IOI_MOVINGMESHTHRESH:
	{
		//parameters from iop
		std::string threshold_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_MOVINGMESHTHRESH, to_text);
		}
	}
	break;

	//Shows electrode box. minorId is the minor Id in STransport::electrode_boxes, auxId is the number of the interactive object in the list (electrode index), textId is the electrode rect as a std::string
	case IOI_ELECTRODERECT:
	{
		//parameters from iop
		int electrodeId_minor = iop.minorId;
		int io_index = iop.auxId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());
			replaceall(to_text, ", ", " ");
			replaceall(to_text, "; ", " ");

			sendCommand_verbose(CMD_SETELECTRODERECT, io_index, to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELELECTRODE, io_index);
	}
	break;

	//Shows electrode potential. minorId is the electrode index, textId is potential value as a std::string
	case IOI_ELECTRODEPOTENTIAL:
	{
		//parameters from iop
		int el_index = iop.minorId;
		std::string potential_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SETELECTRODEPOTENTIAL, el_index, to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELELECTRODE, el_index);
	}
	break;

	//Shows electrode ground setting. minorId is the electrode index, auxId is the setting (0 : not ground, 1 : ground)
	case IOI_ELECTRODEGROUND:
	{
		//parameters from iop
		int el_index = iop.minorId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DESIGNATEGROUND, el_index);
	}
	break;

	//Shows constant current source setting. auxId is the setting.
	case IOI_CONSTANTCURRENTSOURCE:
	{
		//parameters from iop
		bool is_constant_current = (bool)iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) {

			if (is_constant_current) sendCommand_verbose(CMD_SETPOTENTIAL, 0);
			else sendCommand_verbose(CMD_SETCURRENT, 0);
		}
	}
	break;

	//Shows transport solver convergence error. textId is the convergence error value.
	case IOI_TSOLVERCONVERROR:
	{
		//parameters from iop
		double conv_error = ToNum(iop.textId);

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_TSOLVERCONFIG, to_text, SMesh.CallModuleMethod(&STransport::GetConvergenceTimeout));
		}
	}
	break;

	//Shows transport solver timeout iterations. auxId is the timeout value.
	case IOI_TSOLVERTIMEOUT:
	{
		//parameters from iop
		int timeout = iop.auxId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_TSOLVERCONFIG, SMesh.CallModuleMethod(&STransport::GetConvergenceError), to_text);
		}
	}
	break;

	//Shows spin transport solver convergence error. textId is the convergence error value.
	case IOI_SSOLVERCONVERROR:
	{
		//parameters from iop
		double conv_error = ToNum(iop.textId);

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SSOLVERCONFIG, to_text, SMesh.CallModuleMethod(&STransport::GetSConvergenceTimeout));
		}
	}
	break;

	//Shows spin transport solver timeout iterations. auxId is the timeout value.
	case IOI_SSOLVERTIMEOUT:
	{
		//parameters from iop
		int timeout = iop.auxId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SSOLVERCONFIG, SMesh.CallModuleMethod(&STransport::GetSConvergenceError), to_text);
		}
	}
	break;

	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a std::string. (DBL2 since we need different damping values for V and S solvers)
	case IOI_SORDAMPING:
	{
		//parameters from iop
		std::string SOR_damping = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());
			replaceall(to_text, ", ", " ");
			replaceall(to_text, "; ", " ");

			sendCommand_verbose(CMD_SETSORDAMPING, to_text);
		}
	}
	break;

	//Shows tmr type setting. minorId is the unique mesh id number, auxId is the value.
	case IOI_TMRTYPE:
	{
		int meshId = iop.minorId;
		int TMR_type = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_TMRTYPE, SMesh.key_from_meshId(meshId), (TMR_type + 1) % (int)TMR_NUMOPTIONS);
		else if (actionCode == AC_MOUSERIGHTDOWN) {

			if (TMR_type > 0) TMR_type--;
			else TMR_type = (int)TMR_NUMOPTIONS - 1;

			sendCommand_verbose(CMD_TMRTYPE, SMesh.key_from_meshId(meshId), TMR_type);
		}
	}
	break;

	//Static transport solver state. auxId is the value (0/1)
	case IOI_STATICTRANSPORT:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_STATICTRANSPORTSOLVER, !status);
	}
	break;

	//Disabled transport solver state. auxId is the value (0/1)
	case IOI_DISABLEDTRANSPORT:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DISABLETRANSPORTSOLVER, !status);
	}
	break;

	//Shows mesh temperature. minorId is the unique mesh id number, textId is the temperature value
	case IOI_BASETEMPERATURE:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_TEMPERATURE, to_text, SMesh.key_from_meshId(meshId));
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_TEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows ambient temperature for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the temperature value
	case IOI_AMBIENT_TEMPERATURE:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;
		bool enabled = iop.auxId;

		if (enabled) {

			//on double-click make popup edit box to edit the currently displayed value
			if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

			//popup edit box has returned some text - try to set value from it
			else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

				//the actual text returned by the popup edit box
				std::string to_text = trimendspaces(pTO->GetText());

				sendCommand_verbose(CMD_AMBIENTTEMPERATURE, to_text, SMesh.key_from_meshId(meshId));
			}

			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_AMBIENTTEMPERATURE, temp_string, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows alpha value (W/m^2K) for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the value
	case IOI_ROBIN_ALPHA:
	{
		int meshId = iop.minorId;
		std::string alpha_string = iop.textId;
		bool enabled = iop.auxId;

		if (enabled) {

			//on double-click make popup edit box to edit the currently displayed value
			if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

			//popup edit box has returned some text - try to set value from it
			else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

				//the actual text returned by the popup edit box
				std::string to_text = trimendspaces(pTO->GetText());

				sendCommand_verbose(CMD_ROBINALPHA, to_text, SMesh.key_from_meshId(meshId));
			}

			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_ROBINALPHA, alpha_string, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows temperature insulating side setting for heat equation. minorId is the unique mesh id number, auxId is the status (Heat module must be active) : -1 disabled (gray), 0 not insulating (green), 1 insulating (red), textId represents the side : "x", "-x", "y", "-y", "z", "-z"
	case IOI_INSULATINGSIDE:
	{
		int meshId = iop.minorId;
		std::string literal = iop.textId;
		int status = iop.auxId;

		if (status >= 0) {

			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_INSULATINGSIDES, literal, !(bool)status, SMesh.key_from_meshId(meshId));
			if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_INSULATINGSIDES, literal, !(bool)status, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows mesh Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	case IOI_CURIETEMP:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_CURIETEMPERATURE, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows indicative material Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	case IOI_CURIETEMPMATERIAL:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_CURIETEMPERATUREMATERIAL, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows atomic moment multiple of Bohr magneton. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the value
	case IOI_ATOMICMOMENT:
	{
		int meshId = iop.minorId;
		std::string amoment_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_ATOMICMOMENT, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows atomic moment multiple of Bohr magneton for AF meshes. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	case IOI_ATOMICMOMENT_AFM:
	{
		int meshId = iop.minorId;
		std::string amoment_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_ATOMICMOMENT, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows atomic moment multiple of Bohr magneton for AF meshes. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	case IOI_TAU:
	{
		int meshId = iop.minorId;
		std::string tau_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_TAU, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows temperature model type for mesh. minorId is the unique mesh id number, auxId is the model identifier (entry from TMTYPE_ enum)
	case IOI_TMODEL:
	{
		int meshId = iop.minorId;
		int tmodel = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_TMODEL, tmodel % ((int)TMTYPE_NUMMODELS - 1) + 1, SMesh.key_from_meshId(meshId));
		else if (actionCode == AC_MOUSERIGHTDOWN) {

			if (tmodel > 1) tmodel--;
			else tmodel = (int)TMTYPE_NUMMODELS - 1;

			sendCommand_verbose(CMD_TMODEL, tmodel, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows cuda enabled/disabled or n/a state. auxId is enabled (1)/disabled(0)/not available(-1) status.
	case IOI_CUDASTATE:
	{
		int status = iop.auxId;
		if (status >= 0) {

			if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_CUDA, !status);
		}
	}
	break;
	
	//Shows CUDA device information and state. minorId is the device number (from 1 up), auxId is enabled (1)/disabled(0)/not available(-1) status. 
	case IOI_CUDADEVICE:
	{
		int device = iop.minorId;
		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_SELCUDADEV, device);
	}
	break;

	//Shows scale_rects enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCALERECTSSTATUS:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_SCALEMESHRECTS, !status);
	}
	break;

	//Shows coupled_to_dipoles enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_COUPLEDTODIPOLESSTATUS:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_COUPLETODIPOLES, !status);
	}
	break;

	//Shows dipole velocity value. minorId is the unique mesh id number. textId is the value
	case IOI_DIPOLEVELOCITY:
	{
		int meshId = iop.minorId;
		std::string value_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_DIPOLEVELOCITY, SMesh.key_from_meshId(meshId), (DBL3)ToNum(to_text, "m/s"), SMesh[SMesh.key_from_meshId(meshId)]->Get_Dipole_Clipping());
		}
	}
	break;
	
	//Shows diagonal strain set equation. minorId is the unique mesh id number. textId is the equation. auxId is enabled(1)/disabled(0) status.
	case IOI_STRAINEQUATION:
	{
		int meshId = iop.minorId;
		std::string value_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//on right-click clear the equation
		if (actionCode == AC_MOUSERIGHTDOWN) { sendCommand_verbose(CMD_CLEARSTRAINEQUATIONS, SMesh.key_from_meshId(meshId)); }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_STRAINEQUATION, SMesh.key_from_meshId(meshId), to_text);
		}
	}
	break;

	//Shows shear strain set equation. minorId is the unique mesh id number. textId is the equation. auxId is enabled(1)/disabled(0) status.
	case IOI_SHEARSTRAINEQUATION:
	{
		int meshId = iop.minorId;
		std::string value_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//on right-click clear the equation
		if (actionCode == AC_MOUSERIGHTDOWN) { sendCommand_verbose(CMD_CLEARSTRAINEQUATIONS, SMesh.key_from_meshId(meshId)); }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SHEARSTRAINEQUATION, SMesh.key_from_meshId(meshId), to_text);
		}
	}
	break;

	//Shows fixed surface rectangle. minorId is the minor Id in SMElastic::fixed_u_surfaces, auxId is the number of the interactive object in the list (electrode index), textId is the surface rect as a std::string
	case IOI_SURFACEFIX:
	{
		//parameters from iop
		int surfaceId_minor = iop.minorId;
		int io_index = iop.auxId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());
			replaceall(to_text, ", ", " ");
			replaceall(to_text, "; ", " ");

			sendCommand_verbose(CMD_EDITSURFACEFIX, io_index, to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELSURFACEFIX, io_index);
	}
	break;

	//Shows stress surface rectangle. minorId is the minor Id in SMElastic::stress_surfaces_rect, auxId is the number of the interactive object in the list (electrode index), textId is the surface rect as a std::string
	case IOI_SURFACESTRESS:
	{
		//parameters from iop
		int surfaceId_minor = iop.minorId;
		int io_index = iop.auxId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());
			replaceall(to_text, ", ", " ");
			replaceall(to_text, "; ", " ");

			sendCommand_verbose(CMD_EDITSURFACESTRESS, io_index, to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELSURFACESTRESS, io_index);
	}
	break;

	//Shows stress surface equation. minorId is the index in SMElastic::stress_surfaces_equations, auxId is the number of the interactive object in the list (electrode index), textId is the equation
	case IOI_SURFACESTRESSEQ:
	{
		//parameters from iop
		int el_index = iop.minorId;
		std::string equation = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_EDITSURFACESTRESSEQUATION, el_index, to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELSURFACESTRESS, el_index);
	}
	break;

	//Shows dipole shift clipping value. minorId is the unique mesh id number. textId is the value
	case IOI_DIPOLESHIFTCLIP:
	{
		int meshId = iop.minorId;
		std::string value_string = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_DIPOLEVELOCITY, SMesh.key_from_meshId(meshId), SMesh[SMesh.key_from_meshId(meshId)]->Get_Dipole_Velocity(), (DBL3)ToNum(to_text, "m"));
		}
	}
	break;

	//Shows log_error enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_ERRORLOGSTATUS:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_ERRORLOG, !status);
	}
	break;

	//Shows start_check_updates enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_UPDATESTATUSCHECKSTARTUP:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_STARTUPUPDATECHECK, !status);
	}
	break;

	//Shows start_scriptserver enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCRIPTSERVERSTARTUP:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_STARTUPSCRIPTSERVER, !status);
	}
	break;

	//Shows number of threads. auxId is the value.
	case IOI_THREADS:
	{
		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_THREADS, 0);

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_THREADS, to_text);
		}
	}
	break;


	//Shows server port. auxId is the value.
	case IOI_SERVERPORT:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_SERVERPORT, to_text);
		}
	}
	break;

	//Shows server password. textId is the password.
	case IOI_SERVERPWD:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_SERVERPWD, to_text);
		}
	}
	break;

	//Shows server sleep time in ms. auxId is the value.
	case IOI_SERVERSLEEPMS:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_SERVERSLEEPMS, to_text);
		}
	}
	break;

	//Shows neighboring meshes exchange coupling setting for this mesh. minorId is the unique mesh id number, auxId is the status (1/0 : on/off, -1 : not available: must be ferromagnetic mesh)
	case IOI_MESHEXCHCOUPLING:
	{
		bool status = iop.auxId;
		int meshId = iop.minorId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_EXCHANGECOUPLEDMESHES, !status, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, auxId is enabled (1)/disabled(0) status. textId is the value
	case IOI_REFINEROUGHNESS:
	{
		int meshId = iop.minorId;
		std::string refine = iop.textId;
		bool status = iop.auxId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && status) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_REFINEROUGHNESS, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows status of multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_MULTICONV:
	{
		int status = iop.auxId;

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand_verbose(CMD_MULTICONV, !status);
	}
	break;

	//Shows status of gpu kernels demag initialization. auxId is the status (0 : Off, 1 : On)
	case IOI_GPUKERNELS:
	{
		int status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_GPUKERNELS, !status);
	}
	break;

	//Shows status of force 2D multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_2DMULTICONV:
	{
		int status = iop.auxId;

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand_verbose(CMD_2DMULTICONV, (status + 1) % 3);
	}
	break;

	//Shows status of use default n for multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_NCOMMONSTATUS:
	{
		int status = iop.auxId;

		if (status >= 0 && (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN)) sendCommand_verbose(CMD_NCOMMONSTATUS, !status);
	}
	break;

	//Shows n_common for multi-layered convolution. auxId is the status (-1 : N/A, otherwise available). textId is the value as a SZ3.
	case IOI_NCOMMON:
	{
		int status = iop.auxId;
		std::string n_common = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && status >= 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_NCOMMON, to_text);
		}
	}
	break;

	//Shows materials database in use. textId is the name of the database, including the path.
	case IOI_LOCALMDB:
	{
		//parameters from iop
		std::string mdbFile = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_MATERIALSDATABASE, to_text);
		}
	}
	break;

	//Shows relative error fail threshold for ode eval. textId is the value.
	case IOI_ODERELERRFAIL:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_ASTEPCTRL, to_text, SMesh.Get_AStepdTCtrl());
		}
	}
	break;

	//Shows dT increase factor. textId is the value.
	case IOI_ODEDTINCR:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			double y = SMesh.Get_AStepdTCtrl().y;
			double z = SMesh.Get_AStepdTCtrl().z;

			sendCommand_verbose(CMD_ASTEPCTRL, ToString(SMesh.Get_AStepRelErrCtrl()), to_text, ToString(y), ToString(z));
		}
	}
	break;

	//Shows minimum dT value. textId is the value.
	case IOI_ODEDTMIN:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			double x = SMesh.Get_AStepdTCtrl().x;
			double z = SMesh.Get_AStepdTCtrl().z;

			sendCommand_verbose(CMD_ASTEPCTRL, SMesh.Get_AStepRelErrCtrl(), ToString(x), to_text, ToString(z));
		}
	}
	break;

	//Shows maximum dT value. textId is the value.
	case IOI_ODEDTMAX:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			double x = SMesh.Get_AStepdTCtrl().x;
			double y = SMesh.Get_AStepdTCtrl().y;

			sendCommand_verbose(CMD_ASTEPCTRL, SMesh.Get_AStepRelErrCtrl(), ToString(x), ToString(y), to_text);
		}
	}
	break;

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh)
	case IOI_PBC_X:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "x 10");

		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK && images > 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "x", to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "x 0");
	}
	break;

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh)
	case IOI_PBC_Y:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "y 10");
		
		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK && images > 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "y", to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "y 0");
	}
	break;

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh)
	case IOI_PBC_Z:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "z 10");
		
		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK && images > 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "z", to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_PBC, SMesh.key_from_meshId(meshId), "z 0");
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc)
	case IOI_SPBC_X:
	{
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "x 10");

		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK && images > 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "x", to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "x 0");
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc)
	case IOI_SPBC_Y:
	{
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "y 10");

		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK && images > 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "y", to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "y 0");
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc)
	case IOI_SPBC_Z:
	{
		int images = iop.auxId;

		if (actionCode == AC_MOUSELEFTDOWN && images == 0) sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "z 10");
		
		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK && images > 0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "z", to_text);
		}

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_PBC, SMesh.superMeshHandle, "z 0");
	}
	break;

	//Shows individual shape control flag. auxId is the value (0/1)
	case IOI_INDIVIDUALSHAPE:
	{
		bool status = (bool)iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_INDIVIDUALMASKSHAPE, !status);
	}
	break;

	//Shows image cropping settings : textId has the DBL4 value as text
	case IOI_IMAGECROPPING:
	{
		std::string value_text = iop.textId;

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trim(trimendspaces(pTO->GetText()), ",");

			sendCommand_verbose(CMD_IMAGECROPPING, to_text);
		}
	}
	break;

	//Show user constant for text equations : minorId is the index in Simulation::userConstants, auxId is the number of the interactive object in the list as it appears in the console, textId is the constant name and value std::string 
	//Note this entry must always represent the entry in Simulation::userConstants with the index in auxId.
	case IOI_USERCONSTANT:
	{
		//parameters from iop
		int index = iop.minorId;
		int io_index = iop.auxId;

		std::string constant_name = userConstants.get_key_from_index(index);

		if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_DELEQUATIONCONSTANT, constant_name);

		//on double-click make popup edit box to edit the currently displayed value
		else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimspaces(pTO->GetText());

			//After trimming all spaces has the form (last 2 are optional):
			//constant_name:value

			auto fields = split(to_text, ":");
			if (fields.size() == 2) {

				std::string new_constant_name = fields[0];
				double value = ToNum(fields[1]);
				
				if (new_constant_name != constant_name) userConstants.change_key(constant_name, new_constant_name);
				sendCommand_verbose(CMD_EQUATIONCONSTANTS, new_constant_name, value);
			}
		}
	}
	break;

	//Show skypos diameter multiplier : minorId is the unique mesh id number, textId is the multiplier as a std::string
	case IOI_SKYPOSDMUL:
	{
		int meshId = iop.minorId;
		double multiplier = ToNum(iop.textId);

		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK && multiplier > 0.0) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SKYPOSDMUL, to_text, SMesh.key_from_meshId(meshId));
		}
	}
	break;

	//Shows dwpos fitting component. auxId is the value (-1, 0, 1, 2)
	case IOI_DWPOSCOMPONENT:
	{
		//component takes on 0, 1, 2, 3
		int component = iop.auxId + 1;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_DWPOS_COMPONENT, (component + 1) % 4 - 1);
	}
	break;

	//Shows Monte-Carlo computation type (serial/parallel) : minorId is the unique mesh id number, auxId is the status (0 : parallel, 1 : serial, -1 : N/A)
	case IOI_MCCOMPUTATION:
	{
		int meshId = iop.minorId;
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_MCSERIAL, !status, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows Monte-Carlo computefields state flag : auxId is the state (0: disabled, 1: enabled)
	case IOI_MCCOMPUTEFIELDS:
	{
		bool status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_MCCOMPUTEFIELDS, !status);
	}
	break;

	//Shows Monte-Carlo algorithm type : minorId is the unique mesh id number, auxId is the type (0 : classical, 1 : constrained, -1 : N/A), textId is the constrained DBL3 direction.
	case IOI_MCTYPE:
	{
		int meshId = iop.minorId;
		int type = iop.auxId;

		if (type == 0) {

			//Currently in classical mode. Click to toggle to constrained mode with default x axis direction.
			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_MCCONSTRAIN, DBL3(1.0, 0.0, 0.0), SMesh.key_from_meshId(meshId));
		}
		else if (type == 1) {

			//Currently in constrained mode. Right-click to toggle, or double-click to edit direction
			if (actionCode == AC_MOUSERIGHTDOWN) sendCommand_verbose(CMD_MCCONSTRAIN, DBL3(), SMesh.key_from_meshId(meshId));

			//on double-click make popup edit box to edit the currently displayed value
			else if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

			//popup edit box has returned some text - try to set value from it
			else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

				//the actual text returned by the popup edit box
				std::string to_text = trimendspaces(pTO->GetText());

				sendCommand_verbose(CMD_MCCONSTRAIN, to_text, SMesh.key_from_meshId(meshId));
			}
		}
	}
	break;

	//Shows Monte-Carlo disabled/enabled status : minorId is the unique mesh id number, auxId is the status (1 : disabled, 0 : enabled).
	case IOI_MCDISABLED:
	{
		int meshId = iop.minorId;
		int status = iop.auxId;

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_MCDISABLE, !status, SMesh.key_from_meshId(meshId));
	}
	break;

	//Shows shape rotation setting: textId is the value as text (DBL3)
	case IOI_SHAPEROT:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SHAPEMOD_ROT, to_text);
		}
	}
	break;

	//Shows shape repetition setting: textId is the value as text (INT3)
	case IOI_SHAPEREP:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SHAPEMOD_REP, to_text);
		}
	}
	break;

	//Shows shape displacement setting: textId is the value as text (DBL3)
	case IOI_SHAPEDSP:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			sendCommand_verbose(CMD_SHAPEMOD_DISP, to_text);
		}
	}
	break;

	//Shows shape method setting: textId is the value as text (method)
	case IOI_SHAPEMET:
	{
		std::string method = iop.textId;
		int method_id = MeshShape().set_method(method);

		if (actionCode == AC_MOUSERIGHTDOWN || actionCode == AC_MOUSELEFTDOWN) sendCommand_verbose(CMD_SHAPEMOD_METHOD, MeshShape(MSHAPEMETHOD_((method_id + 1) % (int)MSHAPEMETHOD_NUMMETHODS)).method_name());
	}
	break;

	//Shows display render detail level: textId is the value
	case IOI_DISPRENDER_DETAIL:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = pTO->GetText();

			sendCommand_verbose(CMD_DISPLAYDETAILLEVEL, trimspaces(to_text));
		}
	}
	break;

	//Shows display render threshold 1: textId is the value
	case IOI_DISPRENDER_THRESH1:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			INT3 thresholds = BD.GetRenderThresholds();
			thresholds.i = ToNum(to_text);

			sendCommand_verbose(CMD_DISPLAYRENDERTHRESH, ToString(thresholds));
		}
	}
	break;

	//Shows display render threshold 2: textId is the value
	case IOI_DISPRENDER_THRESH2:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			INT3 thresholds = BD.GetRenderThresholds();
			thresholds.j = ToNum(to_text);

			sendCommand_verbose(CMD_DISPLAYRENDERTHRESH, ToString(thresholds));
		}
	}
	break;

	//Shows display render threshold 3: textId is the value
	case IOI_DISPRENDER_THRESH3:
	{
		//on double-click make popup edit box to edit the currently displayed value
		if (actionCode == AC_DOUBLECLICK) { actionOutcome = AO_STARTPOPUPEDITBOX; }

		//popup edit box has returned some text - try to set value from it
		else if (actionCode == AC_POPUPEDITTEXTBOXRETURNEDTEXT) {

			//the actual text returned by the popup edit box
			std::string to_text = trimendspaces(pTO->GetText());

			INT3 thresholds = BD.GetRenderThresholds();
			thresholds.k = ToNum(to_text);

			sendCommand_verbose(CMD_DISPLAYRENDERTHRESH, ToString(thresholds));
		}
	}
	break;


	default:
		break;
	}

	return actionOutcome;
}

#endif