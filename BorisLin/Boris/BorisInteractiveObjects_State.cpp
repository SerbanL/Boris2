#include "stdafx.h"

//Interactive Objects : 
//
//handle user interactions using ConsoleActionHandler
//update their state depending on current program state using ConsoleInteractiveObjectState

//These are not needed in non-graphical mode
#include "CompileFlags.h"
#if GRAPHICS == 1

#include "Simulation.h"

InteractiveObjectStateChange Simulation::ConsoleInteractiveObjectState(InteractiveObjectProperties &iop, TextObject *pTO) 
{
	//!!!IMPORTANT!!!: Do not call for a Refresh in this method, as it is called during a Refresh() : causes infinite loop! 
	//Also, this method was called from within BorisDisplay (through a function pointer), which was thread-safe accessed so the std::mutex is now locked.

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
				//insert a new output data interactive object after this : simMethod_BuildListEntry is a Simulation method which takes an integer argument (the index for which to build the formatted text std::string) and returns a std::string
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
		std::string meshName = iop.textId;

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

				//mesh still exists, it's just the name that has changed (possibly)
				meshName = SMesh().get_key_from_index(meshIdx);
				iop.textId = meshName;
				pTO->set(" " + meshName + " ");
				
				//another possibility is the mesh order has been changed, which will also result in a mismatch in the actual mesh name and that stored in iop.textId
				//in this case also want to rebuild the current mesh entry line as different meshes can have different types of entry lines (e.g. afm and fm meshes have different modules, display options, etc.)
				iop.state = IOS_REPLACINGPARAGRAPH;
				stateChanged.textMessage = CALLFP(this, simMethod_Build_Mesh_ListLine)(meshIdx) + "\n";
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
		std::string meshName = iop.textId;

		//this is the index corresponding to the dataBoxList_idminor - on any mismatch just reconstruct the data box entry to correspond to the element with DataBox_index index in dataBoxList
		int index_in_list = dataBoxList.get_index_from_id(INT2(0, dataBoxList_idminor));

		if (DataBox_index > dataBoxList.last()) {

			//too many fields : get rid of excess fields.
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
			break;
		}

		std::string actualmeshName = dataBoxList[DataBox_index].meshName;

		//if displayed meshname doesn't match the actual mesh name, or if indexes don't match, update Label (the n-th entry in the data box should represent the n-th entry in dataBoxList).
		if ((actualmeshName.length() && meshName != actualmeshName) || index_in_list != DataBox_index) {

			//meshname is set but doesn't match displayes name: update it.
			std::string newObjectText;
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
		MOD_ moduleID = (MOD_)iop.minorId;
		int meshId = iop.auxId;

		//from unique mesh id number get index in pMesh (held in SMesh)
		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0 && SMesh[meshIdx]->IsModuleSet(moduleID)) {

			//the mesh is contained and the module is set : ON color
			pTO->SetBackgroundColor(ONCOLOR);
		}
		else if (meshIdx < 0) {

			//the mesh is not contained : must have been deleted - delete this object
			stateChanged = true;
			iop.state = IOS_DELETING;
		}
		else {

			//the mesh is contained but module not active : OFF color
			pTO->SetBackgroundColor(OFFCOLOR);

			//special treatment of MOD_DEMAG to show if excluded from multilayered demag convolution
			if (moduleID == MOD_DEMAG) {

				if (SMesh.IsSuperMeshModuleSet(MODS_SDEMAG)) {

					if (SMesh[meshIdx]->Get_Demag_Exclusion()) pTO->SetBackgroundColor(ALTONCOLOR);
				}
			}
		}
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

	//Module used for effective field display for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to, textId is the MOD_ value used for display
	case IOI_DISPLAYMODULE:
	{
		//parameters from iop
		MOD_ moduleID = (MOD_)iop.minorId;
		int meshId = iop.auxId;
		int moduleDisplayed = ToNum(iop.textId);

		//from unique mesh id number get index in pMesh (held in SMesh)
		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			//mesh is available
			if (moduleDisplayed != SMesh[meshIdx]->Get_Module_Heff_Display()) {

				//mismatch, so adjust color and store correct value
				iop.textId = ToString(SMesh[meshIdx]->Get_Module_Heff_Display());

				if (SMesh[meshIdx]->Get_Module_Heff_Display() == (int)moduleID) pTO->SetBackgroundColor(ALTONCOLOR);
				else pTO->SetBackgroundColor(OFFCOLOR);
			}
		}
		else if (meshIdx < 0) {

			//the mesh is not contained : must have been deleted - delete this object
			stateChanged = true;
			iop.state = IOS_DELETING;
		}
	}
	break;

	//Available/set ode : minorId is an entry from ODE_ (the equation)
	case IOI_ODE:
	{
		//parameters from iop
		ODE_ odeID = (ODE_)iop.minorId;

		ODE_ actual_odeID;
		SMesh.QueryODE(actual_odeID);

		if (actual_odeID != odeID) {

			pTO->SetBackgroundColor(OFFCOLOR);
		}
		else pTO->SetBackgroundColor(ONCOLOR);
	}
	break;

	//Available/set ode for atomistic meshes: minorId is an entry from ODE_ (the equation)
	case IOI_ATOMODE:
	{
		//parameters from iop
		ODE_ odeID = (ODE_)iop.minorId;

		ODE_ actual_odeID;
		SMesh.QueryAtomODE(actual_odeID);

		if (actual_odeID != odeID) {

			pTO->SetBackgroundColor(OFFCOLOR);
		}
		else pTO->SetBackgroundColor(ONCOLOR);
	}
	break;

	//Set ODE time step: textId is the value
	case IOI_ODEDT:
	{
		//parameters from iop
		std::string dT_string = iop.textId;

		if (dT_string != ToString(SMesh.GetTimeStep(), "s")) {

			iop.textId = ToString(SMesh.GetTimeStep(), "s");

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Set stochastic time-step: textId is the value
	case IOI_STOCHDT:
	{
		//parameters from iop
		std::string dT_string = iop.textId;

		if (dT_string != ToString(SMesh.GetStochTimeStep(), "s")) {

			iop.textId = ToString(SMesh.GetStochTimeStep(), "s");

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Link stochastic time-step to ODE dT flag : auxId is the value
	case IOI_LINKSTOCHDT:
	{
		//parameters from iop
		bool state = (bool)iop.auxId;

		if (state != SMesh.GetLink_dTstoch()) {

			iop.auxId = SMesh.GetLink_dTstoch();

			if (iop.auxId) {

				pTO->set(" On ");
				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				pTO->set(" Off ");
				pTO->SetBackgroundColor(OFFCOLOR);
			}
			
			stateChanged = true;
		}
	}
	break;

	//Set evaluation speedup time-step: textId is the value
	case IOI_SPEEDUPDT:
	{
		//parameters from iop
		std::string dT_string = iop.textId;

		if (dT_string != ToString(SMesh.GetSpeedupTimeStep(), "s")) {

			iop.textId = ToString(SMesh.GetSpeedupTimeStep(), "s");

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Link evaluation speedup time-step to ODE dT flag : auxId is the value
	case IOI_LINKSPEEDUPDT:
	{
		//parameters from iop
		bool state = (bool)iop.auxId;

		if (state != SMesh.GetLink_dTspeedup()) {

			iop.auxId = SMesh.GetLink_dTspeedup();

			if (iop.auxId) {

				pTO->set(" On ");
				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				pTO->set(" Off ");
				pTO->SetBackgroundColor(OFFCOLOR);
			}

			stateChanged = true;
		}
	}
	break;

	//Set heat equation time step: textId is the value
	case IOI_HEATDT:
	{
		//parameters from iop
		std::string dT_string = iop.textId;

		if (SMesh.IsSuperMeshModuleSet(MODS_SHEAT)) {

			std::string actualdT_string = ToString(SMesh.CallModuleMethod(&SHeat::get_heat_dT), "s");

			if (dT_string != actualdT_string) {

				iop.textId = actualdT_string;

				pTO->set(" " + iop.textId + " ");
				if (dT_string == "N/A") pTO->SetBackgroundColor(ONCOLOR);

				stateChanged = true;
			}
		}
		else {

			if (dT_string != "N/A") {

				iop.textId = "N/A";
				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				stateChanged = true;
			}
		}
	}
	break;

	//Set elastodynamics equation time step: textId is the value
	case IOI_ELDT:
	{
		//parameters from iop
		std::string dT_string = iop.textId;

		if (SMesh.IsSuperMeshModuleSet(MODS_SMELASTIC)) {

			std::string actualdT_string = ToString(SMesh.CallModuleMethod(&SMElastic::get_el_dT), "s");

			if (dT_string != actualdT_string) {

				iop.textId = actualdT_string;

				pTO->set(" " + iop.textId + " ");
				if (dT_string == "N/A") pTO->SetBackgroundColor(ONCOLOR);

				stateChanged = true;
			}
		}
		else {

			if (dT_string != "N/A") {

				iop.textId = "N/A";
				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				stateChanged = true;
			}
		}
	}
	break;

	//Set stochastic time-step: textId is the value
	case IOI_LINKELDT:
	{
		//parameters from iop
		bool state = (bool)iop.auxId;

		if (state != SMesh.CallModuleMethod(&SMElastic::get_linked_el_dT)) {

			iop.auxId = SMesh.CallModuleMethod(&SMElastic::get_linked_el_dT);

			if (iop.auxId) {

				pTO->set(" On ");
				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				pTO->set(" Off ");
				pTO->SetBackgroundColor(OFFCOLOR);
			}

			stateChanged = true;
		}
	}
	break;

	//Available/set evaluation method for ode : minorId is an entry from ODE_ as : micromagnetic equation value + 100 * atomistic equation value, auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	case IOI_ODE_EVAL:
	{
		ODE_ actual_odeID, actual_atom_odeID;
		EVAL_ actual_evalID;
		SMesh.QueryODE(actual_odeID, actual_evalID);
		SMesh.QueryAtomODE(actual_atom_odeID);

		//parameters from iop
		int odeCombinedID = iop.minorId;
		EVAL_ evalID = (EVAL_)iop.auxId;

		//check if set ode has changed - if it has we need to update the state of this console object to reflect the set ode properties (e.g. evaluation method might not be available for the set ode)
		if ((int)actual_odeID + 100 * (int)actual_atom_odeID != odeCombinedID) {

			//mismatch in set ODE : update
			iop.minorId = (int)actual_odeID + 100 * (int)actual_atom_odeID;

			//if the newly set ode doesn't have this evaluation method as available then display UNAVAILABLE color
			if (search_vector(odeAllowedEvals(actual_odeID), evalID) < 0 || search_vector(odeAllowedEvals(actual_atom_odeID), evalID) < 0) {

				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				iop.state = IOS_OFF;
			}
			else iop.state = IOS_ON;
		}

		//if evaluation method is available for the set ode then mark it as selected / unselected
		if (iop.state == IOS_ON) {

			if (actual_evalID != evalID) {

				pTO->SetBackgroundColor(OFFCOLOR);
			}
			else pTO->SetBackgroundColor(ONCOLOR);
		}
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMS:
	{
		int meshId = iop.minorId;
		int meshIdx = SMesh.contains_id(meshId);
		std::string meshName = iop.textId;
		
		if (meshIdx >= 0 && meshName != SMesh().get_key_from_index(meshIdx)) {

			//if meshName is mismatched then best to delete this object : mesh positions in list could have been changed and it's problematic to reconstruct the display here
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
		else {

			display_meshIO(&Simulation::Build_MeshParams_Line);
		}
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMSTEMP:
	{
		int meshId = iop.minorId;
		int meshIdx = SMesh.contains_id(meshId);
		std::string meshName = iop.textId;

		if (meshIdx >= 0 && meshName != SMesh().get_key_from_index(meshIdx)) {

			//if meshName is mismatched then best to delete this object : mesh positions in list could have been changed and it's problematic to reconstruct the display here
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
		else {

			display_meshIO(&Simulation::Build_MeshParamsTemp_Text);
		}
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORPARAMSVAR:
	{
		int meshId = iop.minorId;
		int meshIdx = SMesh.contains_id(meshId);
		std::string meshName = iop.textId;

		if (meshIdx >= 0 && meshName != SMesh().get_key_from_index(meshIdx)) {

			//if meshName is mismatched then best to delete this object : mesh positions in list could have been changed and it's problematic to reconstruct the display here
			stateChanged = true;
			iop.state = IOS_DELETINGPARAGRAPH;
		}
		else {

			display_meshIO(&Simulation::Build_MeshParamsVariation_Text);
		}
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORMODULES:
	{
		display_meshIO(&Simulation::Build_Modules_ListLine);
	}
	break;

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently)
	case IOI_MESH_FORDISPLAYMODULES:
	{
		display_meshIO(&Simulation::Build_DisplayModules_ListLine);
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

	case IOI_MESH_FORTMODEL:
	{
		display_meshIO(&Simulation::Build_TemperatureModel_ListLine);
	}
	break;

	case IOI_MESH_FORPBC:
	{
		display_meshIO(&Simulation::Build_PBC_ListLine);
	}
	break;

	case IOI_MESH_FOREXCHCOUPLING:
	{
		display_meshIO(&Simulation::Build_ExchangeCoupledMeshes_ListLine);
	}
	break;

	case IOI_MESH_FORSTOCHASTICITY:
	{
		display_meshIO(&Simulation::Build_Stochasticity_ListLine);
	}
	break;

	case IOI_MESH_FORELASTICITY:
	{
		display_meshIO(&Simulation::Build_Elastodynamics_Equations_ListLine);
	}
	break;

	case IOI_MESH_FORSPEEDUP:
	{
		display_meshIO(&Simulation::Build_Speedup_ListLine);
	}
	break;

	case IOI_MESH_FORSKYPOSDMUL:
	{
		display_meshIO(&Simulation::Build_skypos_dmul_ListLine);
	}
	break;

	case IOI_MESH_FORMC:
	{
		display_meshIO(&Simulation::Build_MCSettings_ListLine);
	}
	break;

	case IOI_MESH_FORDIPOLESHIFT:
	{
		display_meshIO(&Simulation::Build_DipoleShiftingAlgorithm_ListLine);
	}
	break;

	case IOI_MESH_FORTMR:
	{
		display_meshIO(&Simulation::Print_TMRType_ListLine);
	}
	break;

	//Shows mesh rectangle (units m) : minorId is the unique mesh id number, textId is the mesh rectangle
	case IOI_MESHRECTANGLE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		std::string rectValue = iop.textId;

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
		std::string rectValue = iop.textId;

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
		std::string rectValue = iop.textId;

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
		std::string cellsizeValue = iop.textId;

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
		std::string cellsizeValue = iop.textId;

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
		std::string cellsizeValue = iop.textId;

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

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHMCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsizeValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (enabled) {

				if (!SMesh[meshIdx]->MechComputation_Enabled()) {

					iop.textId = "N/A";
					iop.auxId = 0;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(OFFCOLOR);
					stateChanged = true;
				}
				else {
					//update mesh cellsize if not matching
					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshMCellsize();
					if (ToString(meshCellsize, "m") != cellsizeValue) {

						iop.textId = ToString(meshCellsize, "m");
						pTO->set(" " + iop.textId + " ");
						stateChanged = true;
					}
				}
			}
			else {

				if (SMesh[meshIdx]->MechComputation_Enabled()) {

					DBL3 meshCellsize = SMesh[meshIdx]->GetMeshMCellsize();
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

	//Shows stochastic cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	//This is applicable for micromagnetic meshes only
	case IOI_MESHSCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsizeValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (enabled || SMesh[meshIdx]->is_atomistic()) {

				if (!SMesh[meshIdx]->MComputation_Enabled() || SMesh[meshIdx]->is_atomistic()) {

					iop.textId = "N/A";
					iop.auxId = 0;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(OFFCOLOR);
					stateChanged = true;
				}
				else {
					//update mesh cellsize if not matching
					DBL3 meshCellsize = dynamic_cast<Mesh*>(SMesh[meshIdx])->GetMeshSCellsize();
					if (ToString(meshCellsize, "m") != cellsizeValue) {

						iop.textId = ToString(meshCellsize, "m");
						pTO->set(" " + iop.textId + " ");
						stateChanged = true;
					}
				}
			}
			else {

				if (SMesh[meshIdx]->MComputation_Enabled()) {

					DBL3 meshCellsize = dynamic_cast<Mesh*>(SMesh[meshIdx])->GetMeshSCellsize();
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

	//Shows link stochastic flag : minorId is the unique mesh id number, auxId is the value off (0), on (1), N/A (-1)
	//This is applicable for micromagnetic meshes only
	case IOI_LINKSTOCHASTIC:
	{
		//parameters from iop
		int meshId = iop.minorId;
		int status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (status >= 0 && !SMesh[meshIdx]->is_atomistic()) {

				if (status != (int)dynamic_cast<Mesh*>(SMesh[meshIdx])->GetLinkStochastic()) {

					iop.auxId = dynamic_cast<Mesh*>(SMesh[meshIdx])->GetLinkStochastic();

					if (iop.auxId == 0) {

						pTO->set(" Off ");
						pTO->SetBackgroundColor(OFFCOLOR);
					}
					else {

						pTO->set(" On ");
						pTO->SetBackgroundColor(ONCOLOR);
					}

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

	//Shows macrocell size (units m) for atomistic meshes: minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	case IOI_MESHDMCELLSIZE:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool enabled = (bool)iop.auxId;
		std::string cellsizeValue = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (enabled || !SMesh[meshIdx]->is_atomistic()) {

				if (!SMesh[meshIdx]->MComputation_Enabled() || !SMesh[meshIdx]->is_atomistic()) {

					iop.textId = "N/A";
					iop.auxId = 0;
					pTO->set(" " + iop.textId + " ");
					pTO->SetBackgroundColor(OFFCOLOR);
					stateChanged = true;
				}
				else {
					//update mesh cellsize if not matching
					DBL3 meshCellsize = dynamic_cast<Atom_Mesh*>(SMesh[meshIdx])->Get_Demag_Cellsize();
					if (ToString(meshCellsize, "m") != cellsizeValue) {

						iop.textId = ToString(meshCellsize, "m");
						pTO->set(" " + iop.textId + " ");
						stateChanged = true;
					}
				}
			}
			else {

				if (SMesh[meshIdx]->MComputation_Enabled()) {

					DBL3 meshCellsize = dynamic_cast<Atom_Mesh*>(SMesh[meshIdx])->Get_Demag_Cellsize();
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

	//Shows evaluation speedup type: auxId is the type value.
	case IOI_SPEEDUPMODE:
	{
		//parameters from iop
		int option = iop.auxId;

		if (option != SMesh.GetEvaluationSpeedup()) {

			iop.auxId = SMesh.GetEvaluationSpeedup();

			switch (iop.auxId) {

			case 0:
				pTO->set(" None ");
				break;
			case 1:
				pTO->set(" Step ");
				break;
			case 2:
				pTO->set(" Linear ");
				break;
			case 3:
				pTO->set(" Quadratic ");
				break;
			case 4:
				pTO->set(" Cubic ");
				break;
			case 5:
				pTO->set(" Quartic ");
				break;
			case 6:
				pTO->set(" Quintic ");
				break;
			}

			stateChanged = true;
		}
	}
	break;

	//Shows ferromagnetic super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	case IOI_FMSMESHCELLSIZE:
	{
		//parameters from iop
		std::string cellsizeValue = iop.textId;

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
		std::string cellsizeValue = iop.textId;

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
		std::string directory_fromio = iop.textId;

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
		std::string savedataFile_fromiop = iop.textId;

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
		std::string savedataFile_fromiop = iop.textId;

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
		std::string configuredOutData = iop.textId;

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
		std::string configuredSetStage = iop.textId;

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


	//Shows the value to set for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the value as a std::string
	case IOI_SETSTAGEVALUE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		std::string stageValueText = iop.textId;

		//this is the value as a std::string
		std::string actualValuestring = simStages[INT2(0, stageId_minor)].get_value_string();

		if (stageValueText != actualValuestring) {

			iop.textId = actualValuestring;
			pTO->set(" " + actualValuestring + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a std::string
	case IOI_STAGESTOPCONDITION:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		int io_index = iop.auxId;
		std::string stopConditionText = iop.textId;

		if (stopConditionText != Build_SetStages_StopConditionText(io_index)) {

			iop.textId = Build_SetStages_StopConditionText(io_index);
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows the saving condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the DSAVE_ value for this data save type, textId is the save type and value as a std::string
	case IOI_DSAVETYPE:
	{
		//parameters from iop
		int stageId_minor = iop.minorId;
		DSAVE_ dSaveType = (DSAVE_)iop.auxId;
		std::string saveConditionText = iop.textId;

		if (simStages.is_id_set(INT2(0, stageId_minor))) {

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
	}
	break;

	//Shows parameter and value for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter handle and value
	case IOI_MESHPARAM:
	{
		//parameters from iop
		PARAM_ paramId = (PARAM_)iop.minorId;
		int meshId = iop.auxId;
		
		std::string paramText = get_after_match(iop.textId, std::string("\t"));
		std::string meshName = get_before_match(iop.textId, std::string("\t"));

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0 && meshName == SMesh().get_key_from_index(meshIdx)) {

			if (paramText != Build_MeshParams_Text(meshIdx, paramId)) {

				std::string meshParam_text = Build_MeshParams_Text(meshIdx, paramId);
				iop.textId = SMesh().get_key_from_index(meshIdx) + std::string("\t") + meshParam_text;

				pTO->set(" " + meshParam_text + " ");
				stateChanged = true;
			}
		}
		else {

			//this mesh no longer exists, so delete all associated interactive object parameters
			//another possibility is the mesh name has changed; this could mean the mesh positions were swapped in the mesh list, and this is problematic in this case so best to delete this object from console.
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
		
		std::string paramText = get_after_match(iop.textId, std::string("\t"));
		std::string meshName = get_before_match(iop.textId, std::string("\t"));

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0 && meshName == SMesh().get_key_from_index(meshIdx)) {

			if (paramText != SMesh[meshIdx]->get_paraminfo_string(paramId)) {

				std::string meshParam_text = SMesh[meshIdx]->get_paraminfo_string(paramId);
				iop.textId = SMesh().get_key_from_index(meshIdx) + std::string("\t") + meshParam_text;

				pTO->set(" " + meshParam_text + " ");

				if (SMesh[meshIdx]->is_paramtemp_set(paramId)) pTO->SetBackgroundColor(ONCOLOR);
				else pTO->SetBackgroundColor(OFFCOLOR);

				stateChanged = true;
			}

			//red (offcolor) : no temperature dependence in effect, green (oncolor) : temperature dependence in effect
			//this is useful if a text equation has been enetered but it is wrong : no temperature dependence is in effect.
			if (SameColor(pTO->GetFormatSpecifier().bgrndColor, OFFCOLOR) && SMesh[meshIdx]->is_paramtemp_set(paramId)) {

				pTO->SetBackgroundColor(ONCOLOR);
			}
			else if (SameColor(pTO->GetFormatSpecifier().bgrndColor, ONCOLOR) && !SMesh[meshIdx]->is_paramtemp_set(paramId)) {

				pTO->SetBackgroundColor(OFFCOLOR);
			}
		}
		else {

			//this mesh no longer exists, so delete all associated interactive object parameters
			//another possibility is the mesh name has changed; this could mean the mesh positions were swapped in the mesh list, and this is problematic in this case so best to delete this object from console.
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
		
		std::string paramText = get_after_match(iop.textId, std::string("\t"));
		std::string meshName = get_before_match(iop.textId, std::string("\t"));

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0 && meshName == SMesh().get_key_from_index(meshIdx)) {

			if (paramText != SMesh[meshIdx]->get_paramvarinfo_string(paramId)) {

				std::string meshParam_text = SMesh[meshIdx]->get_paramvarinfo_string(paramId);
				iop.textId = SMesh().get_key_from_index(meshIdx) + std::string("\t") + meshParam_text;

				pTO->set(" " + meshParam_text + " ");

				if (SMesh[meshIdx]->is_paramvar_set(paramId)) pTO->SetBackgroundColor(ONCOLOR);
				else pTO->SetBackgroundColor(OFFCOLOR);

				stateChanged = true;
			}
		}
		else {

			//this mesh no longer exists, so delete all associated interactive object parameters
			//another possibility is the mesh name has changed; this could mean the mesh positions were swapped in the mesh list, and this is problematic in this case so best to delete this object from console.
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
		std::string displayHandle = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->GetDisplayedPhysicalQuantity() == displayOption) {

				//this display option enabled
				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				if (SMesh[meshIdx]->GetDisplayedBackgroundPhysicalQuantity() == displayOption && displayOption != MESHDISPLAY_NONE) {

					//this display option enabled
					pTO->SetBackgroundColor(ALTONCOLOR);
				}
				else {

					//this display option disabled
					pTO->SetBackgroundColor(OFFCOLOR);
				}
			}
		}
	}
	break;

	//Shows dual mesh display transparency values : textId is the DBL2 value as a std::string
	case IOI_MESHDISPLAYTRANSPARENCY:
	{
		//parameters from iop
		std::string value = iop.textId;

		if (value != ToString(displayTransparency)) {

			iop.textId = ToString(displayTransparency);
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows mesh display threshold values : textId is the DBL2 value as a std::string
	case IOI_MESHDISPLAYTHRESHOLDS:
	{
		//parameters from iop
		std::string value = iop.textId;

		if (value != ToString(displayThresholds)) {

			iop.textId = ToString(displayThresholds);
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows mesh display threshold trigger type : auxId is the trigger option
	case IOI_MESHDISPLAYTHRESHOLDTRIG:
	{
		//parameters from iop
		int option = iop.auxId;
		
		if (option != displayThresholdTrigger) {

			iop.auxId = displayThresholdTrigger;

			switch (iop.auxId) {

			case 1:
				pTO->set(" X ");
				break;
			case 2:
				pTO->set(" Y ");
				break;
			case 3:
				pTO->set(" Z ");
				break;
			case 5:
				pTO->set(" magnitude ");
				break;
			}

			stateChanged = true;
		}
	}
	break;

	//Shows super-mesh display option : minorId is the MESHDISPLAY_ value, textId is the MESHDISPLAY_ handle
	case IOI_SMESHDISPLAY:
	{
		//parameters from iop
		MESHDISPLAY_ displayOption = (MESHDISPLAY_)iop.minorId;
		std::string displayHandle = iop.textId;

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

	//Shows mesh vectorial quantity display option : minorId is the unique mesh id number, auxId is the display option
	case IOI_MESHVECREP:
	{
		//parameters from iop
		int meshId = iop.minorId;
		int option = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);
		if (meshIdx >= 0) {

			if (option != SMesh[meshIdx]->GetVEC3Rep()) {

				iop.auxId = SMesh[meshIdx]->GetVEC3Rep();

				switch (iop.auxId) {

				case 0:
					pTO->set(" full ");
					break;
				case 1:
					pTO->set(" X ");
					break;
				case 2:
					pTO->set(" Y ");
					break;
				case 3:
					pTO->set(" Z ");
					break;
				case 4:
					pTO->set(" direction ");
					break; 
				case 5:
					pTO->set(" magnitude ");
					break;
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows supermesh vectorial quantity display option : auxId is the display option
	case IOI_SMESHVECREP:
	{
		//parameters from iop
		int option = iop.auxId;

		if (option != SMesh.GetVEC3Rep(SMesh.superMeshHandle)) {

			iop.auxId = SMesh.GetVEC3Rep(SMesh.superMeshHandle);

			switch (iop.auxId) {

			case 0:
				pTO->set(" full ");
				break;
			case 1:
				pTO->set(" X ");
				break;
			case 2:
				pTO->set(" Y ");
				break;
			case 3:
				pTO->set(" Z ");
				break;
			case 4:
				pTO->set(" direction ");
				break;
			case 5:
				pTO->set(" magnitude ");
				break;
			}

			stateChanged = true;
		}
	}
	break;

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	case IOI_MOVINGMESH:
	{
		//parameters from iop
		int meshId = iop.minorId;
		bool moving_mesh = iop.auxId;
		std::string meshName = iop.textId;

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

	//Shows movingmesh threshold : textId is the threshold value as a std::string
	case IOI_MOVINGMESHTHRESH:
	{
		//parameters from iop
		std::string threshold_string = iop.textId;

		if (threshold_string != ToString(SMesh.MoveMeshThreshold())) {

			iop.textId = ToString(SMesh.MoveMeshThreshold());

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows electrode box. minorId is the minor Id in STransport::electrode_boxes, auxId is the number of the interactive object in the list (electrode index), textId is the electrode rect as a std::string
	case IOI_ELECTRODERECT:
	{
		//parameters from iop
		int electrodeId_minor = iop.minorId;
		int io_index = iop.auxId;
		std::string rect_string = iop.textId;

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

	//Shows electrode potential. minorId is the electrode index, textId is potential value as a std::string
	case IOI_ELECTRODEPOTENTIAL:
	{
		//parameters from iop
		int el_index = iop.minorId;
		std::string potential_string = iop.textId;

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

	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a std::string. (DBL2 since we need different damping values for V and S solvers)
	case IOI_SORDAMPING:
	{
		//parameters from iop
		std::string SOR_damping = iop.textId;

		if (SOR_damping != ToString(SMesh.CallModuleMethod(&STransport::GetSORDamping))) {

			iop.textId = ToString(SMesh.CallModuleMethod(&STransport::GetSORDamping));
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;


	//Shows tmr type setting. minorId is the unique mesh id number, auxId is the value.
	case IOI_TMRTYPE:
	{
		int meshId = iop.minorId;
		int TMR_type = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (TMR_type != SMesh[meshIdx]->GetTMRType()) {

				iop.auxId = SMesh[meshIdx]->GetTMRType();

				switch ((TMR_)iop.auxId) {

				case TMR_NONE:
					pTO->set(" N/A ");
					break;

				case TMR_COS:
					pTO->set(" cosine ");
					break;

				case TMR_SLONCZEWSKI:
					pTO->set(" Slonczewski ");
					break;
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Static transport solver state. auxId is the value (0/1)
	case IOI_STATICTRANSPORT:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (status != static_transport_solver) {

			iop.auxId = static_transport_solver;

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
	break;

	//Disabled transport solver state. auxId is the value (0/1)
	case IOI_DISABLEDTRANSPORT:
	{
		//parameters from iop
		bool status = iop.auxId;

		if (status != disabled_transport_solver) {

			iop.auxId = disabled_transport_solver;

			if (iop.auxId == 0) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" Enabled ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Disabled ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows mesh base temperature. minorId is the unique mesh id number, textId is the temperature value
	case IOI_BASETEMPERATURE:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;

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
		std::string temp_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->IsModuleSet(MOD_HEAT) != status) {

				iop.auxId = SMesh[meshIdx]->IsModuleSet(MOD_HEAT);

				if (iop.auxId) {

					pTO->SetBackgroundColor(ONCOLOR);
					iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetAmbientTemperature), "K");
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
			else if (status && ToString(SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetAmbientTemperature), "K") != temp_string) {

				iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetAmbientTemperature), "K");
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
		std::string alpha_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->IsModuleSet(MOD_HEAT) != status) {

				iop.auxId = SMesh[meshIdx]->IsModuleSet(MOD_HEAT);

				if (iop.auxId) {

					pTO->SetBackgroundColor(ONCOLOR);
					iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetAlphaBoundary), "W/m2K");
					pTO->set(" " + iop.textId + " ");
				}
				else {

					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
				}

				stateChanged = true;
			}
			else if (status && ToString(SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetAlphaBoundary), "W/m2K") != alpha_string) {

				iop.textId = ToString(SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetAlphaBoundary), "W/m2K");
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
		std::string literal = iop.textId;
		int status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->IsModuleSet(MOD_HEAT) != bool(status + 1)) {

				bool heat_set = SMesh[meshIdx]->IsModuleSet(MOD_HEAT);
				if (!heat_set) iop.auxId = -1;
				else iop.auxId = SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetInsulatingSide, literal);

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
			else if (status >= 0 && SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetInsulatingSide, literal) != (bool)status) {

				iop.auxId = SMesh[meshIdx]->CallModuleMethod(&HeatBase::GetInsulatingSide, literal);
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
	//This is applicable for micromagnetic meshes only
	case IOI_CURIETEMP:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0 && !SMesh[meshIdx]->is_atomistic()) {

			if (ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetCurieTemperature(), "K") != temp_string) {

				iop.textId = ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetCurieTemperature(), "K");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->Magnetism_Enabled()) != status) {

				iop.auxId = (SMesh[meshIdx]->Magnetism_Enabled());

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
	//This is applicable for micromagnetic meshes only
	case IOI_CURIETEMPMATERIAL:
	{
		int meshId = iop.minorId;
		std::string temp_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0 && !SMesh[meshIdx]->is_atomistic()) {

			if (ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetCurieTemperatureMaterial(), "K") != temp_string) {

				iop.textId = ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetCurieTemperatureMaterial(), "K");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->Magnetism_Enabled()) != status) {

				iop.auxId = (SMesh[meshIdx]->Magnetism_Enabled());

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
	//This is applicable for micromagnetic meshes only
	case IOI_ATOMICMOMENT:
	{
		int meshId = iop.minorId;
		std::string amoment_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0 && !SMesh[meshIdx]->is_atomistic()) {

			if (ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetAtomicMoment(), "uB") != amoment_string) {

				iop.textId = ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetAtomicMoment(), "uB");
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

	//Shows atomic moment multiple of Bohr magneton for AF meshes. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	//This is applicable for micromagnetic meshes only
	case IOI_ATOMICMOMENT_AFM:
	{
		int meshId = iop.minorId;
		std::string amoment_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0 && !SMesh[meshIdx]->is_atomistic()) {

			if (ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetAtomicMoment_AFM(), "uB") != amoment_string) {

				iop.textId = ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetAtomicMoment_AFM(), "uB");
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) != status) {

				iop.auxId = (SMesh[meshIdx]->GetMeshType() == MESH_ANTIFERROMAGNETIC);

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

	//Shows Tc tau couplings. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	//This is applicable for micromagnetic meshes only
	case IOI_TAU:
	{
		int meshId = iop.minorId;
		std::string amoment_string = iop.textId;
		bool status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0 && !SMesh[meshIdx]->is_atomistic()) {

			if (ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetTcCoupling()) != amoment_string) {

				iop.textId = ToString(dynamic_cast<Mesh*>(SMesh[meshIdx])->GetTcCoupling());
				pTO->set(" " + iop.textId + " ");

				stateChanged = true;
			}

			if ((SMesh[meshIdx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) != status) {

				iop.auxId = (SMesh[meshIdx]->GetMeshType() == MESH_ANTIFERROMAGNETIC);

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

	//Shows temperature model type for mesh. minorId is the unique mesh id number, auxId is the model identifier (entry from TMTYPE_ enum)
	case IOI_TMODEL:
	{
		int meshId = iop.minorId;
		int tmodel = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->CallModuleMethod(&HeatBase::Get_TMType) != tmodel) {

				iop.auxId = SMesh[meshIdx]->CallModuleMethod(&HeatBase::Get_TMType);

				switch (iop.auxId) {

				case TMTYPE_NONE:
					pTO->SetBackgroundColor(UNAVAILABLECOLOR);
					pTO->set(" N/A ");
					break;

				case TMTYPE_1TM:
					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" 1TM ");
					break;

				case TMTYPE_2TM:
					pTO->SetBackgroundColor(ONCOLOR);
					pTO->set(" 2TM ");
					break;
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

	//Shows CUDA device information and state. minorId is the device number (from 1 up), auxId is enabled (1)/disabled(0)/not available(-1) status. 
	case IOI_CUDADEVICE:
	{
		int device = iop.minorId;
		int status = iop.auxId;

		//make sure status color matches; if not available then the IOI will be grayed out when created and cannot be selected.
		if ((status == 1 && device != cudaDeviceSelect) || (status == 0 && device == cudaDeviceSelect)) {
			
			if (device == cudaDeviceSelect) {

				iop.auxId = 1;
				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				iop.auxId = 0;
				pTO->SetBackgroundColor(OFFCOLOR);
			}

			stateChanged = true;
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

	//Shows dipole velocity value. minorId is the unique mesh id number. textId is the value
	case IOI_DIPOLEVELOCITY:
	{
		//parameters from iop
		int meshId = iop.minorId;
		std::string value_string = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);

		//is there a value mismatch?
		if (meshIdx >= 0 && value_string != ToString(SMesh[meshIdx]->Get_Dipole_Velocity(), "m/s")) {

			iop.textId = ToString(SMesh[meshIdx]->Get_Dipole_Velocity(), "m/s");
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows diagonal strain set equation. minorId is the unique mesh id number. textId is the equation. auxId is enabled(1)/disabled(0) status.
	case IOI_STRAINEQUATION:
	{
		//parameters from iop
		int meshId = iop.minorId;

		int meshIdx = SMesh.contains_id(meshId);

		//is there a value mismatch?
		if (meshIdx >= 0) {

			std::string text_equation = SMesh[meshIdx]->Get_StrainEquation();
			if (text_equation.length() && iop.textId != text_equation) {

				iop.textId = text_equation;
				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(ONCOLOR);
				stateChanged = true;
			}
			else if (!text_equation.length() && iop.textId != "") {

				pTO->set(" None ");
				pTO->SetBackgroundColor(OFFCOLOR);
				stateChanged = true;
			}
		}
	}
	break;

	//Shows shear strain set equation. minorId is the unique mesh id number. textId is the equation. auxId is enabled(1)/disabled(0) status.
	case IOI_SHEARSTRAINEQUATION:
	{
		//parameters from iop
		int meshId = iop.minorId;

		int meshIdx = SMesh.contains_id(meshId);

		//is there a value mismatch?
		if (meshIdx >= 0) {

			std::string text_equation = SMesh[meshIdx]->Get_ShearStrainEquation();
			if (text_equation.length() && iop.textId != text_equation) {

				iop.textId = text_equation;
				pTO->set(" " + iop.textId + " ");
				pTO->SetBackgroundColor(ONCOLOR);
				stateChanged = true;
			}
			else if (!text_equation.length() && iop.textId != "") {

				pTO->set(" None ");
				pTO->SetBackgroundColor(OFFCOLOR);
				stateChanged = true;
			}
		}
	}
	break;

	//Shows fixed surface rectangle. minorId is the minor Id in SMElastic::fixed_u_surfaces, auxId is the number of the interactive object in the list (electrode index), textId is the surface rect as a std::string
	case IOI_SURFACEFIX:
	{
		//parameters from iop
		int surfaceId_minor = iop.minorId;
		int io_index = iop.auxId;
		std::string rect_string = iop.textId;

		//actual index in surfaces list (should normally be the same as io_index)
		int index_in_list = SMesh.CallModuleMethod(&SMElastic::Get_Fixed_Surface_Index, surfaceId_minor);
		int el_last_index = SMesh.CallModuleMethod(&SMElastic::Get_Num_Fixed_Surfaces) - 1;

		//if there's a mismatch between the object number and the actual index then updating is needed - update the entire object so that it corresponds to the entry at io_index.
		if (io_index <= el_last_index && index_in_list != io_index) {

			//because there are multiple objects on this line, all of them must be replaced. The caller must do this.
			iop.state = IOS_REPLACINGPARAGRAPH;
			stateChanged.textMessage = Build_FixedSurfaces_ListLine(io_index);
			stateChanged = true;
			break;
		}

		//this object is part of a list : make sure this list is updated
		updateList(io_index, el_last_index, &Simulation::Build_FixedSurfaces_ListLine);

		if (rect_string != ToString(SMesh.CallModuleMethod(&SMElastic::Get_Fixed_Surface, io_index), "m")) {

			iop.textId = ToString(SMesh.CallModuleMethod(&SMElastic::Get_Fixed_Surface, io_index), "m");
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows stress surface rectangle. minorId is the minor Id in SMElastic::stress_surfaces_rect, auxId is the number of the interactive object in the list (electrode index), textId is the surface rect as a std::string
	case IOI_SURFACESTRESS:
	{
		//parameters from iop
		int surfaceId_minor = iop.minorId;
		int io_index = iop.auxId;
		std::string rect_string = iop.textId;

		//actual index in surfaces list (should normally be the same as io_index)
		int index_in_list = SMesh.CallModuleMethod(&SMElastic::Get_Stress_Surface_Index, surfaceId_minor);
		int el_last_index = SMesh.CallModuleMethod(&SMElastic::Get_Num_Stress_Surfaces) - 1;

		//if there's a mismatch between the object number and the actual index then updating is needed - update the entire object so that it corresponds to the entry at io_index.
		if (io_index <= el_last_index && index_in_list != io_index) {

			//because there are multiple objects on this line, all of them must be replaced. The caller must do this.
			iop.state = IOS_REPLACINGPARAGRAPH;
			stateChanged.textMessage = Build_StressSurfaces_ListLine(io_index);
			stateChanged = true;
			break;
		}

		//this object is part of a list : make sure this list is updated
		updateList(io_index, el_last_index, &Simulation::Build_StressSurfaces_ListLine);

		if (rect_string != ToString(SMesh.CallModuleMethod(&SMElastic::Get_Stress_Surface, io_index), "m")) {

			iop.textId = ToString(SMesh.CallModuleMethod(&SMElastic::Get_Stress_Surface, io_index), "m");
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows stress surface equation. minorId is the index in SMElastic::stress_surfaces_equations, auxId is the number of the interactive object in the list (electrode index), textId is the equation
	case IOI_SURFACESTRESSEQ:
	{
		//parameters from iop
		int el_index = iop.minorId;
		std::string equation = iop.textId;

		if (SMesh.CallModuleMethod(&SMElastic::Get_Stress_Surface_Equation, el_index) != equation) {

			iop.textId = SMesh.CallModuleMethod(&SMElastic::Get_Stress_Surface_Equation, el_index);
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows dipole shift clipping value. minorId is the unique mesh id number. textId is the value
	case IOI_DIPOLESHIFTCLIP:
	{
		//parameters from iop
		int meshId = iop.minorId;
		std::string value_string = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);

		//is there a value mismatch?
		if (meshIdx >= 0 && value_string != ToString(SMesh[meshIdx]->Get_Dipole_Clipping(), "m")) {

			iop.textId = ToString(SMesh[meshIdx]->Get_Dipole_Clipping(), "m");
			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows log_errors enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_ERRORLOGSTATUS:
	{
		bool status = iop.auxId;

		if (status != log_errors) {

			iop.auxId = log_errors;

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

	//Shows start_check_updates enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_UPDATESTATUSCHECKSTARTUP:
	{
		bool status = iop.auxId;

		if (status != start_check_updates) {

			iop.auxId = start_check_updates;

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

	//Shows start_scriptserver enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCRIPTSERVERSTARTUP:
	{
		bool status = iop.auxId;

		if (status != start_scriptserver) {

			iop.auxId = start_scriptserver;

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

	//Shows number of threads. auxId is the value.
	case IOI_THREADS:
	{
		int threads = iop.auxId;

		if (threads != OmpThreads) {

			iop.auxId = OmpThreads;
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows server port. auxId is the value.
	case IOI_SERVERPORT:
	{
		std::string port = ToString(iop.auxId);
		
		if (port != server_port) {

			iop.auxId = ToNum(server_port);
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows server password. textId is the password.
	case IOI_SERVERPWD:
	{
		std::string password = iop.textId;

		if (password != server_pwd) {

			iop.textId = server_pwd;

			pTO->set(" " + server_pwd + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows server sleep time in ms. auxId is the value.
	case IOI_SERVERSLEEPMS:
	{
		int sleepms = iop.auxId;

		if (sleepms != server_recv_sleep_ms) {

			iop.auxId = server_recv_sleep_ms;
			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows neighboring meshes exchange coupling setting for this mesh. minorId is the unique mesh id number, auxId is the status (1/0 : on/off, -1 : not available: must be ferromagnetic mesh)
	case IOI_MESHEXCHCOUPLING:
	{
		int meshId = iop.minorId;
		int status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (status >= 0 && status != (int)SMesh[meshIdx]->GetMeshExchangeCoupling()) {

				iop.auxId = SMesh[meshIdx]->GetMeshExchangeCoupling();

				if (iop.auxId > 0) {

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

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, auxId is enabled (1)/disabled(0) status. textId is the value
	case IOI_REFINEROUGHNESS:
	{
		int meshId = iop.minorId;
		std::string refine_string = iop.textId;
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

	//Shows status of gpu kernels demag initialization. auxId is the status (0 : Off, 1 : On)
	case IOI_GPUKERNELS:
	{
		int status = iop.auxId;

		if (status != (int)SMesh.Get_Kernel_Initialize_on_GPU()) {

			iop.auxId = (int)SMesh.Get_Kernel_Initialize_on_GPU();

			if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" On ");
			}
			else if (iop.auxId == 0) {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Off ");
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

			if (iop.auxId == 2) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" 2D layered ");

			}
			else if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" 2D meshes ");
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
		std::string common_n = iop.textId;

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
		std::string mdbFile = iop.textId;

		//update name if not matching
		if (mdbFile != mdb.GetDataBaseName()) {

			iop.textId = mdb.GetDataBaseName();
			pTO->set(" " + mdb.GetDataBaseName() + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows relative error fail threshold for ode eval. textId is the value.
	case IOI_ODERELERRFAIL:
	{
		//parameters from iop
		double value = ToNum(iop.textId);

		//update value if not matching
		if (value != SMesh.Get_AStepRelErrCtrl()) {

			iop.textId = ToString(SMesh.Get_AStepRelErrCtrl());
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows dT increase factor. textId is the value.
	case IOI_ODEDTINCR:
	{
		//parameters from iop
		double value = ToNum(iop.textId, "s");

		//update value if not matching
		if (value != SMesh.Get_AStepdTCtrl().i) {

			iop.textId = ToString(SMesh.Get_AStepdTCtrl().i);
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows minimum dT value. textId is the value.
	case IOI_ODEDTMIN:
	{
		//parameters from iop
		double value = ToNum(iop.textId, "s");

		//update value if not matching
		if (value != SMesh.Get_AStepdTCtrl().j) {

			iop.textId = ToString(SMesh.Get_AStepdTCtrl().j);
			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows maximum dT value. textId is the value.
	case IOI_ODEDTMAX:
	{
		//parameters from iop
		double value = ToNum(iop.textId, "s");

		//update value if not matching
		if (value != SMesh.Get_AStepdTCtrl().k) {

			iop.textId = ToString(SMesh.Get_AStepdTCtrl().k);
			pTO->set(" " + iop.textId + " ");

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

	//Shows PBC setting. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh)
	case IOI_PBC_X:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->Is_Demag_Enabled() && SMesh[meshIdx]->CallModuleMethod(&DemagBase::Get_PBC).x != images) {

				iop.auxId = SMesh[meshIdx]->CallModuleMethod(&DemagBase::Get_PBC).x;

				if (iop.auxId != 0) {

					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
				}

				pTO->set(" " + ToString(iop.auxId) + " ");
				stateChanged = true;
			}
			else if (!SMesh[meshIdx]->Is_Demag_Enabled() && SMesh[meshIdx]->Get_Magnetic_PBC().x != images) {

				iop.auxId = SMesh[meshIdx]->Get_Magnetic_PBC().x;

				if (iop.auxId != 0) {

					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
				}

				pTO->set(" " + ToString(iop.auxId) + " ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows PBC setting. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh)
	case IOI_PBC_Y:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->Is_Demag_Enabled() && SMesh[meshIdx]->CallModuleMethod(&DemagBase::Get_PBC).y != images) {

				iop.auxId = SMesh[meshIdx]->CallModuleMethod(&DemagBase::Get_PBC).y;

				if (iop.auxId != 0) {

					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
				}

				pTO->set(" " + ToString(iop.auxId) + " ");
				stateChanged = true;
			}
			else if (!SMesh[meshIdx]->Is_Demag_Enabled() && SMesh[meshIdx]->Get_Magnetic_PBC().y != images) {

				iop.auxId = SMesh[meshIdx]->Get_Magnetic_PBC().y;

				if (iop.auxId != 0) {

					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
				}

				pTO->set(" " + ToString(iop.auxId) + " ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows PBC setting. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh)
	case IOI_PBC_Z:
	{
		int meshId = iop.minorId;
		int images = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (SMesh[meshIdx]->Is_Demag_Enabled() && SMesh[meshIdx]->CallModuleMethod(&DemagBase::Get_PBC).z != images) {

				iop.auxId = SMesh[meshIdx]->CallModuleMethod(&DemagBase::Get_PBC).z;

				if (iop.auxId != 0) {

					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
				}

				pTO->set(" " + ToString(iop.auxId) + " ");
				stateChanged = true;
			}
			else if (!SMesh[meshIdx]->Is_Demag_Enabled() && SMesh[meshIdx]->Get_Magnetic_PBC().z != images) {

				iop.auxId = SMesh[meshIdx]->Get_Magnetic_PBC().z;

				if (iop.auxId != 0) {

					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->SetBackgroundColor(OFFCOLOR);
				}

				pTO->set(" " + ToString(iop.auxId) + " ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc)
	case IOI_SPBC_X:
	{
		int images = iop.auxId;
		int status = iop.minorId;

		if (status && !SMesh.IsSuperMeshModuleSet(MODS_SDEMAG)) {

			pTO->SetBackgroundColor(UNAVAILABLECOLOR);
			pTO->set(" N/A ");
			stateChanged = true;
			iop.auxId = 0;
			iop.minorId = 0;
		}
		else if (SMesh.IsSuperMeshModuleSet(MODS_SDEMAG) && (dynamic_cast<SDemag*>(SMesh.GetSuperMeshModule(MODS_SDEMAG))->Get_PBC().x != images || !status)) {

			iop.minorId = 1;
			iop.auxId = SMesh.CallModuleMethod(&SDemag::Get_PBC).x;

			if (iop.auxId != 0) {

				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
			}

			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc)
	case IOI_SPBC_Y:
	{
		int images = iop.auxId;
		int status = iop.minorId;

		if (status && !SMesh.IsSuperMeshModuleSet(MODS_SDEMAG)) {

			pTO->SetBackgroundColor(UNAVAILABLECOLOR);
			pTO->set(" N/A ");
			stateChanged = true;
			iop.auxId = 0;
			iop.minorId = 0;
		}
		else if (SMesh.IsSuperMeshModuleSet(MODS_SDEMAG) && (dynamic_cast<SDemag*>(SMesh.GetSuperMeshModule(MODS_SDEMAG))->Get_PBC().y != images || !status)) {

			iop.minorId = 1;
			iop.auxId = SMesh.CallModuleMethod(&SDemag::Get_PBC).y;

			if (iop.auxId != 0) {

				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
			}

			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc), minorId is enabled (1)/disabled (0) status
	case IOI_SPBC_Z:
	{
		int images = iop.auxId;
		int status = iop.minorId;

		if (status && !SMesh.IsSuperMeshModuleSet(MODS_SDEMAG)) {

			pTO->SetBackgroundColor(UNAVAILABLECOLOR);
			pTO->set(" N/A ");
			stateChanged = true;
			iop.auxId = 0;
			iop.minorId = 0;
		}
		else if (SMesh.IsSuperMeshModuleSet(MODS_SDEMAG) && (dynamic_cast<SDemag*>(SMesh.GetSuperMeshModule(MODS_SDEMAG))->Get_PBC().z != images || !status)) {

			iop.minorId = 1;
			iop.auxId = SMesh.CallModuleMethod(&SDemag::Get_PBC).z;

			if (iop.auxId != 0) {

				pTO->SetBackgroundColor(ONCOLOR);
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
			}

			pTO->set(" " + ToString(iop.auxId) + " ");
			stateChanged = true;
		}
	}
	break;

	//Shows individual shape control flag. auxId is the value (0/1)
	case IOI_INDIVIDUALSHAPE:
	{
		bool status = (bool)iop.auxId;

		if (status != shape_change_individual) {

			iop.auxId = shape_change_individual;

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
	break;

	//Shows image cropping settings : textId has the DBL4 value as text
	case IOI_IMAGECROPPING:
	{
		std::string value_text = iop.textId;

		DBL4 value = ToNum(value_text);

		if (value != image_cropping) {

			iop.textId = ToString(image_cropping);

			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
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
		std::string userConstant_text = iop.textId;

		if (io_index < userConstants.size() && (index != io_index || userConstant_text != Build_EquationConstants_Text(io_index))) {

			iop.textId = Build_EquationConstants_Text(io_index);

			pTO->set(" " + iop.textId + " ");
			stateChanged = true;
		}

		//this object is part of a list : make sure this list is updated
		updateList(io_index, userConstants.size() - 1, &Simulation::Build_EquationConstants_ListLine);
	}
	break;

	//Show skypos diameter multiplier : minorId is the unique mesh id number, textId is the multiplier as a std::string
	case IOI_SKYPOSDMUL:
	{
		int meshId = iop.minorId;
		double multiplier = ToNum(iop.textId);

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			double set_multiplier = SMesh[meshIdx]->Get_skypos_dmul();

			if (set_multiplier == 0.0) {

				pTO->SetBackgroundColor(UNAVAILABLECOLOR);
				pTO->set(" N/A ");
				stateChanged = true;
				iop.auxId = -1;
			}
			else if (multiplier != set_multiplier) {

				iop.textId = ToString(set_multiplier);

				pTO->set(" " + iop.textId + " ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows dwpos fitting component. auxId is the value (-1, 0, 1, 2)
	case IOI_DWPOSCOMPONENT:
	{
		int component = iop.auxId;

		if (component != SMesh.Get_DWPos_Component()) {

			iop.auxId = SMesh.Get_DWPos_Component();

			if (iop.auxId == -1) {

				pTO->set(" Auto ");
			}
			else if (iop.auxId == 0) {

				pTO->set(" x ");
			}
			else if (iop.auxId == 1) {

				pTO->set(" y ");
			}
			else if (iop.auxId == 2) {

				pTO->set(" z ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows Monte-Carlo computation type (serial/parallel) : minorId is the unique mesh id number, auxId is the status (0 : parallel, 1 : serial, -1 : N/A)
	case IOI_MCCOMPUTATION:
	{
		int meshId = iop.minorId;
		int status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (status >= 0 && (bool)status != SMesh[meshIdx]->Get_MonteCarlo_Serial()) {

				iop.auxId = SMesh[meshIdx]->Get_MonteCarlo_Serial();

				if (iop.auxId == 1) {

					pTO->set(" Serial ");
				}
				else {

					pTO->set(" Parallel ");
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows Monte-Carlo algorithm type : minorId is the unique mesh id number, auxId is the type (0 : classical, 1 : constrained, -1 : N/A), textId is the constrained DBL3 direction.
	case IOI_MCTYPE:
	{
		int meshId = iop.minorId;
		int type = iop.auxId;
		std::string direction = iop.textId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (type >= 0 && 
				((bool)type != SMesh[meshIdx]->Get_MonteCarlo_Constrained() || 
				(type == 1 && (DBL3)ToNum(direction) != SMesh[meshIdx]->Get_MonteCarlo_Constrained_Direction()))) {

				iop.auxId = SMesh[meshIdx]->Get_MonteCarlo_Constrained();
				if (iop.auxId == 0) {

					iop.textId = "Classical";
				}
				else {

					iop.textId = ToString(SMesh[meshIdx]->Get_MonteCarlo_Constrained_Direction());
				}

				pTO->set(" " + iop.textId + " ");
				stateChanged = true;
			}
		}
	}
	break;

	//Shows Monte-Carlo disabled/enabled status : minorId is the unique mesh id number, auxId is the status (1 : disabled, 0 : enabled).
	case IOI_MCDISABLED:
	{
		int meshId = iop.minorId;
		int status = iop.auxId;

		int meshIdx = SMesh.contains_id(meshId);

		if (meshIdx >= 0) {

			if (status >= 0 && (bool)status != SMesh[meshIdx]->Get_MonteCarlo_Disabled()) {

				iop.auxId = SMesh[meshIdx]->Get_MonteCarlo_Disabled();

				if (iop.auxId == 0) {

					pTO->set(" Enabled ");
					pTO->SetBackgroundColor(ONCOLOR);
				}
				else {

					pTO->set(" Disabled ");
					pTO->SetBackgroundColor(OFFCOLOR);
				}

				stateChanged = true;
			}
		}
	}
	break;

	//Shows Monte-Carlo computefields state flag : auxId is the state (0: disabled, 1: enabled)
	case IOI_MCCOMPUTEFIELDS:
	{
		int status = iop.auxId;

		if (status != (int)SMesh.Get_MonteCarlo_ComputeFields()) {

			iop.auxId = SMesh.Get_MonteCarlo_ComputeFields();

			if (iop.auxId == 1) {

				pTO->SetBackgroundColor(ONCOLOR);
				pTO->set(" Enabled ");
			}
			else {

				pTO->SetBackgroundColor(OFFCOLOR);
				pTO->set(" Disabled ");
			}

			stateChanged = true;
		}
	}
	break;

	//Shows shape rotation setting: textId is the value as text (DBL3)
	case IOI_SHAPEROT:
	{
		std::string value_text = iop.textId;

		DBL3 value = ToNum(value_text);

		if (value != shape_rotation) {

			iop.textId = ToString(shape_rotation);

			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows shape repetition setting: textId is the value as text (INT3)
	case IOI_SHAPEREP:
	{
		std::string value_text = iop.textId;

		INT3 value = ToNum(value_text);

		if (value != shape_repetitions) {

			iop.textId = ToString(shape_repetitions);

			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows shape displacement setting: textId is the value as text (DBL3)
	case IOI_SHAPEDSP:
	{
		std::string value_text = iop.textId;

		DBL3 value = ToNum(value_text);

		if (value != shape_displacement) {

			iop.textId = ToString(shape_displacement, "m");

			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows shape method setting: textId is the value as text (method)
	case IOI_SHAPEMET:
	{
		std::string value_text = iop.textId;

		if (value_text != shape_method) {

			iop.textId = shape_method;

			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows display render detail level: textId is the value
	case IOI_DISPRENDER_DETAIL:
	{
		std::string value_text = iop.textId;
		
		if (value_text != ToString(BD.GetDetailLevel(), "m")) {

			iop.textId = ToString(BD.GetDetailLevel(), "m");

			pTO->set(" " + iop.textId + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows display render threshold 1: auxId is the value
	case IOI_DISPRENDER_THRESH1:
	{
		int value = iop.auxId;

		if (value != BD.GetRenderThresholds().i) {

			iop.auxId = BD.GetRenderThresholds().i;

			pTO->set(" " + ToString(iop.auxId) + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows display render threshold 2: auxId is the value
	case IOI_DISPRENDER_THRESH2:
	{
		int value = iop.auxId;

		if (value != BD.GetRenderThresholds().j) {

			iop.auxId = BD.GetRenderThresholds().j;

			pTO->set(" " + ToString(iop.auxId) + " ");

			stateChanged = true;
		}
	}
	break;

	//Shows display render threshold 3: auxId is the value
	case IOI_DISPRENDER_THRESH3:
	{
		int value = iop.auxId;

		if (value != BD.GetRenderThresholds().k) {

			iop.auxId = BD.GetRenderThresholds().k;

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