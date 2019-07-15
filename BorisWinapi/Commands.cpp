#include "stdafx.h"
#include "Simulation.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Scripted communication

void Simulation::Listen_Incoming_Message(void) {

	//Listen for incoming messages - non-blocking call, but this method runs on a thread with an infinite loop which keeps calling Simulation::Listen_Incoming_Message
	string message = commSocket.Listen();

	if (message.length()) {

		//handle command using the usual HandleCommand method, but on a blocking thread call - need to wait for it to finish before responding to client
		set_blocking_thread(THREAD_HANDLEMESSAGE);
		single_call_launch<string>(&Simulation::HandleCommand, message, THREAD_HANDLEMESSAGE);

		//the command was handled using a blocking thread, so now parameters should be available to send
		commSocket.SendDataParams();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Command handler

void Simulation::HandleCommand(string command_string) {

	//cannot change simulation parameters in the middle of an interation. simulationMutex also used by Simulate() method, but there it is non-blocking
	simulationMutex.lock();
	lock_guard<mutex> simlock(simulationMutex, adopt_lock);

	//display messages in console only if this is true. Note, some important messages are displayed regardless (e.g. simulation stopped or simulation started messages.)
	bool verbose = true;

	//if command came from a script client, send back required parameters. Each command must ensure it does this.
	bool script_client_connected = false;

	//commands are always of the form: commandname (param1) (param2) ...
	vector<string> command_fields = split(command_string, " ");
	string command_name = command_fields.front();

	if (!command_name.length()) {

		if (verbose) BD.DisplayConsoleError("ERROR: command not recognized.");
		return;
	}

	//command_fields now contains command parameters (if any)
	command_fields.erase(command_fields.begin());

	//'*' and '>' are special command prefixes used by a script client. They specify 1) a script client issued the command, and 2) if command causes text to be displayed in the console.
	if (command_name[0] == '>') { verbose = true; script_client_connected = true; command_name = command_name.substr(1); }
	if (command_name[0] == '*') { verbose = false; script_client_connected = true; command_name = command_name.substr(1); }
	//'~' is a special command prefix which suppresses console output. Not used by script client but by internal program routines.
	if (command_name[0] == '~') { verbose = false; script_client_connected = false; command_name = command_name.substr(1); }

	//show command usage on ? prefix
	if (command_name[0] == '?') {

		command_name = command_name.substr(1);

		if (commands.has_key(command_name))
			PrintCommandUsage(command_name);
		else if (verbose) BD.DisplayConsoleError("ERROR: command not recognized.");

		return;
	}

	if (commands.has_key(command_name)) {

		CommandSpecifier commandSpec = commands[command_name];

		BError error;

		//valid command name entered : process it
		switch (commandSpec.cmdId) {

		case CMD_RUN:
			RunSimulation();
			break;

		case CMD_STOP:
			StopSimulation();
			break;

		case CMD_RESET:
			ResetSimulation();
			break;

		case CMD_COMPUTEFIELDS:
			ComputeFields();
			break;

		//DEPRECATED
		//This was used by Python scripts to detected when simulation has finished. Still keep it for compatibility with old scripts.
		//Instead you should use the ws.Run() call in Python scripts
		case CMD_ISRUNNING:
			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(is_thread_running(THREAD_LOOP)));
		break;

		case CMD_CENTER:
			UpdateScreen_AutoSet();
			break;

		case CMD_ITERUPDATE:
		{
			int iterations;

			error = commandSpec.GetParameters(command_fields, iterations);

			if (!error) {

				iterUpdate = iterations;
				UpdateScreen();
			}
			else if (verbose) BD.DisplayConsoleListing("iterupdate: " + ToString(iterUpdate));

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(iterUpdate));
		}
		break;

		case CMD_CLEARSCREEN:
		{
			ClearScreen();
		}
		break;

		case CMD_REFRESHSCREEN:
		{
			RefreshScreen();
		}
		break;

		case CMD_UPDATESCREEN:
		{
			UpdateScreen();
		}
		break;

		case CMD_MESHRECT:
		{
			Rect meshRect;

			error = commandSpec.GetParameters(command_fields, meshRect);

			if (!error) {

				StopSimulation();

				std::function<void(string, Rect)> save_data_updater = [&](string meshName, Rect meshRect_old) -> void {

					//update rectangles in saveDataList for the named mesh
					UpdateSaveDataEntries(meshRect_old, SMesh[meshName]->GetMeshRect(), meshName);
				};

				if (!err_hndl.call(&SuperMesh::SetMeshRect, &SMesh, SMesh.GetMeshFocus(), meshRect, save_data_updater)) {

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshRect()));
		}
		break;

		case CMD_SCALEMESHRECTS:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {
			
				SMesh.Set_Scale_Rects(status);

				RefreshScreen();
			}
			else if (verbose) Print_Scale_Rects_Status();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.Get_Scale_Rects()));
		}
		break;

		case CMD_MESH:
		{
			if (verbose) Print_Mesh_List();
		}
		break;

		case CMD_CELLSIZE:
		{
			DBL3 h;
			error = commandSpec.GetParameters(command_fields, h);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&Mesh::SetMeshCellsize, SMesh.active_mesh(), h)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshCellsize()));
		}
		break;

		case CMD_ECELLSIZE:
		{
			DBL3 h_e;
			error = commandSpec.GetParameters(command_fields, h_e);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&Mesh::SetMeshECellsize, SMesh.active_mesh(), h_e)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshECellsize()));
		}
		break;

		case CMD_TCELLSIZE:
		{
			DBL3 h_t;
			error = commandSpec.GetParameters(command_fields, h_t);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&Mesh::SetMeshTCellsize, SMesh.active_mesh(), h_t)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshTCellsize()));
		}
		break;

		case CMD_FMSCELLSIZE:
		{
			DBL3 h;
			error = commandSpec.GetParameters(command_fields, h);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::SetFMSMeshCellsize, &SMesh, h)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetFMSMeshCellsize()));
		}
		break;

		case CMD_ESCELLSIZE:
		{
			DBL3 h;
			error = commandSpec.GetParameters(command_fields, h);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::SetESMeshCellsize, &SMesh, h)) {

					UpdateScreen();
				}
			}
			else if (verbose)  PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetESMeshCellsize()));
		}
		break;

		case CMD_ADDFMESH:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_FERROMAGNETIC, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDMETALMESH:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_METAL, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDINSULATORMESH:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_INSULATOR, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDDIPOLEMESH:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_DIPOLE, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DELMESH:
		{
			string meshName;
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();

				//delete mesh
				if (!err_hndl.qcall(&SuperMesh::DelMesh, &SMesh, meshName)) {

					//delete any affected entries in data box
					DeleteDataBoxFields(meshName);

					//delete any affected entries in saveDataList
					DeleteSaveDataEntries(meshName);

					UpdateScreen_AutoSet();
				}
			}
		}
		break;

		case CMD_RENAMEMESH:
		{
			string oldName = SMesh.GetMeshFocus(), newName;

			error = commandSpec.GetParameters(command_fields, oldName, newName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, newName); oldName = SMesh.GetMeshFocus(); }

			if (!error) {

				//note, mesh name is not allowed to have any spaces - needs to be a single word
				newName = trimspaces(newName);

				if (!err_hndl.qcall(&SuperMesh::RenameMesh, &SMesh, oldName, newName)) {

					StopSimulation();

					//must also change any affected labels in data box
					ChangeDataBoxLabels(oldName, newName);

					//update records in saveDataList
					UpdateSaveDataEntries(oldName, newName);

					//update records in simStage (simulation stages)
					UpdateStageMeshNames(oldName, newName);

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_MESHFOCUS:
		{
			string meshName;
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				string foucusedMesh = SMesh.GetMeshFocus();

				if (!err_hndl.qcall(&SuperMesh::SetMeshFocus, &SMesh, meshName)) {

					if(foucusedMesh != SMesh.GetMeshFocus())
						UpdateScreen_AutoSet();
					else UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetMeshFocus()));
		}
		break;

		case CMD_MESHFOCUS2:
		{
			string meshName;
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				string foucusedMesh = SMesh.GetMeshFocus();

				if (!err_hndl.qcall(&SuperMesh::SetMeshFocus, &SMesh, meshName)) {

					if (foucusedMesh != SMesh.GetMeshFocus())
						UpdateScreen_AutoSet_KeepOrientation();
					else UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetMeshFocus()));
		}
		break;

		case CMD_DELRECT:
		{
			Rect rectangle;
			string meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, rectangle, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, rectangle); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::delrect, &SMesh, meshName, rectangle)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDRECT:
		{
			Rect rectangle;
			string meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, rectangle, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, rectangle); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::setrect, &SMesh, meshName, rectangle)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RESETMESH:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error == BERROR_PARAMMISMATCH) { meshName = SMesh.GetMeshFocus(); error.reset(); }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::resetrect, &SMesh, meshName)) {

					UpdateScreen();
				}
			}
		}
		break;

		case CMD_LOADMASKFILE:
		{
			double zDepth;
			string fileName;

			error = commandSpec.GetParameters(command_fields, zDepth, fileName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, fileName); zDepth = 0; }

			if (!error) {

				StopSimulation();

				if (GetFileTermination(fileName) != ".png")
					fileName += ".png";

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				//define code required to load the bitmap - this will be passed for the mesh to use
				function<vector<BYTE>(string, INT2)> bitmap_loader = [&](string fileName, INT2 n_plane) -> vector<BYTE> {

					vector<BYTE> bitmap;
					BD.BGMethods()->GetBitmapFromImage(fileName, bitmap, n_plane);

					return bitmap;
				};

				error = SMesh[SMesh.GetMeshFocus()]->applymask(zDepth, fileName, bitmap_loader);

				if (error) {

					//If fileName has spaces then the above won't work if a zDepth value was not specified by the user.
					//This happens since after the first space the GetParameters command converts to a zDepth value then gets the rest as a file name.
					//In this case the mask file loading will throw an error (couldn't load file), so try again without getting a zDepth value and see if it works.
					error.reset();
					error = commandSpec.GetParameters(command_fields, fileName);

					if (!error) {

						if (GetFileTermination(fileName) != ".png")
							fileName += ".png";

						if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

						error = SMesh[SMesh.GetMeshFocus()]->applymask(0, fileName, bitmap_loader);
					}
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETANGLE:
		{
			double polar, azim;
			string meshName;

			error = commandSpec.GetParameters(command_fields, polar, azim, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, polar, azim); meshName = ""; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetMagnetisationAngle, &SMesh, meshName, polar, azim)) {

					UpdateScreen();
				}

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_INVERTMAG:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetInvertedMagnetisation, &SMesh, meshName)) {

					UpdateScreen();
				}

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETRECT:
		{
			Rect rectangle;
			double polar, azim;
			string meshName;

			error = commandSpec.GetParameters(command_fields, polar, azim, rectangle, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, polar, azim, rectangle); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetMagnetisationAngle_Rect, &SMesh, meshName, polar, azim, rectangle)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DWALL:
		{
			string longitudinal, transverse;
			double width, position;
			string meshName;

			error = commandSpec.GetParameters(command_fields, longitudinal, transverse, width, position, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, longitudinal, transverse, width, position); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetMagnetisationDomainWall, &SMesh, meshName, longitudinal, transverse, width, position)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SKYRMION:
		{
			int core, chirality;
			double diameter;
			DBL2 position;
			string meshName;

			error = commandSpec.GetParameters(command_fields, core, chirality, diameter, position, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, core, chirality, diameter, position); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetSkyrmion, &SMesh, meshName, core, chirality, diameter, position)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SKYRMIONBLOCH:
		{
			int core, chirality;
			double diameter;
			DBL2 position;
			string meshName;

			error = commandSpec.GetParameters(command_fields, core, chirality, diameter, position, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, core, chirality, diameter, position); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetSkyrmionBloch, &SMesh, meshName, core, chirality, diameter, position)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETFIELD:
		{
			DBL3 field_polar;
			string meshName = "";

			error = commandSpec.GetParameters(command_fields, field_polar, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, field_polar); meshName = ""; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetField, &SMesh, meshName, Polar_to_Cartesian(field_polar))) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&Zeeman::GetField)));
		}
		break;

		case CMD_MODULES:
		{
			if (verbose) Print_Modules_List();
		}
		break;

		case CMD_ADDMODULE:
		{
			string moduleHandle, meshName;

			error = commandSpec.GetParameters(command_fields, meshName, moduleHandle);

			if (!error) {

				StopSimulation();

				MOD_ moduleID = (MOD_)moduleHandles.get_ID_from_value(moduleHandle);

				if(!err_hndl.call(&SuperMesh::AddModule, &SMesh, meshName, moduleID)) {

					RefreshScreen();
				}
			}
			else if (verbose) Print_Modules_List();
		}
		break;

		case CMD_DELMODULE:
		{
			string moduleHandle, meshName;
			
			error = commandSpec.GetParameters(command_fields, meshName, moduleHandle);

			if (!error) {

				StopSimulation();

				MOD_ moduleID = (MOD_)moduleHandles.get_ID_from_value(moduleHandle);

				if (!err_hndl.qcall(&SuperMesh::DelModule, &SMesh, meshName, moduleID)) {

					RefreshScreen();
				}
			}
			else if (verbose) Print_Modules_List();
		}
		break;

		case CMD_MULTICONV:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();

				error = SMesh.CallModuleMethod(&SDemag::Set_Multilayered_Convolution, status);

				RefreshScreen();
			}
			else if (verbose) Print_MultiConvolution_Config();
		}
		break;

		case CMD_2DMULTICONV:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();

				error = SMesh.CallModuleMethod(&SDemag::Set_2D_Multilayered_Convolution, status);

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_NCOMMONSTATUS:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();

				error = SMesh.CallModuleMethod(&SDemag::Set_Default_n_status, status);

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_NCOMMON:
		{
			INT3 n_common;

			error = commandSpec.GetParameters(command_fields, n_common);

			if (!error) {

				StopSimulation();

				error = SMesh.CallModuleMethod(&SDemag::Set_n_common, (SZ3)n_common);

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ODE:
		{
			if (verbose) Print_ODEs();
		}
		break;

		case CMD_SETODE:
		{
			string odeHandle, odeEvalHandle;
			
			error = commandSpec.GetParameters(command_fields, odeHandle, odeEvalHandle);

			if (!error) {

				StopSimulation();

				ODE_ setOde = (ODE_)odeHandles.get_ID_from_value(odeHandle);
				EVAL_ odeEval = (EVAL_)odeEvalHandles.get_ID_from_value(odeEvalHandle);

				if (vector_contains(odeAllowedEvals(setOde), odeEval)) {

					if (!err_hndl.call(&SuperMesh::SetODE, &SMesh, setOde, odeEval)) {

						UpdateScreen();
					}
				}
				else if (verbose) error(BERROR_INCORRECTCONFIG);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETDT:
		{
			double dT;

			error = commandSpec.GetParameters(command_fields, dT);

			if (!error) {

				SMesh.SetTimeStep(dT);
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetTimeStep()));
		}
		break;

		case CMD_SHOWDATA:
		{
			string dataName;
			string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			error = commandSpec.GetParameters(command_fields, dataName, meshName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName, meshName); dataRect = Rect(); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName); meshName = SMesh.GetMeshFocus(); dataRect = Rect(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				if (verbose) {

					string text;

					//mesh name if applicable
					if (!dataDescriptor[dataName].meshless) {

						text += "<b><" + meshName + "> ";
					}

					//box dimensions if applicable
					if (!dataDescriptor[dataName].boxless && !dataRect.IsNull()) {

						text += "<b>(" + ToString(dataRect, "m") + ") ";
					}

					//Label
					text += "<b>" + dataDescriptor[dataName].Label;
					//and the actual value(s)
					text += "</c><i>" + GetDataValueString(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect)) + "</i>";

					BD.DisplayFormattedConsoleMessage(text);
				}

				//for script return the number of returned data fields is variable. This is done by obtaining the a string using GetDataValueString, then splitting it using the usual separators (, or ;)
				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(GetDataValue(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect))));

			}
			else if (verbose) Print_ShowData();
		}
		break;

		case CMD_DATA:
		{
			if (verbose) {

				Print_SetOutputData_List();
				Print_AvailableOutputData();
			}

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(saveDataList.size()));
		}
		break;

		case CMD_DELDATA:
		{
			int index = 0;
			
			error = commandSpec.GetParameters(command_fields, index);

			if (!error) {

				if (index < 0) {

					saveDataList.clear();
					NewSaveDataEntry(DATA_TIME, SMesh.GetMeshFocus(), Rect());
				}
				else if (GoodIdx(saveDataList.last(), index) && saveDataList.size() > 1) {

					saveDataList.erase(index);
				}
				else if (verbose) error(BERROR_INCORRECTACTION);

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDDATA:
		{
			string dataName;
			string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			error = commandSpec.GetParameters(command_fields, dataName, meshName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName, meshName); dataRect = Rect(); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName); meshName = SMesh.GetMeshFocus(); dataRect = Rect(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				//add new save data entry if meshname is correct and box is correct
				Rect meshRect_rel = Rect(SMesh[meshName]->GetMeshDimensions());
				if (meshRect_rel.contains(dataRect)) {

					NewSaveDataEntry((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect);
					RefreshScreen();
				}
				else if (verbose) error(BERROR_PARAMOUTOFBOUNDS);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EDITDATA:
		{
			int index = 0;
			string dataName;
			string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			error = commandSpec.GetParameters(command_fields, index, dataName, meshName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, index, dataName, meshName); dataRect = Rect(); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, index, dataName); meshName = SMesh.GetMeshFocus(); dataRect = Rect(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				//add new save data entry if meshname is correct and box is correct
				Rect meshRect_rel = Rect(SMesh[meshName]->GetMeshDimensions());
				if (meshRect_rel.contains(dataRect)) {

					if (GoodIdx(saveDataList.last(), index)) {

						EditSaveDataEntry(index, (DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect);
						RefreshScreen();
					}
					else if (verbose) error(BERROR_PARAMOUTOFBOUNDS);
				}
				else if (verbose) error(BERROR_PARAMOUTOFBOUNDS);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDPINNEDDATA:
		{
			string dataName, meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, dataName, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName); meshName = SMesh.GetMeshFocus(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				if (dataDescriptor[dataName].meshless)
					NewDataBoxField(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName)));
				else
					NewDataBoxField(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName));

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DELPINNEDDATA:
		{
			int index = 0;

			error = commandSpec.GetParameters(command_fields, index);

			if (!error && GoodIdx(dataBoxList.last(), index)) {

				dataBoxList.erase(index);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CHDIR:
		{
			string directory_;

			error = commandSpec.GetParameters(command_fields, directory_);

			if (!error) {

				directory = directory_;

				if (directory.substr(directory.length() - 1) != "\\" && directory.substr(directory.length() - 1) != "/") {

					directory += "\\";
				}

				RefreshScreen();
			}
			else {

				string text = "[tc1,1,1,1/tc]Current working directory : " + MakeIO(IOI_DIRECTORY, directory);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(directory));
		}
		break;

		case CMD_SAVEDATAFILE:
		{
			string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (!error) {

				string directory_ = ExtractFilenameDirectory(fileName);
				if (directory_.length()) directory = directory_;
				savedataFile = fileName;

				if (!GetFileTermination(savedataFile).length()) savedataFile += ".txt";

				RefreshScreen();
			}
			else {

				string text = "[tc1,1,1,1/tc]Current file for output data : " + MakeIO(IOI_DIRECTORY, directory) + MakeIO(IOI_SAVEDATAFILE, savedataFile);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(savedataFile));
		}
		break;

		case CMD_SAVECOMMENT:
		{
			string fileName,  comment;

			error = commandSpec.GetParameters(command_fields, fileName, comment);

			if (!error) {

				string directory_ = ExtractFilenameDirectory(fileName);
				if (!directory_.length()) fileName = directory + fileName;
				if (!GetFileTermination(fileName).length()) fileName += ".txt";

				ofstream bdout;
				bdout.open(fileName, ios::out | ios::app);
				if (bdout.is_open()) {

					bdout << comment << endl;
					bdout.close();
				}
			}
		}
		break;

		case CMD_SAVEIMAGEFILE:
		{
			string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (!error) {

				string directory_ = ExtractFilenameDirectory(fileName);
				if (directory_.length()) directory = directory_;
				imageSaveFileBase = fileName;

				RefreshScreen();
			}
			else {

				string text = "[tc1,1,1,1/tc]Current file image saving : " + MakeIO(IOI_DIRECTORY, directory) + MakeIO(IOI_SAVEIMAGEFILEBASE, imageSaveFileBase);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(imageSaveFileBase));
		}
		break;

		case CMD_DATASAVEFLAG:
		{
			bool status;
			
			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				saveDataFlag = status;

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(saveDataFlag));
		}
		break;

		case CMD_IMAGESAVEFLAG:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				saveImageFlag = status;

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(saveImageFlag));
		}
		break;

		case CMD_STAGES:
		{
			if (verbose) Print_SetStages_List();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(simStages.size()));
		}
		break;

		case CMD_ADDSTAGE:
		{
			string stageTypeName;
			string meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, stageTypeName, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, stageTypeName); meshName = SMesh.GetMeshFocus(); }

			if (!error && stageDescriptors.has_key(stageTypeName) && (SMesh.contains(meshName) || meshName == SMesh.superMeshHandle)) {

				//add new simulation stage to the schedule with default settings for this stage type (can be edited separately)
				AddGenericStage((SS_)stageDescriptors.get_ID_from_key(stageTypeName), meshName);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DELSTAGE:
		{
			int index = 0;

			error = commandSpec.GetParameters(command_fields, index);

			if (!error) {

				if (index < 0) {

					simStages.clear();
					AddGenericStage(SS_RELAX, SMesh.GetMeshFocus());
				}
				else if (GoodIdx(simStages.last(), index) && simStages.size() > 1) {

					DeleteStage(index);
				}
				else if (verbose) error(BERROR_INCORRECTACTION);

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EDITSTAGE:
		{
			int index;
			string stageTypeName;
			string meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, index, stageTypeName, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, index, stageTypeName); meshName = SMesh.GetMeshFocus(); }

			if (!error && GoodIdx(simStages.last(), index) && stageDescriptors.has_key(stageTypeName) && (SMesh.contains(meshName) || meshName == SMesh.superMeshHandle)) {

				//add new simulation stage to the schedule with default settings for this stage type (can be edited separately)
				EditStageType(index, (SS_)stageDescriptors.get_ID_from_key(stageTypeName), meshName);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EDITSTAGEVALUE:
		{
			int index;
			string value_string;

			error = commandSpec.GetParameters(command_fields, index, value_string);

			if (!error && GoodIdx(simStages.last(), index)) {

				EditStageValue(index, value_string);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EDITSTAGESTOP:
		{
			int index;
			string stageStopName;
			string value_string;

			error = commandSpec.GetParameters(command_fields, index, stageStopName, value_string);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, index, stageStopName); value_string = ""; }

			if (!error && (GoodIdx(simStages.last(), index) || index < 0) && stageStopDescriptors.has_key(stageStopName)) {

				if (index >= 0) EditStageStopCondition(index, (STOP_)stageStopDescriptors.get_ID_from_key(stageStopName), value_string);
				else {

					for (int idx = 0; idx < simStages.size(); idx++)
						EditStageStopCondition(idx, (STOP_)stageStopDescriptors.get_ID_from_key(stageStopName), value_string);
				}

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EDITDATASAVE:
		{
			int index;
			string dSaveType;
			string value_string;

			error = commandSpec.GetParameters(command_fields, index, dSaveType, value_string);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, index, dSaveType); value_string = ""; }

			if (!error && (GoodIdx(simStages.last(), index) || index < 0) && dataSaveDescriptors.has_key(dSaveType)) {

				if (index >= 0) EditDataSaveCondition(index, (DSAVE_)dataSaveDescriptors.get_ID_from_key(dSaveType), value_string);
				else {

					for (int idx = 0; idx < simStages.size(); idx++)
						EditDataSaveCondition(idx, (DSAVE_)dataSaveDescriptors.get_ID_from_key(dSaveType), value_string);
				}

				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PARAMS:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			
			if (error == BERROR_PARAMMISMATCH) {

				error.reset(); 
				meshName = SMesh.GetMeshFocus();
			}

			if (SMesh.contains(meshName) && verbose) Print_MeshParams(meshName);
		}
		break;

		case CMD_SETPARAM:
		{
			string paramName, paramValue, meshName;
			bool set_value = true;

			error = commandSpec.GetParameters(command_fields, meshName, paramName, paramValue);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, paramName, paramValue); meshName = SMesh.GetMeshFocus(); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, paramName); meshName = SMesh.GetMeshFocus(); set_value = false; }

			if (!error && set_value) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::set_meshparam_value, &SMesh, meshName, paramName, paramValue)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) {

				SMesh.get_meshparam_value(meshName, paramName, paramValue);
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(paramValue));
			}
		}
		break;

		case CMD_PARAMSTEMP:
		{
			string meshName;
			
			error = commandSpec.GetParameters(command_fields, meshName);

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				meshName = SMesh.GetMeshFocus();
			}

			if (SMesh.contains(meshName) && verbose) Print_MeshParamsTemperature(meshName);
		}
		break;

		case CMD_CLEARPARAMSTEMP:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				meshName = "";
			}

			StopSimulation();

			if (!err_hndl.qcall(&SuperMesh::clear_meshparam_temp, &SMesh, meshName)) {

				UpdateScreen();
			}
		}
		break;

		case CMD_SETPARAMTEMP:
		{
			string meshName, paramName, formulaName, coeffs_string;

			error = commandSpec.GetParameters(command_fields, meshName, paramName, formulaName, coeffs_string);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName, formulaName); coeffs_string = ""; }

			if (!error) {

				StopSimulation();

				vector<string> coeffs_strings_vec = split(coeffs_string, " ");
				vector<double> coeffs_vec = vec_convert<double, string>(coeffs_strings_vec);

				if (!err_hndl.qcall(&SuperMesh::set_meshparam_formula, &SMesh, meshName, paramName, formulaName, coeffs_vec)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETPARAMTEMPARRAY:
		{
			string paramName;
			string fileName;

			error = commandSpec.GetParameters(command_fields, paramName, fileName);

			//setting parameter with arrays loaded form file
			if (!error) {

				//load from a file directly
				StopSimulation();

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (GetFileTermination(fileName) != ".txt")
					fileName += ".txt";

				vector<vector<double>> load_arrays;
				if (ReadDataColumns(fileName, "\t", load_arrays, { 0, 1 })) {

					error = SMesh.set_meshparam_tscaling_array(SMesh.GetMeshFocus(), paramName, load_arrays[0], load_arrays[1]);

					UpdateScreen();
				}
				else error(BERROR_COULDNOTLOADFILE);
			}

			//if the above failed then try to load array from dp_arrays
			if(error == BERROR_COULDNOTLOADFILE) {

				int dp_T_idx, dp_c_ix;

				error.reset() = commandSpec.GetParameters(command_fields, paramName, dp_T_idx, dp_c_ix);

				if (!error && GoodIdx(dpArr.size(), dp_T_idx, dp_c_ix)) {

					//load from dp arrays
					StopSimulation();

					error = SMesh.set_meshparam_tscaling_array(SMesh.GetMeshFocus(), paramName, dpArr[dp_T_idx], dpArr[dp_c_ix]);
					
					UpdateScreen();
				}
			}
			
			if (verbose && error) PrintCommandUsage(command_name);
			
			if (verbose && !error) BD.DisplayConsoleMessage("Arrays loaded.");

			RefreshScreen();
		}
		break;

		case CMD_COPYPARAMS:
		{
			string meshName_from, meshNames_to_string;

			error = commandSpec.GetParameters(command_fields, meshName_from, meshNames_to_string);

			//multiple mesh names could have been entered
			vector<string> meshNames_to = split(meshNames_to_string, " ");

			if (!error && meshNames_to.size()) {

				StopSimulation();

				for (int idx = 0; idx < meshNames_to.size(); idx++) {

					err_hndl.qcall(&SuperMesh::copy_mesh_parameters, &SMesh, meshName_from, meshNames_to[idx]);
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_COPYMESHDATA:
		{
			string meshName_from, meshNames_to_string;

			error = commandSpec.GetParameters(command_fields, meshName_from, meshNames_to_string);

			//multiple mesh names could have been entered
			vector<string> meshNames_to = split(meshNames_to_string, " ");

			if (!error && meshNames_to.size()) {

				StopSimulation();

				for (int idx = 0; idx < meshNames_to.size(); idx++) {

					err_hndl.qcall(&SuperMesh::copy_mesh_data, &SMesh, meshName_from, meshNames_to[idx]);
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PARAMSVAR:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				meshName = SMesh.GetMeshFocus();
			}

			if (SMesh.contains(meshName) && verbose) Print_MeshParamsVariation(meshName);
		}
		break;

		case CMD_SETDISPLAYEDPARAMSVAR:
		{
			string meshName, paramName;

			error = commandSpec.GetParameters(command_fields, meshName, paramName);

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::set_meshparamvar_display, &SMesh, meshName, paramName)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CLEARPARAMSVAR:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				meshName = "";
			}

			StopSimulation();

			if (!err_hndl.qcall(&SuperMesh::clear_meshparam_variation, &SMesh, meshName)) {

				UpdateScreen();
			}
		}
		break;

		case CMD_CLEARPARAMVAR:
		{
			string meshName, paramName;

			error = commandSpec.GetParameters(command_fields, meshName, paramName);

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::clear_meshparam_variation, &SMesh, meshName, paramName)) {

					UpdateScreen();
				}
			}
		}
		break;

		case CMD_SETPARAMVAR:
		{
			string meshName, paramName, generatorName, generatorArgs;

			error = commandSpec.GetParameters(command_fields, meshName, paramName, generatorName, generatorArgs);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName, generatorName); generatorArgs = ""; }

			if (!error) {

				StopSimulation();

				//used with custom parameter variation generator (set from grayscale png file)
				function<vector<BYTE>(string, INT2)> bitmap_loader = [&](string fileName, INT2 n_plane) -> vector<BYTE> {

					vector<BYTE> bitmap;
					BD.BGMethods()->GetBitmapFromImage(fileName, bitmap, n_plane);

					return bitmap;
				};

				error = SMesh.set_meshparam_var(meshName, paramName, generatorName, generatorArgs, bitmap_loader);

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVESIM:
		{
			string simFileName;

			error = commandSpec.GetParameters(command_fields, simFileName);

			if (error) {

				error.reset();
				simFileName = currentSimulationFile;
			}

			error = SaveSimulation(simFileName);

			if(verbose && !error) BD.DisplayConsoleMessage("Simulation saved : " + simFileName);
		}
		break;

		case CMD_LOADSIM:
		{
			string simFileName;

			error = commandSpec.GetParameters(command_fields, simFileName);

			if (!error) {

				if (!err_hndl.call(&Simulation::LoadSimulation, this, simFileName)) {

					if (verbose) BD.DisplayConsoleMessage("Simulation loaded : " + simFileName);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DEFAULT:
		{
			LoadSimulation(GetDirectory() + std::string("User\\") + "default");
		}
		break;

		case CMD_DISPLAY:
		{
			string name, meshName;

			error = commandSpec.GetParameters(command_fields, name, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, name); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				MESHDISPLAY_ display = (MESHDISPLAY_)displayHandles.get_ID_from_value(name);
				
				if (!err_hndl.qcall(&SuperMesh::SetDisplayedPhysicalQuantity, &SMesh, meshName, (int)display)) {

					UpdateScreen();
				}
			}
			else if(verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_SAVEMESHIMAGE:
		{
			string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				fileName = directory + imageSaveFileBase + ".png";
			}
			else {

				if (GetFileTermination(fileName) != ".png")
					fileName += ".png";

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
			}

			if (BD.SaveMeshImage(fileName)) {

				if (verbose) BD.DisplayConsoleMessage("Saved : " + fileName);
			}
			else if (verbose) error(BERROR_COULDNOTSAVEFILE);
		}
		break;

		case CMD_MAKEVIDEO:
		{
			string fileNameBase;
			double fps;
			int quality;
			
			error = commandSpec.GetParameters(command_fields, fileNameBase, fps, quality);

			if (!error) {

				StopSimulation();

				string sequenceDirectory = directory;
				if (GetFilenameDirectory(fileNameBase).length()) sequenceDirectory = ExtractFilenameDirectory(fileNameBase);

				vector<string> fileNames = GetFilesInDirectory(sequenceDirectory, fileNameBase, ".png");

				BD.DisplayConsoleMessage("Encoding video ... please wait.");

				if (BD.BGMethods()->MakeVideoFromFileSequence(sequenceDirectory + fileNameBase + ".wmv", fileNames, UINT32(fps), 1.0, quality)) {

					if(verbose) BD.DisplayConsoleMessage("Video created.");
				}
				else if (verbose) error(BERROR_COULDNOTSAVEFILE);
			}
			else if(verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_MOVINGMESH:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();

				if (meshName == "0" || meshName == "1") {

					//with "0" : turn off moving mesh trigger. With "1" trigger on first ferromagnetic mesh
					SMesh.SetMoveMeshTrigger(ToNum(meshName));
				}
				else if(SMesh.contains(meshName)) SMesh.SetMoveMeshTrigger(true, meshName);

				UpdateScreen();
			}
			else if (verbose) PrintMovingMeshSettings();
		}
		break;

		case CMD_MOVINGMESHASYM:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();

				SMesh.SetMoveMeshAntisymmetric(status);

				UpdateScreen();
			}
			else if (verbose) PrintMovingMeshSettings();
		}
		break;

		case CMD_MOVINGMESHTHRESH:
		{
			double threshold;

			error = commandSpec.GetParameters(command_fields, threshold);

			if (!error) {

				StopSimulation();

				SMesh.SetMoveMeshThreshold(threshold);

				UpdateScreen();
			}
			else if (verbose) PrintMovingMeshSettings();
		}
		break;

		case CMD_PREPAREMOVINGMESH:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				err_hndl.call(&SuperMesh::PrepareMovingMesh, &SMesh, meshName);				
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PREPAREMOVINGBLOCHMESH:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				err_hndl.call(&SuperMesh::PrepareMovingMesh_Bloch, &SMesh, meshName);
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PREPAREMOVINGNEELMESH:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				err_hndl.call(&SuperMesh::PrepareMovingMesh_Neel, &SMesh, meshName);
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PREPAREMOVINGSKYRMIONMESH:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				err_hndl.call(&SuperMesh::PrepareMovingMesh_Skyrmion, &SMesh, meshName);
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_COUPLETODIPOLES:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();

				SMesh.Set_Coupled_To_Dipoles(status);

				UpdateScreen();
			}
			else if (verbose) Print_CoupledToDipoles_Settings();
		}
		break;

		case CMD_ADDELECTRODE:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				Rect electrode_rect;

				error = commandSpec.GetParameters(command_fields, electrode_rect);

				if (!error) {

					StopSimulation();

					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, electrode_rect);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_DELELECTRODE:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				int index;

				error = commandSpec.GetParameters(command_fields, index);

				if (!error) {

					StopSimulation();

					if (!SMesh.CallModuleMethod(&STransport::DelElectrode, index))
						if (verbose) error(BERROR_PARAMOUTOFBOUNDS);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_CLEARELECTRODES:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				StopSimulation();

				SMesh.CallModuleMethod(&STransport::ClearElectrodes);

				UpdateScreen();
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETDEFAULTELECTRODES:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				StopSimulation();

				//first clear all electrodes
				SMesh.CallModuleMethod(&STransport::ClearElectrodes);

				Rect smeshRect = SMesh.GetESMeshRect();

				//left hand side electrode
				SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(smeshRect.s, DBL3(smeshRect.s.x, smeshRect.e.y, smeshRect.e.z)));
				SMesh.CallModuleMethod(&STransport::DesignateGroundElectrode, 0);

				//right hand side electrode
				SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(DBL3(smeshRect.e.x, smeshRect.s.y, smeshRect.s.z), smeshRect.e));

				SMesh.CallModuleMethod(&STransport::SetPotential, 0.0);

				UpdateScreen();
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETPOTENTIAL:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				double potential;

				error = commandSpec.GetParameters(command_fields, potential);

				if (!error) {

					SMesh.CallModuleMethod(&STransport::SetPotential, potential);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&STransport::GetPotential)));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETCURRENT:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				double current;

				error = commandSpec.GetParameters(command_fields, current);

				if (!error) {

					SMesh.CallModuleMethod(&STransport::SetCurrent, current);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&STransport::GetCurrent)));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_ELECTRODES:
		{
			Print_Electrodes_List();
		}
		break;

		case CMD_SETELECTRODERECT:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				int el_index;
				Rect electrode_rect;

				error = commandSpec.GetParameters(command_fields, el_index, electrode_rect);

				if (!error) {

					StopSimulation();

					if (!SMesh.CallModuleMethod(&STransport::SetElectrodeRect, el_index, electrode_rect))
						if (verbose) error(BERROR_PARAMOUTOFBOUNDS);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETELECTRODEPOTENTIAL:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				double electrode_potential;
				int el_index = 0;

				error = commandSpec.GetParameters(command_fields, el_index, electrode_potential);

				if (!error) {

					if (!SMesh.CallModuleMethod(&STransport::SetElectrodePotential, el_index, electrode_potential))
						if (verbose) error(BERROR_PARAMOUTOFBOUNDS);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, el_index).second));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_DESIGNATEGROUND:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				int el_index;

				error = commandSpec.GetParameters(command_fields, el_index);

				if (!error) {

					SMesh.CallModuleMethod(&STransport::DesignateGroundElectrode, el_index);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_TSOLVERCONFIG:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				double conv_error;
				int iters_timeout;

				error = commandSpec.GetParameters(command_fields, conv_error, iters_timeout);
				if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, conv_error); iters_timeout = 0; }

				if (!error) {

					SMesh.CallModuleMethod(&STransport::SetConvergenceError, conv_error, iters_timeout);

					UpdateScreen();
				}
				else if (verbose) PrintTransportSolverConfig();

				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(
						DBL2(SMesh.CallModuleMethod(&STransport::GetConvergenceError),
						SMesh.CallModuleMethod(&STransport::GetConvergenceTimeout))));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SSOLVERCONFIG:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				double conv_error;
				int iters_timeout;

				error = commandSpec.GetParameters(command_fields, conv_error, iters_timeout);
				if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, conv_error); iters_timeout = 0; }

				if (!error) {

					SMesh.CallModuleMethod(&STransport::SetSConvergenceError, conv_error, iters_timeout);

					UpdateScreen();
				}
				else if (verbose) PrintTransportSolverConfig();

				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(
						DBL2(SMesh.CallModuleMethod(&STransport::GetSConvergenceError),
							SMesh.CallModuleMethod(&STransport::GetSConvergenceTimeout))));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETFIXEDSOR:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				bool fixed_SOR_damping;

				error = commandSpec.GetParameters(command_fields, fixed_SOR_damping);

				if (!error) {

					SMesh.CallModuleMethod(&STransport::SetSORType, fixed_SOR_damping);

					UpdateScreen();
				}
				else if (verbose) PrintTransportSolverConfig();

				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&STransport::IsFixedSORdamping)));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETSORDAMPING:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				DBL2 SOR_damping;

				error = commandSpec.GetParameters(command_fields, SOR_damping);

				if (!error) {

					SMesh.CallModuleMethod(&STransport::SetSORDamping, SOR_damping);

					UpdateScreen();
				}
				else if (verbose) PrintTransportSolverConfig();

				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&STransport::GetSORDamping)));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_TEMPERATURE:
		{
			double Temperature;
			string meshName;

			error = commandSpec.GetParameters(command_fields, Temperature, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, Temperature); meshName = SMesh.superMeshHandle; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetBaseTemperature, &SMesh, meshName, Temperature)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_MeshTemperature_List();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetBaseTemperature()));
		}
		break;

		case CMD_SETHEATDT:
		{
			double dT;

			error = commandSpec.GetParameters(command_fields, dT);

			if (!error) {

				StopSimulation();

				SMesh.CallModuleMethod(&SHeat::set_heat_dT, dT);

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&SHeat::get_heat_dT)));
		}
		break;

		case CMD_AMBIENTTEMPERATURE:
		{
			double T_ambient;
			string meshName;

			error = commandSpec.GetParameters(command_fields, T_ambient, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, T_ambient); meshName = SMesh.superMeshHandle; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetAmbientTemperature, &SMesh, meshName, T_ambient)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_HeatBoundaries_List();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&Heat::GetAmbientTemperature)));
		}
		break;

		case CMD_ROBINALPHA:
		{
			double alpha_boundary;
			string meshName;

			error = commandSpec.GetParameters(command_fields, alpha_boundary, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, alpha_boundary); meshName = SMesh.superMeshHandle; }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::SetAlphaHeatBoundary, &SMesh, meshName, alpha_boundary)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_HeatBoundaries_List();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&Heat::GetAlphaBoundary)));
		}
		break;

		case CMD_INSULATINGSIDES:
		{
			string literal;
			bool status;
			string meshName;

			error = commandSpec.GetParameters(command_fields, literal, status, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, literal, status); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::SetInsulatingSides, &SMesh, meshName, literal, status)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_HeatBoundaries_List();

			if (script_client_connected) {
				if (SMesh.active_mesh()->IsModuleSet(MOD_HEAT)) {

					vector<bool> insulating = SMesh.active_mesh()->CallModuleMethod(&Heat::GetInsulatingSides);
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(insulating[0], insulating[1], insulating[2], insulating[3], insulating[4], insulating[5]));
				}
			}
		}
		break;

		case CMD_CURIETEMPERATURE:
		{
			double T_Curie;
			string meshName;

			error = commandSpec.GetParameters(command_fields, T_Curie, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, T_Curie); meshName = SMesh.superMeshHandle; }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::SetCurieTemperature, &SMesh, meshName, T_Curie)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetCurieTemperature()));
		}
		break;

		case CMD_CURIETEMPERATUREMATERIAL:
		{
			double T_Curie_material;
			string meshName;

			error = commandSpec.GetParameters(command_fields, T_Curie_material, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, T_Curie_material); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::SetCurieTemperatureMaterial, &SMesh, meshName, T_Curie_material)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetCurieTemperatureMaterial()));
		}
		break;

		case CMD_ATOMICMOMENT:
		{
			double atomic_moment;
			string meshName;

			error = commandSpec.GetParameters(command_fields, atomic_moment, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, atomic_moment); meshName = SMesh.superMeshHandle; }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::SetAtomicMagneticMoment, &SMesh, meshName, atomic_moment)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetAtomicMoment()));
		}
		break;

		case CMD_CUDA:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);
			
			if (!error) {

				if (cudaAvailable) {

					if (status != cudaEnabled) {

						StopSimulation();
						
						if (!err_hndl.qcall(&SuperMesh::SwitchCUDAState, &SMesh, status)) {

							cudaEnabled = status;
						}
						else {
							
							SMesh.SwitchCUDAState(false);
							cudaEnabled = false;
						}
						
						UpdateScreen();
					}
				}
				else error(BERROR_NOTAVAILABLE);
			}
			else if (verbose) Print_CUDAStatus();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(cudaEnabled));
		}
		break;

		case CMD_MEMORY:
		{
			if (verbose) Print_MemoryInfo();
		}
		break;

		case CMD_OPENMANUAL:
		{
			string directory = GetDirectory() + "Manual\\";
			string fileName = "BorisManual-v" + ToString(Program_Version) + ".pdf";

			open_file(directory + fileName);
		}
		break;

		case CMD_REFINEROUGHNESS:
		{
			string meshName;
			INT3 refine;

			error = commandSpec.GetParameters(command_fields, refine, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, refine); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetMeshRoughnessRefinement, &SMesh, meshName, refine)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_MeshRoughnessRefinement(SMesh.GetMeshFocus());
		}
		break;

		case CMD_CLEARROUGHNESS:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::ClearMeshRoughness, &SMesh, meshName)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ROUGHENMESH:
		{
			double depth;
			string axis;
			int seed;

			error = commandSpec.GetParameters(command_fields, depth, axis, seed);
			if (error == BERROR_PARAMMISMATCH) { 
				
				error.reset() = commandSpec.GetParameters(command_fields, depth, axis); 
				seed = 1;

				if (error == BERROR_PARAMMISMATCH) {

					error.reset() = commandSpec.GetParameters(command_fields, depth);
					axis = "z";
					seed = 1;
				}
			}

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::RoughenMeshSides, &SMesh, SMesh.GetMeshFocus(), axis, depth, seed)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;
		
		case CMD_SURFROUGHENJAGGED:
		{
			double depth, spacing;
			string sides;
			int seed;

			error = commandSpec.GetParameters(command_fields, depth, spacing, seed, sides);
			if (error == BERROR_PARAMMISMATCH) {

				error.reset() = commandSpec.GetParameters(command_fields, depth, spacing, seed);

				if (error == BERROR_PARAMMISMATCH) {

					error.reset() = commandSpec.GetParameters(command_fields, depth, spacing);
					seed = 1;
				}
			}

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::RoughenMeshSurfaces_Jagged, &SMesh, SMesh.GetMeshFocus(), depth, spacing, seed, sides)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GENERATE2DGRAINS:
		{
			double spacing;
			int seed;

			error = commandSpec.GetParameters(command_fields, spacing, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, spacing); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::GenerateGrains2D, &SMesh, SMesh.GetMeshFocus(), spacing, seed)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GENERATE3DGRAINS:
		{
			double spacing;
			int seed;

			error = commandSpec.GetParameters(command_fields, spacing, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, spacing); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::GenerateGrains3D, &SMesh, SMesh.GetMeshFocus(), spacing, seed)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_BENCHTIME:
		{
			unsigned int benchtime_ms = 0;

			if (sim_end_ms > sim_start_ms) benchtime_ms = sim_end_ms - sim_start_ms;

			if (verbose) BD.DisplayConsoleListing("Simulation duration : " + ToString(benchtime_ms) + " ms.");

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(benchtime_ms));
		}
		break;

		case CMD_MATERIALSDATABASE:
		{
			string mdbName;

			error = commandSpec.GetParameters(command_fields, mdbName);

			if (!error) {

				if (!err_hndl.call(&MaterialsDB::SwitchDataBase, &mdb, mdbName)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_MaterialsDatabase();
		}
		break;

		case CMD_ADDMDBENTRY:
		{
			string meshName, materialName;

			error = commandSpec.GetParameters(command_fields, meshName, materialName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName); materialName = meshName; }

			if (!error) {

				if (SMesh.contains(meshName)) {

					error = mdb.AddMDBEntry(materialName, *SMesh[meshName], SMesh[meshName]->GetMeshType());

					if (!error) BD.DisplayConsoleMessage("Material added to local database.");
				}
				else BD.DisplayConsoleError("Mesh name doesn't exist.");

				UpdateScreen();

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DELMDBENTRY:
		{
			string materialName;

			error = commandSpec.GetParameters(command_fields, materialName);

			if (!error) {

				err_hndl.call(&MaterialsDB::DelMDBEntry, &mdb, materialName);

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_REFRESHMDB:
		{
			err_hndl.call(&MaterialsDB::RefreshMDB, &mdb);
		}
		break;

		case CMD_ADDMATERIAL:
		{
			string materialName;
			Rect meshRect;

			error = commandSpec.GetParameters(command_fields, materialName, meshRect);

			if (!error) {

				int meshType;

				StopSimulation();

				if (!err_hndl.call(&MaterialsDB::LoadMaterial, &mdb, materialName, &meshType)) {

					//material loaded in mdb and meshType available.
					
					string meshName = materialName;

					//does this mesh name already exist?
					while (SMesh.contains(meshName)) {

						//if meshName exists then add _#, where # = 1 to start
						//if meshName already contains a termination of the form _# then increase # until meshName available
						size_t pos = meshName.find_last_of('_');

						if (pos == std::string::npos) {

							meshName += string("_1");
						}
						else {

							string termination = meshName.substr(pos + 1);
							if (has_digits_only(termination)) {

								int number = ToNum(termination);
								number++;
								termination = ToString(number);

								meshName = meshName.substr(0, pos + 1) + termination;
							}
							else meshName += string("_1");
						}
					}

					//first make the mesh with the correct rectangle and mesh type
					if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, (MESH_)meshType, meshRect)) {

						//mesh created, now copy parameter values
						mdb.copy_parameters(*SMesh[meshName]);
					}
				}

				UpdateScreen();

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETMATERIAL:
		{
			string materialName;

			error = commandSpec.GetParameters(command_fields, materialName);

			if (!error) {

				if (!err_hndl.call(&MaterialsDB::LoadMaterial, &mdb, materialName, (int*)nullptr)) {

					//material loaded in mdb available.
	
					//now copy parameter values - this is done even if mesh types are not the same. e.g. you might want to copy from a ferromagnetic mesh to a dipole mesh
					mdb.copy_parameters(*SMesh.active_mesh());
				}

				UpdateScreen();

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_REQMDBSYNC:
		{
			string materialName, emailContact;

			error = commandSpec.GetParameters(command_fields, materialName, emailContact);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, materialName); emailContact = "N/A"; }

			if (!error) {

				string returnMessage;

				if (!err_hndl.call(&MaterialsDB::RequestMDBSync, &mdb, materialName, domain_name, mdb_entry_handler, emailContact, &returnMessage)) {

					BD.DisplayConsoleMessage(returnMessage);
				}
				else {

					BD.DisplayConsoleError(returnMessage);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_UPDATEMDB:
		{
			string returnMessage;

			if (!err_hndl.call(&MaterialsDB::UpdateMDB, &mdb,domain_name, mdb_update_handler, &returnMessage)) {

				BD.DisplayConsoleMessage("Materials database successfully updated.");
			}
			else BD.DisplayConsoleError(returnMessage);
		}
		break;

		case CMD_SHOWLENGHTS:
		{
			if (verbose) {

				if (SMesh.active_mesh()->Magnetisation_Enabled()) {
					
					double A = SMesh.active_mesh()->A;
					double Ms = SMesh.active_mesh()->Ms;
					double Ku = SMesh.active_mesh()->K1;
					double D = SMesh.active_mesh()->D;

					string l_ex, l_Bloch("N/A"), l_sky("N/A");

					l_ex = ToString(sqrt(2 * A / (MU0 * Ms*Ms)), "m");

					if (IsNZ(Ku)) l_Bloch = ToString(sqrt(A / Ku), "m");
					if (IsNZ(Ku) && IsNZ(D)) l_sky = ToString(sqrt(PI * D / (4 * Ku)), "m");

					BD.DisplayConsoleMessage("l_ex = " + l_ex + ", l_Bloch = " + l_Bloch + ", l_sky = " + l_sky);
				}
				else BD.DisplayConsoleError("Focused mesh is not a ferromagnetic mesh.");
			}
		}
		break;

		case CMD_SHOWMCELLS:
		{
			if (SMesh.active_mesh()->Magnetisation_Enabled()) {

				if (verbose) BD.DisplayConsoleMessage("M discretisation cells : " + ToString(SMesh.active_mesh()->n));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->n));
			}
			else if (verbose) BD.DisplayConsoleError("Focused mesh is not a ferromagnetic mesh.");
		}
		break;

		//---------------- CMD_DP_ commands here

		case CMD_DP_CLEARALL:
		{
			dpArr.clear_all();
		}
		break;

		case CMD_DP_CLEAR:
		{
			unsigned arr_idx;

			error = commandSpec.GetParameters(command_fields, arr_idx);

			if (!error) {

				dpArr.clear(arr_idx);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;
		
		case CMD_DP_SHOWSIZES:
		{
			int dp_arr_idx;

			error = commandSpec.GetParameters(command_fields, dp_arr_idx);
			if (error) { error.reset(); dp_arr_idx = -1; }

			if (verbose) {

				if (dp_arr_idx < 0) {

					for (int arr_idx = 0; arr_idx < dpArr.size(); arr_idx++) {

						if (dpArr[arr_idx].size()) BD.DisplayConsoleListing("dp array " + ToString(arr_idx) + " has size " + ToString((int)dpArr[arr_idx].size()));
					}
				}
				else if (dpArr[dp_arr_idx].size()) BD.DisplayConsoleListing("dp array " + ToString(dp_arr_idx) + " has size " + ToString((int)dpArr[dp_arr_idx].size()));
			}

			if (script_client_connected && dp_arr_idx >= 0)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters((int)dpArr[dp_arr_idx].size()));
		}
		break;

		case CMD_DP_GET:
		{
			int dp_arr, index;

			error = commandSpec.GetParameters(command_fields, dp_arr, index);

			if (!error && index < dpArr[dp_arr].size()) {

				BD.DisplayConsoleListing("dpArr[" + ToString(dp_arr) + "][" + ToString(index) + "] = " + ToString(dpArr[dp_arr][index]));
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected && index < dpArr[dp_arr].size())
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(dpArr[dp_arr][index]));
		}
		break;

		case CMD_DP_SET:
		{
			int dp_arr, index;
			double value;

			error = commandSpec.GetParameters(command_fields, dp_arr, index, value);

			if (!error && index < dpArr[dp_arr].size()) {

				dpArr[dp_arr][index] = value;
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_LOAD:
		{
			string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric string (the file path and file name) then a numeric string.
			vector<string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				string indexes_string = entries.back();

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (GetFileTermination(fileName) != ".txt")
					fileName += ".txt";

				int rows_read;

				error = dpArr.load_arrays(fileName, vec_convert<int, string>(split(indexes_string, " ")), &rows_read);

				if (verbose && !error) BD.DisplayConsoleMessage("Number of rows read : " + ToString(rows_read));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVE:
		{
			string fileName;
			string indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName, indexes_string);

			if (!error) {

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (GetFileTermination(fileName) != ".txt")
					fileName += ".txt";

				error = dpArr.save_arrays(fileName, vec_convert<int, string>(split(indexes_string, " ")));

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_GETPROFILE:
		{
			DBL3 start, end;
			int arr_idx;

			error = commandSpec.GetParameters(command_fields, start, end, arr_idx);

			if (!error) {

				error = dpArr.get_profile(start, end, &SMesh, arr_idx);

				if (verbose && !error) BD.DisplayConsoleMessage("Extracted profile between : " + ToString(start, "m") + " and " + ToString(end, "m") + ".");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_APPEND:
		{
			int dp_original, dp_new;

			error = commandSpec.GetParameters(command_fields, dp_original, dp_new);

			if (!error) {

				error = dpArr.append_array(dp_original, dp_new);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SEQUENCE:
		{
			int dp_idx, points;
			double start_value, increment;

			error = commandSpec.GetParameters(command_fields, dp_idx, start_value, increment, points);

			if (!error) {

				error = dpArr.generate_sequence(dp_idx, start_value, increment, points);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_RAREFY:
		{
			int dp_in, dp_out, skip;

			error = commandSpec.GetParameters(command_fields, dp_in, dp_out, skip);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in, dp_out); skip = 1; }

			if (!error) {

				error = dpArr.rarefy(dp_in, dp_out, skip);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_EXTRACT:
		{
			int dp_in, dp_out, start_index, length;

			error = commandSpec.GetParameters(command_fields, dp_in, dp_out, start_index, length);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in, dp_out, start_index); length = -1; }

			if (!error) {

				error = dpArr.extract(dp_in, dp_out, start_index, length);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_ERASE:
		{
			int dp_index, start_index, length;

			error = commandSpec.GetParameters(command_fields, dp_index, start_index, length);

			if (!error) {

				error = dpArr.erase(dp_index, start_index, length);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_AVERAGEMESHRECT:
		{
			Rect rect;

			error = commandSpec.GetParameters(command_fields, rect);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); rect = SMesh.active_mesh()->GetMeshRect(); }

			if (!error) {
				
				if (verbose) BD.DisplayConsoleMessage("Average value = " + dpArr.get_meshaverage(&SMesh, SMesh.GetMeshFocus(), rect));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_TOPOCHARGE:
		{
			double x, y, radius;

			error = commandSpec.GetParameters(command_fields, x, y, radius);
			if (error == BERROR_PARAMMISMATCH) { 
				
				error.reset(); 
				x = 0.0; y = 0.0;
				DBL3 meshDims = SMesh.active_mesh()->GetMeshDimensions();
				radius = sqrt(meshDims.x * meshDims.x + meshDims.y * meshDims.y);
			}

			if (!error) {

				if (SMesh.active_mesh()->Magnetisation_Enabled()) {

					double Q = 0.0;

					error = dpArr.get_topological_charge(SMesh.active_mesh()->Get_M(), x, y, radius, &Q);

					if (!error) {

						if (verbose) BD.DisplayConsoleMessage("Q = " + ToString(Q));

						if (script_client_connected)
							commSocket.SetSendData(commandSpec.PrepareReturnParameters(Q));
					}
				}
				else {

					if (verbose) BD.DisplayConsoleError("Focused mesh is not a ferromagnetic mesh.");
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_ADD:
		{
			int dp_source, dp_dest;
			double value;

			error = commandSpec.GetParameters(command_fields, dp_source, value, dp_dest);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_source, value); dp_dest = dp_source; }

			if (!error) {

				error = dpArr.add(dp_source, dp_dest, value);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SUB:
		{
			int dp_source, dp_dest;
			double value;

			error = commandSpec.GetParameters(command_fields, dp_source, value, dp_dest);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_source, value); dp_dest = dp_source; }

			if (!error) {

				error = dpArr.subtract(dp_source, dp_dest, value);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_MUL:
		{
			int dp_source, dp_dest;
			double value;

			error = commandSpec.GetParameters(command_fields, dp_source, value, dp_dest);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_source, value); dp_dest = dp_source; }

			if (!error) {

				error = dpArr.multiply(dp_source, dp_dest, value);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_DIV:
		{
			int dp_source, dp_dest;
			double value;

			error = commandSpec.GetParameters(command_fields, dp_source, value, dp_dest);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_source, value); dp_dest = dp_source; }

			if (!error) {

				error = dpArr.divide(dp_source, dp_dest, value);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_DOTPROD:
		{
			int dp_vector, dp_out;
			DBL3 u_vector;

			error = commandSpec.GetParameters(command_fields, dp_vector, u_vector, dp_out);

			if (!error) {

				error = dpArr.dotproduct(dp_vector, dp_vector + 1, dp_vector + 2, u_vector, dp_out);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_ADDDP:
		{
			int dp_x1, dp_x2, dp_dest;

			error = commandSpec.GetParameters(command_fields, dp_x1, dp_x2, dp_dest);

			if (!error) {

				error = dpArr.add_dparr(dp_x1, dp_x2, dp_dest);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SUBDP:
		{
			int dp_x1, dp_x2, dp_dest;

			error = commandSpec.GetParameters(command_fields, dp_x1, dp_x2, dp_dest);

			if (!error) {

				error = dpArr.subtract_dparr(dp_x1, dp_x2, dp_dest);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_MULDP:
		{
			int dp_x1, dp_x2, dp_dest;

			error = commandSpec.GetParameters(command_fields, dp_x1, dp_x2, dp_dest);

			if (!error) {

				error = dpArr.multiply_dparr(dp_x1, dp_x2, dp_dest);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_DIVDP:
		{
			int dp_x1, dp_x2, dp_dest;

			error = commandSpec.GetParameters(command_fields, dp_x1, dp_x2, dp_dest);

			if (!error) {

				error = dpArr.divide_dparr(dp_x1, dp_x2, dp_dest);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_DOTPRODDP:
		{
			int dp_x1, dp_x2;

			error = commandSpec.GetParameters(command_fields, dp_x1, dp_x2);

			if (!error) {

				double value;

				error = dpArr.dotproduct_dparr(dp_x1, dp_x2, &value);

				if (!error & verbose) BD.DisplayConsoleMessage("dot product = " + ToString(value));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_MINMAX:
		{
			int dp_source;

			error = commandSpec.GetParameters(command_fields, dp_source);

			if (!error) {

				DBL2 min_max;
				INT2 min_max_indexes;

				error = dpArr.get_min_max(dp_source, &min_max, &min_max_indexes);

				if(verbose && !error) BD.DisplayConsoleMessage("Minimum = " + ToString(min_max.x) + " at index " + ToString(min_max_indexes.x) + ". Maximum = " + ToString(min_max.y) + " at index " + ToString(min_max_indexes.y));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(min_max.i, min_max_indexes.i, min_max.j, min_max_indexes.j));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_MEAN:
		{
			int dp_source;

			error = commandSpec.GetParameters(command_fields, dp_source);

			if (!error) {

				DBL2 mean_stdev;

				error  = dpArr.get_mean(dp_source, &mean_stdev);

				if (verbose && !error) BD.DisplayConsoleMessage("Mean = " + ToString(mean_stdev.x) + " +/- " + ToString(mean_stdev.y));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(mean_stdev));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_GETAMPLI:
		{
			int dp_source, pointsPeriod;

			error = commandSpec.GetParameters(command_fields, dp_source, pointsPeriod);

			if (!error) {

				double amplitude;

				error = dpArr.get_amplitude(dp_source, pointsPeriod, &amplitude);

				if (verbose && !error) BD.DisplayConsoleMessage("Max Amplitude = " + ToString(amplitude));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(amplitude));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;


		case CMD_DP_LINREG:
		{
			int dp_x, dp_y;

			int dp_z, dp_out;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y, dp_z, dp_out);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_x, dp_y); dp_z = -1; dp_out = -1; }

			if (!error) {

				DBL2 gradient, intercept;

				error = dpArr.linreg(dp_x, dp_y, dp_z, dp_out, &gradient, &intercept);

				if (verbose && !error) BD.DisplayConsoleMessage("y = gx + c, where g = " + ToString(gradient.x) + " +/- " + ToString(gradient.y) + " and c = " + ToString(intercept.x) + " +/- " + ToString(intercept.y));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(gradient, intercept));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_COERCIVITY:
		{
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				DBL3 Hc_up, Hc_dn;
				error = dpArr.get_coercivity(dp_x, dp_y, &Hc_up, &Hc_dn);

				if (verbose && !error) BD.DisplayConsoleMessage("Hc_up = " + ToString(Hc_up.x) + " -/+ (" + ToString(Hc_up.y) + ", "  + ToString(Hc_up.z) + ") A/m. " +
																"Hc_dn = " + ToString(Hc_dn.x) + " -/+ (" + ToString(Hc_dn.y) + ", " + ToString(Hc_dn.z) + ") A/m. ");

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Hc_up, Hc_dn));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_REMANENCE:
		{
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				double Mr_up, Mr_dn;
				error = dpArr.get_remanence(dp_x, dp_y, &Mr_up, &Mr_dn);

				if (verbose && !error) BD.DisplayConsoleMessage("Mr_up = " + ToString(Mr_up) + " A/m. " + "Mr_dn = " + ToString(Mr_dn) + " A/m.");

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Mr_up, Mr_dn));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_COMPLETEHYSTERLOOP:
		{
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				error = dpArr.complete_hysteresis(dp_x, dp_y);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_DUMPTDEP:
		{
			string meshName, paramName;
			int dp_arr;
			double max_temperature;

			error = commandSpec.GetParameters(command_fields, meshName, paramName, max_temperature, dp_arr);

			if (!error) {

				error = dpArr.dump_tdep(&SMesh, meshName, paramName, max_temperature, dp_arr);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_FITLORENTZ:
		{
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				DBL2 S, H0, dH, y0;

				error = dpArr.fit_lorentz(dp_x, dp_y, &S, &H0, &dH, &y0);

				if (verbose && !error) {

					BD.DisplayConsoleMessage(
						"S = " + ToString(S.major) + " +/- " + ToString(S.minor) + ", " +
						"H0 = " + ToString(H0.major) + " +/- " + ToString(H0.minor) + ", " +
						"dH = " + ToString(dH.major) + " +/- " + ToString(dH.minor) + ", " +
						"y0 = " + ToString(y0.major) + " +/- " + ToString(y0.minor)
						);
				}

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(S.major, H0.major, dH.major, y0.major, S.minor, H0.minor, dH.minor, y0.minor));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_FITSKYRMION:
		{
			int dp_x, dp_y;
			double w_val = 0.0;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y, w_val);
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, dp_x, dp_y); w_val = 0.0; }

			if (!error) {

				DBL2 R, Ms, w;

				if (IsNZ(w_val)) w.major = w_val;

				error = dpArr.fit_skyrmion(dp_x, dp_y, &R, &Ms, &w);

				if (verbose && !error) {

					BD.DisplayConsoleMessage(
						"R = " + ToString(R.major) + " +/- " + ToString(R.minor) + ", " +
						"Ms = " + ToString(Ms.major) + " +/- " + ToString(Ms.minor) + ", " +
						"w = " + ToString(w.major) + " +/- " + ToString(w.minor)
					);
				}

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(R.major, Ms.major, w.major, R.minor, Ms.minor, w.minor));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_REPLACEREPEATS:
		{
			int dp_in, dp_out;

			error = commandSpec.GetParameters(command_fields, dp_in, dp_out);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in); dp_out = dp_in; }

			if (!error) {

				error = dpArr.replace_repeats(dp_in, dp_out);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_REMOVEOFFSET:
		{
			int dp_in, dp_out;

			error = commandSpec.GetParameters(command_fields, dp_in, dp_out);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in); dp_out = dp_in; }

			if (!error) {

				error = dpArr.remove_offset(dp_in, dp_out);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CARTESIANTOPOLAR:
		{
			int dp_in_x, dp_in_y, dp_out_r, dp_out_theta;

			error = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_r, dp_out_theta);
			if (error == BERROR_PARAMMISMATCH) { 
				
				error.reset() = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y); 
				dp_out_r = dp_in_x; 
				dp_out_theta = dp_in_y;
			}

			if (!error) {

				error = dpArr.Cartesian_to_Polar(dp_in_x, dp_in_y, dp_out_r, dp_out_theta);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SMOOTH:
		{
			int dp_in, dp_out, window_size;

			error = commandSpec.GetParameters(command_fields, dp_in, dp_out, window_size);

			if (!error) {

				error = dpArr.adjacent_averaging(dp_in, dp_out, window_size);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		//---------------- finished CMD_DP_ commands

		case CMD_TEST:
		{
			/*
			vector<string> commands_output;
			commands_output.resize(commands.size());

			for (int idx = 0; idx < commands.size(); idx++) {

				commands_output[idx] = commands.get_key_from_index(idx);
			}

			quicksort(commands_output);

			
			string commands_description;

			for (int idx = 0; idx < commands.size(); idx++) {

				commands_description += commands_output[idx] + "\n\n";

				commands_description += commands[commands_output[idx]].usage + "\n\n";
				commands_description += commands[commands_output[idx]].descr + "\n";
				
				if (commands[commands_output[idx]].return_descr.length()) {

					commands_description += commands[commands_output[idx]].return_descr + "\n\n";
				}
				else commands_description += "\n";
			}

			commands_description = trim(commands_description, "[tc0,0.5,0,1/tc]");
			commands_description = trim(commands_description, "[tc0,0.5,0.5,1/tc]");
			commands_description = trim(commands_description, "<b>");
			commands_description = trim(commands_description, "</b>");
			commands_description = trim(commands_description, "<i>");
			commands_description = trim(commands_description, "</i>");

			SaveTextToFile("c:/commands.txt", commands_description);
			*/

			BD.DisplayConsoleMessage(ToString(SMesh[SMesh.GetMeshFocus()]->n));
			
		}
		break;

		default:
			break;
		}

		if(error) err_hndl.handle_error(error);
	}
	else if (verbose) BD.DisplayConsoleError("ERROR : command not recognized.");
}