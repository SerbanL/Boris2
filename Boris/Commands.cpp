#include "stdafx.h"
#include "Simulation.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Scripted communication

//disable / enable script server
void Simulation::Script_Server_Control(bool status)
{
	if (status && !is_thread_running(THREAD_NETWORK)) {

		//start window sockets thread to listen for incoming messages
		infinite_loop_launch(&Simulation::Listen_Incoming_Message, THREAD_NETWORK);

		BD.DisplayConsoleMessage("Script server running.");
	}
	else if (!status && is_thread_running(THREAD_NETWORK)) {

		stop_thread(THREAD_NETWORK);

		BD.DisplayConsoleMessage("Script server stopped.");
	}
}

void Simulation::Listen_Incoming_Message(void) 
{
	//Listen for incoming messages - non-blocking call, but this method runs on a thread with an infinite loop which keeps calling Simulation::Listen_Incoming_Message
	string message = commSocket.Listen();

	if (message.length()) {

		//handle command using the usual HandleCommand method, but on a blocking thread call - need to wait for it to finish before responding to client
		set_blocking_thread(THREAD_HANDLEMESSAGE);
		//Keep trying to launch the command handler : if THREAD_HANDLEMESSAGE is busy or has not yet been marked as not active (this is possible even if the previous call to THREAD_HANDLEMESSAGE has finished and we returned parameters to client - it can take much longer than usual sometimes),
		//we receive THREAD_GENERATENEW response back until THREAD_HANDLEMESSAGE is available (marked as not active).
		//If the command handler was launched and completed successfully, the single_call_launch method will respond with THREAD_HANDLEMESSAGE after it returns.
		while (single_call_launch<string>(&Simulation::HandleCommand, message, THREAD_HANDLEMESSAGE) != THREAD_HANDLEMESSAGE);

		//the command was handled using a blocking thread, so now parameters should be available to send
		commSocket.SendDataParams();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Command handler

void Simulation::HandleCommand(string command_string) 
{
	//cannot change simulation parameters in the middle of an interation. simulationMutex also used by Simulate() method, but there it is non-blocking
	simulationMutex.lock();
	lock_guard<mutex> simlock(simulationMutex, adopt_lock);

	//ensure we do not handle a command during a display refresh
	//This method locks and unlocks the display mutex (cannot leave it locked and unlock it at the end of HandleCommand, since HandleCommand can call other BD methods)
	//It is important to do this since HandleCommand can cause changes to data which are used during a display refresh and could result in a crash.
	//It is theoretically possible a BD method is called right after WaitForDisplayEnd but before HandleCommand returns - this could happen in response to a user action only.
	//This is virtually impossible however, as the user would have to issue commands quicker than it takes for the program to make changes to critical data.
	//Not ideal solution, the only foolproof solution is to lock the display mutex and unlock it before calling any BD methods here (or before returning), but after changes to data have been made.
	//This is messy so I'd rather not go down that route since it only guards against an extremely unlikely scenario.
	BD.WaitForDisplayEnd();

	//display messages in console only if this is true. Note, some important messages are displayed regardless (e.g. simulation stopped or simulation started messages.)
	bool verbose = true;

	//if command came from a script client, send back required parameters. Each command must ensure it does this.
	bool script_client_connected = false;

	//commands are always of the form: commandname (param1) (param2) ...
	vector<string> command_fields = split(command_string, " ");
	string command_name = command_fields.front();

	if (!command_name.length()) {

		err_hndl.show_error(BERROR_COMMAND_NOTRECOGNIZED, verbose);
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
		else err_hndl.show_error(BERROR_COMMAND_NOTRECOGNIZED, verbose);

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
			else if (verbose) {

				PrintCommandUsage(command_name);
				BD.DisplayConsoleListing("iterupdate: " + ToString(iterUpdate));
			}

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
					//update rectangles in dataBoxList for the named mesh
					UpdateDataBoxEntries(meshRect_old, SMesh[meshName]->GetMeshRect(), meshName);
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

				if (!err_hndl.call(&MeshBase::SetMeshCellsize, SMesh.active_mesh(), h)) {

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

				if (!err_hndl.call(&MeshBase::SetMeshECellsize, SMesh.active_mesh(), h_e)) {

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

				if (!err_hndl.call(&MeshBase::SetMeshTCellsize, SMesh.active_mesh(), h_t)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshTCellsize()));
		}
		break;

		case CMD_SCELLSIZE:
		{
			DBL3 h_s;
			error = commandSpec.GetParameters(command_fields, h_s);

			if (!error && !SMesh.active_mesh()->is_atomistic()) {

				StopSimulation();

				if (!err_hndl.call(&Mesh::SetMeshSCellsize, dynamic_cast<Mesh*>(SMesh.active_mesh()), h_s, false)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected && !SMesh.active_mesh()->is_atomistic()) {

				commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh.active_mesh())->GetMeshSCellsize()));
			}
		}
		break;

		case CMD_MCELLSIZE:
		{
			DBL3 h_m;
			error = commandSpec.GetParameters(command_fields, h_m);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&MeshBase::SetMeshMCellsize, SMesh.active_mesh(), h_m)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshMCellsize()));
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

		case CMD_ATOMDMCELLSIZE:
		{
			DBL3 h_dm;
			error = commandSpec.GetParameters(command_fields, h_dm);

			if (!error) {

				StopSimulation();

				Atom_Mesh *pMesh = dynamic_cast<Atom_Mesh*>(SMesh.active_mesh());

				if (pMesh) {

					if (!err_hndl.call(&Atom_Mesh::Set_Demag_Cellsize, pMesh, h_dm)) {

						UpdateScreen();
					}
				}
				else if (verbose) error(BERROR_NOTATOMISTIC);
			}
			else if (verbose) Print_Speedup_List();

			if (script_client_connected) {

				Atom_Mesh *pMesh = dynamic_cast<Atom_Mesh*>(SMesh.active_mesh());

				if (pMesh) {

					commSocket.SetSendData(commandSpec.PrepareReturnParameters(pMesh->Get_Demag_Cellsize()));
				}
			}
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

		case CMD_ADDAFMESH:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_ANTIFERROMAGNETIC, meshRect)) {

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

		case CMD_ADDDIAMAGNETICMESH:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_DIAMAGNETIC, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDAMESHCUBIC:
		{
			string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, meshName, MESH_ATOM_CUBIC, meshRect)) {

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
				function<vector<unsigned char>(string, INT2)> bitmap_loader = [&](string fileName, INT2 n_plane) -> vector<unsigned char> {

					vector<unsigned char> bitmap;
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

		case CMD_INDIVIDUALMASKSHAPE:
		{
			int status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				shape_change_individual = (bool)status;

				UpdateScreen();
			}
			else if (verbose) Print_IndividualShapeStatus();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(shape_change_individual));
		}
		break;

		case CMD_SETANGLE:
		{
			double polar, azim;
			string meshName;

			error = commandSpec.GetParameters(command_fields, polar, azim, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, polar, azim); meshName = ""; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetMagAngle, &SMesh, meshName, polar, azim)) {

					UpdateScreen();
				}

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RANDOM:
		{
			string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetRandomMag, &SMesh, meshName)) {

					UpdateScreen();
				}

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_INVERTMAG:
		{
			string meshName;
			string inputString;
			bool x = true, y = true, z = true;

			error = commandSpec.GetParameters(command_fields, inputString);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (inputString.length()) {

					vector<string> fields = split(inputString, " ");

					int num_options;

					if (fields.back() != "x" && fields.back() != "y" && fields.back() != "z") {

						meshName = fields.back();
						num_options = fields.size() - 1;
					}
					else {

						meshName = SMesh.GetMeshFocus();
						num_options = fields.size();
					}

					if (num_options > 0) {

						x = false; 
						y = false;
						z = false;

						for (int idx = 0; idx < num_options; idx++) {

							if (fields[idx] == "x") x = true;
							if (fields[idx] == "y") y = true;
							if (fields[idx] == "z") z = true;
						}
					}
				}
				else meshName = SMesh.GetMeshFocus();
				
				if (!err_hndl.qcall(&SuperMesh::SetInvertedMag, &SMesh, meshName, x, y, z)) {

					UpdateScreen();
				}

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_MIRRORMAG:
		{
			string meshName;
			string axis;

			error = commandSpec.GetParameters(command_fields, axis, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, axis); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetMirroredMag, &SMesh, meshName, axis)) {

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

				if (!err_hndl.qcall(&SuperMesh::SetMagAngle_Rect, &SMesh, meshName, polar, azim, rectangle)) {

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

				if (!err_hndl.qcall(&SuperMesh::SetMagDomainWall, &SMesh, meshName, longitudinal, transverse, width, position)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_VORTEX:
		{
			int longitudinal, rotation, core;
			Rect rect;
			string meshName;

			error = commandSpec.GetParameters(command_fields, longitudinal, rotation, core, rect, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, longitudinal, rotation, core, rect); meshName = SMesh.GetMeshFocus(); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, longitudinal, rotation, core); rect = Rect(); meshName = SMesh.GetMeshFocus(); }

			if (!error && longitudinal && rotation && core) {

				StopSimulation();

				string vortexFile;
				bool invertMag = false;

				//l = +1, r = -1, c = -1 : vortex_hh_clk_dn
				//l = +1, r = -1, c = +1 : vortex_hh_clk_up
				//l = +1, r = +1, c = -1 : vortex_hh_cclk_dn
				//l = +1, r = +1, c = +1 : vortex_hh_cclk_up

				//l = -1, r = -1, c = -1 : vortex_hh_cclk_up, then invert magnetisation.
				//l = -1, r = -1, c = +1 : vortex_hh_cclk_dn, then invert magnetisation.
				//l = -1, r = +1, c = -1 : vortex_hh_clk_up, then invert magnetisation.
				//l = -1, r = +1, c = +1 : vortex_hh_clk_dn, then invert magnetisation.

				if (longitudinal > 0) {

					if (rotation < 0 && core < 0) vortexFile = "vortex_hh_clk_dn.ovf";
					else if (rotation < 0 && core > 0) vortexFile = "vortex_hh_clk_up.ovf";
					else if (rotation > 0 && core < 0) vortexFile = "vortex_hh_cclk_dn.ovf";
					else vortexFile = "vortex_hh_cclk_up.ovf";
				}
				else {

					invertMag = true;

					if (rotation < 0 && core < 0) vortexFile = "vortex_hh_cclk_up.ovf";
					else if (rotation < 0 && core > 0) vortexFile = "vortex_hh_cclk_dn.ovf";
					else if (rotation > 0 && core < 0) vortexFile = "vortex_hh_clk_up.ovf";
					else vortexFile = "vortex_hh_clk_dn.ovf";
				}

				VEC<DBL3> data;

				OVF2 ovf2;
				error = ovf2.Read_OVF2_VEC(GetExeDirectory() + vortexFile, data);

				if (!error) {

					//data loaded correctly, so resize currently focused mesh (if ferromagnetic) then copy magnetisation data to it.
					if (SMesh.active_mesh()->is_atomistic()) {

						data.renormalize(dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->mu_s.get0());
					}
					else {

						data.renormalize(dynamic_cast<Mesh*>(SMesh.active_mesh())->Ms.get0());
					}

					if (invertMag) data *= -1.0;

					if (SMesh.active_mesh()->Magnetism_Enabled()) {

						SMesh.active_mesh()->SetMagFromData(data, rect);
						UpdateScreen();
					}
					else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
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

		case CMD_PBC:
		{
			string meshName;
			string flag;
			int images;

			error = commandSpec.GetParameters(command_fields, meshName, flag, images);

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::Set_PBC, &SMesh, meshName, flag, images)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_PBC();
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

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&ZeemanBase::GetField)));
		}
		break;

		case CMD_SETSTRESS:
		{
			DBL3 stress_polar;
			string meshName = "";

			error = commandSpec.GetParameters(command_fields, stress_polar, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, stress_polar); meshName = ""; }

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetUniformStress, &SMesh, meshName, Polar_to_Cartesian(stress_polar))) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&MElastic::GetUniformStress)));
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
			int status;

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

		case CMD_EXCLUDEMULTICONVDEMAG:
		{
			bool status;
			string meshName;

			error = commandSpec.GetParameters(command_fields, status, meshName);

			if (!error) {

				StopSimulation();

				error = SMesh.Set_Demag_Exclusion(status, meshName);

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

				ODE_ atom_odeID;
				SMesh.QueryAtomODE(atom_odeID);

				if (setOde != ODE_ERROR && odeEval != EVAL_ERROR && vector_contains(odeAllowedEvals(setOde), odeEval) && vector_contains(odeAllowedEvals(atom_odeID), odeEval)) {

					if (!err_hndl.call(&SuperMesh::SetODE, &SMesh, setOde, odeEval)) {

						UpdateScreen();
					}
				}
				else if (verbose) error(BERROR_INCORRECTCONFIG);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETATOMODE:
		{
			string odeHandle, odeEvalHandle;

			error = commandSpec.GetParameters(command_fields, odeHandle, odeEvalHandle);

			if (!error) {

				StopSimulation();

				ODE_ setatom_Ode = (ODE_)atom_odeHandles.get_ID_from_value(odeHandle);
				EVAL_ odeEval = (EVAL_)odeEvalHandles.get_ID_from_value(odeEvalHandle);

				ODE_ odeID;
				SMesh.QueryODE(odeID);

				if (setatom_Ode != ODE_ERROR && odeEval != EVAL_ERROR && vector_contains(odeAllowedEvals(setatom_Ode), odeEval) && vector_contains(odeAllowedEvals(odeID), odeEval)) {

					if (!err_hndl.call(&SuperMesh::SetAtomisticODE, &SMesh, setatom_Ode, odeEval)) {

						UpdateScreen();
					}
				}
				else if (verbose) error(BERROR_INCORRECTCONFIG);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETODEEVAL:
		{
			string odeEvalHandle;

			error = commandSpec.GetParameters(command_fields, odeEvalHandle);

			if (!error) {

				StopSimulation();

				EVAL_ odeEval = (EVAL_)odeEvalHandles.get_ID_from_value(odeEvalHandle);

				ODE_ odeID, atom_odeID;
				SMesh.QueryODE(odeID);
				SMesh.QueryAtomODE(atom_odeID);

				if (odeEval != EVAL_ERROR && vector_contains(odeAllowedEvals(odeID), odeEval) && vector_contains(odeAllowedEvals(atom_odeID), odeEval)) {

					if (!err_hndl.call(&SuperMesh::SetODEEval, &SMesh, odeEval)) {

						UpdateScreen();
					}
				}
				else if (verbose) error(BERROR_INCORRECTCONFIG);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EVALSPEEDUP:
		{
			int status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();

				SMesh.SetEvaluationSpeedup(status);

				UpdateScreen();
			}
			else if (verbose) Print_Speedup_List();
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

		case CMD_ASTEPCTRL:
		{
			double err_fail, err_high, err_low, dT_incr, dT_min, dT_max;

			error = commandSpec.GetParameters(command_fields, err_fail, err_high, err_low, dT_incr, dT_min, dT_max);

			if (!error) {

				SMesh.SetAdaptiveTimeStepCtrl(err_fail, err_high, err_low, dT_incr, dT_min, dT_max);
				UpdateScreen();
			}
			else if (verbose) Print_AStepCtrl();
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

		case CMD_SETDATA:
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

					//edit first data entry (there's always at least one)
					EditSaveDataEntry(0, (DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect);
					//get rid of all other data entries
					while(saveDataList.size() > 1) saveDataList.erase(1);

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
			Rect dataRect;

			error = commandSpec.GetParameters(command_fields, dataName, meshName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName, meshName); dataRect = Rect(); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dataName); meshName = SMesh.GetMeshFocus(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				if (dataDescriptor[dataName].meshless) {

					NewDataBoxField(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName)));
				}
				else {

					NewDataBoxField(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect));
				}

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

				if (directory.substr(directory.length() - 1) != "/" && directory.substr(directory.length() - 1) != "/") {

					directory += "/";
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

		case CMD_SETSTAGE:
		{
			string stageTypeName;
			string meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, stageTypeName, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, stageTypeName); meshName = SMesh.GetMeshFocus(); }

			if (!error && stageDescriptors.has_key(stageTypeName) && (SMesh.contains(meshName) || meshName == SMesh.superMeshHandle)) {

				//edit stage 0 to new settings (there's always at least one stage)
				EditStageType(0, (SS_)stageDescriptors.get_ID_from_key(stageTypeName), meshName);
				//get rid of all other stages
				while(simStages.size() > 1) DeleteStage(1);

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
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName); set_value = false; }

			if (!error && set_value) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::set_meshparam_value, &SMesh, meshName, paramName, paramValue)) {

					UpdateScreen();
				}
			}
			else if (verbose && set_value) PrintCommandUsage(command_name);

			if (script_client_connected && !error && !set_value) {

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
			string meshName, paramName;

			error = commandSpec.GetParameters(command_fields, meshName, paramName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName); paramName = ""; }

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				meshName = "";
			}

			StopSimulation();

			if (!err_hndl.qcall(&SuperMesh::clear_meshparam_temp, &SMesh, meshName, paramName)) {

				UpdateScreen();
			}
		}
		break;

		case CMD_SETPARAMTEMPEQUATION:
		{
			string meshName, paramName, equationText;

			error = commandSpec.GetParameters(command_fields, meshName, paramName, equationText);

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::set_meshparam_t_equation, &SMesh, meshName, paramName, equationText)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETPARAMTEMPARRAY:
		{
			string meshName, paramName, fileName;

			error = commandSpec.GetParameters(command_fields, meshName, paramName, fileName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName); fileName = ""; }

			//setting parameter with arrays loaded form file
			if (!error) {

				//load from a file directly if given
				StopSimulation();

				if (fileName.length()) {

					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

					if (GetFileTermination(fileName) != ".txt")
						fileName += ".txt";

					vector<vector<double>> load_arrays;
					if (ReadDataColumns(fileName, "\t", load_arrays, { 0, 1, 2, 3 })) {

						vector<double> empty;

						if (load_arrays.size() == 4 && load_arrays[3].size()) error = SMesh.set_meshparam_tscaling_array(meshName, paramName, load_arrays[0], load_arrays[1], load_arrays[2], load_arrays[3], fileName);
						else if (load_arrays.size() == 3 && load_arrays[2].size()) error = SMesh.set_meshparam_tscaling_array(meshName, paramName, load_arrays[0], load_arrays[1], load_arrays[2], empty, fileName);
						else if (load_arrays.size() == 2 && load_arrays[1].size()) error = SMesh.set_meshparam_tscaling_array(meshName, paramName, load_arrays[0], load_arrays[1], empty, empty, fileName);
						else error(BERROR_COULDNOTLOADFILE);

						UpdateScreen();
					}
					else error(BERROR_COULDNOTLOADFILE);
				}
				else {
					
					//default array
					vector<double> temp(2), scaling(2), empty;
					temp = { 0.0, 1.0 };
					scaling = { 1.0, 1.0 };
					error = SMesh.set_meshparam_tscaling_array(meshName, paramName, temp, scaling, empty, empty);
					
					UpdateScreen();
				}
			}

			//if the above failed then try to load array from dp_arrays
			if(error == BERROR_COULDNOTLOADFILE) {

				int dp_T_idx, dp_c_ix, dp_c_iy, dp_c_iz;

				error.reset() = commandSpec.GetParameters(command_fields, paramName, dp_T_idx, dp_c_ix, dp_c_iy, dp_c_iz);
				if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_T_idx, dp_c_ix, dp_c_iy); dp_c_iz = -1; }
				if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_T_idx, dp_c_ix); dp_c_iy = -1; }

				if (!error && GoodIdx(dpArr.size(), dp_T_idx, dp_c_ix)) {

					//load from dp arrays
					StopSimulation();

					vector<double> empty;

					if (GoodIdx(dpArr.size(), dp_c_iy, dp_c_iz)) error = SMesh.set_meshparam_tscaling_array(SMesh.GetMeshFocus(), paramName, dpArr[dp_T_idx], dpArr[dp_c_ix], dpArr[dp_c_iy], dpArr[dp_c_iz]);
					else if (GoodIdx(dpArr.size(), dp_c_iy)) error = SMesh.set_meshparam_tscaling_array(SMesh.GetMeshFocus(), paramName, dpArr[dp_T_idx], dpArr[dp_c_ix], dpArr[dp_c_iy], empty);
					else error = SMesh.set_meshparam_tscaling_array(SMesh.GetMeshFocus(), paramName, dpArr[dp_T_idx], dpArr[dp_c_ix], empty, empty);

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

				//generatorName can be specified as equation, but typically this would not be specified and simply generatorName has the equation text
				//in this case generatorName is not a key in vargenerator_descriptor - this is how we recognize this type of input
				if (vargenerator_descriptor.get_ID_from_key(generatorName) == MATPVAR_EQUATION || !vargenerator_descriptor.has_key(generatorName)) {

					string equationText;

					if (!vargenerator_descriptor.has_key(generatorName)) {

						//generatorName not a key in vargenerator_descriptor : this is the equation text, or part of it (it's possible the equation text is split between generatorName and generatorArgs if there was a space character
						equationText = generatorName + " " + generatorArgs;
					}
					else {

						//default equation if none specified
						if (!generatorArgs.length()) generatorArgs = "1";

						//generatorName specifically has the key corresponding to MATPVAR_EQUATION, thus generatorArgs holds the equation text.
						equationText = generatorArgs;
					}

					if (!err_hndl.qcall(&SuperMesh::set_meshparam_s_equation, &SMesh, meshName, paramName, equationText)) {

						UpdateScreen();
					}
				}
				else {

					//used with custom parameter variation generator (set from grayscale png file)
					function<vector<unsigned char>(string, INT2)> bitmap_loader = [&](string fileName, INT2 n_plane) -> vector<unsigned char> {

						vector<unsigned char> bitmap;
						BD.BGMethods()->GetBitmapFromImage(fileName, bitmap, n_plane);

						return bitmap;
					};

					error = SMesh.set_meshparam_var(meshName, paramName, generatorName, generatorArgs, bitmap_loader);
				}

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

			if (verbose) BD.DisplayConsoleMessage("Saving simulation ... please wait.");

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
			LoadSimulation(GetUserDocumentsPath() + boris_data_directory + boris_simulations_directory + "default");
			BD.DisplayConsoleMessage("Default state restored.");
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

		case CMD_DISPLAYBACKGROUND:
		{
			string name, meshName;

			error = commandSpec.GetParameters(command_fields, name, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, name); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				MESHDISPLAY_ display = (MESHDISPLAY_)displayHandles.get_ID_from_value(name);

				if (!err_hndl.qcall(&SuperMesh::SetDisplayedBackgroundPhysicalQuantity, &SMesh, meshName, (int)display)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_VECREP:
		{
			string meshName;
			int vecreptype;

			error = commandSpec.GetParameters(command_fields, meshName, vecreptype);

			if (!error) {

				if (!err_hndl.qcall(&SuperMesh::SetVEC3Rep, &SMesh, meshName, (int)vecreptype)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_DISPLAYTRANSPARENCY:
		{
			double foreground, background;

			error = commandSpec.GetParameters(command_fields, foreground, background);

			if (!error) {

				displayTransparency = DBL2(foreground, background);

				UpdateScreen();
			}
			else if (verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_DISPLAYTHRESHOLDS:
		{
			double minimum, maximum;

			error = commandSpec.GetParameters(command_fields, minimum, maximum);

			if (!error) {

				displayThresholds = DBL2(minimum, maximum);

				UpdateScreen();
			}
			else if (verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_DISPLAYTHRESHOLDTRIGGER:
		{
			int trigtype;

			error = commandSpec.GetParameters(command_fields, trigtype);

			if (!error) {

				displayThresholdTrigger = trigtype;

				UpdateScreen();
			}
			else if (verbose) Print_MeshDisplay_List();
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

			if (BD.SaveMeshImage(fileName, image_cropping)) {

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

				if (BD.BGMethods()->MakeVideoFromFileSequence(sequenceDirectory + fileNameBase + ".wmv", fileNames, unsigned(fps), 1.0, quality)) {

					if(verbose) BD.DisplayConsoleMessage("Video created.");
				}
				else if (verbose) error(BERROR_COULDNOTSAVEFILE);
			}
			else if(verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_IMAGECROPPING:
		{
			double left, bottom, right, top;

			error = commandSpec.GetParameters(command_fields, left, bottom, right, top);

			if (!error) {

				image_cropping = DBL4(left, bottom, right, top);

				UpdateScreen();
			}
			else if (verbose) {

				string text = "[tc1,1,1,1/tc]Image cropping settings : " + MakeIO(IOI_IMAGECROPPING);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}
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

		case CMD_CLEARMOVINGMESH:
		{
			StopSimulation();
			SMesh.ClearMovingMesh();
			UpdateScreen();
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

		case CMD_EXCHANGECOUPLEDMESHES:
		{
			bool status;
			string meshName;

			error = commandSpec.GetParameters(command_fields, status, meshName);
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, status); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				SMesh.Set_ExchangeCoupledMeshes(status, meshName);

				UpdateScreen();
			}
			else if (verbose) Print_ExchangeCoupledMeshes_List();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->GetMeshExchangeCoupling()));
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

				SMesh.CallModuleMethod(&STransport::SetPotential, 0.0, true);

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

					SMesh.CallModuleMethod(&STransport::SetPotential, potential, true);

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

					SMesh.CallModuleMethod(&STransport::SetCurrent, current, true);

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

		case CMD_STATICTRANSPORTSOLVER:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				bool status;

				error = commandSpec.GetParameters(command_fields, status);

				if (!error) {

					static_transport_solver = status;

					UpdateScreen();
				}
				else if (verbose) PrintTransportSolverConfig();

				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(static_transport_solver));
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
			else if (verbose) {

				string heatdT = "[tc1,1,1,1/tc]Heat Equation Time Step: " + MakeIO(IOI_HEATDT) + "</c>";
				BD.DisplayFormattedConsoleMessage(heatdT);
			}

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
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&HeatBase::GetAmbientTemperature)));
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
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->CallModuleMethod(&HeatBase::GetAlphaBoundary)));
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

					vector<bool> insulating = SMesh.active_mesh()->CallModuleMethod(&HeatBase::GetInsulatingSides);
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

			if (script_client_connected && !SMesh.active_mesh()->is_atomistic())
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh.active_mesh())->GetCurieTemperature()));
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

			if (script_client_connected && !SMesh.active_mesh()->is_atomistic())
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh.active_mesh())->GetCurieTemperatureMaterial()));
		}
		break;

		case CMD_ATOMICMOMENT:
		{
			double atomic_moment;
			DBL2 atomic_moment_AFM;
			string meshName;

			error = commandSpec.GetParameters(command_fields, atomic_moment_AFM, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, atomic_moment, meshName); atomic_moment_AFM = DBL2(atomic_moment); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, atomic_moment_AFM); meshName = SMesh.superMeshHandle; }

			if (!error) {

				StopSimulation();

				if (!err_hndl.qcall(&SuperMesh::SetAtomicMagneticMoment, &SMesh, meshName, atomic_moment_AFM)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected) {

				if (SMesh.active_mesh()->GetMeshType() == MESH_ANTIFERROMAGNETIC) 
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh.active_mesh())->GetAtomicMoment_AFM()));
				else if (!SMesh.active_mesh()->is_atomistic()) 
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh.active_mesh())->GetAtomicMoment()));
			}
		}
		break;

		case CMD_TAU:
		{
			DBL2 tau_ii;
			string entry;
			DBL2 tau_ij;
			string meshName;

			error = commandSpec.GetParameters(command_fields, tau_ii, entry);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, tau_ii); meshName = SMesh.superMeshHandle; tau_ij = DBL2(-1.0); }

			if (!error) {

				vector<string> fields = split(entry, " ");

				if (fields.size() == 3) {

					//both tau_ij and meshname specified (at least can assume to be attempted)
					meshName = fields[2];
					tau_ij = ToNum(fields[0] + " " + fields[1]);
				}
				else if (fields.size() == 2) {

					//only attempted tau_ij entry
					meshName = SMesh.superMeshHandle;
					tau_ij = ToNum(fields[0] + " " + fields[1]);
				}
				else if (fields.size() == 1) {

					//only attempted meshName entry
					meshName = fields[0];
					tau_ij = DBL2(-1.0);
				}

				StopSimulation();

				if (tau_ij >= DBL2(0.0)) {

					//set both intra and inter terms
					if (!err_hndl.qcall(&SuperMesh::SetTcCoupling, &SMesh, meshName, tau_ii, tau_ij)) {

						UpdateScreen();
					}
				}

				//only intra terms
				else if (!err_hndl.qcall(&SuperMesh::SetTcCoupling_Intra, &SMesh, meshName, tau_ii)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected && !SMesh.active_mesh()->is_atomistic()) {

				commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh.active_mesh())->GetTcCoupling()));
			}
		}
		break;

		case CMD_TMODEL:
		{
			int tmodeltype;
			string meshName;

			error = commandSpec.GetParameters(command_fields, tmodeltype, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, tmodeltype); meshName = SMesh.GetMeshFocus(); }

			if (!error) {

				StopSimulation();

				SMesh.SetTemperatureModel(meshName, tmodeltype);
				
				UpdateScreen();
			}
			else if (verbose) Print_TemperatureModel_List();
		}
		break;

		case CMD_STOCHASTIC:
		{
			if (verbose) Print_Stochasticity_List();
		}
		break;

		case CMD_LINKSTOCHASTIC:
		{
			bool flag;
			string meshName;

			error = commandSpec.GetParameters(command_fields, flag, meshName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, flag); meshName = SMesh.superMeshHandle; }

			if (!error) {

				StopSimulation();

				SMesh.SetLinkStochastic(flag, meshName);

				UpdateScreen();
			}
			else if (verbose) Print_Stochasticity_List();
		}
		break;

		case CMD_SETDTSTOCH:
		{
			double dTstoch;

			error = commandSpec.GetParameters(command_fields, dTstoch);

			if (!error) {

				StopSimulation();

				SMesh.SetStochTimeStep(dTstoch);
				UpdateScreen();
			}
			else if (verbose) Print_Stochasticity_List();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetStochTimeStep()));
		}
		break;

		case CMD_LINKDTSTOCHASTIC:
		{
			bool flag;

			error = commandSpec.GetParameters(command_fields, flag);

			if (!error) {

				StopSimulation();

				SMesh.SetLink_dTstoch(flag);
				UpdateScreen();
			}
			else if (verbose) Print_Stochasticity_List();
		}
		break;

		case CMD_SETDTSPEEDUP:
		{
			double dTspeedup;

			error = commandSpec.GetParameters(command_fields, dTspeedup);

			if (!error) {

				StopSimulation();

				SMesh.SetSpeedupTimeStep(dTspeedup);
				UpdateScreen();
			}
			else if (verbose) Print_Speedup_List();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetSpeedupTimeStep()));
		}
		break;
			
		case CMD_LINKDTSPEEDUP:
		{
			bool flag;

			error = commandSpec.GetParameters(command_fields, flag);

			if (!error) {

				StopSimulation();

				SMesh.SetLink_dTspeedup(flag);
				UpdateScreen();
			}
			else if (verbose) Print_Speedup_List();
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
			string directory = GetUserDocumentsPath() + boris_data_directory;
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

				if (SMesh.contains(meshName) && !SMesh[meshName]->is_atomistic()) {

					error = mdb.AddMDBEntry(materialName, *dynamic_cast<Mesh*>(SMesh[meshName]), SMesh[meshName]->GetMeshType());

					if (!error) BD.DisplayConsoleMessage("Material added to local database.");
				}
				else {

					if (!SMesh.contains(meshName)) err_hndl.show_error(BERROR_MESHNAMEINEXISTENT, verbose);
					else err_hndl.show_error(BERROR_INCORRECTNAME, verbose);
				}

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
						mdb.copy_parameters(*dynamic_cast<Mesh*>(SMesh[meshName]));
					}

					if (script_client_connected)
						commSocket.SetSendData({ meshName });
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
	
					if (!SMesh.active_mesh()->is_atomistic()) {

						//now copy parameter values - this is done even if mesh types are not the same. e.g. you might want to copy from a ferromagnetic mesh to a dipole mesh
						//this only works for micromagnetics meshes
						mdb.copy_parameters(*dynamic_cast<Mesh*>(SMesh.active_mesh()));
					}
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

			if (!err_hndl.call(&MaterialsDB::UpdateMDB, &mdb, domain_name, mdb_update_handler, &returnMessage)) {

				BD.DisplayConsoleMessage("Materials database successfully updated.");
			}
			else BD.DisplayConsoleError(returnMessage);
		}
		break;

		case CMD_SHOWLENGHTS:
		{
			if (verbose) {

				if (SMesh.active_mesh()->Magnetism_Enabled() && !SMesh.active_mesh()->is_atomistic()) {
					
					double A = dynamic_cast<Mesh*>(SMesh.active_mesh())->A;
					double Ms = dynamic_cast<Mesh*>(SMesh.active_mesh())->Ms;
					double Ku = dynamic_cast<Mesh*>(SMesh.active_mesh())->K1;
					double D = dynamic_cast<Mesh*>(SMesh.active_mesh())->D;

					string l_ex, l_Bloch("N/A"), l_sky("N/A");

					l_ex = ToString(sqrt(2 * A / (MU0 * Ms*Ms)), "m");

					if (IsNZ(Ku)) l_Bloch = ToString(sqrt(A / Ku), "m");
					if (IsNZ(Ku) && IsNZ(D)) l_sky = ToString(sqrt(PI * D / (4 * Ku)), "m");

					BD.DisplayConsoleMessage("l_ex = " + l_ex + ", l_Bloch = " + l_Bloch + ", l_sky = " + l_sky);
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
		}
		break;

		case CMD_SHOWMCELLS:
		{
			if (SMesh.active_mesh()->Magnetism_Enabled()) {

				if (verbose) BD.DisplayConsoleMessage("M discretisation cells : " + ToString(SMesh.active_mesh()->n));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.active_mesh()->n));
			}
			else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
		}
		break;

		case CMD_LOADOVF2MESH:
		{
			string params_string;

			error = commandSpec.GetParameters(command_fields, params_string);
			if (!error) {

				StopSimulation();

				double renormalize_value = 0.0;
				string fileName;

				//using split_numeric approach since the file name path can contain spaces.
				vector<string> entries = split_numeric(params_string);

				if (entries.size() == 2) {

					renormalize_value = ToNum(entries[0]);
					fileName = entries[1];
				}
				else fileName = params_string;

				if (GetFileTermination(fileName) != ".ovf")
					fileName += ".ovf";

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				VEC<DBL3> data;

				OVF2 ovf2;
				error = ovf2.Read_OVF2_VEC(fileName, data);

				if (!error) {

					//data loaded correctly, so make a ferromagnetic mesh and set M values from data

					if (IsNZ(renormalize_value)) data.renormalize(renormalize_value);

					string new_mesh_name = "permalloy";

					int meshnum = 0;

					while (SMesh.contains(new_mesh_name)) {

						meshnum++;
						new_mesh_name = string("permalloy_") + ToString(meshnum);
					}

					if (!err_hndl.call(&SuperMesh::AddMesh, &SMesh, new_mesh_name, MESH_FERROMAGNETIC, data.rect)) {

						if (!err_hndl.call(&MeshBase::SetMeshCellsize, SMesh[new_mesh_name], data.h)) {

							SMesh[new_mesh_name]->SetMagFromData(data);
						}
					}
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2MAG:
		{
			string params_string;

			error = commandSpec.GetParameters(command_fields, params_string);
			if (!error) {

				StopSimulation();

				double renormalize_value = 0.0;
				string fileName;

				//using split_numeric approach since the file name path can contain spaces.
				vector<string> entries = split_numeric(params_string);

				if (entries.size() == 2) {

					renormalize_value = ToNum(entries[0]);
					fileName = entries[1];
				}
				else fileName = params_string;

				if (GetFileTermination(fileName) != ".ovf")
					fileName += ".ovf";

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				VEC<DBL3> data;

				OVF2 ovf2;
				error = ovf2.Read_OVF2_VEC(fileName, data);

				if (!error) {

					//data loaded correctly, so resize currently focused mesh (if ferromagnetic) then copy magnetisation data to it.

					if (IsNZ(renormalize_value)) data.renormalize(renormalize_value);

					if (SMesh.active_mesh()->Magnetism_Enabled()) {

						SMesh.active_mesh()->SetMagFromData(data);
						UpdateScreen();
					}
					else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVEOVF2MAG:
		{
			string parameters;

			error = commandSpec.GetParameters(command_fields, parameters);

			if (!error) {

				if (SMesh.active_mesh()->Magnetism_Enabled() && !SMesh.active_mesh()->is_atomistic()) {

					bool normalize = false;
					string data_type = "bin8";
					string fileName;

					vector<string> fields = split(parameters, " ");

					int oparams = 0;

					if (fields[0] == "n" || fields[0] == "bin4" || fields[0] == "bin8" || fields[0] == "text") {

						oparams++;

						if (fields[0] == "n") normalize = true;
						else data_type = fields[0];
					}

					if (fields.size() > 1 && (fields[1] == "bin4" || fields[1] == "bin8" || fields[1] == "text")) {

						oparams++;
						data_type = fields[1];
					}

					fileName = combine(subvec(fields, oparams), " ");

					if (GetFileTermination(fileName) != ".ovf")
						fileName += ".ovf";

					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

					double Ms0 = dynamic_cast<Mesh*>(SMesh.active_mesh())->Ms.get0();
					if (!normalize) Ms0 = 1.0;

					OVF2 ovf2;
					error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Mesh*>(SMesh.active_mesh())->M, data_type, Ms0);
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVEOVF2:
		{
			string parameters;

			error = commandSpec.GetParameters(command_fields, parameters);

			if (!error) {
				
				string data_type = "bin8";

				vector<string> fields = split(parameters, " ");

				int oparams = 0;

				if (fields[0] == "bin4" || fields[0] == "bin8" || fields[0] == "text") {

					oparams++;
					data_type = fields[0];
				}
				
				if (fields.size() > oparams) {

					//get filename
					string fileName = combine(subvec(fields, oparams), " ");

					if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

					if (!err_hndl.call(&SuperMesh::SaveOnScreenPhysicalQuantity, &SMesh, fileName, data_type)) {

						BD.DisplayConsoleMessage("Data saved : " + fileName);
					}
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVEOVF2PARAMVAR:
		{
			string parameters;

			//expecting (data_type) (meshname) paramname (directory/)filename
			error = commandSpec.GetParameters(command_fields, parameters);

			if (!error) {

				string data_type = "bin8";
				string fileName, meshname, paramname;

				vector<string> fields = split(parameters, " ");

				int oparams = 0;

				if (fields[0] == "bin4" || fields[0] == "bin8" || fields[0] == "text") {

					//data type specified (if not, default stands)
					oparams++;
					data_type = fields[0];
				}

				if (fields.size() > oparams) {

					meshname = fields[oparams];
					if (SMesh.contains(meshname)) {

						//meshname specified
						oparams++;
					}
					else {

						//meshname not specified : use focused mesh
						meshname = SMesh.GetMeshFocus();
					}
				}

				if (fields.size() > oparams + 1) {

					//get parameter name
					paramname = fields[oparams];

					if (!SMesh[meshname]->contains_param(paramname)) {

						err_hndl.show_error(BERROR_INCORRECTNAME, verbose);
						break;
					}

					//get filename
					fileName = combine(subvec(fields, oparams + 1), " ");
				}
				//something wrong : not enough parameters specified
				else {

					err_hndl.show_error(BERROR_PARAMMISMATCH_SHOW, verbose);
					break;
				}

				if (GetFileTermination(fileName) != ".ovf")
					fileName += ".ovf";

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				OVF2 ovf2;
				PARAM_ paramID = (PARAM_)SMesh[meshname]->get_meshparam_id(paramname);

				if (SMesh[meshname]->is_paramvarequation_set(paramID)) {

					if (SMesh[meshname]->is_paramvar_scalar(paramID)) {

						VEC<double> s_scaling(SMesh[meshname]->get_paramtype_cellsize(paramID), SMesh[meshname]->meshRect);
						SMesh[meshname]->calculate_meshparam_s_scaling(paramID, s_scaling, SMesh.GetStageTime());

						error = ovf2.Write_OVF2_SCA(fileName, s_scaling, data_type);
					}
					else {

						VEC<DBL3> s_scaling(SMesh[meshname]->get_paramtype_cellsize(paramID), SMesh[meshname]->meshRect);
						SMesh[meshname]->calculate_meshparam_s_scaling(paramID, s_scaling, SMesh.GetStageTime());

						error = ovf2.Write_OVF2_VEC(fileName, s_scaling, data_type);
					}
				}
				else {

					void* s_scaling = SMesh[meshname]->get_meshparam_s_scaling(paramID);

					if (SMesh[meshname]->is_paramvar_scalar(paramID)) error = ovf2.Write_OVF2_SCA(fileName, *reinterpret_cast<VEC<double>*>(s_scaling), data_type);
					else error = ovf2.Write_OVF2_VEC(fileName, *reinterpret_cast<VEC<DBL3>*>(s_scaling), data_type);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2DISP:
		{
			string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);
			if (!error) {

				StopSimulation();

				if (SMesh.active_mesh()->Magnetism_Enabled() && SMesh.active_mesh()->IsModuleSet(MOD_MELASTIC)) {

					if (GetFileTermination(fileName) != ".ovf")
						fileName += ".ovf";

					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

					if (!error) {

						error = SMesh.active_mesh()->CallModuleMethod(&MElastic::Load_Displacement_OVF2, fileName);
						UpdateScreen();
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2STRAIN:
		{
			string parameters;

			error = commandSpec.GetParameters(command_fields, parameters);
			if (!error) {

				StopSimulation();

				if (SMesh.active_mesh()->Magnetism_Enabled() && SMesh.active_mesh()->IsModuleSet(MOD_MELASTIC)) {

					vector<string> fields = split(parameters, ".ovf");

					if (fields.size() >= 2) {

						string fileName_diag, fileName_odiag;

						fileName_diag = fields[0] + ".ovf";
						fileName_odiag = fields[1] + ".ovf";

						if (!GetFilenameDirectory(fileName_diag).length()) fileName_diag = directory + fileName_diag;
						if (!GetFilenameDirectory(fileName_odiag).length()) fileName_odiag = directory + fileName_odiag;

						if (!error) {

							error = SMesh.active_mesh()->CallModuleMethod(&MElastic::Load_Strain_OVF2, fileName_diag, fileName_odiag);
							UpdateScreen();
						}
					}
					else err_hndl.show_error(BERROR_COULDNOTOPENFILE, verbose);
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SCRIPTSERVER:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				Script_Server_Control(status);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CHECKUPDATES:
		{
			//Check with "www.boris-spintronics.uk" if program version is up to date, and get latest update time for materials database
			single_call_launch(&Simulation::CheckUpdate, THREAD_HANDLEMESSAGE2);
		}
		break;

		case CMD_EQUATIONCONSTANTS:
		{
			string userconstant_name;
			double value;

			error = commandSpec.GetParameters(command_fields, userconstant_name, value);
			
			if (error) {

				error.reset() = commandSpec.GetParameters(command_fields, userconstant_name);

				if (!error && verbose) {

					if (userConstants.has_key(userconstant_name))
						BD.DisplayConsoleListing(userconstant_name + " = " + ToString(userConstants[userconstant_name]));
					else err_hndl.show_error(BERROR_NOTDEFINED, verbose);
				}
				else if (verbose) Print_EquationConstants();
			}
			else {

				StopSimulation();

				if (!userConstants.has_key(userconstant_name)) userConstants.push_back(value, userconstant_name);
				else userConstants[userconstant_name] = value;

				//now update all objects which have a text equation
				SMesh.UpdateConfiguration_Values(UPDATECONFIG_TEQUATION_CONSTANTS);

				UpdateScreen();
			}
		}
		break;

		case CMD_DELEQUATIONCONSTANT:
		{
			string userconstant_name;

			error = commandSpec.GetParameters(command_fields, userconstant_name);

			if (!error) {

				StopSimulation();

				if (userConstants.has_key(userconstant_name)) userConstants.erase(userconstant_name);

				//now update all objects which have a text equation
				SMesh.UpdateConfiguration_Values(UPDATECONFIG_TEQUATION_CONSTANTS);

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CLEAREQUATIONCONSTANTS:
		{
			StopSimulation();

			userConstants.clear();

			//now update all objects which have a text equation
			SMesh.UpdateConfiguration_Values(UPDATECONFIG_TEQUATION_CONSTANTS);

			UpdateScreen();
		}
		break;

		case CMD_FLUSHERRORLOG:
		{
			ofstream bdout;
			bdout.open(errorlog_fileName.c_str(), ios::out);
			bdout.close();
		}
		break;

		case CMD_ERRORLOG:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				log_errors = status;
				Save_Startup_Flags();
			}
			else if (verbose) Print_ErrorLogStatus();
		}
		break;

		case CMD_STARTUPUPDATECHECK:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				start_check_updates = status;
				Save_Startup_Flags();
			}
			else if (verbose) Print_StartupUpdateCheckStatus();
		}
		break;

		case CMD_STARTUPSCRIPTSERVER:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				start_scriptserver = status;
				Save_Startup_Flags();
			}
			else if (verbose) Print_StartupScriptServerStatus();
		}
		break;

		case CMD_SHOWTC:
		{
			if (SMesh.active_mesh()->Magnetism_Enabled() && SMesh.active_mesh()->is_atomistic()) {

				double Tc = dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->Show_Transition_Temperature();

				if (verbose) BD.DisplayConsoleMessage("Tc = " + ToString(Tc));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Tc));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SHOWMS:
		{
			if (SMesh.active_mesh()->Magnetism_Enabled() && SMesh.active_mesh()->is_atomistic()) {

				double Ms = dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->Show_Ms();

				if (verbose) BD.DisplayConsoleMessage("Ms = " + ToString(Ms));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Ms));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SHOWA:
		{
			if (SMesh.active_mesh()->Magnetism_Enabled() && SMesh.active_mesh()->is_atomistic()) {

				double A = dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->Show_A();

				if (verbose) BD.DisplayConsoleMessage("A = " + ToString(A));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(A));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SHOWK:
		{
			if (SMesh.active_mesh()->Magnetism_Enabled() && SMesh.active_mesh()->is_atomistic()) {

				double Ku = dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->Show_Ku();

				if (verbose) BD.DisplayConsoleMessage("Ku = " + ToString(Ku));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Ku));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
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

					bool found_nonempty = false;

					for (int arr_idx = 0; arr_idx < dpArr.size(); arr_idx++) {

						found_nonempty |= (bool)dpArr[arr_idx].size();

						if (dpArr[arr_idx].size()) BD.DisplayConsoleListing("dp array " + ToString(arr_idx) + " has size " + ToString((int)dpArr[arr_idx].size()));
					}

					if (!found_nonempty) BD.DisplayConsoleMessage("All arrays are empty.");
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

				error = dpArr.load_arrays(fileName, vec_convert<int, string>(split(indexes_string)), &rows_read);

				if (verbose && !error) BD.DisplayConsoleMessage("Number of rows read : " + ToString(rows_read));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVE:
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

				error = dpArr.save_arrays(fileName, vec_convert<int, string>(split(indexes_string)));

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVEASROW:
		{
			string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric string (the file path and file name) then a numeric string.
			vector<string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				string index_string = entries.back();

				int dp_index = ToNum(index_string);

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (GetFileTermination(fileName) != ".txt")
					fileName += ".txt";

				error = dpArr.save_array_transposed(fileName, dp_index);

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVEAPPENDASROW:
		{
			string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric string (the file path and file name) then a numeric string.
			vector<string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				string index_string = entries.back();

				int dp_index = ToNum(index_string);

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (GetFileTermination(fileName) != ".txt")
					fileName += ".txt";

				error = dpArr.save_array_transposed(fileName, dp_index, true);

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_NEWFILE:
		{
			string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (!error) {

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (GetFileTermination(fileName) != ".txt")
					fileName += ".txt";

				std::ofstream bdout;
				bdout.open(fileName, std::ios::out);
				bdout.close();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVEAPPEND:
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

				error = dpArr.save_arrays(fileName, vec_convert<int, string>(split(indexes_string)), true);

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

		case CMD_DP_GETEXACTPROFILE:
		{
			DBL3 start, end;
			double step;
			int arr_idx;

			error = commandSpec.GetParameters(command_fields, start, end, step, arr_idx);

			if (!error) {

				//prepare data displayed on screen ready to be read out below
				SMesh.PrepareDisplayedMeshValue();

				int num_points = round((end - start).norm() / step) + 1;

				vector<double> data_x(num_points), data_y(num_points), data_z(num_points);

				#pragma omp parallel for
				for (int idx = 0; idx < num_points; idx++) {

					DBL3 pos = start + idx * (end - start) / (num_points - 1);

					DBL3 value = SMesh.GetDisplayedMeshValue(pos);

					data_x[idx] = value.x;
					data_y[idx] = value.y;
					data_z[idx] = value.z;
				}

				dpArr.set_array(arr_idx, data_x);
				dpArr.set_array(arr_idx + 1, data_y);
				dpArr.set_array(arr_idx + 2, data_z);

				if (verbose) BD.DisplayConsoleMessage("Path extracted.");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_GETPATH:
		{
			int arr_idx_path, arr_idx_data;

			error = commandSpec.GetParameters(command_fields, arr_idx_path, arr_idx_data);

			if (!error) {

				vector<double> x, y, z;
				x = dpArr[arr_idx_path];
				y = dpArr[arr_idx_path + 1];
				z = dpArr[arr_idx_path + 2];

				int pathsize = minimum(x.size(), y.size(), z.size());

				vector<double> data_x(pathsize), data_y(pathsize), data_z(pathsize);

				//prepare data displayed on screen ready to be read out below
				SMesh.PrepareDisplayedMeshValue();

				#pragma omp parallel for
				for (int idx = 0; idx < pathsize; idx++) {

					DBL3 value = SMesh.GetDisplayedMeshValue(DBL3(x[idx], y[idx], z[idx]));

					data_x[idx] = value.x;
					data_y[idx] = value.y;
					data_z[idx] = value.z;
				}

				dpArr.set_array(arr_idx_data, data_x);
				dpArr.set_array(arr_idx_data + 1, data_y);
				dpArr.set_array(arr_idx_data + 2, data_z);

				if (verbose) BD.DisplayConsoleMessage("Path extracted.");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GETVALUE:
		{
			DBL3 abs_pos;

			error = commandSpec.GetParameters(command_fields, abs_pos);

			if (!error) {

				//prepare data displayed on screen ready to be read out below
				SMesh.PrepareDisplayedMeshValue();
				Any value = SMesh.GetDisplayedMeshValue(abs_pos);
				if (verbose) BD.DisplayConsoleMessage(value.convert_to_string());

				if (script_client_connected) {

					//Longer version so it compiles with C++14
					if (value.is_type(btype_info<double>())) {
						
						double value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL2>())) {
						
						DBL2 value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL3>())) {
						
						DBL3 value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
				}
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

		case CMD_AVERAGEMESHRECT:
		{
			Rect rect;

			error = commandSpec.GetParameters(command_fields, rect);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); rect = Rect(); }

			if (!error) {
				
				Any value = SMesh.GetAverageDisplayedMeshValue(rect);
				if (verbose) BD.DisplayConsoleMessage("Average value = " + value.convert_to_string());

				if (script_client_connected) {

					//Longer version so it compiles with C++14
					if (value.is_type(btype_info<double>())) {
						
						double value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL2>())) {
						
						DBL2 value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL3>())) {
						
						DBL3 value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
				}
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

				if (SMesh.active_mesh()->Magnetism_Enabled()) {

					double Q = 0.0;

					if (!SMesh.active_mesh()->is_atomistic()) {

						error = dpArr.get_topological_charge(dynamic_cast<Mesh*>(SMesh.active_mesh())->Get_M(), x, y, radius, &Q);
					}
					else {

						error = dpArr.get_topological_charge(dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->Get_M1(), x, y, radius, &Q);
					}

					if (!error) {

						if (verbose) BD.DisplayConsoleMessage("Q = " + ToString(Q));

						if (script_client_connected)
							commSocket.SetSendData(commandSpec.PrepareReturnParameters(Q));
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_COUNTSKYRMIONS:
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

				if (SMesh.active_mesh()->Magnetism_Enabled()) {

					double Q = 0.0;

					if (!SMesh.active_mesh()->is_atomistic()) {

						error = dpArr.count_skyrmions(dynamic_cast<Mesh*>(SMesh.active_mesh())->Get_M(), x, y, radius, &Q);
					}
					else {

						error = dpArr.count_skyrmions(dynamic_cast<Atom_Mesh*>(SMesh.active_mesh())->Get_M1(), x, y, radius, &Q);
					}

					if (!error) {

						if (verbose) BD.DisplayConsoleMessage("Qmag = " + ToString(Q));

						if (script_client_connected)
							commSocket.SetSendData(commandSpec.PrepareReturnParameters(Q));
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_HISTOGRAM:
		{
			double bin = 0.0, min = 0.0, max = 0.0;
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y, bin, min, max);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_x, dp_y); bin = 0.0; min = 0.0; max = 0.0; }

			if (!error) {

				if (SMesh.active_mesh()->Magnetism_Enabled() && !SMesh.active_mesh()->is_atomistic()) {

					error = dpArr.calculate_histogram(dynamic_cast<Mesh*>(SMesh.active_mesh())->Get_M(), dp_x, dp_y, bin, min, max);

					if (!error) {

						if (verbose) BD.DisplayConsoleMessage("Histogram calculated.");
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
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

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				DBL2 R, x0, Ms, w;

				error = dpArr.fit_skyrmion(dp_x, dp_y, &R, &x0, &Ms, &w);

				if (verbose && !error) {

					BD.DisplayConsoleMessage(
						"R = " + ToString(R.major) + " +/- " + ToString(R.minor) + ", " +
						"x0 = " + ToString(x0.major) + " +/- " + ToString(x0.minor) + ", " +
						"Ms = " + ToString(Ms.major) + " +/- " + ToString(Ms.minor) + ", " +
						"w = " + ToString(w.major) + " +/- " + ToString(w.minor)
					);
				}

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(R.major, x0.major, Ms.major, w.major, R.minor, x0.minor, Ms.minor, w.minor));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_FITSTT:
		{
			Rect rectangle;

			error = commandSpec.GetParameters(command_fields, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); rectangle = SMesh.active_mesh()->GetMeshRect(); }

			if (!error) {

				string meshName = SMesh.GetMeshFocus();

				if (SMesh[meshName]->GetMeshType() == MESH_FERROMAGNETIC && SMesh[meshName]->IsModuleSet(MOD_TRANSPORT) && SMesh.SolveSpinCurrent() && 
					(IsNZ(dynamic_cast<Mesh*>(SMesh[meshName])->ts_eff.get0()) || IsNZ(dynamic_cast<Mesh*>(SMesh[meshName])->tsi_eff.get0()))) {

					VEC_VC<DBL3>& M = dynamic_cast<Mesh*>(SMesh[meshName])->Get_M();
					VEC_VC<DBL3>& J = dynamic_cast<Mesh*>(SMesh[meshName])->Get_Jc();
					VEC<DBL3> T(M.h, M.rect);

					if (IsZ(dynamic_cast<Mesh*>(SMesh[meshName])->ts_eff.get0())) {

						//interfacial spin torque only
						T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque());

						if (!T.linear_size()) {

							err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
							break;
						}
					}
					else if (IsZ(dynamic_cast<Mesh*>(SMesh[meshName])->tsi_eff.get0())) {

						//bulk spin torque only
						T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_SpinTorque());

						if (!T.linear_size()) {

							err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
							break;
						}
					}
					else {

						//both bulk and interfacial spin torque
						T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque());

						if (!T.linear_size()) {

							err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
							break;
						}

						T.add_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_SpinTorque());
					}		

					DBL2 P, beta;
					double Rsq = 0.0;
					error = dpArr.fit_stt(T, M, J, rectangle, &P, &beta, &Rsq);

					if (verbose && !error) {

						BD.DisplayConsoleMessage(
							"P = " + ToString(P.major) + " +/- " + ToString(P.minor) + ", " +
							"beta = " + ToString(beta.major) + " +/- " + ToString(beta.minor) +
							", R^2 = " + ToString(Rsq));
					}

					if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(P.major, beta.major, P.minor, beta.minor, Rsq));
				}
				else err_hndl.show_error(BERROR_SPINSOLVER_FIT, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_FITADIABATIC:
		case CMD_DP_FITNONADIABATIC:
		{
			string meshName = SMesh.GetMeshFocus();

			double abs_error_threshold, Rsq_threshold, T_ratio_threshold;
			int stencil = 3;
			bool user_thresholds = true;

			error = commandSpec.GetParameters(command_fields, abs_error_threshold, Rsq_threshold, T_ratio_threshold, stencil);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, abs_error_threshold, Rsq_threshold, T_ratio_threshold);  stencil = 3; }
			if (error == BERROR_PARAMMISMATCH) { error.reset(); user_thresholds = false; }

			if (SMesh[meshName]->GetMeshType() == MESH_FERROMAGNETIC && SMesh[meshName]->IsModuleSet(MOD_TRANSPORT) && SMesh.SolveSpinCurrent() &&
				(IsNZ(dynamic_cast<Mesh*>(SMesh[meshName])->ts_eff.get0()) || IsNZ(dynamic_cast<Mesh*>(SMesh[meshName])->tsi_eff.get0()))) {

				VEC_VC<DBL3>& M = dynamic_cast<Mesh*>(SMesh[meshName])->Get_M();
				VEC_VC<DBL3>& J = dynamic_cast<Mesh*>(SMesh[meshName])->Get_Jc();
				VEC<DBL3> T(M.h, M.rect);

				if (IsZ(dynamic_cast<Mesh*>(SMesh[meshName])->ts_eff.get0())) {

					//interfacial spin torque only
					T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque());

					if (!T.linear_size()) {

						err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
						break;
					}
				}
				else if (IsZ(dynamic_cast<Mesh*>(SMesh[meshName])->tsi_eff.get0())) {

					//bulk spin torque only
					T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_SpinTorque());

					if (!T.linear_size()) {

						err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
						break;
					}
				}
				else {

					//both bulk and interfacial spin torque
					T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque());

					if (!T.linear_size()) {

						err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
						break;
					}

					T.add_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_SpinTorque());
				}

				if (verbose) BD.DisplayConsoleMessage("Fitting... Please wait.");

				if (user_thresholds) {

					error = dpArr.fit_stt_variation(T, M, J, SMesh[meshName]->displayVEC_SCA, commandSpec.cmdId == CMD_DP_FITADIABATIC, abs_error_threshold, Rsq_threshold, T_ratio_threshold, stencil);
				}
				else {

					error = dpArr.fit_stt_variation(T, M, J, SMesh[meshName]->displayVEC_SCA, commandSpec.cmdId == CMD_DP_FITADIABATIC);
				}

				if (!error) {

					if (commandSpec.cmdId == CMD_DP_FITADIABATIC) {

						if (verbose) BD.DisplayConsoleMessage("Finished fitting P - see results in Cust_S display.");
					}
					else {

						if (verbose) BD.DisplayConsoleMessage("Finished fitting beta - see results in Cust_S display.");
					}

					UpdateScreen();
				}
			}
			else err_hndl.show_error(BERROR_SPINSOLVER_FIT2, verbose);
		}
		break;

		case CMD_DP_FITSOT:
		{
			Rect rectangle;
			string hm_mesh;

			error = commandSpec.GetParameters(command_fields, hm_mesh, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, hm_mesh); rectangle = SMesh.active_mesh()->GetMeshRect(); }

			if (!error) {

				string meshName = SMesh.GetMeshFocus();

				if (SMesh[meshName]->GetMeshType() == MESH_FERROMAGNETIC && SMesh[meshName]->IsModuleSet(MOD_TRANSPORT) && SMesh.SolveSpinCurrent()
					&& SMesh.contains(hm_mesh) && SMesh[hm_mesh]->GetMeshType() == MESH_METAL && SMesh[hm_mesh]->IsModuleSet(MOD_TRANSPORT)) {

					VEC<DBL3>& T = dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque();
					VEC_VC<DBL3>& M = dynamic_cast<Mesh*>(SMesh[meshName])->Get_M();
					VEC_VC<DBL3>& J = dynamic_cast<Mesh*>(SMesh[hm_mesh])->Get_Jc();

					DBL2 SHAeff, flST;
					double Rsq = 0.0;
					error = dpArr.fit_sot(T, M, J, rectangle, &SHAeff, &flST, &Rsq);

					if (verbose && !error) {

						BD.DisplayConsoleMessage(
							"SHAeff = " + ToString(SHAeff.major) + " +/- " + ToString(SHAeff.minor) + ", " +
							"flST = " + ToString(flST.major) + " +/- " + ToString(flST.minor) +
							", R^2 = " + ToString(Rsq));
					}

					if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SHAeff.major, flST.major, SHAeff.minor, flST.minor, Rsq));
				}
				else err_hndl.show_error(BERROR_SPINSOLVER_FIT3, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_FITSOTSTT:
		{
			Rect rectangle;
			string hm_mesh;

			error = commandSpec.GetParameters(command_fields, hm_mesh, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, hm_mesh); rectangle = SMesh.active_mesh()->GetMeshRect(); }

			if (!error) {

				string meshName = SMesh.GetMeshFocus();

				if (SMesh[meshName]->GetMeshType() == MESH_FERROMAGNETIC && SMesh[meshName]->IsModuleSet(MOD_TRANSPORT) && SMesh.SolveSpinCurrent()
					&& SMesh.contains(hm_mesh) && SMesh[hm_mesh]->GetMeshType() == MESH_METAL && SMesh[hm_mesh]->IsModuleSet(MOD_TRANSPORT)) {

					VEC_VC<DBL3>& M = dynamic_cast<Mesh*>(SMesh[meshName])->Get_M();
					VEC_VC<DBL3>& J_hm = dynamic_cast<Mesh*>(SMesh[hm_mesh])->Get_Jc();
					VEC_VC<DBL3>& J_fm = dynamic_cast<Mesh*>(SMesh[meshName])->Get_Jc();
					VEC<DBL3> T(M.h, M.rect);

					if (IsZ(dynamic_cast<Mesh*>(SMesh[meshName])->ts_eff.get0())) {

						//interfacial spin torque only
						T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque());

						if (!T.linear_size()) {

							err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
							break;
						}
					}
					else if (IsZ(dynamic_cast<Mesh*>(SMesh[meshName])->tsi_eff.get0())) {

						//bulk spin torque only
						T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_SpinTorque());

						if (!T.linear_size()) {

							err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
							break;
						}
					}
					else {

						//both bulk and interfacial spin torque
						T.copy_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_InterfacialSpinTorque());

						if (!T.linear_size()) {

							err_hndl.show_error(BERROR_NOTCOMPUTED, verbose);
							break;
						}

						T.add_values(dynamic_cast<Mesh*>(SMesh[meshName])->Get_SpinTorque());
					}

					DBL2 SHAeff, flST, P, beta;
					double Rsq = 0.0;
					error = dpArr.fit_sotstt(T, M, J_hm, J_fm, rectangle, &SHAeff, &flST, &P, &beta, &Rsq);

					if (verbose && !error) {

						BD.DisplayConsoleMessage(
							"SHAeff = " + ToString(SHAeff.major) + " +/- " + ToString(SHAeff.minor) + ", " +
							"flST = " + ToString(flST.major) + " +/- " + ToString(flST.minor) + ", " +
							"P = " + ToString(P.major) + " +/- " + ToString(P.minor) + ", " +
							"beta = " + ToString(beta.major) + " +/- " + ToString(beta.minor) +
							", R^2 = " + ToString(Rsq));
					}

					if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SHAeff.major, flST.major, P.major, beta.major, SHAeff.minor, flST.minor, P.minor, beta.minor, Rsq));
				}
				else err_hndl.show_error(BERROR_SPINSOLVER_FIT3, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_CALCSOT:
		{
			string hm_mesh, fm_mesh;

			error = commandSpec.GetParameters(command_fields, hm_mesh, fm_mesh);

			if (!error) {

				if (SMesh.contains(fm_mesh) && SMesh[fm_mesh]->GetMeshType() == MESH_FERROMAGNETIC && 
					SMesh.contains(hm_mesh) && SMesh[hm_mesh]->GetMeshType() == MESH_METAL) {

					double SHAeff, flST;

					double SHA = dynamic_cast<Mesh*>(SMesh[hm_mesh])->SHA.get0();
					double d_N = dynamic_cast<Mesh*>(SMesh[hm_mesh])->GetMeshDimensions().z;
					double lsf_N = dynamic_cast<Mesh*>(SMesh[hm_mesh])->l_sf.get0();
					double sigma_N = dynamic_cast<Mesh*>(SMesh[hm_mesh])->elecCond.get0();
					DBL2 G;

					if (SMesh[hm_mesh]->GetOrigin().z > SMesh[fm_mesh]->GetOrigin().z) {

						//hm mesh on top of fm mesh : hm mesh sets Gmix
						G = dynamic_cast<Mesh*>(SMesh[hm_mesh])->Gmix.get0();
					}
					else {

						//fm mesh on top of hm mesh : fm mesh sets Gmix
						G = dynamic_cast<Mesh*>(SMesh[fm_mesh])->Gmix.get0();
					}

					//this is G tilda
					G *= 2 / sigma_N;

					double Nlam = tanh(d_N / lsf_N) / lsf_N;

					SHAeff = SHA * (1.0 - 1.0 / cosh(d_N / lsf_N)) * (Nlam*G.i + G.i*G.i - G.j*G.j) / ((Nlam + G.i)*(Nlam + G.i) + G.j*G.j);
					flST = (Nlam*G.j + 2 * G.i*G.j) / (G.i*G.i - G.j*G.j + Nlam*G.i);

					if (verbose) {

						BD.DisplayConsoleMessage("SHAeff = " + ToString(SHAeff) + ", " + "flST = " + ToString(flST));
					}

					if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SHAeff, flST));
				}
				else err_hndl.show_error(BERROR_SPINSOLVER_FIT4, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_CALCEXCHANGE:
		{
			string meshName = SMesh.GetMeshFocus();

			if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->ExchangeComputation_Enabled()) {

				SMesh[meshName]->Compute_Exchange();
				UpdateScreen();
			}
			else err_hndl.show_error(BERROR_INCORRECTCONFIG, verbose);
		}
		break;

		case CMD_DP_CALCTOPOCHARGEDENSITY:
		{
			string meshName = SMesh.GetMeshFocus();

			if (SMesh[meshName]->Magnetism_Enabled()) {

				SMesh[meshName]->Compute_TopoChargeDensity();
				UpdateScreen();
			}
			else err_hndl.show_error(BERROR_INCORRECTCONFIG, verbose);
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

		case CMD_DP_MONOTONIC:
		{
			int dp_in_x, dp_in_y, dp_out_x, dp_out_y;

			error = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y);

			if (!error) {

				error = dpArr.extract_monotonic(dp_in_x, dp_in_y, dp_out_x, dp_out_y);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_CROSSINGSHISTOGRAM:
		{
			int dp_in_x, dp_in_y, dp_out_x, dp_out_y, steps;

			error = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y, steps);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y); steps = 100; }

			if (!error) {

				error = dpArr.crossings_histogram(dp_in_x, dp_in_y, dp_out_x, dp_out_y, steps);

				if (verbose) BD.DisplayConsoleMessage("Histogram computed.");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_CROSSINGSFREQUENCY:
		{
			int dp_in_x, dp_in_y, dp_out_x, dp_out_y, dp_out_y2, steps;

			error = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y, dp_out_y2, steps);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y, dp_out_y2); steps = 100; }

			if (!error) {

				error = dpArr.crossings_frequencies(dp_in_x, dp_in_y, dp_out_x, dp_out_y, dp_out_y2, steps);

				if (verbose) BD.DisplayConsoleMessage("Histograms computed.");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_PEAKSFREQUENCY:
		{
			int dp_in_x, dp_in_y, dp_out_x, dp_out_y, steps;

			error = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y, steps);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_in_x, dp_in_y, dp_out_x, dp_out_y); steps = 100; }

			if (!error) {

				error = dpArr.peaks_frequencies(dp_in_x, dp_in_y, dp_out_x, dp_out_y, steps);

				if (verbose) BD.DisplayConsoleMessage("Histogram computed.");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		//---------------- finished CMD_DP_ commands

		case CMD_TEST:
		{
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
		}
		break;

		default:
			break;
		}

		if (error) err_hndl.show_error(error, verbose);
	}
	else err_hndl.show_error(BERROR_COMMAND_NOTRECOGNIZED, verbose);
}