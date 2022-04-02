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
	std::string message = commSocket.Listen();

	if (message.length()) {

		//must make sure there is no clash between direct user input and scripted input, otherwise crashes can result
		userInputMutex.lock();

		//handle command using the usual HandleCommand method, but on a blocking thread call - need to wait for it to finish before responding to client
		set_blocking_thread(THREAD_HANDLEMESSAGE);
		//Keep trying to launch the command handler : if THREAD_HANDLEMESSAGE is busy or has not yet been marked as not active (this is possible even if the previous call to THREAD_HANDLEMESSAGE has finished and we returned parameters to client - it can take much longer than usual sometimes),
		//we receive THREAD_GENERATENEW response back until THREAD_HANDLEMESSAGE is available (marked as not active).
		//If the command handler was launched and completed successfully, the single_call_launch method will respond with THREAD_HANDLEMESSAGE after it returns.
		while (single_call_launch<std::string>(&Simulation::HandleCommand, message, THREAD_HANDLEMESSAGE) != THREAD_HANDLEMESSAGE);

		//the command was handled using a blocking thread, so now parameters should be available to send
		commSocket.SendDataParams();

		userInputMutex.unlock();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Command handler

void Simulation::HandleCommand(std::string command_string) 
{
	////////////////////////////////////////////////////////////////////////////
	//
	// AUXILIARY LAMBDAS

	//adjust new_meshName such that it doesn't clash with any existing mesh names; return adjusted name
	auto adjust_meshname = [&](std::string new_meshName) -> std::string {

		//does this mesh name already exist? If yes then modify to a name which doesn't exist
		while (SMesh.contains(new_meshName)) {

			//if meshName exists then add _#, where # = 1 to start
			//if meshName already contains a termination of the form _# then increase # until meshName available
			size_t pos = new_meshName.find_last_of('_');

			if (pos == std::string::npos) {

				new_meshName += std::string("_1");
			}
			else {

				std::string termination = new_meshName.substr(pos + 1);
				if (has_digits_only(termination)) {

					int number = ToNum(termination);
					number++;
					termination = ToString(number);

					new_meshName = new_meshName.substr(0, pos + 1) + termination;
				}
				else new_meshName += std::string("_1");
			}
		}

		return new_meshName;
	};

	//delete all meshes, except existingmeshName (which must exist)
	auto delete_all_meshes_except = [&](std::string existingmeshName) {

		if (!SMesh.contains(existingmeshName)) return;

		std::vector<std::string> meshNames;
		for (int idx = 0; idx < SMesh.size(); idx++) {

			meshNames.push_back(SMesh.key_from_meshIdx(idx));
		}

		for (auto& meshName : meshNames) {

			if (meshName != existingmeshName) {

				SMesh.DelMesh(meshName);
				//delete any affected entries in data box
				DeleteDataBoxFields(meshName);
				//delete any affected entries in saveDataList
				DeleteSaveDataEntries(meshName);
			}
		}
	};

	//invoked by commands where the command has meshName as an optional parameter, with default value of focused mesh
	//adjust command fields so the meshName appears at the start explicitly
	//This lambda guarantees the first entry in command_fields will be an existing meshname
	//If check_supermesh is true, then include supermesh handle in possible mesh names
	auto optional_meshname_check_focusedmeshdefault = [&](std::vector<std::string>& command_fields, bool check_supermesh = false) -> bool {

		for (int idx = 0; idx < command_fields.size(); idx++) {

			if (SMesh.contains(command_fields[idx]) || (check_supermesh && command_fields[idx] == SMesh.superMeshHandle)) {

				if (idx) {

					std::string meshName = command_fields[idx];
					command_fields.erase(command_fields.begin() + idx);
					command_fields.insert(command_fields.begin(), meshName);
				}

				//return true to indicate mesh name was included already
				return true;
			}

			if (idx == command_fields.size() - 1 && !SMesh.contains(command_fields[idx]) && !(check_supermesh && command_fields[idx] == SMesh.superMeshHandle)) {

				command_fields.insert(command_fields.begin(), SMesh.GetMeshFocus());
				
				//return false to indicate mesh name was not previously included
				return false;
			}
		}

		if (!command_fields.size()) {

			command_fields.push_back(SMesh.GetMeshFocus());
			return false;
		}

		//will never get here but compiler doesn't know that
		return true;
	};

	//check if the optional display quantity name has been given
	auto optional_meshname_check_focusedmeshdefault_quantityname = [&](std::vector<std::string>& command_fields, bool check_supermesh = false) -> bool {

		//first fix meshname
		optional_meshname_check_focusedmeshdefault(command_fields, check_supermesh);

		for (int idx = 0; idx < command_fields.size(); idx++) {
			//a quantity name found, so return true to indicate it is available
			if (displayHandles.has_value(command_fields[idx])) return true;
		}

		//quantity name not found so return false
		//Also insert "Nothing" at index 1 : this will go after meshname
		command_fields.insert(command_fields.begin() + 1, displayHandles(MESHDISPLAY_NONE));
		return false;
	};

	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	//cannot change simulation parameters in the middle of an interation. simulationMutex also used by Simulate() method, but there it is non-blocking
	simulationMutex.lock();
	std::lock_guard<std::mutex> simlock(simulationMutex, std::adopt_lock);

	//ensure we do not handle a command during a display refresh
	//This method locks and unlocks the display std::mutex (cannot leave it locked and unlock it at the end of HandleCommand, since HandleCommand can call other BD methods)
	//It is important to do this since HandleCommand can cause changes to data which are used during a display refresh and could result in a crash.
	//It is theoretically possible a BD method is called right after WaitForDisplayEnd but before HandleCommand returns - this could happen in response to a user action only.
	//This is virtually impossible however, as the user would have to issue commands quicker than it takes for the program to make changes to critical data.
	//Not ideal solution, the only foolproof solution is to lock the display std::mutex and unlock it before calling any BD methods here (or before returning), but after changes to data have been made.
	//This is messy so I'd rather not go down that route since it only guards against an extremely unlikely scenario.
	BD.WaitForDisplayEnd();

	//display messages in console only if this is true. Note, some important messages are displayed regardless (e.g. simulation stopped or simulation started messages.)
	bool verbose = true;

	//if command came from a script client, send back required parameters. Each command must ensure it does this.
	bool script_client_connected = false;

	//commands are always of the form: commandname (param1) (param2) ...
	std::vector<std::string> command_fields = split(command_string, " ");
	std::string command_name = command_fields.front();

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

#if COMPILECUDA == 1
		//Commands are executed on newly spawned threads, so if cuda is on and we are not using device 0 (default device) we must switch to required device, otherwise 0 will be used
		if (cudaEnabled && cudaDeviceSelect != 0) cudaSetDevice(cudaDeviceSelect);
#endif

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

		case CMD_RUNSTAGE:
		{
			int stage;
			error = commandSpec.GetParameters(command_fields, stage);

			if (!error && stage < simStages.size()) {

				stage_step = INT2(stage, 0);
				SetSimulationStageValue();
				single_stage_run = true;
				RunSimulation();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
			break;

		case CMD_SETSCHEDULESTAGE:
		{
			int stage;
			error = commandSpec.GetParameters(command_fields, stage);

			if (!error && stage < simStages.size()) {

				stage_step = INT2(stage, 0);
				SetSimulationStageValue();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(stage_step.major));
		}
		break;

		case CMD_NEXTSTAGE:
		{
			AdvanceSimulationSchedule();
		}
		break;

		case CMD_COMPUTEFIELDS:
			ComputeFields();
			break;
			
		//DEPRECATED
		//This was used by Python scripts to detect when simulation has finished. Still keep it for compatibility with old scripts.
		//Instead you should use the ns.Run() call in Python scripts
		case CMD_ISRUNNING:
			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(is_thread_running(THREAD_LOOP)));
			break;

		case CMD_RUNSCRIPT:
		{
#if PYTHON_EMBEDDING == 1
			std::string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (GetFileTermination(fileName) != ".py") fileName += ".py";
			if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

			if (!error) {

				//launch script on dedicated non-blocking thread to avoid blocking simulationMutex
				single_call_launch<std::string>(&Simulation::RunPythonScript, fileName, THREAD_PYTHONSCRIPT);
				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(std::string("Running script in embedded mode : " + fileName)));
			}
#endif
		}
			break;

		case CMD_STOPSCRIPT:
		{
#if PYTHON_EMBEDDING == 1
			if (is_thread_running(THREAD_PYTHONSCRIPT)) {

				python_script_running = false;

				while (is_thread_running(THREAD_LOOP)) {
					StopSimulation();
					std::this_thread::sleep_for(std::chrono::milliseconds(100));
				}
			}
			else BD.DisplayConsoleError("No Python scripts currently running.");
#endif
		}
		break;

		case CMD_SCRIPTPRINT:
		{
			std::string message = combine(command_fields, " ");
			BD.DisplayConsoleMessage(message);
		}
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
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::function<void(std::string, Rect)> save_data_updater = [&](std::string meshName, Rect meshRect_old) -> void {

					//update rectangles in saveDataList for the named mesh
					UpdateSaveDataEntries(meshRect_old, SMesh[meshName]->GetMeshRect(), meshName);
					//update rectangles in dataBoxList for the named mesh
					UpdateDataBoxEntries(meshRect_old, SMesh[meshName]->GetMeshRect(), meshName);
				};

				if (!err_hndl.call(error, &SuperMesh::SetMeshRect, &SMesh, meshName, meshRect, save_data_updater)) UpdateScreen_AutoSet();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetMeshRect()));
		}
		break;

		case CMD_SHIFTDIPOLE:
		{
			DBL3 shift;
			std::string meshName;

			error = commandSpec.GetParameters(command_fields, meshName, shift);

			if (!error) {

				if (SMesh.contains(meshName)) SMesh[meshName]->Shift_Dipole(shift);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DIPOLEVELOCITY:
		{
			DBL3 velocity, clipping;
			std::string meshName;

			error = commandSpec.GetParameters(command_fields, meshName, velocity, clipping);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, velocity); clipping = DBL3(DIPOLESHIFTCLIP); }

			if (!error) {

				if (SMesh.contains(meshName)) SMesh[meshName]->Set_Dipole_Velocity(velocity, clipping);
				UpdateScreen();
			}
			else if (verbose) Print_DipoleShiftingAlgorithm_List();

			if (script_client_connected) {

				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->Get_Dipole_Velocity(), SMesh[meshName]->Get_Dipole_Clipping()));
			}
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
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, h);

			if (!error) {

				StopSimulation();
				if (!err_hndl.call(error, &MeshBase::SetMeshCellsize, SMesh[meshName], h)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetMeshCellsize()));
		}
		break;

		case CMD_ECELLSIZE:
		{
			DBL3 h_e;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, h_e);

			if (!error) {

				StopSimulation();
				if (!err_hndl.call(error, &MeshBase::SetMeshECellsize, SMesh[meshName], h_e)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetMeshECellsize()));
		}
		break;

		case CMD_TCELLSIZE:
		{
			DBL3 h_t;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, h_t);

			if (!error) {

				StopSimulation();
				if (!err_hndl.call(error, &MeshBase::SetMeshTCellsize, SMesh[meshName], h_t)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetMeshTCellsize()));
		}
		break;

		case CMD_SCELLSIZE:
		{
			DBL3 h_s;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, h_s);

			if (!error && !SMesh[meshName]->is_atomistic()) {

				StopSimulation();
				if (!err_hndl.call(error, &Mesh::SetMeshSCellsize, dynamic_cast<Mesh*>(SMesh[meshName]), h_s, false)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected && !SMesh[meshName]->is_atomistic()) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh[meshName])->GetMeshSCellsize()));
		}
		break;

		case CMD_MCELLSIZE:
		{
			DBL3 h_m;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, h_m);

			if (!error) {

				StopSimulation();
				if (!err_hndl.call(error, &MeshBase::SetMeshMCellsize, SMesh[meshName], h_m)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetMeshMCellsize()));
		}
		break;

		case CMD_FMSCELLSIZE:
		{
			DBL3 h;
			error = commandSpec.GetParameters(command_fields, h);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::SetFMSMeshCellsize, &SMesh, h)) {

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

				if (!err_hndl.call(error, &SuperMesh::SetESMeshCellsize, &SMesh, h)) {

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
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, h_dm);

			if (!error && SMesh[meshName]->is_atomistic()) {

				StopSimulation();
				if (!err_hndl.call(error, &Atom_Mesh::Set_Demag_Cellsize, dynamic_cast<Atom_Mesh*>(SMesh[meshName]), h_dm)) UpdateScreen();
				else if (verbose) error(BERROR_NOTATOMISTIC);
			}
			else if (verbose) Print_Speedup_List();

			if (script_client_connected && SMesh[meshName]->is_atomistic()) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Get_Demag_Cellsize()));
		}
		break;

		case CMD_ADDFMESH:
		{
			std::string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, MESH_FERROMAGNETIC, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETFMESH:
		{
			std::string new_meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, new_meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::string new_meshName_actual = new_meshName;
				new_meshName = adjust_meshname(new_meshName);

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, MESH_FERROMAGNETIC, meshRect)) {

					SMesh.SetMeshFocus(new_meshName);
					UpdateScreen_AutoSet_KeepOrientation();

					delete_all_meshes_except(new_meshName);

					//if we had to change name, now rename the mesh to the actual required mesh name
					if (new_meshName != new_meshName_actual) SMesh.RenameMesh(new_meshName, new_meshName_actual);

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDAFMESH:
		{
			std::string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, MESH_ANTIFERROMAGNETIC, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETAFMESH:
		{
			std::string new_meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, new_meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::string new_meshName_actual = new_meshName;
				new_meshName = adjust_meshname(new_meshName);

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, MESH_ANTIFERROMAGNETIC, meshRect)) {

					SMesh.SetMeshFocus(new_meshName);
					UpdateScreen_AutoSet_KeepOrientation();

					delete_all_meshes_except(new_meshName);

					//if we had to change name, now rename the mesh to the actual required mesh name
					if (new_meshName != new_meshName_actual) SMesh.RenameMesh(new_meshName, new_meshName_actual);

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDMETALMESH:
		{
			std::string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, MESH_METAL, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETMETALMESH:
		{
			std::string new_meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, new_meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::string new_meshName_actual = new_meshName;
				new_meshName = adjust_meshname(new_meshName);

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, MESH_METAL, meshRect)) {

					SMesh.SetMeshFocus(new_meshName);
					UpdateScreen_AutoSet_KeepOrientation();

					delete_all_meshes_except(new_meshName);

					//if we had to change name, now rename the mesh to the actual required mesh name
					if (new_meshName != new_meshName_actual) SMesh.RenameMesh(new_meshName, new_meshName_actual);

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDINSULATORMESH:
		{
			std::string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, MESH_INSULATOR, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETINSULATORMESH:
		{
			std::string new_meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, new_meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::string new_meshName_actual = new_meshName;
				new_meshName = adjust_meshname(new_meshName);

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, MESH_INSULATOR, meshRect)) {

					SMesh.SetMeshFocus(new_meshName);
					UpdateScreen_AutoSet_KeepOrientation();

					delete_all_meshes_except(new_meshName);

					//if we had to change name, now rename the mesh to the actual required mesh name
					if (new_meshName != new_meshName_actual) SMesh.RenameMesh(new_meshName, new_meshName_actual);

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDDIPOLEMESH:
		{
			std::string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, MESH_DIPOLE, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETDIPOLEMESH:
		{
			std::string new_meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, new_meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::string new_meshName_actual = new_meshName;
				new_meshName = adjust_meshname(new_meshName);

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, MESH_DIPOLE, meshRect)) {

					SMesh.SetMeshFocus(new_meshName);
					UpdateScreen_AutoSet_KeepOrientation();

					delete_all_meshes_except(new_meshName);

					//if we had to change name, now rename the mesh to the actual required mesh name
					if (new_meshName != new_meshName_actual) SMesh.RenameMesh(new_meshName, new_meshName_actual);

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDAMESHCUBIC:
		{
			std::string meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, meshName, meshRect);

			if (!error) {

				StopSimulation();

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, MESH_ATOM_CUBIC, meshRect)) {

					UpdateScreen();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETAMESHCUBIC:
		{
			std::string new_meshName;
			Rect meshRect;

			//note, mesh name is not allowed to have any spaces - needs to be a single word
			error = commandSpec.GetParameters(command_fields, new_meshName, meshRect);

			if (!error) {

				StopSimulation();

				std::string new_meshName_actual = new_meshName;
				new_meshName = adjust_meshname(new_meshName);

				if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, MESH_ATOM_CUBIC, meshRect)) {

					SMesh.SetMeshFocus(new_meshName);
					UpdateScreen_AutoSet_KeepOrientation();

					delete_all_meshes_except(new_meshName);

					//if we had to change name, now rename the mesh to the actual required mesh name
					if (new_meshName != new_meshName_actual) SMesh.RenameMesh(new_meshName, new_meshName_actual);

					UpdateScreen_AutoSet();
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DELMESH:
		{
			std::string meshName;
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();

				//delete mesh
				if (!err_hndl.qcall(error, &SuperMesh::DelMesh, &SMesh, meshName)) {

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
			std::string oldName = SMesh.GetMeshFocus(), newName;

			error = commandSpec.GetParameters(command_fields, oldName, newName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, newName); oldName = SMesh.GetMeshFocus(); }

			if (!error) {

				//note, mesh name is not allowed to have any spaces - needs to be a single word
				newName = trimspaces(newName);

				if (!err_hndl.qcall(error, &SuperMesh::RenameMesh, &SMesh, oldName, newName)) {

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
			std::string meshName;
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				std::string foucusedMesh = SMesh.GetMeshFocus();

				if (!err_hndl.qcall(error, &SuperMesh::SetMeshFocus, &SMesh, meshName)) {

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
			std::string meshName;
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				std::string foucusedMesh = SMesh.GetMeshFocus();

				if (!err_hndl.qcall(error, &SuperMesh::SetMeshFocus, &SMesh, meshName)) {

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
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); rectangle = Rect(DBL3(), SMesh[meshName]->GetMeshDimensions()); }

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::delrect, &SMesh, meshName, rectangle)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADDRECT:
		{
			Rect rectangle;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, rectangle);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::setrect, &SMesh, meshName, rectangle)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RESETMESH:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::resetrect, &SMesh, meshName)) UpdateScreen();
			}
		}
		break;

		case CMD_LOADMASKFILE:
		{
			std::string meshName, params_string;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, params_string);

			if (!error) {

				StopSimulation();

				double zDepth;
				std::string fileName;

				//using split_numeric approach since the file name path can contain spaces.
				std::vector<std::string> entries = split_numeric(params_string);

				if (entries.size() == 2) {

					zDepth = ToNum(entries[0]);
					fileName = entries[1];
				}
				else fileName = params_string;

				if (GetFileTermination(fileName) != ".png") fileName += ".png";
				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				ScanFileNameData(fileName);

				//define code required to load the bitmap - this will be passed for the mesh to use
				std::function<std::vector<unsigned char>(std::string, INT2)> bitmap_loader = [&](std::string fileName, INT2 n_plane) -> std::vector<unsigned char> {

					std::vector<unsigned char> bitmap;
					BD.BGMethods()->GetBitmapFromImage(fileName, bitmap, n_plane);

					return bitmap;
				};

				error = SMesh[meshName]->applymask(zDepth, fileName, bitmap_loader);
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
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, polar, azim);
			if (!meshName_specified) meshName = "";

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagAngle, &SMesh, meshName, polar, azim)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_FLOWER:
		{
			int direction;
			DBL3 centre;
			double radius, thickness;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, direction, radius, thickness, centre);
			if (error == BERROR_PARAMMISMATCH) { 

				error.reset() = commandSpec.GetParameters(command_fields, meshName, direction);
				if (error == BERROR_PARAMMISMATCH) { error.reset(); direction = 1; }

				//entire mesh
				radius = 0.0;
				//entire thickness
				thickness = 0.0;
				//centered
				centre = SMesh[meshName]->GetMeshDimensions() / 2;
			}

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagFlower, &SMesh, meshName, direction, centre, radius, thickness)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ONION:
		{
			int direction;
			DBL3 centre;
			double radius1, radius2, thickness;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, direction, radius1, radius2, thickness, centre);
			if (error == BERROR_PARAMMISMATCH) {

				error.reset() = commandSpec.GetParameters(command_fields, meshName, direction);
				if (error == BERROR_PARAMMISMATCH) { error.reset(); direction = 1; }
				
				//entire mesh
				radius1 = 0.0; radius2 = 0.0;
				//entire thickness
				thickness = 0.0;
				//centered
				centre = SMesh[meshName]->GetMeshDimensions() / 2;
			}

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagOnion, &SMesh, meshName, direction, centre, radius1, radius2, thickness)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CROSSTIE:
		{
			int direction;
			DBL3 centre;
			double radius, thickness;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, direction, radius, thickness, centre);
			if (error == BERROR_PARAMMISMATCH) {

				error.reset() = commandSpec.GetParameters(command_fields, meshName, direction);
				if (error == BERROR_PARAMMISMATCH) { error.reset(); direction = 1; }

				//entire mesh
				radius = 0.0;
				//entire thickness
				thickness = 0.0;
				//centered
				centre = SMesh[meshName]->GetMeshDimensions() / 2;
			}

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagCrosstie, &SMesh, meshName, direction, centre, radius, thickness)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETOBJECTANGLE:
		{
			double polar, azim;
			DBL3 position;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, polar, azim, position);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagAngle_Object, &SMesh, meshName, polar, azim, position)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RANDOM:
		{
			int seed;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetRandomMag, &SMesh, meshName, seed)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RANDOMXY:
		{
			int seed;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetRandomXYMag, &SMesh, meshName, seed)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_INVERTMAG:
		{
			std::string meshName = SMesh.GetMeshFocus();
			std::string components;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, components);
			if (error == BERROR_PARAMMISMATCH) error.reset();

			if (!error) {

				bool x = true, y = true, z = true;

				if (components.length()) {

					std::vector<std::string> fields = split(components, " ");

					x = false; y = false; z = false;
					for (int idx = 0; idx < fields.size(); idx++) {

						if (fields[idx] == "x") x = true;
						if (fields[idx] == "y") y = true;
						if (fields[idx] == "z") z = true;
					}
				}
				
				if (!err_hndl.qcall(error, &SuperMesh::SetInvertedMag, &SMesh, meshName, x, y, z)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_MIRRORMAG:
		{
			std::string meshName = SMesh.GetMeshFocus();
			std::string axis;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, axis);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMirroredMag, &SMesh, meshName, axis)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETRECT:
		{
			Rect rectangle;
			double polar, azim;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, polar, azim, rectangle);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagAngle_Rect, &SMesh, meshName, polar, azim, rectangle)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DWALL:
		{
			std::string longitudinal, transverse;
			double width, position;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, longitudinal, transverse, width, position);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMagDomainWall, &SMesh, meshName, longitudinal, transverse, width, position)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_VORTEX:
		{
			int longitudinal, rotation, core;
			Rect rect;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, longitudinal, rotation, core, rect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, longitudinal, rotation, core); rect = Rect(); }

			if (!error && longitudinal && rotation && core) {

				StopSimulation();

				std::string vortexFile;
				bool invertMag = false;

				//l = +1, r = -1, c = -1 : vortex_hh_clk_dn
				//l = +1, r = -1, c = +1 : vortex_hh_clk_up
				//l = +1, r = +1, c = -1 : vortex_hh_cclk_dn
				//l = +1, r = +1, c = +1 : vortex_hh_cclk_up

				//l = -1, r = -1, c = -1 : vortex_hh_cclk_up, then invert magnetization.
				//l = -1, r = -1, c = +1 : vortex_hh_cclk_dn, then invert magnetization.
				//l = -1, r = +1, c = -1 : vortex_hh_clk_up, then invert magnetization.
				//l = -1, r = +1, c = +1 : vortex_hh_clk_dn, then invert magnetization.

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

					//data loaded correctly, so resize currently focused mesh (if ferromagnetic) then copy magnetization data to it.
					if (SMesh[meshName]->is_atomistic()) {

						data.renormalize(dynamic_cast<Atom_Mesh*>(SMesh[meshName])->mu_s.get0());
					}
					else {

						data.renormalize(dynamic_cast<Mesh*>(SMesh[meshName])->Ms.get0());
					}

					if (invertMag) data *= -1.0;

					if (SMesh[meshName]->Magnetism_Enabled()) {

						SMesh[meshName]->SetMagFromData(data, rect);
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
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, core, chirality, diameter, position);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetSkyrmion, &SMesh, meshName, core, chirality, diameter, position)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SKYRMIONBLOCH:
		{
			int core, chirality;
			double diameter;
			DBL2 position;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, core, chirality, diameter, position);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetSkyrmionBloch, &SMesh, meshName, core, chirality, diameter, position)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PBC:
		{
			std::string meshName;
			std::string flag;
			int images;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, flag, images);
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::Set_PBC, &SMesh, meshName, flag, images)) UpdateScreen();
			}
			else if (verbose) Print_PBC();
		}
		break;

		case CMD_SETFIELD:
		{
			DBL3 field_polar;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, field_polar);
			if (!meshName_specified) meshName = "";

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetField, &SMesh, meshName, Polar_to_Cartesian(field_polar))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) {

				if (!meshName_specified) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->CallModuleMethod(&ZeemanBase::GetField)));
			}
		}
		break;

		case CMD_SETSTRESS:
		{
			DBL3 stress_polar;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, stress_polar);
			if (!meshName_specified) meshName = "";

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetUniformStress, &SMesh, meshName, Polar_to_Cartesian(stress_polar))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) {

				if (!meshName_specified) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->CallModuleMethod(&MElastic::GetUniformStress)));
			}
		}
		break;

		case CMD_MODULES:
		{
			std::string meshName;
			std::string modules;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			if (meshName_specified) {

				error = commandSpec.GetParameters(command_fields, meshName, modules);
				if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName); };
			}

			if (meshName_specified) {

				std::vector<std::string> set_modules_handles;

				auto get_modules_handles = [&]() -> std::vector<std::string> {

					std::vector<std::string> modules_handles;

					std::vector<MOD_> modules_IDs = SMesh[meshName]->GetModulesIDs();
					for (int idx = 0; idx < modules_IDs.size(); idx++) {

						modules_handles.push_back(moduleHandles(modules_IDs[idx]));
					}

					return modules_handles;
				};

				//get currently set modules
				set_modules_handles = get_modules_handles();

				if (modules.length()) {

					//handles of modules to set
					std::vector<std::string> toset_modules_handles = split(modules, " ");

					//check input module handles - if any name errors found then need to issue error to user otherwise configuration will not be as expected
					for (int idx = 0; idx < toset_modules_handles.size(); idx++) {
						
						if (!moduleHandles.has_value(toset_modules_handles[idx])) {

							error(BERROR_INCORRECTNAME);
							break;
						}
					}

					for (int idx = 0; idx < set_modules_handles.size(); idx++) {

						//if any module currently set is not in the list of modules to set, then delete it (we need to delete those first before setting required modules, since some are exclusive)
						if (!vector_contains(toset_modules_handles, set_modules_handles[idx]))
							SMesh.DelModule(meshName, (MOD_)moduleHandles.get_ID_from_value(set_modules_handles[idx]));
					}

					//now add all required modules if not set already
					for (int idx = 0; idx < toset_modules_handles.size(); idx++) {

						MOD_ modID = (MOD_)moduleHandles.get_ID_from_value(toset_modules_handles[idx]);
						if (!SMesh[meshName]->IsModuleSet(modID)) SMesh.AddModule(meshName, modID);
					}

					//update list of set modules
					set_modules_handles = get_modules_handles();
				}
							
				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(combine(set_modules_handles, " ")));
			}
			else if (verbose) Print_Modules_List();
		}
		break;

		case CMD_ADDMODULE:
		{
			std::string moduleHandle, meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, moduleHandle);
			if (!meshName_specified) meshName = "";

			if (!error) {
				
				StopSimulation();
				if(!err_hndl.call(error, &SuperMesh::AddModule, &SMesh, meshName, (MOD_)moduleHandles.get_ID_from_value(moduleHandle))) RefreshScreen();
			}
			else if (verbose) Print_Modules_List();
		}
		break;

		case CMD_DELMODULE:
		{
			std::string moduleHandle, meshName;
			
			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, moduleHandle);
			if (!meshName_specified) meshName = "";

			if (!error) {
				
				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::DelModule, &SMesh, meshName, (MOD_)moduleHandles.get_ID_from_value(moduleHandle))) RefreshScreen();
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
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, status);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::Set_Demag_Exclusion, &SMesh, status, meshName)) RefreshScreen();
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GPUKERNELS:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::Set_Kernel_Initialize_on_GPU, &SMesh, status)) RefreshScreen();
				RefreshScreen();
			}
			else if (verbose) Print_GPUKernels_Config();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.Get_Kernel_Initialize_on_GPU()));
		}
		break;

		case CMD_ODE:
		{
			if (verbose) Print_ODEs();
		}
		break;

		case CMD_SETODE:
		{
			std::string odeHandle, odeEvalHandle;
			
			error = commandSpec.GetParameters(command_fields, odeHandle, odeEvalHandle);

			if (!error) {

				StopSimulation();

				ODE_ setOde = (ODE_)odeHandles.get_ID_from_value(odeHandle);
				EVAL_ odeEval = (EVAL_)odeEvalHandles.get_ID_from_value(odeEvalHandle);

				if (setOde != ODE_ERROR && odeEval != EVAL_ERROR && vector_contains(odeAllowedEvals(setOde), odeEval)) {

					if (!err_hndl.call(error, &SuperMesh::SetODE, &SMesh, setOde, odeEval)) {

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
			std::string odeHandle, odeEvalHandle;

			error = commandSpec.GetParameters(command_fields, odeHandle, odeEvalHandle);

			if (!error) {

				StopSimulation();

				ODE_ setatom_Ode = (ODE_)atom_odeHandles.get_ID_from_value(odeHandle);
				EVAL_ odeEval = (EVAL_)odeEvalHandles.get_ID_from_value(odeEvalHandle);

				ODE_ odeID;
				SMesh.QueryODE(odeID);

				if (setatom_Ode != ODE_ERROR && odeEval != EVAL_ERROR && vector_contains(odeAllowedEvals(setatom_Ode), odeEval) && vector_contains(odeAllowedEvals(odeID), odeEval)) {

					if (!err_hndl.call(error, &SuperMesh::SetAtomisticODE, &SMesh, setatom_Ode, odeEval)) {

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
			std::string odeEvalHandle;

			error = commandSpec.GetParameters(command_fields, odeEvalHandle);

			if (!error) {

				StopSimulation();

				EVAL_ odeEval = (EVAL_)odeEvalHandles.get_ID_from_value(odeEvalHandle);

				ODE_ odeID, atom_odeID;
				SMesh.QueryODE(odeID);
				SMesh.QueryAtomODE(atom_odeID);

				if (odeEval != EVAL_ERROR && vector_contains(odeAllowedEvals(odeID), odeEval) && vector_contains(odeAllowedEvals(atom_odeID), odeEval)) {

					if (!err_hndl.call(error, &SuperMesh::SetODEEval, &SMesh, odeEval)) {

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

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.GetEvaluationSpeedup()));
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
			double err_fail, dT_incr, dT_min, dT_max;

			error = commandSpec.GetParameters(command_fields, err_fail, dT_incr, dT_min, dT_max);

			if (!error) {

				SMesh.SetAdaptiveTimeStepCtrl(err_fail, dT_incr, dT_min, dT_max);
				UpdateScreen();
			}
			else if (verbose) Print_AStepCtrl();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(DBL4(SMesh.Get_AStepRelErrCtrl(), SMesh.Get_AStepdTCtrl().i, SMesh.Get_AStepdTCtrl().j, SMesh.Get_AStepdTCtrl().k)));
		}
		break;

		case CMD_SHOWDATA:
		{
			std::string dataName;
			std::string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dataName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dataName); dataRect = Rect(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				if (verbose) {

					std::string text;

					//mesh name if applicable
					if (!dataDescriptor[dataName].meshless) text += "<b><" + meshName + "> ";

					//box dimensions if applicable
					if (!dataDescriptor[dataName].boxless && !dataRect.IsNull()) text += "<b>(" + ToString(dataRect, "m") + ") ";

					//Label
					text += "<b>" + dataDescriptor[dataName].Label;
					//and the actual value(s)
					text += "</c><i>" + GetDataValueString(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect)) + "</i>";

					BD.DisplayFormattedConsoleMessage(text);
				}

				//for script return the number of returned data fields is variable. This is done by obtaining the a std::string using GetDataValueString, then splitting it using the usual separators (, or ;)
				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(GetDataValue(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect))));

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
			std::string dataName;
			std::string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dataName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dataName); dataRect = Rect(); }

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
			std::string dataName;
			std::string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dataName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dataName); dataRect = Rect(); }

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
			std::string dataName;
			std::string meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, index, dataName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, index, dataName); dataRect = Rect(); }

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

		case CMD_SAVEDATA:
		{
			int append_option;
			error = commandSpec.GetParameters(command_fields, append_option);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); append_option = -1; }

			//0 : make new file and save to it immediately
			if (append_option == 0) {

				//first reset buffer positions
				savedata_diskbuffer_position = 0;

				//Get the data to save (into buffer)
				SaveData();

				//and immediately flush buffer, making sure to create new file
				appendToDataFile = false;
				if (savedata_diskbuffer_position) SaveData_DiskBufferFlush(&savedata_diskbuffer, &savedata_diskbuffer_position);
			}
			else {

				if (append_option >= 1) appendToDataFile = true;

				//if simulation is currently running then save data in the usual way (append it to buffer, which will eventually be emptied when full or simulation stops)
				if (is_thread_running(THREAD_LOOP)) SaveData();
				else {

					//if simulation is not running then we want to save immediately - so SaveData() then empty buffer
					SaveData();
					if (savedata_diskbuffer_position) SaveData_DiskBufferFlush(&savedata_diskbuffer, &savedata_diskbuffer_position);
				}
			}
		}
		break;

		case CMD_ADDPINNEDDATA:
		{
			std::string dataName, meshName = SMesh.GetMeshFocus();
			Rect dataRect;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dataName, dataRect);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dataName); dataRect = Rect(); }

			if (!error && dataDescriptor.has_key(dataName) && SMesh.contains(meshName)) {

				if (dataDescriptor[dataName].meshless) NewDataBoxField(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName)));
				else NewDataBoxField(DatumConfig((DATA_)dataDescriptor.get_ID_from_key(dataName), meshName, dataRect));

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
			std::string directory_;

			error = commandSpec.GetParameters(command_fields, directory_);

			if (!error) {

				directory = FixedDirectorySlashes(directory_);

				if (directory.substr(directory.length() - 1) != "/" && directory.substr(directory.length() - 1) != "/") {

					directory += "/";
				}

				RefreshScreen();
			}
			else {

				std::string text = "[tc1,1,1,1/tc]Current working directory : " + MakeIO(IOI_DIRECTORY, directory);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(directory));
		}
		break;

		case CMD_SAVEDATAFILE:
		{
			std::string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (!error) {

				std::string directory_ = ExtractFilenameDirectory(fileName);
				if (directory_.length()) directory = directory_;
				savedataFile = fileName;

				if (!GetFileTermination(savedataFile).length()) savedataFile += ".txt";

				RefreshScreen();
			}
			else {

				std::string text = "[tc1,1,1,1/tc]Current file for output data : " + MakeIO(IOI_DIRECTORY, directory) + MakeIO(IOI_SAVEDATAFILE, savedataFile);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(savedataFile));
		}
		break;

		case CMD_SAVECOMMENT:
		{
			std::string fileName,  comment;

			error = commandSpec.GetParameters(command_fields, fileName, comment);

			if (!error) {

				std::string directory_ = ExtractFilenameDirectory(fileName);
				if (!directory_.length()) fileName = directory + fileName;
				if (!GetFileTermination(fileName).length()) fileName += ".txt";
				
				std::ofstream bdout;
				bdout.open(ScanFileNameData(fileName), std::ios::out | std::ios::app);
				if (bdout.is_open()) {

					bdout << comment << std::endl;
					bdout.close();
				}
			}
		}
		break;

		case CMD_SAVEIMAGEFILE:
		{
			std::string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (!error) {

				std::string directory_ = ExtractFilenameDirectory(fileName);
				if (directory_.length()) directory = directory_;
				imageSaveFileBase = fileName;

				RefreshScreen();
			}
			else {

				std::string text = "[tc1,1,1,1/tc]Current file image saving : " + MakeIO(IOI_DIRECTORY, directory) + MakeIO(IOI_SAVEIMAGEFILEBASE, imageSaveFileBase);
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

		case CMD_DATAPRECISION:
		{
			int precision;

			error = commandSpec.GetParameters(command_fields, precision);

			if (!error) {

				StopSimulation();

				dataprecision = precision;
				Conversion::tostringconversion::precision() = dataprecision;

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dataprecision));
		}
		break;

		case CMD_DISKBUFFERLINES:
		{
			int bufferlines;

			error = commandSpec.GetParameters(command_fields, bufferlines);

			if (!error) {

				StopSimulation();

				savedata_diskbuffer_size = bufferlines;
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(savedata_diskbuffer_size));
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
			std::string stageTypeName;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, stageTypeName);

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
			std::string stageTypeName;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, stageTypeName);

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
			std::string stageTypeName;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, index, stageTypeName);

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
			std::string value_string;

			error = commandSpec.GetParameters(command_fields, index, value_string);

			if (!error && GoodIdx(simStages.last(), index)) {

				EditStageValue(index, value_string);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETSTAGEVALUE:
		{
			int index;

			error = commandSpec.GetParameters(command_fields, index);

			if (!error && GoodIdx(simStages.last(), index)) {

				SetSimulationStageValue(index);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_EDITSTAGESTOP:
		{
			int index;
			std::string stageStopName;
			std::string value_string;

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
			std::string dSaveType;
			std::string value_string;

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
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh.contains(meshName) && verbose) Print_MeshParams(meshName);
		}
		break;

		case CMD_SETPARAM:
		{
			std::string paramName, paramValue, meshName;
			bool set_value = true;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName, paramValue);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName); set_value = false; }

			if (!error && set_value) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::set_meshparam_value, &SMesh, meshName, paramName, paramValue)) UpdateScreen();
			}
			else if (verbose && set_value) PrintCommandUsage(command_name);

			if (script_client_connected && !set_value) {

				SMesh.get_meshparam_value(meshName, paramName, paramValue);
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(paramValue));
			}
		}
		break;

		case CMD_SETKTENS:
		{
			std::string kterms;
			std::string meshName = SMesh.GetMeshFocus();

			error = commandSpec.GetParameters(command_fields, kterms);

			if (!error) {
		
				//input terms should of the form dxn1yn2zn3
				//Here d is a number (floating point), n1, n2, n3 are integers. x, y, z are string literals and we expect to find them
				std::vector<std::string> ktens = split(kterms, " ");

				//extract d n1 n2 n3 values
				std::vector<DBL4> Kt;

				for (int idx = 0; idx < ktens.size(); idx++) {

					double d = 1.0;
					int n1 = 0, n2 = 0, n3 = 0;
					
					ktens[idx] += " ";
					size_t pos = ktens[idx].find_first_of("xyz");
					if (pos != std::string::npos) {

						if (pos) if (ktens[idx].substr(0, pos) != "-") d = ToNum(ktens[idx].substr(0, pos)); else d = -1;
						std::vector<std::string> entry_fields = split(ktens[idx].substr(pos), { "x", "y", "z" });
						
						for (int fidx = 1; fidx < entry_fields.size(); fidx++) {

							std::string component = ktens[idx].substr(pos, 1);
							if (component == "x") if (entry_fields[fidx].length() && entry_fields[fidx] != " ") n1 = ToNum(entry_fields[fidx]); else n1 = 1;
							if (component == "y") if (entry_fields[fidx].length() && entry_fields[fidx] != " ") n2 = ToNum(entry_fields[fidx]); else n2 = 1;
							if (component == "z") if (entry_fields[fidx].length() && entry_fields[fidx] != " ") n3 = ToNum(entry_fields[fidx]); else n3 = 1;
							pos += entry_fields[fidx].length() + 1;
						}
					}
					else d = ToNum(ktens[idx]);
					
					if (d != 0.0 && n1 >= 0 && n2 >= 0 && n3 >= 0) Kt.push_back(DBL4(d, n1, n2, n3));
				}

				error = SMesh[meshName]->set_tensorial_anisotropy(Kt);
			}
			else if (verbose) BD.DisplayConsoleMessage(SMesh[meshName]->get_tensorial_anisotropy_string());
		}
		break;

		case CMD_PARAMSTEMP:
		{
			std::string meshName;
			
			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh.contains(meshName) && verbose) Print_MeshParamsTemperature(meshName);
		}
		break;

		case CMD_CLEARPARAMSTEMP:
		{
			std::string meshName, paramName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName);
			if (error == BERROR_PARAMMISMATCH && meshName_specified) { error.reset(); paramName = ""; }
			if (!meshName_specified) { meshName = ""; }

			StopSimulation();
			if (!err_hndl.qcall(error, &SuperMesh::clear_meshparam_temp, &SMesh, meshName, paramName)) UpdateScreen();
		}
		break;

		case CMD_SETPARAMTEMPEQUATION:
		{
			std::string meshName, paramName, equationText;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName, equationText);
			if (!meshName_specified) { meshName = ""; }

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::set_meshparam_t_equation, &SMesh, meshName, paramName, equationText)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETPARAMTEMPARRAY:
		{
			std::string meshName, paramName, fileName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName, fileName);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName); fileName = ""; }
			if (!meshName_specified) { meshName = ""; }

			//setting parameter with arrays loaded form file
			if (!error) {

				//load from a file directly if given
				StopSimulation();

				if (fileName.length()) {

					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
					if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

					ScanFileNameData(fileName);

					std::vector<std::vector<double>> load_arrays;
					if (ReadDataColumns(fileName, "\t", load_arrays, { 0, 1 })) {

						std::vector<double> empty;

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
					std::vector<double> temp(2), scaling(2), empty;
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

					std::vector<double> empty;

					if (GoodIdx(dpArr.size(), dp_c_iy, dp_c_iz)) error = SMesh.set_meshparam_tscaling_array(meshName, paramName, dpArr[dp_T_idx], dpArr[dp_c_ix], dpArr[dp_c_iy], dpArr[dp_c_iz]);
					else if (GoodIdx(dpArr.size(), dp_c_iy)) error = SMesh.set_meshparam_tscaling_array(meshName, paramName, dpArr[dp_T_idx], dpArr[dp_c_ix], dpArr[dp_c_iy], empty);
					else error = SMesh.set_meshparam_tscaling_array(meshName, paramName, dpArr[dp_T_idx], dpArr[dp_c_ix], empty, empty);

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
			std::string meshName_from, meshNames_to_string;

			error = commandSpec.GetParameters(command_fields, meshName_from, meshNames_to_string);

			//multiple mesh names could have been entered
			std::vector<std::string> meshNames_to = split(meshNames_to_string, " ");

			if (!error && meshNames_to.size()) {

				StopSimulation();

				for (int idx = 0; idx < meshNames_to.size(); idx++) {

					err_hndl.qcall(error, &SuperMesh::copy_mesh_parameters, &SMesh, meshName_from, meshNames_to[idx]);
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_COPYMESHDATA:
		{
			std::string meshName_from, meshNames_to_string;

			error = commandSpec.GetParameters(command_fields, meshName_from, meshNames_to_string);

			//multiple mesh names could have been entered
			std::vector<std::string> meshNames_to = split(meshNames_to_string, " ");

			if (!error && meshNames_to.size()) {

				StopSimulation();

				for (int idx = 0; idx < meshNames_to.size(); idx++) {

					err_hndl.qcall(error, &SuperMesh::copy_mesh_data, &SMesh, meshName_from, meshNames_to[idx]);
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PARAMSVAR:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh.contains(meshName) && verbose) Print_MeshParamsVariation(meshName);
		}
		break;

		case CMD_SETDISPLAYEDPARAMSVAR:
		{
			std::string meshName, paramName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::set_meshparamvar_display, &SMesh, meshName, paramName)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_CLEARPARAMSVAR:
		{
			std::string meshName, paramName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName);
			if (error == BERROR_PARAMMISMATCH && meshName_specified) { error.reset(); paramName = ""; }
			if (!meshName_specified) { meshName = ""; }

			StopSimulation();
			if (!err_hndl.qcall(error, &SuperMesh::clear_meshparam_variation, &SMesh, meshName, paramName)) UpdateScreen();
		}
		break;

		case CMD_SETPARAMVAR:
		{
			std::string meshName, paramName, generatorName, generatorArgs;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, paramName, generatorName, generatorArgs);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, paramName, generatorName); generatorArgs = ""; }
			if (!meshName_specified) meshName = "";

			if (!error) {

				StopSimulation();

				//generatorName can be specified as equation, but typically this would not be specified and simply generatorName has the equation text
				//in this case generatorName is not a key in vargenerator_descriptor - this is how we recognize this type of input
				if (vargenerator_descriptor.get_ID_from_key(generatorName) == MATPVAR_EQUATION || !vargenerator_descriptor.has_key(generatorName)) {

					std::string equationText;

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

					if (!err_hndl.qcall(error, &SuperMesh::set_meshparam_s_equation, &SMesh, meshName, paramName, equationText)) {

						SMesh.set_meshparamvar_display(meshName, paramName);
						UpdateScreen();
					}
				}
				else {

					//used with custom parameter variation generator (set from grayscale png file)
					std::function<std::vector<unsigned char>(std::string, INT2)> bitmap_loader = [&](std::string fileName, INT2 n_plane) -> std::vector<unsigned char> {

						std::vector<unsigned char> bitmap;
						BD.BGMethods()->GetBitmapFromImage(fileName, bitmap, n_plane);

						return bitmap;
					};

					error = SMesh.set_meshparam_var(meshName, paramName, generatorName, generatorArgs, bitmap_loader);
					if (!error) SMesh.set_meshparamvar_display(meshName, paramName);
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GETMESHTYPE:
		{
			std::string meshName;

			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				if (SMesh.contains(meshName)) {

					MESH_ meshType = SMesh[meshName]->GetMeshType();

					std::string meshtypeHandle = meshtypeHandles(meshType);

					if (verbose) BD.DisplayConsoleListing(meshtypeHandle);

					if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(meshtypeHandle));
				}
				else error(BERROR_INCORRECTNAME);
			}
		}
		break;

		case CMD_SAVESIM:
		{
			std::string simFileName;

			error = commandSpec.GetParameters(command_fields, simFileName);

			ScanFileNameData(simFileName);

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
			std::string simFileName;

			error = commandSpec.GetParameters(command_fields, simFileName);

			ScanFileNameData(simFileName);

			if (!error) {

				if (!err_hndl.call(error, &Simulation::LoadSimulation, this, simFileName)) {

					if (verbose) BD.DisplayConsoleMessage("Simulation loaded : " + simFileName);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DEFAULT:
		{
			int current_cudaDeviceSelect = cudaDeviceSelect;
			LoadSimulation(GetUserDocumentsPath() + boris_data_directory + boris_simulations_directory + "default");
			BD.DisplayConsoleMessage("Default state restored.");
			cudaDeviceSelect = current_cudaDeviceSelect;
			error = SMesh.SwitchCUDAState(false, current_cudaDeviceSelect);
		}
		break;

		case CMD_DISPLAY:
		{
			std::string name, meshName;

			optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, name);

			if (!error) {

				MESHDISPLAY_ display = (MESHDISPLAY_)displayHandles.get_ID_from_value(name);
				if (!err_hndl.qcall(error, &SuperMesh::SetDisplayedPhysicalQuantity, &SMesh, meshName, (int)display)) UpdateScreen();
			}
			else if(verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_DISPLAYMODULE:
		{
			std::string moduleName, meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, moduleName);

			if (!error) {

				if (SMesh.contains(meshName)) {

					//need to stop simulation so on initialization memory is allocated correctly
					StopSimulation();

					if (moduleHandles.has_value(moduleName)) SMesh[meshName]->Set_Module_Display(moduleHandles.get_ID_from_value(moduleName));
					else SMesh[meshName]->Set_Module_Display(MOD_ALL);

					UpdateScreen();
				}
			}
			else if (verbose) Print_DisplayModules_List();
		}
		break;

		case CMD_ROTCAMABOUTORIGIN:
		{
			double dAzim, dPolar;

			error = commandSpec.GetParameters(command_fields, dAzim, dPolar);

			if (!error) {

				BD.RotateCameraAboutOrigin(dAzim, dPolar);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ROTCAMABOUTAXIS:
		{
			double dAngle;

			error = commandSpec.GetParameters(command_fields, dAngle);

			if (!error) {

				BD.RotateCameraView(dAngle);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ADJUSTCAMDISTANCE:
		{
			double dZ;

			error = commandSpec.GetParameters(command_fields, dZ);

			if (!error) {

				BD.AdjustCameraDistanceFromOrigin(dZ);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SHIFTCAMORIGIN:
		{
			double dX, dY;

			error = commandSpec.GetParameters(command_fields, dX, dY);

			if (!error) {

				BD.Shift3DOriginPixelPosition(dX, dY);
				RefreshScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DISPLAYDETAILLEVEL:
		{
			double detail_level;

			error = commandSpec.GetParameters(command_fields, detail_level);

			if (!error) {

				BD.SetDetailLevel(detail_level);
				UpdateScreen();
			}
			else if (verbose) Print_DisplayRenderSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(BD.GetDetailLevel()));
		}
		break;

		case CMD_DISPLAYRENDERTHRESH:
		{
			INT3 renderthresholds;

			error = commandSpec.GetParameters(command_fields, renderthresholds);

			if (!error) {

				BD.SetRenderThresholds(renderthresholds);
				UpdateScreen();
			}
			else if (verbose) Print_DisplayRenderSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(BD.GetRenderThresholds()));
		}
		break;

		case CMD_DISPLAYBACKGROUND:
		{
			std::string name, meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, name);

			if (!error) {

				MESHDISPLAY_ display = (MESHDISPLAY_)displayHandles.get_ID_from_value(name);
				if (!err_hndl.qcall(error, &SuperMesh::SetDisplayedBackgroundPhysicalQuantity, &SMesh, meshName, (int)display)) UpdateScreen();
			}
			else if (verbose) Print_MeshDisplay_List();
		}
		break;

		case CMD_VECREP:
		{
			std::string meshName;
			int vecreptype;

			optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, vecreptype);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetVEC3Rep, &SMesh, meshName, (int)vecreptype)) UpdateScreen();
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
			std::string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			ScanFileNameData(fileName);

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

		case CMD_SAVEIMAGE:
		{
			std::string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			ScanFileNameData(fileName);

			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				fileName = directory + imageSaveFileBase + ".png";
			}
			else {

				if (GetFileTermination(fileName) != ".png")
					fileName += ".png";

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
			}

			if (BD.SaveImage(fileName, SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel()))) {

				if (verbose) BD.DisplayConsoleMessage("Saved : " + fileName);
				RefreshScreen();
			}
			else if (verbose) error(BERROR_COULDNOTSAVEFILE);
		}
		break;

		case CMD_MAKEVIDEO:
		{
			std::string fileNameBase;
			double fps;
			int quality;
			
			error = commandSpec.GetParameters(command_fields, fileNameBase, fps, quality);

			if (!error) {

				StopSimulation();

				std::string sequenceDirectory = directory;
				if (GetFilenameDirectory(fileNameBase).length()) sequenceDirectory = ExtractFilenameDirectory(fileNameBase);

				std::vector<std::string> fileNames = GetFilesInDirectory(sequenceDirectory, fileNameBase, ".png");

				BD.DisplayConsoleMessage("Encoding video ... please wait.");

				if (BD.BGMethods()->MakeVideoFromFileSequence(sequenceDirectory + fileNameBase + ".wmv", fileNames, unsigned(fps), 1.0, quality)) {

					if(verbose) BD.DisplayConsoleMessage("Video created.");
				}
				else if (verbose) error(BERROR_COULDNOTSAVEFILE);
			}
			else if (verbose) PrintCommandUsage(command_name);
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

				std::string text = "[tc1,1,1,1/tc]Image cropping settings : " + MakeIO(IOI_IMAGECROPPING);
				if (verbose) BD.DisplayFormattedConsoleMessage(text);
			}
		}
		break;

		case CMD_MOVINGMESH:
		{
			std::string meshName;

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
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();
				err_hndl.call(error, &SuperMesh::PrepareMovingMesh, &SMesh, meshName);				
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PREPAREMOVINGBLOCHMESH:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();
				err_hndl.call(error, &SuperMesh::PrepareMovingMesh_Bloch, &SMesh, meshName);
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PREPAREMOVINGNEELMESH:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();
				err_hndl.call(error, &SuperMesh::PrepareMovingMesh_Neel, &SMesh, meshName);
				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_PREPAREMOVINGSKYRMIONMESH:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				StopSimulation();
				err_hndl.call(error, &SuperMesh::PrepareMovingMesh_Skyrmion, &SMesh, meshName);
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
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, status);

			if (!error) {

				StopSimulation();
				SMesh.Set_ExchangeCoupledMeshes(status, meshName);
				UpdateScreen();
			}
			else if (verbose) Print_ExchangeCoupledMeshes_List();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetMeshExchangeCoupling()));
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
				
				std::string sides;

				error = commandSpec.GetParameters(command_fields, sides);
				if (error) { error.reset(); sides = "x"; }

				StopSimulation();

				//first clear all electrodes
				SMesh.CallModuleMethod(&STransport::ClearElectrodes);

				Rect smeshRect = SMesh.GetESMeshRect();

				if (sides == "z") {
				
					//left hand side electrode
					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(smeshRect.s, DBL3(smeshRect.e.x, smeshRect.e.y, smeshRect.s.z)));
					SMesh.CallModuleMethod(&STransport::DesignateGroundElectrode, 0);

					//right hand side electrode
					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(DBL3(smeshRect.s.x, smeshRect.s.y, smeshRect.e.z), smeshRect.e));
				}
				else if (sides == "y") {

					//left hand side electrode
					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(smeshRect.s, DBL3(smeshRect.e.x, smeshRect.s.y, smeshRect.e.z)));
					SMesh.CallModuleMethod(&STransport::DesignateGroundElectrode, 0);

					//right hand side electrode
					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(DBL3(smeshRect.s.x, smeshRect.e.y, smeshRect.s.z), smeshRect.e));
				}
				else {

					//left hand side electrode
					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(smeshRect.s, DBL3(smeshRect.s.x, smeshRect.e.y, smeshRect.e.z)));
					SMesh.CallModuleMethod(&STransport::DesignateGroundElectrode, 0);

					//right hand side electrode
					SMesh.CallModuleMethod(&STransport::AddElectrode, 0.0, Rect(DBL3(smeshRect.e.x, smeshRect.s.y, smeshRect.s.z), smeshRect.e));
				}

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
					disabled_transport_solver = false;
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
					disabled_transport_solver = false;
					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&STransport::GetCurrent)));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_SETCURRENTDENSITY:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				DBL3 Jcvalue;
				std::string meshName;

				optional_meshname_check_focusedmeshdefault(command_fields);
				error = commandSpec.GetParameters(command_fields, meshName, Jcvalue);

				if (!error) {

					if (!err_hndl.call(error, &SuperMesh::SetEFromJcValue, &SMesh, Jcvalue, meshName)) {

						disabled_transport_solver = true;
						UpdateScreen();
					}
				}
				else if (verbose) PrintCommandUsage(command_name);
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

		case CMD_DISABLETRANSPORTSOLVER:
		{
			if (SMesh.IsSuperMeshModuleSet(MODS_STRANSPORT)) {

				bool status;

				error = commandSpec.GetParameters(command_fields, status);

				if (!error) {

					StopSimulation();
					disabled_transport_solver = status;
					SMesh.UpdateConfiguration(UPDATECONFIG_TRANSPORT);
					UpdateScreen();
				}
				else if (verbose) PrintTransportSolverConfig();

				if (script_client_connected)
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(disabled_transport_solver));
			}
			else error(BERROR_INCORRECTACTION);
		}
		break;

		case CMD_TMRTYPE:
		{
			int setting;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, setting);
			if (!meshName_specified) meshName = "";

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetTMRType, &SMesh, meshName, (TMR_)setting)) UpdateScreen();
			}
			else if (verbose) Print_TMRType_List();

			if (script_client_connected) {

				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters((int)SMesh[meshName]->GetTMRType()));
			}
		}
		break;

		case CMD_RAPBIAS_EQUATION:
		{
			std::string meshName, text_equation;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, text_equation);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); text_equation = ""; }

			if (meshName_specified) {

				if (!err_hndl.qcall(error, &SuperMesh::SetTMR_BiasEquationParallel, &SMesh, meshName, text_equation)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RAAPBIAS_EQUATION:
		{
			std::string meshName, text_equation;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, text_equation);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); text_equation = ""; }

			if (meshName_specified) {

				if (!err_hndl.qcall(error, &SuperMesh::SetTMR_BiasEquationAntiParallel, &SMesh, meshName, text_equation)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_TEMPERATURE:
		{
			double Temperature;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, Temperature);
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetBaseTemperature, &SMesh, meshName, Temperature)) UpdateScreen();
			}
			else if (verbose) Print_MeshTemperature_List();

			if (script_client_connected) {

				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetBaseTemperature()));
			}
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

				std::string heatdT = "[tc1,1,1,1/tc]Heat Equation Time Step: " + MakeIO(IOI_HEATDT) + "</c>";
				BD.DisplayFormattedConsoleMessage(heatdT);
			}

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.CallModuleMethod(&SHeat::get_heat_dT)));
		}
		break;

		case CMD_AMBIENTTEMPERATURE:
		{
			double T_ambient;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, T_ambient);
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetAmbientTemperature, &SMesh, meshName, T_ambient)) UpdateScreen();
			}
			else if (verbose) Print_HeatBoundaries_List();

			if (script_client_connected) {

				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->CallModuleMethod(&HeatBase::GetAmbientTemperature)));
			}
		}
		break;

		case CMD_ROBINALPHA:
		{
			double alpha_boundary;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, alpha_boundary);
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::SetAlphaHeatBoundary, &SMesh, meshName, alpha_boundary)) UpdateScreen();
			}
			else if (verbose) Print_HeatBoundaries_List();

			if (script_client_connected) {

				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->CallModuleMethod(&HeatBase::GetAlphaBoundary)));
			}
		}
		break;

		case CMD_INSULATINGSIDES:
		{
			std::string literal;
			bool status;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, literal, status);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::SetInsulatingSides, &SMesh, meshName, literal, status)) UpdateScreen();
			}
			else if (verbose) Print_HeatBoundaries_List();

			if (script_client_connected) {
				if (SMesh[meshName]->IsModuleSet(MOD_HEAT)) {

					std::vector<bool> insulating = SMesh[meshName]->CallModuleMethod(&HeatBase::GetInsulatingSides);
					commSocket.SetSendData(commandSpec.PrepareReturnParameters(insulating[0], insulating[1], insulating[2], insulating[3], insulating[4], insulating[5]));
				}
			}
		}
		break;

		case CMD_CURIETEMPERATURE:
		{
			double T_Curie;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, T_Curie);
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::SetCurieTemperature, &SMesh, meshName, T_Curie)) UpdateScreen();
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected) {

				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				if (!SMesh[meshName]->is_atomistic()) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh[meshName])->GetCurieTemperature()));
			}
		}
		break;

		case CMD_CURIETEMPERATUREMATERIAL:
		{
			double T_Curie_material;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, T_Curie_material);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::SetCurieTemperatureMaterial, &SMesh, meshName, T_Curie_material)) UpdateScreen();
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected && !SMesh[meshName]->is_atomistic()) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh[meshName])->GetCurieTemperatureMaterial()));
		}
		break;

		case CMD_ATOMICMOMENT:
		{
			double atomic_moment;
			DBL2 atomic_moment_AFM;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, atomic_moment_AFM);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, atomic_moment); atomic_moment_AFM = DBL2(atomic_moment); }
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &SuperMesh::SetAtomicMagneticMoment, &SMesh, meshName, atomic_moment_AFM)) UpdateScreen();
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected) {
				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				if (SMesh[meshName]->GetMeshType() == MESH_ANTIFERROMAGNETIC) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh[meshName])->GetAtomicMoment_AFM()));
				else if (!SMesh[meshName]->is_atomistic()) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh[meshName])->GetAtomicMoment()));
			}
		}
		break;

		case CMD_TAU:
		{
			DBL2 tau_ii, tau_ij;
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, tau_ii, tau_ij);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, tau_ii); tau_ij = DBL2(-1.0); }
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

			if (!error) {

				StopSimulation();

				if (tau_ij >= DBL2(0.0)) {

					//set both intra and inter terms
					if (!err_hndl.qcall(error, &SuperMesh::SetTcCoupling, &SMesh, meshName, tau_ii, tau_ij)) UpdateScreen();
				}

				//only intra terms
				else if (!err_hndl.qcall(error, &SuperMesh::SetTcCoupling_Intra, &SMesh, meshName, tau_ii)) UpdateScreen();
			}
			else if (verbose) Print_CurieandMoment_List();

			if (script_client_connected) {
				if (!SMesh.contains(meshName)) meshName = SMesh.GetMeshFocus();
				if (!SMesh[meshName]->is_atomistic()) commSocket.SetSendData(commandSpec.PrepareReturnParameters(dynamic_cast<Mesh*>(SMesh[meshName])->GetTcCoupling()));
			}
		}
		break;

		case CMD_TMODEL:
		{
			int tmodeltype;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, tmodeltype);

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
			std::string meshName;

			bool meshName_specified = optional_meshname_check_focusedmeshdefault(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, flag);
			if (!meshName_specified) meshName = SMesh.superMeshHandle;

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
						
						if (!err_hndl.qcall(error, &SuperMesh::SwitchCUDAState, &SMesh, status, cudaDeviceSelect)) {

							cudaEnabled = status;
						}
						else {
							
							SMesh.SwitchCUDAState(false, cudaDeviceSelect);
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

		case CMD_SELCUDADEV:
		{
			int device;

			error = commandSpec.GetParameters(command_fields, device);

			if (!error && device < cudaDeviceVersions.size()) {

#if COMPILECUDA == 1
				if (cudaAvailable) {

					if (cudaDeviceVersions[device].first != __CUDA_ARCH__) error(BERROR_CUDAVERSIONMISMATCH_NCRIT);
					else {

						if (cudaDeviceSelect != device) {

							StopSimulation();

							//if CUDA on, switch it off then back on so everything switches over to the new GPU selection (memory transfer from old GPU to CPU, then to new GPU seamlessly).
							if (cudaEnabled) {

								//first switch off with current device still selected
								err_hndl.qcall(error, &SuperMesh::SwitchCUDAState, &SMesh, false, cudaDeviceSelect);

								//now select new device and turn cuda back on
								cudaDeviceSelect = device;
								err_hndl.qcall(error, &SuperMesh::SwitchCUDAState, &SMesh, true, cudaDeviceSelect);
							}
							else {

								cudaDeviceSelect = device;
							}

							UpdateScreen();
						}
					}
				}
				else error(BERROR_NOTAVAILABLE);
#endif
			}
			else if (verbose) Print_CUDAStatus();

			if (script_client_connected)
				commSocket.SetSendData(commandSpec.PrepareReturnParameters(cudaDeviceSelect));
		}
		break;

		case CMD_OPENMANUAL:
		{
			std::string directory = GetUserDocumentsPath() + boris_data_directory;
			std::string fileName = "BorisManual-v" + ToString(Program_Version) + ".pdf";

			open_file(directory + fileName);
		}
		break;

		case CMD_REFINEROUGHNESS:
		{
			std::string meshName;
			INT3 refine;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, refine);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::SetMeshRoughnessRefinement, &SMesh, meshName, refine)) UpdateScreen();
			}
			else if (verbose) Print_MeshRoughnessRefinement(meshName);
		}
		break;

		case CMD_CLEARROUGHNESS:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::ClearMeshRoughness, &SMesh, meshName)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_ROUGHENMESH:
		{
			std::string meshName;
			double depth;
			std::string side;
			int seed;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, depth, side, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, depth, side); seed = 1; }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, depth); side = "z"; seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::RoughenMeshSides, &SMesh, meshName, side, depth, seed)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;
		
		case CMD_SURFROUGHENJAGGED:
		{
			std::string meshName;
			double depth, spacing;
			std::string sides;
			int seed;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, depth, spacing, seed, sides);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, depth, spacing, seed); }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, depth, spacing); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::RoughenMeshSurfaces_Jagged, &SMesh, meshName, depth, spacing, seed, sides)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GENERATE2DGRAINS:
		{
			std::string meshName;
			double spacing;
			int seed;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, spacing, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, spacing); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::GenerateGrains2D, &SMesh, meshName, spacing, seed)) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GENERATE3DGRAINS:
		{
			std::string meshName;
			double spacing;
			int seed;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, spacing, seed);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, spacing); seed = 1; }

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::GenerateGrains3D, &SMesh, meshName, spacing, seed)) UpdateScreen();
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
			std::string mdbName;

			error = commandSpec.GetParameters(command_fields, mdbName);

			if (!error) {

				if (!err_hndl.call(error, &MaterialsDB::SwitchDataBase, &mdb, mdbName)) {

					UpdateScreen();
				}
			}
			else if (verbose) Print_MaterialsDatabase();
		}
		break;

		case CMD_ADDMDBENTRY:
		{
			std::string meshName, materialName;

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
			std::string materialName;

			error = commandSpec.GetParameters(command_fields, materialName);

			if (!error) {

				err_hndl.call(error, &MaterialsDB::DelMDBEntry, &mdb, materialName);

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_REFRESHMDB:
		{
			err_hndl.call(error, &MaterialsDB::RefreshMDB, &mdb);
		}
		break;

		case CMD_ADDMATERIAL:
		{
			std::string materialName;
			Rect meshRect;

			error = commandSpec.GetParameters(command_fields, materialName, meshRect);

			if (!error) {

				int meshType;

				StopSimulation();

				if (!err_hndl.call(error, &MaterialsDB::LoadMaterial, &mdb, materialName, &meshType)) {

					//material loaded in mdb and meshType available.
					
					std::string meshName = materialName;

					//does this mesh name already exist?
					while (SMesh.contains(meshName)) {

						//if meshName exists then add _#, where # = 1 to start
						//if meshName already contains a termination of the form _# then increase # until meshName available
						size_t pos = meshName.find_last_of('_');

						if (pos == std::string::npos) {

							meshName += std::string("_1");
						}
						else {

							std::string termination = meshName.substr(pos + 1);
							if (has_digits_only(termination)) {

								int number = ToNum(termination);
								number++;
								termination = ToString(number);

								meshName = meshName.substr(0, pos + 1) + termination;
							}
							else meshName += std::string("_1");
						}
					}

					//first make the mesh with the correct rectangle and mesh type
					if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, meshName, (MESH_)meshType, meshRect)) {

						//mesh created, now copy parameter values
						mdb.copy_parameters(*dynamic_cast<Mesh*>(SMesh[meshName]));
					}

					if (script_client_connected) commSocket.SetSendData({ meshName });
				}

				UpdateScreen();

			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETMATERIAL:
		{
			std::string materialName;
			Rect meshRect;

			error = commandSpec.GetParameters(command_fields, materialName, meshRect);

			if (!error) {

				int meshType;

				StopSimulation();

				if (!err_hndl.call(error, &MaterialsDB::LoadMaterial, &mdb, materialName, &meshType)) {

					//material loaded in mdb and meshType available.

					//does this mesh name already exist?
					std::string new_meshName = adjust_meshname(materialName);

					//first make the mesh with the correct rectangle and mesh type
					if (!err_hndl.call(error, &SuperMesh::AddMesh, &SMesh, new_meshName, (MESH_)meshType, meshRect)) {

						//mesh created, now copy parameter values
						mdb.copy_parameters(*dynamic_cast<Mesh*>(SMesh[new_meshName]));

						SMesh.SetMeshFocus(new_meshName);
						UpdateScreen_AutoSet_KeepOrientation();

						delete_all_meshes_except(new_meshName);

						UpdateScreen_AutoSet();

						if (script_client_connected)
							commSocket.SetSendData({ new_meshName });
					}
				}

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_REQMDBSYNC:
		{
			std::string materialName, emailContact;

			error = commandSpec.GetParameters(command_fields, materialName, emailContact);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, materialName); emailContact = "N/A"; }

			if (!error) {

				std::string returnMessage;

				if (!err_hndl.call(error, &MaterialsDB::RequestMDBSync, &mdb, materialName, domain_name, mdb_entry_handler, emailContact, &returnMessage)) {

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
			std::string returnMessage;

			if (!err_hndl.call(error, &MaterialsDB::UpdateMDB, &mdb, domain_name, mdb_update_handler, &returnMessage)) {

				BD.DisplayConsoleMessage("Materials database successfully updated.");
			}
			else BD.DisplayConsoleError(returnMessage);
		}
		break;

		case CMD_SHOWLENGHTS:
		{
			std::string meshName;
			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (verbose) {

				if (SMesh[meshName]->Magnetism_Enabled() && !SMesh[meshName]->is_atomistic()) {
					
					double A = dynamic_cast<Mesh*>(SMesh[meshName])->A;
					double Ms = dynamic_cast<Mesh*>(SMesh[meshName])->Ms;
					double Ku = dynamic_cast<Mesh*>(SMesh[meshName])->K1;
					double D = dynamic_cast<Mesh*>(SMesh[meshName])->D;

					std::string l_ex, l_Bloch("N/A"), l_sky("N/A");

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
			std::string meshName;
			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh[meshName]->Magnetism_Enabled()) {

				if (verbose) BD.DisplayConsoleMessage("M discretisation cells : " + ToString(SMesh[meshName]->n) + " Total cells: " + ToString(SMesh[meshName]->n.dim()));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->n));
			}
			else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
		}
		break;

		case CMD_LOADOVF2MAG:
		{
			std::string meshName, params_string;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, params_string);

			if (!error) {

				StopSimulation();

				double renormalize_value = 0.0;
				std::string fileName;

				//using split_numeric approach since the file name path can contain spaces.
				std::vector<std::string> entries = split_numeric(params_string);

				if (entries.size() == 2) {

					renormalize_value = ToNum(entries[0]);
					fileName = entries[1];
				}
				else fileName = params_string;

				if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				VEC<DBL3> data;

				OVF2 ovf2;
				error = ovf2.Read_OVF2_VEC(ScanFileNameData(fileName), data);

				if (!error) {

					//data loaded correctly, so copy magnetization data to focused mesh.

					if (IsNZ(renormalize_value)) data.renormalize(renormalize_value);

					if (SMesh[meshName]->Magnetism_Enabled()) {

						SMesh[meshName]->SetMagFromData(data);
						UpdateScreen();
					}
					else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2FIELD:
		{
			std::string meshName, fileName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, fileName);
			
			if (!error) {

				StopSimulation();

				if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				if (SMesh[meshName]->Magnetism_Enabled()) {

					if (SMesh[meshName]->IsModuleSet(MOD_ZEEMAN)) {

						error = SMesh[meshName]->CallModuleMethod(&ZeemanBase::SetFieldVEC_FromOVF2, ScanFileNameData(fileName));
					}
					
					UpdateScreen();
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2TEMP:
		{
			std::string meshName, fileName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, fileName);
			
			if (!error) {

				StopSimulation();

				if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				VEC<double> data;

				OVF2 ovf2;
				error = ovf2.Read_OVF2_SCA(ScanFileNameData(fileName), data);

				if (!error) {

					//data loaded correctly, so set temperature from it
					if (SMesh[meshName]->TComputation_Enabled()) {

						SMesh[meshName]->SetTempFromData(data);
						UpdateScreen();
					}
					else err_hndl.show_error(BERROR_NOHEAT, verbose);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2CURR:
		{
			std::string meshName, fileName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, fileName);

			if (!error) {

				StopSimulation();

				if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				VEC<DBL3> data;

				OVF2 ovf2;
				error = ovf2.Read_OVF2_VEC(ScanFileNameData(fileName), data);

				if (!error) {

					//data loaded correctly, so set electric field from it (we've loaded current density so divide by conductivity)

					if (SMesh[meshName]->EComputation_Enabled()) {

						SMesh[meshName]->SetEFromJcData(data);
						UpdateScreen();
					}
					else err_hndl.show_error(BERROR_NOTRANSPORT, verbose);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVEOVF2MAG:
		{
			std::string meshName, parameters;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, parameters);

			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled() && !SMesh[meshName]->is_atomistic()) {

					bool normalize = false;
					std::string data_type = "bin8";
					std::string fileName;

					std::vector<std::string> fields = split(parameters, " ");

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

					double Ms0 = dynamic_cast<Mesh*>(SMesh[meshName])->Ms.get0();
					if (!normalize) Ms0 = 1.0;

					OVF2 ovf2;
					error = ovf2.Write_OVF2_VEC(ScanFileNameData(fileName), dynamic_cast<Mesh*>(SMesh[meshName])->Get_M(), data_type, Ms0);
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVEOVF2:
		{
			std::string meshName, quantity, parameters;

			optional_meshname_check_focusedmeshdefault_quantityname(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, quantity, parameters);

			if (!error) {
				
				std::string data_type = "bin8";

				std::vector<std::string> fields = split(parameters, " ");

				int oparams = 0;

				if (fields[0] == "bin4" || fields[0] == "bin8" || fields[0] == "text") {

					oparams++;
					data_type = fields[0];
				}
				
				if (fields.size() > oparams) {

					//get filename
					std::string fileName = combine(subvec(fields, oparams), " ");

					if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

					if (!err_hndl.call(error, &SuperMesh::SaveOnScreenPhysicalQuantity, &SMesh, 
						meshName, ScanFileNameData(fileName), data_type,
						(MESHDISPLAY_)displayHandles.get_ID_from_value(quantity))) {

						BD.DisplayConsoleMessage("Data saved : " + fileName);
					}
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SAVEOVF2PARAMVAR:
		{
			std::string meshName, parameters;

			optional_meshname_check_focusedmeshdefault(command_fields);
			//expecting meshname (data_type) paramname (directory/)filename after meshname
			error = commandSpec.GetParameters(command_fields, meshName, parameters);

			if (!error) {

				std::string data_type = "bin8";
				std::string fileName, paramname;

				std::vector<std::string> fields = split(parameters, " ");

				int oparams = 0;
				if (fields.size() && (fields[0] == "bin4" || fields[0] == "bin8" || fields[0] == "text")) {

					//data type specified (if not, default stands)
					oparams++;
					data_type = fields[0];
				}

				if (fields.size() > oparams) {

					//get parameter name
					paramname = fields[oparams];

					if (!SMesh[meshName]->contains_param(paramname)) {

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

				if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

				ScanFileNameData(fileName);

				OVF2 ovf2;
				PARAM_ paramID = (PARAM_)SMesh[meshName]->get_meshparam_id(paramname);

				if (SMesh[meshName]->is_paramvarequation_set(paramID)) {

					if (SMesh[meshName]->is_paramvar_scalar(paramID)) {

						VEC<double> s_scaling(SMesh[meshName]->get_paramtype_cellsize(paramID), SMesh[meshName]->meshRect);
						SMesh[meshName]->calculate_meshparam_s_scaling(paramID, s_scaling, SMesh.GetStageTime());

						error = ovf2.Write_OVF2_SCA(fileName, s_scaling, data_type);
					}
					else {

						VEC<DBL3> s_scaling(SMesh[meshName]->get_paramtype_cellsize(paramID), SMesh[meshName]->meshRect);
						SMesh[meshName]->calculate_meshparam_s_scaling(paramID, s_scaling, SMesh.GetStageTime());

						error = ovf2.Write_OVF2_VEC(fileName, s_scaling, data_type);
					}
				}
				else {

					void* s_scaling = SMesh[meshName]->get_meshparam_s_scaling(paramID);

					if (SMesh[meshName]->is_paramvar_scalar(paramID)) error = ovf2.Write_OVF2_SCA(fileName, *reinterpret_cast<VEC<double>*>(s_scaling), data_type);
					else error = ovf2.Write_OVF2_VEC(fileName, *reinterpret_cast<VEC<DBL3>*>(s_scaling), data_type);
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_LOADOVF2DISP:
		{
			std::string meshName, fileName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, fileName);
			if (!error) {

				StopSimulation();

				if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->IsModuleSet(MOD_MELASTIC)) {

					if (GetFileTermination(fileName) != ".ovf") fileName += ".ovf";
					if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;

					ScanFileNameData(fileName);

					if (!error) {

						error = SMesh[meshName]->CallModuleMethod(&MElastic::Load_Displacement_OVF2, fileName);
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
			std::string meshName, parameters;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, parameters);
			if (!error) {

				StopSimulation();

				if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->IsModuleSet(MOD_MELASTIC)) {

					std::vector<std::string> fields = split(parameters, ".ovf");

					if (fields.size() >= 2) {

						std::string fileName_diag, fileName_odiag;

						fileName_diag = fields[0] + ".ovf";
						fileName_odiag = fields[1] + ".ovf";

						if (!GetFilenameDirectory(fileName_diag).length()) fileName_diag = directory + fileName_diag;
						if (!GetFilenameDirectory(fileName_odiag).length()) fileName_odiag = directory + fileName_odiag;

						if (!error) {

							error = SMesh[meshName]->CallModuleMethod(&MElastic::Load_Strain_OVF2, ScanFileNameData(fileName_diag), ScanFileNameData(fileName_odiag));
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
			std::string userconstant_name;
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
			std::string userconstant_name;
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
			std::ofstream bdout;
			bdout.open(errorlog_fileName.c_str(), std::ios::out);
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

		case CMD_THREADS:
		{
			int threads;

			error = commandSpec.GetParameters(command_fields, threads);

			if (!error) {

				StopSimulation();

				OmpThreads = threads;
				if (!OmpThreads || OmpThreads > omp_get_num_procs()) OmpThreads = omp_get_num_procs();
				Save_Startup_Flags();

				UpdateScreen();
			}
			else if (verbose && error == BERROR_PARAMOUTOFBOUNDS) PrintCommandUsage(command_name);
			else if (verbose) Print_Threads();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(OmpThreads));
		}
		break;

		case CMD_SERVERPORT:
		{
			int port;

			error = commandSpec.GetParameters(command_fields, port);

			if (!error) {

				StopSimulation();

				server_port = ToString(port);
				commSocket.Change_Port(server_port);
				Save_Startup_Flags();

				UpdateScreen();
			}
			else if (verbose && error == BERROR_PARAMOUTOFBOUNDS) PrintCommandUsage(command_name);
			else if (verbose) Print_ServerInfo();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(server_port));
		}
		break;

		case CMD_SERVERPWD:
		{
			std::string password;

			error = commandSpec.GetParameters(command_fields, password);

			if (!error) {

				StopSimulation();

				server_pwd = password;
				commSocket.Change_Password(server_pwd);
				Save_Startup_Flags();

				UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;
			
		case CMD_SERVERSLEEPMS:
		{
			int sleep_ms;

			error = commandSpec.GetParameters(command_fields, sleep_ms);

			if (!error) {

				StopSimulation();

				server_recv_sleep_ms = sleep_ms;
				commSocket.Change_RecvSleep(server_recv_sleep_ms);
				Save_Startup_Flags();

				UpdateScreen();
			}
			else if (verbose && error == BERROR_PARAMOUTOFBOUNDS) PrintCommandUsage(command_name);
			else if (verbose) Print_ServerInfo();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(server_recv_sleep_ms));
		}
		break;

		case CMD_NEWINSTANCE:
		{
			int port;
			int cudaDevice;
			std::string password;

			error = commandSpec.GetParameters(command_fields, port, cudaDevice, password);
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, port, cudaDevice); password = ""; }
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, port); cudaDevice = -1; password = ""; }

			if (!error) {

				std::string boris_path = GetExeDirectory();

				if (password.length())
					open_file(boris_path + progName, "-p " + ToString(port) + " " + password + " -g " + ToString(cudaDevice) + " -w back");
				else
					open_file(boris_path + progName, "-p " + ToString(port) + " -g " + ToString(cudaDevice) + " -w back");
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		//-------------------------------------------VERSION UPDATE-------------------------------------------

		case CMD_VERSIONUPDATE:
		{
			std::string action;
			int targetver = 0;

			error = commandSpec.GetParameters(command_fields, action, targetver);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, action); targetver = 0; }

			if (!error) {

				if (action == "check") CheckDeltaUpdate();
				if (action == "download") DownloadDeltaUpdate(targetver);
				if (action == "install") InstallDeltaUpdate(targetver);
				if (action == "rollback") RollbackUpdate(targetver);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SHOWTC:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->is_atomistic()) {

				double Tc = dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Show_Transition_Temperature();
				if (verbose) BD.DisplayConsoleMessage("Tc = " + ToString(Tc));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Tc));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SHOWMS:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->is_atomistic()) {

				double Ms = dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Show_Ms();
				if (verbose) BD.DisplayConsoleMessage("Ms = " + ToString(Ms));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Ms));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SHOWA:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->is_atomistic()) {

				double A = dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Show_A();
				if (verbose) BD.DisplayConsoleMessage("A = " + ToString(A));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(A));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SHOWK:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName);

			if (SMesh[meshName]->Magnetism_Enabled() && SMesh[meshName]->is_atomistic()) {

				double Ku = dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Show_Ku();
				if (verbose) BD.DisplayConsoleMessage("Ku = " + ToString(Ku));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Ku));
			}
			else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
		}
		break;

		case CMD_SKYPOSDMUL:
		{
			double multiplier;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, multiplier);

			if (!error) {

				SMesh[meshName]->Set_skypos_dmul(multiplier);
				UpdateScreen();
			}
			else if (verbose) Print_skypos_dmul();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->Get_skypos_dmul()));
		}
		break;

		case CMD_DWPOS_COMPONENT:
		{
			int component;

			error = commandSpec.GetParameters(command_fields, component);

			if (!error) {

				SMesh.Set_DWPos_Component(component);
				UpdateScreen();
			}
			else if (verbose) Print_DWPos_Component();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.Get_DWPos_Component()));
		}
		break;

		case CMD_MCSERIAL:
		{
			int status;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, status);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::Set_MonteCarlo_Serial, &SMesh, (bool)status, meshName)) UpdateScreen();
			}
			else if (verbose) Print_MCSettings();
		}
		break;

		case CMD_MCCONSTRAIN:
		{
			DBL3 value;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, value);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::Set_MonteCarlo_Constrained, &SMesh, !value.IsNull(), value, meshName)) UpdateScreen();
			}
			else if (verbose) Print_MCSettings();
		}
		break;

		case CMD_MCDISABLE:
		{
			bool status;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, status);

			if (!error) {

				if (!err_hndl.qcall(error, &SuperMesh::Set_MonteCarlo_Disabled, &SMesh, status, meshName)) UpdateScreen();
			}
			else if (verbose) Print_MCSettings();
		}
		break;

		case CMD_MCCOMPUTEFIELDS:
		{
			bool status;

			error = commandSpec.GetParameters(command_fields, status);

			if (!error) {

				SMesh.Set_MonteCarlo_ComputeFields(status);

				UpdateScreen();
			}
			else if (verbose) Print_MCSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.Get_MonteCarlo_ComputeFields()));
		}
		break;

		case CMD_MCCONEANGLELIMITS:
		{
			DBL2 cone_angles;

			error = commandSpec.GetParameters(command_fields, cone_angles);

			if (!error) {

				SMesh.Set_MonteCarlo_ConeAngleLimits(cone_angles);

				UpdateScreen();
			}
			else if (verbose) Print_MCSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh.Get_MonteCarlo_ConeAngleLimits()));
		}
		break;

		case CMD_SHAPEMOD_ROT:
		{
			DBL3 rotation;

			error = commandSpec.GetParameters(command_fields, rotation);

			if (!error) {

				shape_rotation = rotation;
				UpdateScreen();
			}
			else if (verbose) Print_ShapeSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(shape_rotation));
		}
		break;

		case CMD_SHAPEMOD_REP:
		{
			INT3 repetitions;

			error = commandSpec.GetParameters(command_fields, repetitions);

			if (!error) {

				shape_repetitions = repetitions;
				UpdateScreen();
			}
			else if (verbose) Print_ShapeSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(shape_repetitions));
		}
		break;

		case CMD_SHAPEMOD_DISP:
		{
			DBL3 displacement;

			error = commandSpec.GetParameters(command_fields, displacement);

			if (!error) {

				shape_displacement = displacement;
				UpdateScreen();
			}
			else if (verbose) Print_ShapeSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(shape_displacement));
		}
		break;

		case CMD_SHAPEMOD_METHOD: 
		{
			std::string method;

			error = commandSpec.GetParameters(command_fields, method);

			if (!error) {
				
				//make sure method name is correct
				if (MeshShape().set_method(method) != MSHAPEMETHOD_NONE) {

					shape_method = method;
					UpdateScreen();
				}
				else error(BERROR_INCORRECTSTRING);
			}
			else if (verbose) Print_ShapeSettings();

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(shape_method));
		}
		break;

		case CMD_SHAPE_DISK:
		{
			DBL2 diameters;
			DBL2 position;
			DBL2 z_coords;
			std::string meshName;
			DBL3 dimensions, centre_pos;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, diameters, position, z_coords);
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, meshName, diameters, position); z_coords = DBL2(0, SMesh[meshName]->GetMeshDimensions().z); }

			if (!error) {

				dimensions = DBL3(diameters.x, diameters.y, z_coords.j - z_coords.i);
				centre_pos = DBL3(position.x, position.y, (z_coords.j + z_coords.i) / 2);

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_disk, SMesh[meshName], MeshShape(MSHAPE_DISK, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_DISK).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_RECT:
		{
			DBL2 lengths;
			DBL2 position;
			DBL2 z_coords;
			std::string meshName = SMesh.GetMeshFocus();
			DBL3 dimensions, centre_pos;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, lengths, position, z_coords);
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, meshName, lengths, position); z_coords = DBL2(0, SMesh[meshName]->GetMeshDimensions().z); }

			if (!error) {

				dimensions = DBL3(lengths.x, lengths.y, z_coords.j - z_coords.i);
				centre_pos = DBL3(position.x, position.y, (z_coords.j + z_coords.i) / 2);

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_rect, SMesh[meshName], MeshShape(MSHAPE_DISK, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_RECT).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_TRIANGLE:
		{
			DBL2 lengths;
			DBL2 position;
			DBL2 z_coords;
			std::string meshName = SMesh.GetMeshFocus();
			DBL3 dimensions, centre_pos;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, lengths, position, z_coords);
			if (error) { error.reset() = commandSpec.GetParameters(command_fields, meshName, lengths, position); z_coords = DBL2(0, SMesh[meshName]->GetMeshDimensions().z); }

			if (!error) {

				dimensions = DBL3(lengths.x, lengths.y, z_coords.j - z_coords.i);
				centre_pos = DBL3(position.x, position.y, (z_coords.j + z_coords.i) / 2);

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_triangle, SMesh[meshName], MeshShape(MSHAPE_TRIANGLE, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_TRIANGLE).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_ELLIPSOID:
		{
			DBL3 dimensions;
			DBL3 centre_pos;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dimensions, centre_pos);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_ellipsoid, SMesh[meshName], MeshShape(MSHAPE_ELLIPSOID, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_ELLIPSOID).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_PYRAMID:
		{
			DBL3 dimensions;
			DBL3 centre_pos;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dimensions, centre_pos);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_pyramid, SMesh[meshName], MeshShape(MSHAPE_PYRAMID, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_PYRAMID).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_TETRAHEDRON:
		{
			DBL3 dimensions;
			DBL3 centre_pos;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dimensions, centre_pos);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_tetrahedron, SMesh[meshName], MeshShape(MSHAPE_TETRAHEDRON, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_TETRAHEDRON).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_CONE:
		{
			DBL3 dimensions;
			DBL3 centre_pos;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dimensions, centre_pos);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_cone, SMesh[meshName], MeshShape(MSHAPE_CONE, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_CONE).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_TORUS:
		{
			DBL3 dimensions;
			DBL3 centre_pos;
			std::string meshName = SMesh.GetMeshFocus();

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dimensions, centre_pos);

			if (!error) {

				StopSimulation();
				if (!err_hndl.qcall(error, &MeshBase::shape_torus, SMesh[meshName], MeshShape(MSHAPE_TORUS, dimensions, centre_pos, shape_rotation * PI / 180, shape_repetitions, shape_displacement, MeshShape().set_method(shape_method)))) UpdateScreen();
			}
			else if (verbose) PrintCommandUsage(command_name);

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(MeshShape(MSHAPE_TORUS).name(), dimensions, centre_pos, shape_rotation, shape_repetitions, shape_displacement, shape_method));
		}
		break;

		case CMD_SHAPE_SET:
		{
			//first get meshName, then remove it
			std::string meshName;
			optional_meshname_check_focusedmeshdefault(command_fields);
			meshName = command_fields[0];
			command_fields.erase(command_fields.begin());

			//name : 1, dimensions: 3, centre_pos: 3, shape_rotation: 3, shape_repetitions: 3, shape_displacement: 3, shape_method: 1. Total space-sparated fields per shape: 17
			int num_fields_per_shape = 17;
			int num_shapes = command_fields.size() / num_fields_per_shape;
			
			if (num_shapes >= 1) {
				
				StopSimulation();

				std::string shape_name;
				DBL3 dimensions;
				DBL3 centre_pos;
				DBL3 rotation;
				INT3 repetitions;
				DBL3 displacement;
				std::string method;

				std::vector<MeshShape> shapes;
				for (int shape_idx = 0; shape_idx < num_shapes; shape_idx++) {

					std::vector<std::string> shape_fields = subvec(command_fields, shape_idx * num_fields_per_shape, (shape_idx + 1) * num_fields_per_shape);
					error = commandSpec.GetParameters(shape_fields, shape_name, dimensions, centre_pos, rotation, repetitions, displacement, method);

					MSHAPE_ shape_id = MeshShape().set_shape(shape_name);
					MSHAPEMETHOD_ method_id = MeshShape().set_method(method);

					//make sure shape and method names are correct
					if (shape_id == MSHAPE_NONE || method_id == MSHAPEMETHOD_NONE) {

						error(BERROR_INCORRECTNAME);
						break;
					}

					if (!error) shapes.push_back(MeshShape(shape_id, dimensions, centre_pos, rotation * PI / 180, repetitions, displacement, method_id));
					else break;
				}

				if (!error) {
				
					if (!err_hndl.qcall(error, &MeshBase::shape_set, SMesh[meshName], shapes)) UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SHAPE_GET:
		{
			//first get meshName, then remove it
			std::string meshName, quantity;
			optional_meshname_check_focusedmeshdefault_quantityname(command_fields, true);
			meshName = command_fields[0];
			command_fields.erase(command_fields.begin());
			quantity = command_fields[0];
			command_fields.erase(command_fields.begin());

			//name : 1, dimensions: 3, centre_pos: 3, shape_rotation: 3, shape_repetitions: 3, shape_displacement: 3, shape_method: 1. Total space-sparated fields per shape: 17
			int num_fields_per_shape = 17;
			int num_shapes = command_fields.size() / num_fields_per_shape;

			if (num_shapes >= 1) {

				StopSimulation();

				std::string shape_name;
				DBL3 dimensions;
				DBL3 centre_pos;
				DBL3 rotation;
				INT3 repetitions;
				DBL3 displacement;
				std::string method;

				std::vector<MeshShape> shapes;
				for (int shape_idx = 0; shape_idx < num_shapes; shape_idx++) {

					std::vector<std::string> shape_fields = subvec(command_fields, shape_idx * num_fields_per_shape, (shape_idx + 1) * num_fields_per_shape);
					error = commandSpec.GetParameters(shape_fields, shape_name, dimensions, centre_pos, rotation, repetitions, displacement, method);

					MSHAPE_ shape_id = MeshShape().set_shape(shape_name);
					MSHAPEMETHOD_ method_id = MeshShape().set_method(method);

					//make sure shape and method names are correct
					if (shape_id == MSHAPE_NONE || method_id == MSHAPEMETHOD_NONE) {

						error(BERROR_INCORRECTNAME);
						break;
					}

					if (!error) shapes.push_back(MeshShape(shape_id, dimensions, centre_pos, rotation * PI / 180, repetitions, displacement, method_id));
					else break;
				}

				if (!error) {

					Any value = SMesh.GetAverageDisplayedMeshValue(
						meshName, Rect(), shapes,
						(MESHDISPLAY_)displayHandles.get_ID_from_value(quantity));

					if (script_client_connected) {

						//Longer version so it compiles with C++14
						if (value.is_type(btype_info<double>()) || value.is_type(btype_info<float>())) {

							double value_converted = value;
							commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
						}
						else if (value.is_type(btype_info<DBL2>()) || value.is_type(btype_info<FLT2>())) {

							DBL2 value_converted = value;
							commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
						}
						else if (value.is_type(btype_info<DBL3>()) || value.is_type(btype_info<FLT3>())) {

							DBL3 value_converted = value;
							commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
						}
					}
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SETSHAPEANGLE:
		{
			//first get meshName, then remove it
			std::string meshName;
			optional_meshname_check_focusedmeshdefault(command_fields);
			meshName = command_fields[0];
			command_fields.erase(command_fields.begin());

			//name : 1, dimensions: 3, centre_pos: 3, shape_rotation: 3, shape_repetitions: 3, shape_displacement: 3, shape_method: 1. Total space-sparated fields per shape: 17
			//two extra parameters for angle values
			int num_fields_per_shape = 17;
			int num_shapes = ((int)command_fields.size() - 2) / num_fields_per_shape;

			if (num_shapes >= 1) {

				StopSimulation();

				std::string shape_name;
				DBL3 dimensions;
				DBL3 centre_pos;
				DBL3 rotation;
				INT3 repetitions;
				DBL3 displacement;
				std::string method;
				DBL2 angle;

				std::vector<MeshShape> shapes;
				for (int shape_idx = 0; shape_idx < num_shapes; shape_idx++) {

					std::vector<std::string> shape_fields = subvec(command_fields, shape_idx * num_fields_per_shape, (shape_idx + 1) * num_fields_per_shape);
					JoinToVector(shape_fields, subvec(command_fields, (int)command_fields.size() - 2));
					error = commandSpec.GetParameters(shape_fields, shape_name, dimensions, centre_pos, rotation, repetitions, displacement, method, angle);

					MSHAPE_ shape_id = MeshShape().set_shape(shape_name);
					MSHAPEMETHOD_ method_id = MeshShape().set_method(method);
					
					//make sure shape and method names are correct
					if (shape_id == MSHAPE_NONE || method_id == MSHAPEMETHOD_NONE) {

						error(BERROR_INCORRECTNAME);
						break;
					}
					
					if (!error) shapes.push_back(MeshShape(shape_id, dimensions, centre_pos, rotation * PI / 180, repetitions, displacement, method_id));
					else break;
				}

				if (!error) {

					if (!err_hndl.qcall(error, &SuperMesh::SetMagAngle_Shape, &SMesh, meshName, angle.i, angle.j, shapes)) {

						UpdateScreen();
					}
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_SHAPE_SETPARAM:
		{
			//first get meshName, then remove it
			std::string meshName;
			optional_meshname_check_focusedmeshdefault(command_fields);
			meshName = command_fields[0];
			command_fields.erase(command_fields.begin());

			//name : 1, dimensions: 3, centre_pos: 3, shape_rotation: 3, shape_repetitions: 3, shape_displacement: 3, shape_method: 1. Total space-sparated fields per shape: 17
			//extra parameter for paramName and also for value
			int num_fields_per_shape = 17;
			int num_shapes = ((int)command_fields.size() - 1) / num_fields_per_shape;
			int num_fields_value = (int)command_fields.size() - num_shapes * num_fields_per_shape - 1;

			if (num_shapes >= 1 && num_fields_value >= 1) {

				StopSimulation();

				std::string paramName;
				std::string shape_name;
				DBL3 dimensions;
				DBL3 centre_pos;
				DBL3 rotation;
				INT3 repetitions;
				DBL3 displacement;
				std::string method;
				std::string value_text;

				std::vector<MeshShape> shapes;
				for (int shape_idx = 0; shape_idx < num_shapes; shape_idx++) {

					std::vector<std::string> shape_fields = subvec(command_fields, 0, 1);
					JoinToVector(shape_fields, subvec(command_fields, shape_idx * num_fields_per_shape + 1, (shape_idx + 1) * num_fields_per_shape + 1));
					JoinToVector(shape_fields, subvec(command_fields, (int)command_fields.size() - num_fields_value));
					error = commandSpec.GetParameters(shape_fields, paramName, shape_name, dimensions, centre_pos, rotation, repetitions, displacement, method, value_text);

					MSHAPE_ shape_id = MeshShape().set_shape(shape_name);
					MSHAPEMETHOD_ method_id = MeshShape().set_method(method);

					//make sure shape and method names are correct
					if (shape_id == MSHAPE_NONE || method_id == MSHAPEMETHOD_NONE) {

						error(BERROR_INCORRECTNAME);
						break;
					}

					if (!error) shapes.push_back(MeshShape(shape_id, dimensions, centre_pos, rotation * PI / 180, repetitions, displacement, method_id));
					else break;
				}

				if (!error) {

					error = SMesh.set_meshparam_shape(meshName, paramName, shapes, value_text);
					if (!error) SMesh.set_meshparamvar_display(meshName, paramName);

					UpdateScreen();
				}
				else if (verbose) PrintCommandUsage(command_name);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_RUNCOMMBUFFER:
			while (single_call_launch(&Simulation::RunCommandBuffer, THREAD_HANDLEMESSAGE3) != THREAD_HANDLEMESSAGE3);
			break;

		case CMD_CLEARCOMMBUFFER:
			command_buffer.clear();
			break;

		case CMD_BUFFERCOMMAND:
		{
			std::string command_with_params;

			error = commandSpec.GetParameters(command_fields, command_with_params);

			if (!error) {

				//buffer the command with current verbosity level if recognized as a command
				//NOTE : must stop CMD_RUNCOMMBUFFER, CMD_CLEARCOMMBUFFER and CMD_BUFFERCOMMAND from themselves being buffered, to avoid causing program hangups if for whatever reason a user decides buffering these is a good idea!

				std::vector<std::string> command_fields = split(command_with_params, " ");
				std::string command_name = command_fields.front();

				if (commands.has_key(command_name)) {

					CommandSpecifier commandSpec = commands[command_name];

					if (commandSpec.cmdId != CMD_RUNCOMMBUFFER && commandSpec.cmdId != CMD_CLEARCOMMBUFFER && commandSpec.cmdId != CMD_BUFFERCOMMAND) {

						if (verbose) command_buffer.push_back(command_with_params);
						else command_buffer.push_back("~" + command_with_params);
					}
				}
			}
			else if (verbose) {

				if (!command_buffer.size()) BD.DisplayConsoleMessage("Command buffer is empty.");

				for (int idx = 0; idx < command_buffer.size(); idx++) {

					std::vector<std::string> fields = split(command_buffer[idx]);
				
					std::string message;
					for (int fidx = 0; fidx < fields.size(); fidx++) {

						if (fidx == 0) message += "[tc1,1,1,1/tc]<b>" + fields[0] + "</b>";
						else message += " <i>" + fields[fidx] + "</i>";
					}

					BD.DisplayFormattedConsoleMessage(message);
				}
			}
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
			std::string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric std::string (the file path and file name) then a numeric std::string.
			std::vector<std::string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				std::string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				std::string indexes_string = entries.back();

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
				if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

				int rows_read;

				error = dpArr.load_arrays(ScanFileNameData(fileName), vec_convert<int, std::string>(split(indexes_string)), &rows_read);

				if (verbose && !error) BD.DisplayConsoleMessage("Number of rows read : " + ToString(rows_read));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVE:
		{
			std::string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric std::string (the file path and file name) then a numeric std::string.
			std::vector<std::string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				std::string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				std::string indexes_string = entries.back();

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
				if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

				error = dpArr.save_arrays(ScanFileNameData(fileName), vec_convert<int, std::string>(split(indexes_string)));

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVEASROW:
		{
			std::string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric std::string (the file path and file name) then a numeric std::string.
			std::vector<std::string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				std::string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				std::string index_string = entries.back();

				int dp_index = ToNum(index_string);

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
				if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

				error = dpArr.save_array_transposed(ScanFileNameData(fileName), dp_index);

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVEAPPENDASROW:
		{
			std::string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric std::string (the file path and file name) then a numeric std::string.
			std::vector<std::string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				std::string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				std::string index_string = entries.back();

				int dp_index = ToNum(index_string);

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
				if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

				error = dpArr.save_array_transposed(ScanFileNameData(fileName), dp_index, true);

				if (verbose && !error) BD.DisplayConsoleMessage("Data saved in : " + fileName);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_NEWFILE:
		{
			std::string fileName;

			error = commandSpec.GetParameters(command_fields, fileName);

			if (!error) {

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
				if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

				std::ofstream bdout;
				bdout.open(ScanFileNameData(fileName), std::ios::out);
				bdout.close();
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SAVEAPPEND:
		{
			std::string fileName_indexes_string;

			error = commandSpec.GetParameters(command_fields, fileName_indexes_string);

			//using split_numeric approach since the file name path can contain spaces. If parameters correct then we'll have first a non-numeric std::string (the file path and file name) then a numeric std::string.
			std::vector<std::string> entries = split_numeric(fileName_indexes_string);
			if (entries.size() < 2) error(BERROR_PARAMMISMATCH);

			if (!error) {

				std::string fileName = combine(subvec(entries, 0, entries.size() - 1), " ");
				std::string indexes_string = entries.back();

				if (!GetFilenameDirectory(fileName).length()) fileName = directory + fileName;
				if (GetFileTermination(fileName) != ".txt") fileName += ".txt";

				error = dpArr.save_arrays(ScanFileNameData(fileName), vec_convert<int, std::string>(split(indexes_string)), true);

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
			DBL3 start, end, stencil;
			double step;
			int arr_idx;
			std::string meshName, quantity;

			optional_meshname_check_focusedmeshdefault_quantityname(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, quantity, start, end, step, arr_idx, stencil);
			if (error) { 
				error.reset() = commandSpec.GetParameters(command_fields, meshName, quantity, start, end, step, arr_idx);
				stencil = DBL3(); 
			}

			if (!error && start != end) {

				int num_points = round((end - start).norm() / step) + 1;

				std::vector<double> position(num_points);
				std::vector<double>* pprofile_dbl = nullptr;
				std::vector<DBL3>* pprofile_dbl3 = nullptr;

				SMesh.GetPhysicalQuantityProfile(
					start, end, step, stencil, 
					pprofile_dbl3, pprofile_dbl, 
					meshName, arr_idx < 0, false,
					(MESHDISPLAY_)displayHandles.get_ID_from_value(quantity));

				if (arr_idx >= 0) {

#pragma omp parallel for
					for (int idx = 0; idx < num_points; idx++) position[idx] = idx * step;

					dpArr.set_array(arr_idx, position);
					if (pprofile_dbl) { if (!dpArr.set_array(arr_idx + 1, *pprofile_dbl)) error(BERROR_INCORRECTARRAYS); }
					else if (pprofile_dbl3) { error = dpArr.set_arrays(arr_idx + 1, *pprofile_dbl3); }
					if (verbose && !error && (pprofile_dbl || pprofile_dbl3)) BD.DisplayConsoleMessage("Path extracted.");
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_GETAVERAGEDPROFILE:
		{
			DBL3 start, end;
			double step;
			int arr_idx;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, start, end, step, arr_idx);

			if (!error && start != end) {

				int num_points = round((end - start).norm() / step) + 1;

				std::vector<double> position(num_points);
				std::vector<double>* pprofile_dbl = nullptr;
				std::vector<DBL3>* pprofile_dbl3 = nullptr;

				SMesh.GetPhysicalQuantityProfile(start, end, step, DBL3(), pprofile_dbl3, pprofile_dbl, meshName, false, true);

#pragma omp parallel for
				for (int idx = 0; idx < num_points; idx++) position[idx] = idx * step;

				dpArr.set_array(arr_idx, position);
				if (pprofile_dbl) { if (!dpArr.set_array(arr_idx + 1, *pprofile_dbl)) error(BERROR_INCORRECTARRAYS); }
				else if (pprofile_dbl3) { error = dpArr.set_arrays(arr_idx + 1, *pprofile_dbl3); }
				if (verbose && !error && (pprofile_dbl || pprofile_dbl3)) BD.DisplayConsoleMessage("Average path extracted. Averaging reset.");
			}
			else if (verbose) BD.DisplayConsoleMessage("Current number of averages : " + ToString(SMesh[meshName]->GetProfileAverages()));

			if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(SMesh[meshName]->GetProfileAverages()));
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
			std::string meshName, quantity;
			Rect rect;
			int dp_index;

			optional_meshname_check_focusedmeshdefault_quantityname(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, quantity, rect, dp_index);
			if (error == BERROR_PARAMMISMATCH) { 
				error.reset() = commandSpec.GetParameters(command_fields, meshName, quantity, rect);
				dp_index = -1;
			}
			if (error == BERROR_PARAMMISMATCH) {
				error.reset() = commandSpec.GetParameters(command_fields, meshName, quantity);
				rect = Rect(); dp_index = -1; }

			if (!error) {
				
				Any value = SMesh.GetAverageDisplayedMeshValue(
					meshName, rect, {},
					(MESHDISPLAY_)displayHandles.get_ID_from_value(quantity));

				if (verbose) BD.DisplayConsoleMessage("Average value = " + value.convert_to_string());

				if (script_client_connected || dp_index >= 0) {

					//Longer version so it compiles with C++14
					if (value.is_type(btype_info<double>()) || value.is_type(btype_info<float>())) {
						
						double value_converted = value;

						if (dp_index >= 0) dpArr.push_value(dp_index, value_converted);
						if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL2>()) || value.is_type(btype_info<FLT2>())) {
						
						DBL2 value_converted = value;

						if (dp_index >= 0) dpArr.push_value(dp_index, value_converted);
						if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL3>()) || value.is_type(btype_info<FLT3>())) {
						
						DBL3 value_converted = value;

						if (dp_index >= 0) dpArr.push_value(dp_index, value_converted);
						if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
				}
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_GETVALUE:
		{
			std::string meshName, quantity;
			DBL3 position;

			optional_meshname_check_focusedmeshdefault_quantityname(command_fields, true);
			error = commandSpec.GetParameters(command_fields, meshName, quantity, position);

			if (!error) {

				Any value = SMesh.GetAverageDisplayedMeshValue(
					meshName, Rect(position, position), {}, 
					(MESHDISPLAY_)displayHandles.get_ID_from_value(quantity));

				if (verbose) BD.DisplayConsoleMessage("Value = " + value.convert_to_string());

				if (script_client_connected) {

					//Longer version so it compiles with C++14
					if (value.is_type(btype_info<double>()) || value.is_type(btype_info<float>())) {

						double value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL2>()) || value.is_type(btype_info<FLT2>())) {

						DBL2 value_converted = value;
						commSocket.SetSendData(commandSpec.PrepareReturnParameters(value_converted));
					}
					else if (value.is_type(btype_info<DBL3>()) || value.is_type(btype_info<FLT3>())) {

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
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, x, y, radius);
			if (error == BERROR_PARAMMISMATCH) { 
				
				error.reset(); 
				x = 0.0; y = 0.0;
				DBL3 meshDims = SMesh[meshName]->GetMeshDimensions();
				radius = sqrt(meshDims.x * meshDims.x + meshDims.y * meshDims.y);
			}

			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled()) {

					double Q = 0.0;

					if (!SMesh[meshName]->is_atomistic()) {

						error = dpArr.get_topological_charge(dynamic_cast<Mesh*>(SMesh[meshName])->Get_M(), x, y, radius, &Q);
					}
					else {

						error = dpArr.get_topological_charge(dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Get_M1(), x, y, radius, &Q);
					}

					if (!error) {

						if (verbose) BD.DisplayConsoleMessage("Q = " + ToString(Q));

						if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Q));
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
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, x, y, radius);
			if (error == BERROR_PARAMMISMATCH) {

				error.reset();
				x = 0.0; y = 0.0;
				DBL3 meshDims = SMesh[meshName]->GetMeshDimensions();
				radius = sqrt(meshDims.x * meshDims.x + meshDims.y * meshDims.y);
			}

			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled()) {

					double Q = 0.0;

					if (!SMesh[meshName]->is_atomistic()) {

						error = dpArr.count_skyrmions(dynamic_cast<Mesh*>(SMesh[meshName])->Get_M(), x, y, radius, &Q);
					}
					else {

						error = dpArr.count_skyrmions(dynamic_cast<Atom_Mesh*>(SMesh[meshName])->Get_M1(), x, y, radius, &Q);
					}

					if (!error) {

						if (verbose) BD.DisplayConsoleMessage("Qmag = " + ToString(Q));

						if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(Q));
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_HISTOGRAM:
		{
			int num_bins = 0;
			INT3 macrocell_dims = INT3(1);
			double min = 0.0, max = 0.0;
			int dp_x, dp_y;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims, num_bins, min, max);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims); num_bins = 0; min = 0.0; max = 0.0; }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y); macrocell_dims = INT3(1); num_bins = 0; min = 0.0; max = 0.0; }
			
			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled()) {

					if (!dpArr.GoodArrays_Unique(dp_x, dp_y)) error(BERROR_INCORRECTARRAYS);
					else {

						if (SMesh[meshName]->Get_Histogram(dpArr[dp_x], dpArr[dp_y], num_bins, min, max, macrocell_dims))
							if (verbose) BD.DisplayConsoleMessage("Histogram calculated.");
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_THAVHISTOGRAM:
		{
			int num_bins = 0;
			INT3 macrocell_dims = INT3(1);
			double min = 0.0, max = 0.0;
			int dp_x, dp_y;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims, num_bins, min, max);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims); num_bins = 0; min = 0.0; max = 0.0; }

			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled()) {

					if (!dpArr.GoodArrays_Unique(dp_x, dp_y)) error(BERROR_INCORRECTARRAYS);
					else {

						if (SMesh[meshName]->Get_ThAvHistogram(dpArr[dp_x], dpArr[dp_y], num_bins, min, max, macrocell_dims))
							if (verbose) BD.DisplayConsoleMessage("Histogram calculated.");
					}
				}
				else err_hndl.show_error(BERROR_NOTMAGNETIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_ANGHISTOGRAM:
		{
			int num_bins = 0;
			INT3 macrocell_dims = INT3(1);
			DBL3 ndir;
			double min = 0.0, max = 0.0;
			int dp_x, dp_y;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims, ndir, num_bins, min, max);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims, ndir); num_bins = 0; min = 0.0; max = 0.0; }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims); ndir = DBL3(); num_bins = 0; min = 0.0; max = 0.0; }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y); macrocell_dims = INT3(1); ndir = DBL3(); num_bins = 0; min = 0.0; max = 0.0; }

			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled()) {

					if (!dpArr.GoodArrays_Unique(dp_x, dp_y)) error(BERROR_INCORRECTARRAYS);
					else {

						if (SMesh[meshName]->Get_AngHistogram(dpArr[dp_x], dpArr[dp_y], num_bins, min, max, macrocell_dims, ndir))
							if (verbose) BD.DisplayConsoleMessage("Angular deviation histogram calculated.");
					}
				}
				else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_THAVANGHISTOGRAM:
		{
			int num_bins = 0;
			INT3 macrocell_dims = INT3(1);
			DBL3 ndir;
			double min = 0.0, max = 0.0;
			int dp_x, dp_y;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims, ndir, num_bins, min, max);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims, ndir); num_bins = 0; min = 0.0; max = 0.0; }
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, macrocell_dims); ndir = DBL3(); num_bins = 0; min = 0.0; max = 0.0; }

			if (!error) {

				if (SMesh[meshName]->Magnetism_Enabled()) {

					if (!dpArr.GoodArrays_Unique(dp_x, dp_y)) error(BERROR_INCORRECTARRAYS);
					else {

						if (SMesh[meshName]->Get_ThAvAngHistogram(dpArr[dp_x], dpArr[dp_y], num_bins, min, max, macrocell_dims, ndir))
							if (verbose) BD.DisplayConsoleMessage("Angular deviation histogram calculated.");
					}
				}
				else err_hndl.show_error(BERROR_NOTATOMISTIC, verbose);
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_HISTOGRAM2:
		{
			int num_bins = 0;
			double min = 0.0, max = 0.0;
			double M2 = 0.0, deltaM2 = 0.0;
			int dp_x, dp_y;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y, num_bins, min, max, M2, deltaM2);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, dp_x, dp_y); num_bins = 0; min = 0.0; max = 0.0; M2 = 0.0; deltaM2 = 0.0; }

			if (!error) {
				
				if (SMesh[meshName]->GetMeshType() == MESH_ANTIFERROMAGNETIC && !SMesh[meshName]->is_atomistic()) {

					if (IsZ(M2) && IsZ(deltaM2)) {

						Mesh* pMesh = dynamic_cast<Mesh*>(SMesh[meshName]);
						if (pMesh) {

							DBL2 Ms_AFM = pMesh->Ms_AFM;
							//equilibrium magnetization on sub-lattice B at set mesh base temperature
							M2 = Ms_AFM.j;
							deltaM2 = M2 * 0.01;
						}
					}

					error = dpArr.calculate_histogram2(
						dynamic_cast<Mesh*>(SMesh[meshName])->Get_M(),
						dynamic_cast<Mesh*>(SMesh[meshName])->Get_M2(),
						dp_x, dp_y, 
						num_bins, min, max, M2, deltaM2);

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

		case CMD_DP_POW:
		{
			int dp_source, dp_dest;
			double exponent;

			error = commandSpec.GetParameters(command_fields, dp_source, exponent, dp_dest);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_source, exponent); dp_dest = dp_source; }

			if (!error) {

				error = dpArr.exponentiate(dp_source, dp_dest, exponent);
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
			double exclusion_ratio = 0.0;

			error = commandSpec.GetParameters(command_fields, dp_source, exclusion_ratio);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, dp_source); exclusion_ratio = 0.0; }

			if (!error) {

				DBL2 mean_stdev;

				error  = dpArr.get_mean(dp_source, &mean_stdev, exclusion_ratio);

				if (verbose && !error) BD.DisplayConsoleMessage("Mean = " + ToString(mean_stdev.x) + " +/- " + ToString(mean_stdev.y));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(mean_stdev));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_SUM:
		{
			int dp_source;

			error = commandSpec.GetParameters(command_fields, dp_source);

			if (!error) {

				double sum = 0.0;
				error = dpArr.get_sum(dp_source, &sum);

				if (verbose && !error) BD.DisplayConsoleMessage("Sum = " + ToString(sum));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(sum));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_CHUNKEDSTD:
		{
			int dp_source, chunk;

			error = commandSpec.GetParameters(command_fields, dp_source, chunk);

			if (!error) {

				double chunked_stdev;

				error = dpArr.get_chunkedstd(dp_source, chunk, &chunked_stdev);

				if (verbose && !error) BD.DisplayConsoleMessage("Chunked std = " + ToString(chunked_stdev));

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(chunked_stdev));
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
			std::string meshName, paramName;
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

		case CMD_DP_FITLORENTZ2:
		{
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				DBL2 S, A, H0, dH, y0;

				error = dpArr.fit_lorentz2(dp_x, dp_y, &S, &A, &H0, &dH, &y0);

				if (verbose && !error) {

					BD.DisplayConsoleMessage(
						"S = " + ToString(S.major) + " +/- " + ToString(S.minor) + ", " +
						"A = " + ToString(A.major) + " +/- " + ToString(A.minor) + ", " +
						"H0 = " + ToString(H0.major) + " +/- " + ToString(H0.minor) + ", " +
						"dH = " + ToString(dH.major) + " +/- " + ToString(dH.minor) + ", " +
						"y0 = " + ToString(y0.major) + " +/- " + ToString(y0.minor)
					);
				}

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(S.major, A.major, H0.major, dH.major, y0.major, S.minor, A.minor, H0.minor, dH.minor, y0.minor));
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

		case CMD_DP_FITDW:
		{
			int dp_x, dp_y;

			error = commandSpec.GetParameters(command_fields, dp_x, dp_y);

			if (!error) {

				DBL2 A, x0, D;

				error = dpArr.fit_domainwall(dp_x, dp_y, &A, &x0, &D);

				if (verbose && !error) {

					BD.DisplayConsoleMessage(
						"D = " + ToString(D.major) + " +/- " + ToString(D.minor) + ", " +
						"x0 = " + ToString(x0.major) + " +/- " + ToString(x0.minor) + ", " +
						"A = " + ToString(A.major) + " +/- " + ToString(A.minor)
					);
				}

				if (script_client_connected) commSocket.SetSendData(commandSpec.PrepareReturnParameters(D.major, x0.major, A.major, D.minor, x0.minor, A.minor));
			}
			else if (verbose) PrintCommandUsage(command_name);
		}
		break;

		case CMD_DP_FITSTT:
		{
			Rect rectangle;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset(); rectangle = SMesh[meshName]->GetMeshRect() - SMesh[meshName]->GetOrigin(); }

			if (!error) {

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
			double abs_error_threshold, Rsq_threshold, T_ratio_threshold;
			int stencil = 3;
			bool user_thresholds = true;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, abs_error_threshold, Rsq_threshold, T_ratio_threshold, stencil);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, abs_error_threshold, Rsq_threshold, T_ratio_threshold);  stencil = 3; }
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
			std::string hm_mesh;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, hm_mesh, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, hm_mesh); rectangle = SMesh[meshName]->GetMeshRect() - SMesh[meshName]->GetOrigin(); }

			if (!error) {

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
			std::string hm_mesh;
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			error = commandSpec.GetParameters(command_fields, meshName, hm_mesh, rectangle);
			if (error == BERROR_PARAMMISMATCH) { error.reset() = commandSpec.GetParameters(command_fields, meshName, hm_mesh); rectangle = SMesh[meshName]->GetMeshRect() - SMesh[meshName]->GetOrigin(); }

			if (!error) {

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
			std::string hm_mesh, fm_mesh;

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

		case CMD_DP_CALCTOPOCHARGEDENSITY:
		{
			std::string meshName;

			optional_meshname_check_focusedmeshdefault(command_fields);
			meshName = command_fields[0];

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
			/*
			//DUMP ALL COMMANDS
			std::vector<std::string> commands_output;
			commands_output.resize(commands.size());

			for (int idx = 0; idx < commands.size(); idx++) {

				commands_output[idx] = commands.get_key_from_index(idx);
			}

			std::sort(commands_output.begin(), commands_output.end());

			std::string commands_description;

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

			SaveTextToFile("C:/Boris/PythonScripting/commands.txt", commands_description);

			//DUMP ALL PARAMETER NAMES FOR EACH MESH TYPE
			auto dump_params = [&](MESH_ meshType, std::string meshtypename) -> void {

				std::string meshname_dump = "dump_params";
				std::string params_text = "This goes in class Parameters_" + meshtypename + "\n\n";

				SMesh.AddMesh(meshname_dump, meshType, Rect(DBL3(5e-9, 5e-9, 5e-9)));

				std::vector<PARAM_>& params = params_for_meshtype(meshType);
				for (PARAM_ param : params) {

					std::string param_handle = SMesh[meshname_dump]->get_meshparam_handle(param);
					params_text += param_handle + " = ''\n";
				}

				params_text += "\n\nThis goes in class Parameters_" + meshtypename + " constructor\n\n";

				for (PARAM_ param : params) {

					std::string param_handle = SMesh[meshname_dump]->get_meshparam_handle(param);
					params_text += "self." + param_handle + " = ns.Param(ns, '" + param_handle + "', meshname)\n";
				}

				SMesh.DelMesh(meshname_dump);
				SaveTextToFile("C:/Boris/PythonScripting/" + meshtypename + "_params.txt", params_text);
			};
			
			dump_params(MESH_FERROMAGNETIC, "Ferromagnet");
			dump_params(MESH_ANTIFERROMAGNETIC, "AntiFerromagnet");
			dump_params(MESH_DIPOLE, "Dipole");
			dump_params(MESH_METAL, "Metal");
			dump_params(MESH_INSULATOR, "Insulator");
			dump_params(MESH_ATOM_CUBIC, "Atomistic");

			//DUMP ALL QUANTITY NAMES FOR EACH MESH TYPE
			auto dump_quantities = [&](MESH_ meshType, std::string meshtypename) -> void {

				std::string quantities_text = "This goes in class Quantities_" + meshtypename + "\n\n";
				std::vector<MESHDISPLAY_>& quantities = meshAllowedDisplay(meshType);
				for (MESHDISPLAY_ quantity : quantities) {
					
					std::string quantity_handle = displayHandles(quantity);
					quantities_text += quantity_handle + " = ''\n";
				}

				quantities_text += "\n\nThis goes in class Quantities_" + meshtypename + " constructor\n\n";

				for (MESHDISPLAY_ quantity : quantities) {

					std::string quantity_handle = displayHandles(quantity);
					quantities_text += "self." + quantity_handle + " = ns.Quantity(ns, '" + quantity_handle + "', meshname)\n";
				}

				SaveTextToFile("C:/Boris/PythonScripting/" + meshtypename + "_quantities.txt", quantities_text);
			};

			dump_quantities(MESH_FERROMAGNETIC, meshtypeHandles(MESH_FERROMAGNETIC));
			dump_quantities(MESH_ANTIFERROMAGNETIC, meshtypeHandles(MESH_ANTIFERROMAGNETIC));
			dump_quantities(MESH_DIPOLE, meshtypeHandles(MESH_DIPOLE));
			dump_quantities(MESH_METAL, meshtypeHandles(MESH_METAL));
			dump_quantities(MESH_INSULATOR, meshtypeHandles(MESH_INSULATOR));
			dump_quantities(MESH_ATOM_CUBIC, meshtypeHandles(MESH_ATOM_CUBIC));
			*/
		}
		break;

		default:
			break;
		}

		if (error || error.warning_set()) {

			err_hndl.show_error(error, verbose);

			//show error in Python console through return parameter, but not if the error is due to a parameter mismatch as this is also caused if the command is issued without all parameter in order to retrieve data
			if (script_client_connected && error != BERROR_PARAMMISMATCH && !err_hndl.is_silent_error(error)) commSocket.SetSendData({err_hndl.get_error_text(error)});
		}
	}
	else err_hndl.show_error(BERROR_COMMAND_NOTRECOGNIZED, verbose);
}