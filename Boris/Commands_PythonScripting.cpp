#include "stdafx.h"
#include "Simulation.h"

#if PYTHON_EMBEDDING == 1

#include "PythonScripting.h"

//Receive and execute a message from embedded Python script, and return response.
std::string Simulation::NewPythonMessage(std::string message)
{
	if (!is_thread_running(THREAD_LOOP)) {

		set_blocking_thread(THREAD_HANDLEMESSAGE);
		while (single_call_launch<std::string>(&Simulation::HandleCommand, message, THREAD_HANDLEMESSAGE) != THREAD_HANDLEMESSAGE);
		return commSocket.GetDataParams_String();
	}
	else {

		HandleCommand(message);
		return commSocket.GetDataParams_String();
	}
}

//Received run command from embedded Python script - special handling required.
void Simulation::NewPythonMessage_Run(void)
{
	//This is similar to RunSimulation() but, instead keep THREAD_LOOP busy with a dummy routine, and run the real thing in this thread
	PrepareRunSimulation();
	infinite_loop_launch(&Simulation::Simulate_Dummy, &Simulation::SetupRunSimulation, THREAD_LOOP);

	//infinite loop on main thread : blocks embedded Python code execution until done (also set required number of Omp threads on this thread)
	SetupRunSimulation();
	while (is_thread_running(THREAD_LOOP)) { Simulate(); }

	//THREAD_LOOP_STOP was used to stop THREAD_LOOP, so make sure it's done before returning.
	while (is_thread_running(THREAD_LOOP_STOP)) {}
}

void Simulation::NewPythonMessage_RunStage(int stage)
{
	stage_step = INT2(stage, 0);
	SetSimulationStageValue();
	single_stage_run = true;

	//This is similar to RunSimulation() but, instead keep THREAD_LOOP busy with a dummy routine, and run the real thing in this thread
	PrepareRunSimulation();
	infinite_loop_launch(&Simulation::Simulate_Dummy, &Simulation::SetupRunSimulation, THREAD_LOOP);

	//infinite loop on main thread : blocks embedded Python code execution until done (also set required number of Omp threads on this thread)
	SetupRunSimulation();
	while (is_thread_running(THREAD_LOOP)) { Simulate(); }

	//THREAD_LOOP_STOP was used to stop THREAD_LOOP, so make sure it's done before returning.
	while (is_thread_running(THREAD_LOOP_STOP)) {}
}

//Received run command from embedded Python script, which also required further Python code to be executed every simulation iteration - special handling required.
//return 1 if simulation is still running, 0 if not
int Simulation::NewPythonMessage_RunWithCode(void)
{
	//Here do not call Simulate on an infinite loop, but just once - the loop is moved to the embedded Python code! (very little overhead incurred so a good solution)
	//This allows further Python code to be executed after each iteration in the same scope.
	Simulate();

	//is simulation still required to run? 
	//If THREAD_LOOP_STOP was used to stop THREAD_LOOP, then make sure it's done before returning state
	while (is_thread_running(THREAD_LOOP_STOP)) {}
	//now return true state of THREAD_LOOP
	return is_thread_running(THREAD_LOOP);
}

//before we can use NewPythonMessage_RunWithCode (in a loop), must call NewPythonMessage_PrepareRunWithCode
void Simulation::NewPythonMessage_PrepareRunWithCode(void)
{
	//This is similar to RunSimulation() but, instead keep THREAD_LOOP busy with a dummy routine, and run the real thing in this thread
	PrepareRunSimulation();
	infinite_loop_launch(&Simulation::Simulate_Dummy, &Simulation::SetupRunSimulation, THREAD_LOOP);
	
	//Also set required number of Omp threads on this thread
	SetupRunSimulation();
}

//This method loads and runs the Python script - when called from HandleCommand must be launched on new thread to avoid mutex blocking
void Simulation::RunPythonScript(std::string fileName)
{
	//1. Scan Python file into a vector of strings, which we'll execute using PyRun_SimpleString after modifying it as needed
	std::vector<std::string> PythonScript_lines;

	//If using pragma omp parallel for, additional scripts will be generated here to be executed in parallel
	std::vector<std::vector<std::string>> mGPU_PythonScript_lines;

	std::ifstream bdin;
	bdin.open(fileName.c_str(), std::ios::in);

	if (bdin.is_open()) {

		char line[FILEROWCHARS];
		while (bdin.getline(line, FILEROWCHARS)) {

			PythonScript_lines.push_back(std::string(line));
		}

		bdin.close();
	}
	else BD.DisplayConsoleError("Could not load Python script (" + fileName + ").");
	
	if (!PythonScript_lines.size()) return;
	
	//2. Check for #pragma parallel for directive if we have more than 1 GPU available, or stride is greater than 1 (python_script_parallel)
	
	if (python_script_mGPU) {

		int stride = 1;
		std::vector<int> offsets;
		if (python_script_parallel.size() >= 2) {

			stride = python_script_parallel[0];
			offsets = subvec(python_script_parallel, 1);
		}

		if (cudaDeviceVersions.size() > 1 || stride > 1) {

			if (stride == 1) {

				//Stride and offsets not specified, but we have multiple GPUs : set stride to number of GPUs, and include all offsets
				stride = cudaDeviceVersions.size();
				for (int idx = 0; idx < stride; idx++) offsets.push_back(idx);
			}

			//now scan for directive and generate scripts using stride and offsets information (one for each offset)
			for (int idx = 0; idx < PythonScript_lines.size(); idx++) {

				//look for pragma directive
				size_t pos = PythonScript_lines[idx].find("#pragma parallel for");
				if (pos != std::string::npos && idx + 1 < PythonScript_lines.size()) {

					std::string indent_pragma = std::string(PythonScript_lines[idx]).substr(0, pos);

					//next line must be a for loop
					if (idx + 2 < PythonScript_lines.size()) {

						pos = PythonScript_lines[idx + 1].find("for");
						if (pos != std::string::npos) {

							pos = PythonScript_lines[idx + 2].find_first_not_of(" ");
							std::string indent_for = std::string(PythonScript_lines[idx + 2]).substr(0, pos);

							//since we do have a pragma directive specified, then make space for new scripts
							if (!mGPU_PythonScript_lines.size()) {

								for (int offset_idx = 1; offset_idx < offsets.size(); offset_idx++) mGPU_PythonScript_lines.push_back(PythonScript_lines);
							}

							//1
							PythonScript_lines[idx] = indent_pragma + "pragma_mGPU_parallel_for_index = -1";
							for (int mGPU_idx = 0; mGPU_idx < mGPU_PythonScript_lines.size(); mGPU_idx++)
								mGPU_PythonScript_lines[mGPU_idx][idx] = PythonScript_lines[idx];

							//2
							PythonScript_lines.insert(PythonScript_lines.begin() + idx + 2, indent_for + "pragma_mGPU_parallel_for_index += 1");
							for (int mGPU_idx = 0; mGPU_idx < mGPU_PythonScript_lines.size(); mGPU_idx++)
								mGPU_PythonScript_lines[mGPU_idx].insert(mGPU_PythonScript_lines[mGPU_idx].begin() + idx + 2, indent_for + "pragma_mGPU_parallel_for_index += 1");

							//3
							PythonScript_lines.insert(PythonScript_lines.begin() + idx + 3, indent_for + "if (pragma_mGPU_parallel_for_index + " + ToString(offsets[0]) + ") % " + ToString(stride) + " != 0: continue");
							for (int mGPU_idx = 0; mGPU_idx < mGPU_PythonScript_lines.size(); mGPU_idx++)
								mGPU_PythonScript_lines[mGPU_idx].insert(mGPU_PythonScript_lines[mGPU_idx].begin() + idx + 3, indent_for + "if (pragma_mGPU_parallel_for_index + " + ToString(offsets[mGPU_idx + 1]) + ") % " + ToString(stride) + " != 0: continue");
						}
					}
				}
			}
		}

		if (mGPU_PythonScript_lines.size()) {

			//for loop split on this machine, so try to assign gpus in sequential order (there might be no GPUs on this machine - in which case user shouldn't be assigning more than 1 offset here; in this case all of them will be executed on CPU anyway)
			//start from next one up after cudaDeviceSelect
			int next_gpu = cudaDeviceSelect + 1;
			if (cudaDeviceVersions.size() >= 1) next_gpu = next_gpu % cudaDeviceVersions.size();
			else next_gpu = 0;

			for (int mGPU_idx = 0; mGPU_idx < mGPU_PythonScript_lines.size(); mGPU_idx++) {

				std::string fileName_mGPU = fileName + "__" + ToString(stride) + "_" + ToString(offsets[mGPU_idx + 1]) + ".py";

				std::ofstream bdout;
				bdout.open(ScanFileNameData(fileName_mGPU), std::ios::out);
				for (int idx = 0; idx < PythonScript_lines.size(); idx++) bdout << mGPU_PythonScript_lines[mGPU_idx][idx] << std::endl;
				bdout.close();

				int script_gpu = offsets[mGPU_idx + 1];

				open_file(GetExeDirectory() + progName, "-g " + ToString(next_gpu) + " -s " + fileName_mGPU + " -st 1 -sm 0 -sd 1");
				if (cudaDeviceVersions.size() >= 1) next_gpu = (next_gpu + 1) % cudaDeviceVersions.size();
			}
		}
	}

	//3. Add wrapper and polling function so we can stop the script before completion

	//this will be the string to execute, after we modify it as needed
	std::string WrappedPythonScript_as_string;

	std::string embeddedPythonWrapper_fileNamewithPath = GetExeDirectory() + embeddedPythonWraper_fileName;
	bdin.open(embeddedPythonWrapper_fileNamewithPath.c_str(), std::ios::in);
	if (bdin.is_open()) {

		char line[FILEROWCHARS];
		while (bdin.getline(line, FILEROWCHARS)) {

			WrappedPythonScript_as_string += std::string(line) + std::string("\n");

			size_t pos = std::string(line).find("#user script goes here in its entirety - do not modify this line");

			if (pos != std::string::npos) {

				WrappedPythonScript_as_string += std::string("\n");

				std::string indent = std::string(line).substr(0, pos);

				for (std::string scriptline : PythonScript_lines) {

					WrappedPythonScript_as_string += indent + scriptline + std::string("\n");
				}
			}
		}

		bdin.close();
	}
	else BD.DisplayConsoleError("Could not load embedded Python wrapper file : " + embeddedPythonWrapper_fileNamewithPath);

	//4. Make sure current script directory is in the Python path (it's not when we use embedding)
	PyObject* sysPath = PySys_GetObject((char*)"path");
	PyList_Append(sysPath, (PyUnicode_FromString(fileName.c_str())));

	//set working directory same as script file directory
	directory = GetFilenameDirectory(fileName);

	//5. Delete source script?
	if (python_script_deletesource) {

#if OPERATING_SYSTEM == OS_WIN
		DeleteFileA(fileName.c_str());
#elif OPERATING_SYSTEM == OS_LIN
		remove(fileName.c_str());
#endif
	}

	//6. Scan script for any Mesh objects, and insert mesh name - capturing object name using traceback doesn't work with PyRun_SimpleString
	
	auto find_end_bracket = [&](size_t pos_start, std::string& text) -> size_t {

		//pos_start must be the index after an opening bracket, (
		size_t pos_end = pos_start;
		int num_open_brackets = 1;
		while (pos_end < text.length() && num_open_brackets > 0) {

			if (text[pos_end] == '(') num_open_brackets++;
			if (text[pos_end] == ')') num_open_brackets--;
			pos_end++;
		}

		//return index after closing bracket, or overflow text length if none could be found
		return pos_end;
	};

	auto find_meshname = [&](size_t pos_start, std::string& text) -> std::string {

		//we are looking for something of the form meshname = ns.Ferromagnet(...
		//pos_start must be the index before the object to the right of '=' - get the object name on the left of '='

		auto is_var_char = [](char c) -> bool { return isalnum(c) || c == '_'; };

		size_t pos_left = 0, pos_right = 0;
		bool assignment_found = false;
		while (pos_start > 0) {

			if (assignment_found) {

				if (is_var_char(text[pos_start])) {

					if (pos_right == 0) pos_right = pos_start;
				}
				else {

					if (!(text[pos_start] == ' ' && pos_right == 0)) {

						pos_left = pos_start + 1;
						break;
					}
				}
			}
			else if (text[pos_start] == '=') assignment_found = true;

			pos_start--;
		}

		return text.substr(pos_left, pos_right + 1 - pos_left);
	};

	for(int mesh_idx = 0; mesh_idx < meshtypeHandles.size() + 1; mesh_idx++) {

		std::string meshType;
		if (mesh_idx < meshtypeHandles.size()) meshType = "." + meshtypeHandles[mesh_idx] + "(";
		else meshType = ".Material(";
	
		size_t pos = 0;
		while (true) {

			pos = WrappedPythonScript_as_string.find(meshType, pos);

			if (pos != std::string::npos) {

				std::string meshname = find_meshname(pos - 1, WrappedPythonScript_as_string);
				size_t pos_end = find_end_bracket(pos + meshType.length() + 1, WrappedPythonScript_as_string);
				if (pos_end < WrappedPythonScript_as_string.length()) {

					//is meshname already given? If not then insert it, otherwise leave it alone
					int num_commas = 0;
					int num_open_brackets = 0;
					for (int idx = pos; idx < pos_end; idx++) {
						
						if (WrappedPythonScript_as_string[idx] == ',' && !num_open_brackets) num_commas++;
						if (WrappedPythonScript_as_string[idx] == '[') num_open_brackets++;
						if (WrappedPythonScript_as_string[idx] == ']') num_open_brackets--;
					}

					//if there are fewer than 3 parameters passed (2 commas) then meshname was not given (or not commas for Dipole mesh)
					if (mesh_idx < meshtypeHandles.size()) {
						if (meshtypeHandles.get_ID_from_index(mesh_idx) == MESH_DIPOLE) {
							if (!num_commas) WrappedPythonScript_as_string.insert(pos_end - 1, ", meshname = '" + meshname + "'");
						}
						else if (num_commas < 2)
							WrappedPythonScript_as_string.insert(pos_end - 1, ", meshname = '" + meshname + "'");
					}
					else if (num_commas < 3)
						WrappedPythonScript_as_string.insert(pos_end - 1, ", meshname = '" + meshname + "'");
				}

				pos = pos_end;
			}
			else break;
		}
	}

	//7. Finally execute the processed Python script
	BD.DisplayConsoleMessage("Started Python script : " + fileName);
	python_script_running = true;
	PyRun_SimpleString(WrappedPythonScript_as_string.c_str());
	python_script_running = false;
	BD.DisplayConsoleMessage("Finished Python script : " + fileName);

	//8. Do we have to terminate program after script execution? Only use this if script launched from command line
	if (python_script_terminate && python_script.length()) {

#if OPERATING_SYSTEM == OS_WIN && GRAPHICS == 1
		PostMessage(hWnd, WM_CLOSE, 0, GetCurrentTime());
#elif OPERATING_SYSTEM == OS_LIN
		exit(3);
#endif
	}
}

#endif