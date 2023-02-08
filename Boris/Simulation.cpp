#include "stdafx.h"
#include "Simulation.h"

#if GRAPHICS == 1

Simulation::Simulation(
	HWND hWnd, 
	int Program_Version, std::string progName_,
	std::string server_port_, std::string server_pwd_,
	int cudaDevice,
	std::string python_script_, int python_script_mGPU_, std::vector<int> python_script_parallel_, int python_script_terminate_, int python_script_deletesource_) :
	err_hndl(this),
	BD(hWnd, new SimTOFunct(this, &Simulation::ConsoleActionHandler, &Simulation::ConsoleInteractiveObjectState)),
	SimulationSharedData(true),
	mdb(params_for_meshtype(MESH_SUPERMESH)),
	ProgramStateNames(this,
		{
			VINFO(BD),
			VINFO(directory), VINFO(savedataFile), VINFO(imageSaveFileBase), VINFO(currentSimulationFile), VINFO(appendToDataFile), VINFO(saveDataFlag), VINFO(saveImageFlag), VINFO(dataprecision), VINFO(savedata_diskbuffer_size),
			VINFO(saveDataList), VINFO(dataBoxList),
			VINFO(stage_step),
			VINFO(simStages), VINFO(iterUpdate), VINFO(autocomplete),
			VINFO(SMesh),
			VINFO(cudaEnabled), VINFO(cudaDeviceSelect),
			VINFO(shape_change_individual),
			VINFO(static_transport_solver), VINFO(disabled_transport_solver),
			VINFO(image_cropping), VINFO(displayTransparency), VINFO(displayThresholds), VINFO(displayThresholdTrigger),
			VINFO(shape_rotation), VINFO(shape_repetitions), VINFO(shape_displacement), VINFO(shape_method),
			VINFO(userConstants),
			VINFO(command_buffer)
		}, {})

#else

Simulation::Simulation(
	int Program_Version, std::string progName_,
	std::string server_port_, std::string server_pwd_,
	int cudaDevice,
	std::string python_script_, int python_script_mGPU_, std::vector<int> python_script_parallel_, int python_script_terminate_, int python_script_deletesource_) :
	err_hndl(this),
	BD(),
	SimulationSharedData(true),
	mdb(params_for_meshtype(MESH_SUPERMESH)),
	ProgramStateNames(this,
		{
			VINFO(BD),
			VINFO(directory), VINFO(savedataFile), VINFO(imageSaveFileBase), VINFO(currentSimulationFile), VINFO(appendToDataFile), VINFO(saveDataFlag), VINFO(saveImageFlag), VINFO(dataprecision), VINFO(savedata_diskbuffer_size),
			VINFO(saveDataList), VINFO(dataBoxList),
			VINFO(stage_step),
			VINFO(simStages), VINFO(iterUpdate), VINFO(autocomplete),
			VINFO(SMesh),
			VINFO(cudaEnabled), VINFO(cudaDeviceSelect),
			VINFO(shape_change_individual),
			VINFO(static_transport_solver), VINFO(disabled_transport_solver),
			VINFO(image_cropping), VINFO(displayTransparency), VINFO(displayThresholds), VINFO(displayThresholdTrigger),
			VINFO(shape_rotation), VINFO(shape_repetitions), VINFO(shape_displacement), VINFO(shape_method),
			VINFO(userConstants),
			VINFO(command_buffer)
		}, {})

#endif

{
#if GRAPHICS == 1 && OPERATING_SYSTEM == OS_WIN
	//Save window handle if in graphics mode on windows
	this->hWnd = hWnd;
#endif

	progName = progName_;
	ExtractFilenameDirectory(progName);

	MakeIOInfo();

	//---------------------------------------------------------------- SETTINGS

	this->Program_Version = Program_Version;

	//default directory
	directory = GetUserDocumentsPath() + boris_data_directory + boris_simulations_directory;

	//---------------------------------------------------------------- CUDA

#if COMPILECUDA == 1

	int deviceCount;
	cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);

	if (cudaResultCode != cudaSuccess) cudaAvailable = false;
	else cudaAvailable = true;

	if (cudaAvailable) {

		for (int idx = 0; idx < deviceCount; idx++)
		{
			cudaDeviceProp deviceProperties;
			cudaGetDeviceProperties(&deviceProperties, idx);

			std::string device_info = deviceProperties.name + std::string("; Compute: ") + ToString(deviceProperties.major);

			cudaDeviceVersions.push_back(std::pair<int, std::string>(deviceProperties.major * 100, device_info));
		}

		//with the default parameter cudaDevice = -1 we need to determine a value for cudaDeviceSelect
		if (cudaDevice < 0) {

			//Select first device with matches current CUDA version
			for (int idx = 0; idx < deviceCount; idx++) {

				if (cudaDeviceVersions[idx].first == __CUDA_ARCH__) {

					cudaDeviceSelect = idx;
					break;
				}
			}
		}
		else if (cudaDevice < deviceCount) {

			if (cudaDeviceVersions[cudaDevice].first == __CUDA_ARCH__) cudaDeviceSelect = cudaDevice;
		}

		cudaSetDevice(cudaDeviceSelect);
		cudaDeviceReset();
	}
#endif

	//---------------------------------------------------------------- COMMANDS

	commands.insert(CMD_RUN, CommandSpecifier(CMD_RUN), "run");
	commands[CMD_RUN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>run</b>";
	commands[CMD_RUN].descr = "[tc0,0.5,0.5,1/tc]Run simulation from current state.";

	commands.insert(CMD_STOP, CommandSpecifier(CMD_STOP), "stop");
	commands[CMD_STOP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>stop</b>";
	commands[CMD_STOP].descr = "[tc0,0.5,0.5,1/tc]Stop simulation without resetting it.";

	commands.insert(CMD_RESET, CommandSpecifier(CMD_RESET), "reset");
	commands[CMD_RESET].usage = "[tc0,0.5,0,1/tc]USAGE : <b>reset</b>";
	commands[CMD_RESET].descr = "[tc0,0.5,0.5,1/tc]Reset simulation state to the starting state.";

	commands.insert(CMD_RUNSTAGE, CommandSpecifier(CMD_RUNSTAGE), "runstage");
	commands[CMD_RUNSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>runstage</b> <i>stage</i>";
	commands[CMD_RUNSTAGE].limits = { { int(0), Any() } };
	commands[CMD_RUNSTAGE].descr = "[tc0,0.5,0.5,1/tc]Run given simulation stage only.";

	commands.insert(CMD_SETSCHEDULESTAGE, CommandSpecifier(CMD_SETSCHEDULESTAGE), "setschedulestage");
	commands[CMD_SETSCHEDULESTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setschedulestage</b> <i>stage</i>";
	commands[CMD_SETSCHEDULESTAGE].limits = { { int(0), Any() } };
	commands[CMD_SETSCHEDULESTAGE].descr = "[tc0,0.5,0.5,1/tc]Set stage value, but do not run simulation. Must be a valid stage number.";
	commands[CMD_SETSCHEDULESTAGE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>stage</i>";
	
	commands.insert(CMD_NEXTSTAGE, CommandSpecifier(CMD_NEXTSTAGE), "nextstage");
	commands[CMD_NEXTSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>nextstage</b>";
	commands[CMD_NEXTSTAGE].descr = "[tc0,0.5,0.5,1/tc]Increment stage value by one.";

	commands.insert(CMD_COMPUTEFIELDS, CommandSpecifier(CMD_COMPUTEFIELDS), "computefields");
	commands[CMD_COMPUTEFIELDS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>computefields</b>";
	commands[CMD_COMPUTEFIELDS].descr = "[tc0,0.5,0.5,1/tc]Run simulation from current state for a single iteration without advancing the simulation time.";

	commands.insert(CMD_ISRUNNING, CommandSpecifier(CMD_ISRUNNING), "isrunning");
	commands[CMD_ISRUNNING].usage = "[tc0,0.5,0,1/tc]USAGE : <b>isrunning</b>";
	commands[CMD_ISRUNNING].descr = "[tc0,0.5,0.5,1/tc]Checks if the simulation is running and sends state value to the calling script.";
	commands[CMD_ISRUNNING].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>state</i> - return simulation running state.";

	commands.insert(CMD_RUNSCRIPT, CommandSpecifier(CMD_RUNSCRIPT), "runscript");
	commands[CMD_RUNSCRIPT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>runscript</b> <i>(directory/)filename</i>";
	commands[CMD_RUNSCRIPT].descr = "[tc0,0.5,0.5,1/tc]Execute a Python-scripted simulation (filename must be a Python script). If directory not specified then default directory is used.";

	commands.insert(CMD_STOPSCRIPT, CommandSpecifier(CMD_STOPSCRIPT), "stopscript");
	commands[CMD_STOPSCRIPT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>stopscript</b>";
	commands[CMD_STOPSCRIPT].descr = "[tc0,0.5,0.5,1/tc]Stop any currently running Python script.";
	
	commands.insert(CMD_SCRIPTPRINT, CommandSpecifier(CMD_SCRIPTPRINT), "bprint");
	commands[CMD_SCRIPTPRINT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>bprint</b> <i>message</i>";
	commands[CMD_SCRIPTPRINT].descr = "[tc0,0.5,0.5,1/tc]Print message in console.";

	commands.insert(CMD_CENTER, CommandSpecifier(CMD_CENTER), "center");
	commands[CMD_CENTER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>center</b>";
	commands[CMD_CENTER].descr = "[tc0,0.5,0.5,1/tc]Center mesh view and scale to fit window size.";

	commands.insert(CMD_ITERUPDATE, CommandSpecifier(CMD_ITERUPDATE), "iterupdate");
	commands[CMD_ITERUPDATE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>iterupdate</b> <i>iterations</i>";
	commands[CMD_ITERUPDATE].limits = { { int(0), Any() } };
	commands[CMD_ITERUPDATE].descr = "[tc0,0.5,0.5,1/tc]Update mesh display every given number of iterations during a simulation.";
	commands[CMD_ITERUPDATE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>iterations</i> - return number of iterations for display update.";

	commands.insert(CMD_CLEARSCREEN, CommandSpecifier(CMD_CLEARSCREEN), "clearscreen");
	commands[CMD_CLEARSCREEN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearscreen</b>";
	commands[CMD_CLEARSCREEN].descr = "[tc0,0.5,0.5,1/tc]Clear all console text.";

	commands.insert(CMD_REFRESHSCREEN, CommandSpecifier(CMD_REFRESHSCREEN), "refreshscreen");
	commands[CMD_REFRESHSCREEN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>refreshscreen</b>";
	commands[CMD_REFRESHSCREEN].descr = "[tc0,0.5,0.5,1/tc]Refreshes entire screen.";

	commands.insert(CMD_UPDATESCREEN, CommandSpecifier(CMD_UPDATESCREEN), "updatescreen");
	commands[CMD_UPDATESCREEN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>updatescreen</b>";
	commands[CMD_UPDATESCREEN].descr = "[tc0,0.5,0.5,1/tc]Updates all displayed values on screen and also refreshes.";

	commands.insert(CMD_ADDFMESH, CommandSpecifier(CMD_ADDFMESH), "addmesh");
	commands[CMD_ADDFMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addmesh</b> <i>name rectangle</i>";
	commands[CMD_ADDFMESH].limits = { { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDFMESH].descr = "[tc0,0.5,0.5,1/tc]Add a ferromagnetic mesh with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDFMESH].unit = "m";

	commands.insert(CMD_SETFMESH, CommandSpecifier(CMD_SETFMESH), "setmesh");
	commands[CMD_SETFMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setmesh</b> <i>name rectangle</i>";
	commands[CMD_SETFMESH].limits = { { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETFMESH].descr = "[tc0,0.5,0.5,1/tc]Set a single ferromagnetic mesh (deleting all other meshes) with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETFMESH].unit = "m";

	commands.insert(CMD_ADDAFMESH, CommandSpecifier(CMD_ADDAFMESH), "addafmesh");
	commands[CMD_ADDAFMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addafmesh</b> <i>name rectangle</i>";
	commands[CMD_ADDAFMESH].limits = { { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDAFMESH].descr = "[tc0,0.5,0.5,1/tc]Add antiferromagnetic mesh with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDAFMESH].unit = "m";

	commands.insert(CMD_SETAFMESH, CommandSpecifier(CMD_SETAFMESH), "setafmesh");
	commands[CMD_SETAFMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setafmesh</b> <i>name rectangle</i>";
	commands[CMD_SETAFMESH].limits = { { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETAFMESH].descr = "[tc0,0.5,0.5,1/tc]Set a single antiferromagnetic mesh (deleting all other meshes) with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETAFMESH].unit = "m";

	commands.insert(CMD_ADDMETALMESH, CommandSpecifier(CMD_ADDMETALMESH), "addconductor");
	commands[CMD_ADDMETALMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addconductor</b> <i>name rectangle</i>";
	commands[CMD_ADDMETALMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDMETALMESH].descr = "[tc0,0.5,0.5,1/tc]Add a normal metal mesh with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDMETALMESH].unit = "m";

	commands.insert(CMD_SETMETALMESH, CommandSpecifier(CMD_SETMETALMESH), "setconductor");
	commands[CMD_SETMETALMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setconductor</b> <i>name rectangle</i>";
	commands[CMD_SETMETALMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETMETALMESH].descr = "[tc0,0.5,0.5,1/tc]Set a single metal mesh (deleting all other meshes) with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETMETALMESH].unit = "m";

	commands.insert(CMD_ADDDIPOLEMESH, CommandSpecifier(CMD_ADDDIPOLEMESH), "adddipole");
	commands[CMD_ADDDIPOLEMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDDIPOLEMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>adddipole</b> <i>name rectangle</i>";
	commands[CMD_ADDDIPOLEMESH].descr = "[tc0,0.5,0.5,1/tc]Add a rectangular dipole with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDDIPOLEMESH].unit = "m";

	commands.insert(CMD_SETDIPOLEMESH, CommandSpecifier(CMD_SETDIPOLEMESH), "setdipole");
	commands[CMD_SETDIPOLEMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETDIPOLEMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdipole</b> <i>name rectangle</i>";
	commands[CMD_SETDIPOLEMESH].descr = "[tc0,0.5,0.5,1/tc]Set a single dipole (deleting all other meshes) with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETDIPOLEMESH].unit = "m";

	commands.insert(CMD_ADDINSULATORMESH, CommandSpecifier(CMD_ADDINSULATORMESH), "addinsulator");
	commands[CMD_ADDINSULATORMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addinsulator</b> <i>name rectangle</i>";
	commands[CMD_ADDINSULATORMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDINSULATORMESH].descr = "[tc0,0.5,0.5,1/tc]Add an insulator mesh with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDINSULATORMESH].unit = "m";

	commands.insert(CMD_SETINSULATORMESH, CommandSpecifier(CMD_SETINSULATORMESH), "setinsulator");
	commands[CMD_SETINSULATORMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setinsulator</b> <i>name rectangle</i>";
	commands[CMD_SETINSULATORMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETINSULATORMESH].descr = "[tc0,0.5,0.5,1/tc]Set a single insulator mesh (deleting all other meshes) with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETINSULATORMESH].unit = "m";

	commands.insert(CMD_ADDAMESHCUBIC, CommandSpecifier(CMD_ADDAMESHCUBIC), "addameshcubic");
	commands[CMD_ADDAMESHCUBIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addameshcubic</b> <i>name rectangle</i>";
	commands[CMD_ADDAMESHCUBIC].limits = { { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDAMESHCUBIC].descr = "[tc0,0.5,0.5,1/tc]Add an atomistic mesh with simple cubic structure, with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDAMESHCUBIC].unit = "m";

	commands.insert(CMD_SETAMESHCUBIC, CommandSpecifier(CMD_SETAMESHCUBIC), "setameshcubic");
	commands[CMD_SETAMESHCUBIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setameshcubic</b> <i>name rectangle</i>";
	commands[CMD_SETAMESHCUBIC].limits = { { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETAMESHCUBIC].descr = "[tc0,0.5,0.5,1/tc]Set a single atomistic mesh (deleting all other meshes) with simple cubic structure, with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETAMESHCUBIC].unit = "m";

	commands.insert(CMD_DELMESH, CommandSpecifier(CMD_DELMESH), "delmesh");
	commands[CMD_DELMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delmesh</b> <i>meshname</i>";
	commands[CMD_DELMESH].descr = "[tc0,0.5,0.5,1/tc]Delete mesh with given name.";

	commands.insert(CMD_RENAMEMESH, CommandSpecifier(CMD_RENAMEMESH), "renamemesh");
	commands[CMD_RENAMEMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>renamemesh</b> <i>(old_name) new_name</i>";
	commands[CMD_RENAMEMESH].descr = "[tc0,0.5,0.5,1/tc]Rename mesh. If old_name not specified then the focused mesh is renamed.";

	commands.insert(CMD_MESHFOCUS, CommandSpecifier(CMD_MESHFOCUS), "meshfocus");
	commands[CMD_MESHFOCUS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshfocus</b> <i>meshname</i>";
	commands[CMD_MESHFOCUS].descr = "[tc0,0.5,0.5,1/tc]Change mesh focus to given mesh name.";
	commands[CMD_MESHFOCUS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshname</i> - return name of focused mesh.";

	commands.insert(CMD_MESHFOCUS2, CommandSpecifier(CMD_MESHFOCUS2), "meshfocus2");
	commands[CMD_MESHFOCUS2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshfocus2</b> <i>meshname</i>";
	commands[CMD_MESHFOCUS2].descr = "[tc0,0.5,0.5,1/tc]Change mesh focus to given mesh name but do not change camera orientation.";
	commands[CMD_MESHFOCUS2].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshname</i> - return name of focused mesh.";

	commands.insert(CMD_MESHRECT, CommandSpecifier(CMD_MESHRECT), "meshrect");
	commands[CMD_MESHRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshrect</b> <i>(meshname) rectangle</i>";
	commands[CMD_MESHRECT].limits = { {Any(), Any()}, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_MESHRECT].descr = "[tc0,0.5,0.5,1/tc]Change rectangle of meshname (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin. If meshname not given use the focused mesh.";
	commands[CMD_MESHRECT].unit = "m";
	commands[CMD_MESHRECT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>rectangle</i> - return rectangle of focused mesh, or of meshname if given.";

	commands.insert(CMD_SHIFTDIPOLE, CommandSpecifier(CMD_SHIFTDIPOLE), "shiftdipole");
	commands[CMD_SHIFTDIPOLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshrect</b> <i>meshname shift</i>";
	commands[CMD_SHIFTDIPOLE].limits = { {Any(), Any()}, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) } };
	commands[CMD_SHIFTDIPOLE].descr = "[tc0,0.5,0.5,1/tc]Shift the dipole mesh by the given amount (will only take effect when called on a dipole mesh). shift is specified using x y z coordinates.";
	commands[CMD_SHIFTDIPOLE].unit = "m";

	commands.insert(CMD_DIPOLEVELOCITY, CommandSpecifier(CMD_DIPOLEVELOCITY), "dipolevelocity");
	commands[CMD_DIPOLEVELOCITY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dipolevelocity</b> <i>meshname velocity (clipping)</i>";
	commands[CMD_DIPOLEVELOCITY].limits = { {Any(), Any()}, { DBL3(-DIPOLEMAXVELOCITY), DBL3(DIPOLEMAXVELOCITY) }, {DBL3(), DBL3(MAXSIMSPACE)} };
	commands[CMD_DIPOLEVELOCITY].descr = "[tc0,0.5,0.5,1/tc]Set velocity for dipole shifting algorithm (x, y, z, components). Optionally specify a clipping distance (x, y, z components) - i.e. dipole shift clipped. Default is 0.5 nm; to disable clipping set to zero.";
	commands[CMD_DIPOLEVELOCITY].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>velocity clipping</i>";

	commands.insert(CMD_SCALEMESHRECTS, CommandSpecifier(CMD_SCALEMESHRECTS), "scalemeshrects");
	commands[CMD_SCALEMESHRECTS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>scalemeshrects</b> <i>status</i>";
	commands[CMD_SCALEMESHRECTS].descr = "[tc0,0.5,0.5,1/tc]When changing a mesh rectangle scale and shift all other mesh rectangles in proportion if status set.";
	commands[CMD_SCALEMESHRECTS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_MESH, CommandSpecifier(CMD_MESH), "mesh");
	commands[CMD_MESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mesh</b>";
	commands[CMD_MESH].descr = "[tc0,0.5,0.5,1/tc]Display information for all meshes.";

	commands.insert(CMD_CELLSIZE, CommandSpecifier(CMD_CELLSIZE), "cellsize");
	commands[CMD_CELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>cellsize</b> <i>(meshname) value</i>";
	commands[CMD_CELLSIZE].limits = { {Any(), Any()} , { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_CELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of meshname (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>. If meshname not given use the focused mesh.";
	commands[CMD_CELLSIZE].unit = "m";
	commands[CMD_CELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return cellsize of focused mesh, or of meshname if given.";

	commands.insert(CMD_ECELLSIZE, CommandSpecifier(CMD_ECELLSIZE), "ecellsize");
	commands[CMD_ECELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ecellsize</b> <i>(meshname) value</i>";
	commands[CMD_ECELLSIZE].limits = { {Any(), Any()} , { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_ECELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of meshname for electrical conduction (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>. If meshname not given use the focused mesh.";
	commands[CMD_ECELLSIZE].unit = "m";
	commands[CMD_ECELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return electrical conduction cellsize of focused mesh, or of meshname if given.";

	commands.insert(CMD_TCELLSIZE, CommandSpecifier(CMD_TCELLSIZE), "tcellsize");
	commands[CMD_TCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tcellsize</b> <i>(meshname) value</i>";
	commands[CMD_TCELLSIZE].limits = { {Any(), Any()} , { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_TCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of meshname for thermal conduction (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>. If meshname not given use the focused mesh.";
	commands[CMD_TCELLSIZE].unit = "m";
	commands[CMD_TCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return thermal conduction cellsize of focused mesh, or of meshname if given.";

	commands.insert(CMD_SCELLSIZE, CommandSpecifier(CMD_SCELLSIZE), "scellsize");
	commands[CMD_SCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>scellsize</b> <i>(meshname) value</i>";
	commands[CMD_SCELLSIZE].limits = { {Any(), Any()} , { DBL3(MINMESHSPACE / 2), Any() } };
	commands[CMD_SCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of meshname for stochastic properties (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>. If meshname not given use the focused mesh.";
	commands[CMD_SCELLSIZE].unit = "m";
	commands[CMD_SCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return stochastic properties cellsize of focused mesh, or of meshname if given.";

	commands.insert(CMD_MCELLSIZE, CommandSpecifier(CMD_MCELLSIZE), "mcellsize");
	commands[CMD_MCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcellsize</b> <i>(meshname) value</i>";
	commands[CMD_MCELLSIZE].limits = { {Any(), Any()} , { DBL3(MINMESHSPACE / 2), Any() } };
	commands[CMD_MCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of meshname for mechanical solver (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>. If meshname not given use the focused mesh.";
	commands[CMD_MCELLSIZE].unit = "m";
	commands[CMD_MCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return mechanical cellsize of focused mesh, or of meshname if given.";

	commands.insert(CMD_FMSCELLSIZE, CommandSpecifier(CMD_FMSCELLSIZE), "fmscellsize");
	commands[CMD_FMSCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>fmscellsize</b> <i>value</i>";
	commands[CMD_FMSCELLSIZE].limits = { { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_FMSCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize for ferromagnetic super-mesh (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_FMSCELLSIZE].unit = "m";
	commands[CMD_FMSCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return cellsize for ferromagnetic super-mesh.";

	commands.insert(CMD_ESCELLSIZE, CommandSpecifier(CMD_ESCELLSIZE), "escellsize");
	commands[CMD_ESCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>escellsize</b> <i>value</i>";
	commands[CMD_ESCELLSIZE].limits = { { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_ESCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize for electric super-mesh (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_ESCELLSIZE].unit = "m";
	commands[CMD_ESCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return cellsize for electric super-mesh.";

	commands.insert(CMD_ATOMDMCELLSIZE, CommandSpecifier(CMD_ATOMDMCELLSIZE), "dmcellsize");
	commands[CMD_ATOMDMCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dmcellsize</b> <i>(meshname) value</i>";
	commands[CMD_ATOMDMCELLSIZE].limits = { {Any(), Any()} , { DBL3(MINMESHSPACE / 2), Any() } };
	commands[CMD_ATOMDMCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change demagnetizing field macrocell size of meshname, for atomistic meshes (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>. If meshname not given use the focused mesh.";
	commands[CMD_ATOMDMCELLSIZE].unit = "m";
	commands[CMD_ATOMDMCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return demagnetizing field macrocell size of focused mesh, or of meshname if given.";

	commands.insert(CMD_DELRECT, CommandSpecifier(CMD_DELRECT), "delrect");
	commands[CMD_DELRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delrect</b> <i>((meshname) (rectangle))</i>";
	commands[CMD_DELRECT].limits = { { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_DELRECT].descr = "[tc0,0.5,0.5,1/tc]Void rectangle (m) within given mesh (focused mesh if not specified). The rectangle coordinates are relative to mesh. If rectangle not given void entire mesh.";
	commands[CMD_DELRECT].unit = "m";

	commands.insert(CMD_ADDRECT, CommandSpecifier(CMD_ADDRECT), "addrect");
	commands[CMD_ADDRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addrect</b> <i>(meshname) rectangle</i>";
	commands[CMD_ADDRECT].limits = { { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_ADDRECT].descr = "[tc0,0.5,0.5,1/tc]Fill rectangle (m) within given mesh (focused mesh if not specified). The rectangle coordinates are relative to mesh.";
	commands[CMD_ADDRECT].unit = "m";

	commands.insert(CMD_RESETMESH, CommandSpecifier(CMD_RESETMESH), "resetmesh");
	commands[CMD_RESETMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>resetmesh</b> <i>(meshname)</i>";
	commands[CMD_RESETMESH].limits = { { Any(), Any() } };
	commands[CMD_RESETMESH].descr = "[tc0,0.5,0.5,1/tc]Reset to constant magnetization in given mesh (focused mesh if not specified).";

	commands.insert(CMD_LOADMASKFILE, CommandSpecifier(CMD_LOADMASKFILE), "loadmaskfile");
	commands[CMD_LOADMASKFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadmaskfile</b> <i>(meshname) (z_depth) (directory/)filename</i>";
	commands[CMD_LOADMASKFILE].descr = "[tc0,0.5,0.5,1/tc]Apply .png mask file to magnetization in given mesh (focused mesh if not specified); i.e. transfer shape from .png file to mesh - white means empty cells. If image is in grayscale then void cells up to given depth top down (z_depth > 0) or down up (z_depth < 0). If z-depth = 0 then void top down up to all z cells.";
	commands[CMD_LOADMASKFILE].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_INDIVIDUALMASKSHAPE, CommandSpecifier(CMD_INDIVIDUALMASKSHAPE), "individualshape");
	commands[CMD_INDIVIDUALMASKSHAPE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>individualshape</b> <i>status</i>";
	commands[CMD_INDIVIDUALMASKSHAPE].limits = { { int(0), int(1) } };
	commands[CMD_INDIVIDUALMASKSHAPE].descr = "[tc0,0.5,0.5,1/tc]When changing the shape inside a mesh set this flag to true so the shape is applied only to the primary displayed physical quantity. If set to false then all relevant physical quantities are shaped.";
	commands[CMD_INDIVIDUALMASKSHAPE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_SETANGLE, CommandSpecifier(CMD_SETANGLE), "setangle");
	commands[CMD_SETANGLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setangle</b> <i>(meshname) polar azimuthal</i>";
	commands[CMD_SETANGLE].limits = { { Any(), Any() }, { double(-360.0), double(360.0) },{ double(-360.0), double(360.0) } };
	commands[CMD_SETANGLE].descr = "[tc0,0.5,0.5,1/tc]Set magnetization angle (deg.) in mesh uniformly using polar coordinates (all magnetic meshes if name not given).";

	commands.insert(CMD_FLOWER, CommandSpecifier(CMD_FLOWER), "flower");
	commands[CMD_FLOWER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>flower</b> <i>(meshname) (direction (radius thickness centre))</i>";
	commands[CMD_FLOWER].limits = { { Any(), Any() }, {int(-1), int(+1)}, { double(0.0), Any() }, { double(0.0), Any() }, { DBL3(), Any() } };
	commands[CMD_FLOWER].descr = "[tc0,0.5,0.5,1/tc]Set magnetization in a flower state with given direction (-1 or +1) starting at given centre, over given radius and thickness. If not specified then entire mesh is used, centred. If mesh name not specified, the focused mesh is used.";
	commands[CMD_FLOWER].unit = "m";

	commands.insert(CMD_ONION, CommandSpecifier(CMD_ONION), "onion");
	commands[CMD_ONION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>onion</b> <i>(meshname) (direction (radius1 radius2 thickness centre))</i>";
	commands[CMD_ONION].limits = { { Any(), Any() }, {int(-1), int(+1)}, { double(0.0), Any() }, { double(0.0), Any() }, { double(0.0), Any() }, { DBL3(), Any() } };
	commands[CMD_ONION].descr = "[tc0,0.5,0.5,1/tc]Set magnetization in an onion state with given direction (-1 clockwise, +1 anti-clockwise) starting at given centre, from radius1 to radius2 and thickness. If not specified then entire mesh is used, centred. If mesh name not specified, the focused mesh is used.";
	commands[CMD_ONION].unit = "m";

	commands.insert(CMD_CROSSTIE, CommandSpecifier(CMD_CROSSTIE), "crosstie");
	commands[CMD_CROSSTIE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>crosstie</b> <i>(meshname) (direction (radius thickness centre))</i>";
	commands[CMD_CROSSTIE].limits = { { Any(), Any() }, {int(-1), int(+1)}, { double(0.0), Any() }, { double(0.0), Any() }, { DBL3(), Any() } };
	commands[CMD_CROSSTIE].descr = "[tc0,0.5,0.5,1/tc]Set magnetization in a crosstie state with given direction (-1 or +1) starting at given centre, over given radius and thickness. If not specified then entire mesh is used, centred. If mesh name not specified, the focused mesh is used.";
	commands[CMD_CROSSTIE].unit = "m";

	commands.insert(CMD_SETOBJECTANGLE, CommandSpecifier(CMD_SETOBJECTANGLE), "setobjectangle");
	commands[CMD_SETOBJECTANGLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setobjectangle</b> <i>(meshname) polar azimuthal position</i>";
	commands[CMD_SETOBJECTANGLE].limits = { { Any(), Any() }, { double(-360.0), double(360.0) },{ double(-360.0), double(360.0) }, { DBL3(), Any() } };
	commands[CMD_SETOBJECTANGLE].descr = "[tc0,0.5,0.5,1/tc]Set magnetization angle in solid object only containing given relative position (x y z) uniformly using polar coordinates. If mesh name not specified, the focused mesh is used.";
	commands[CMD_SETOBJECTANGLE].unit = "m";

	commands.insert(CMD_RANDOM, CommandSpecifier(CMD_RANDOM), "random");
	commands[CMD_RANDOM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>random</b> <i>(meshname (seed))</i>";
	commands[CMD_RANDOM].limits = { { Any(), Any() }, { int(1), Any() } };
	commands[CMD_RANDOM].descr = "[tc0,0.5,0.5,1/tc]Set random magnetization distribution in mesh, with pseud-random number generator seed (default 1 if not specified). If mesh name not specified, the focused mesh is used.";

	commands.insert(CMD_RANDOMXY, CommandSpecifier(CMD_RANDOMXY), "randomxy");
	commands[CMD_RANDOMXY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>randomxy</b> <i>(meshname (seed))</i>";
	commands[CMD_RANDOMXY].limits = { { Any(), Any() }, { int(1), Any() } };
	commands[CMD_RANDOMXY].descr = "[tc0,0.5,0.5,1/tc]Set random magnetization distribution in the xy plane in mesh, with pseud-random number generator seed (default 1 if not specified). If mesh name not specified, the focused mesh is used.";

	commands.insert(CMD_SETRECT, CommandSpecifier(CMD_SETRECT), "setrect");
	commands[CMD_SETRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setrect</b> <i>(meshname) polar azimuthal rectangle</i>";
	commands[CMD_SETRECT].limits = { { Any(), Any() }, { double(-360.0), double(360.0) },{ double(-360.0), double(360.0) }, { Rect(), Any() } };
	commands[CMD_SETRECT].descr = "[tc0,0.5,0.5,1/tc]Set magnetization angle in given rectangle of mesh (relative coordinates) uniformly using polar coordinates. If mesh name not specified, the focused mesh is used.";
	commands[CMD_SETRECT].unit = "m";

	commands.insert(CMD_INVERTMAG, CommandSpecifier(CMD_INVERTMAG), "invertmag");
	commands[CMD_INVERTMAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>invertmag</b> <i>(meshname) (components)</i>";
	commands[CMD_INVERTMAG].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_INVERTMAG].descr = "[tc0,0.5,0.5,1/tc]Invert magnetization direction. If mesh name not specified, the focused mesh is used. You can choose to invert just one or two components instead of the entire vector: specify components as x y z, e.g. invertmag x";

	commands.insert(CMD_MIRRORMAG, CommandSpecifier(CMD_MIRRORMAG), "mirrormag");
	commands[CMD_MIRRORMAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mirrormag</b> <i>(meshname) axis</i>";
	commands[CMD_MIRRORMAG].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_MIRRORMAG].descr = "[tc0,0.5,0.5,1/tc]Mirror magnetization in a given axis, specified as x, y, z, e.g. mirrormag x. If mesh name not specified, the focused mesh is used.";

	commands.insert(CMD_DWALL, CommandSpecifier(CMD_DWALL), "dwall");
	commands[CMD_DWALL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dwall</b> <i>(meshname) longitudinal transverse width position</i>";
	commands[CMD_DWALL].limits = { { Any(), Any() }, { Any(), Any() } , { Any(), Any() }, { double(0.0), Any() }, { double(0.0), Any() } };
	commands[CMD_DWALL].unit = "m";
	commands[CMD_DWALL].descr = "[tc0,0.5,0.5,1/tc]Create an idealised domain wall (tanh profile for longitudinal component, 1/cosh profile for transverse component) along the x-axis direction in the given mesh (focused mesh if not specified). For longitudinal and transverse specify the components of magnetization as x, -x, y, -y, z, -z, i.e. specify using these std::string literals. For width and position use metric units.";

	commands.insert(CMD_VORTEX, CommandSpecifier(CMD_VORTEX), "vortex");
	commands[CMD_VORTEX].usage = "[tc0,0.5,0,1/tc]USAGE : <b>vortex</b> <i>(meshname) longitudinal rotation core (rectangle)</i>";
	commands[CMD_VORTEX].limits = { { Any(), Any() }, { int(-1), int(+1) } , { int(-1), int(+1) }, { int(-1), int(+1) }, { Rect(), Any() } };
	commands[CMD_VORTEX].unit = "m";
	commands[CMD_VORTEX].descr = "[tc0,0.5,0.5,1/tc]Create a vortex domain wall with settings: longitudinal (-1: tail-to-tail, 1: head-to-head), rotation (-1: clockwise, 1: counter-clockwise), core (-1: down, 1: up). The vortex may be set in the given rectangle (entire mesh if not given), in the given mesh (focused mesh if not specified).";

	commands.insert(CMD_SKYRMION, CommandSpecifier(CMD_SKYRMION), "skyrmion");
	commands[CMD_SKYRMION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyrmion</b> <i>(meshname) core chirality diameter position</i>";
	commands[CMD_SKYRMION].limits = { { Any(), Any() }, { int(-1), int(1) }, { int(-1), int(1) }, { double(1e-9), Any() }, { DBL2(), Any() } };
	commands[CMD_SKYRMION].unit = "m";
	commands[CMD_SKYRMION].descr = "[tc0,0.5,0.5,1/tc]Create an idealised Neel-type skyrmion with given diameter and centre position in the x-y plane (2 relative coordinates needed only) of the given mesh (focused mesh if not specified). Core specifies the skyrmion core direction: -1 for down, 1 for up. Chirality specifies the radial direction rotation: 1 for towards core, -1 away from core. For diameter and position use metric units.";

	commands.insert(CMD_SKYRMIONBLOCH, CommandSpecifier(CMD_SKYRMIONBLOCH), "skyrmionbloch");
	commands[CMD_SKYRMIONBLOCH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyrmionbloch</b> <i>(meshname) core chirality diameter position</i>";
	commands[CMD_SKYRMIONBLOCH].limits = { { Any(), Any() }, { int(-1), int(1) }, { int(-1), int(1) }, { double(1e-9), Any() }, { DBL2(), Any() } };
	commands[CMD_SKYRMIONBLOCH].unit = "m";
	commands[CMD_SKYRMIONBLOCH].descr = "[tc0,0.5,0.5,1/tc]Create an idealised Bloch-type skyrmion with given diameter and centre position in the x-y plane (2 relative coordinates needed only) of the given mesh (focused mesh if not specified). Core specifies the skyrmion core direction: -1 for down, 1 for up. Chirality specifies the radial direction rotation: 1 for clockwise, -1 for anti-clockwise. For diameter and position use metric units.";

	commands.insert(CMD_PBC, CommandSpecifier(CMD_PBC), "pbc");
	commands[CMD_PBC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>pbc</b> <i>(meshname) flag images</i>";
	commands[CMD_PBC].descr = "[tc0,0.5,0.5,1/tc]Set periodic boundary conditions for magnetization in given mesh (must be magnetic, or supermesh; supermesh/all meshes if name not given). Flags specify types of perodic boundary conditions: x, y, or z; images specify the number of mesh images to use either side for the given direction when calculating the demagnetising kernel - a value of zero disables pbc. e.g. ""pbc x 10"" sets x periodic boundary conditions with 10 images either side for all meshes; ""pbc x 0"" clears pbc for the x axis.";
	commands[CMD_PBC].limits = { { Any(), Any() }, { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETFIELD, CommandSpecifier(CMD_SETFIELD), "setfield");
	commands[CMD_SETFIELD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setfield</b> <i>(meshname) magnitude polar azimuthal</i>";
	commands[CMD_SETFIELD].limits = { { Any(), Any() }, { DBL3(-MAXFIELD, -360.0, -360.0), DBL3(MAXFIELD, 360.0, 360.0) } };
	commands[CMD_SETFIELD].descr = "[tc0,0.5,0.5,1/tc]Set uniform magnetic field (A/m) using polar coordinates. If mesh name not specified, this is set for all magnetic meshes - must have Zeeman module added.";
	commands[CMD_SETFIELD].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i><Ha_x, Ha_y, Ha_z></i> - applied field in Cartesian coordinates for, or of meshname if given.";
	commands[CMD_SETFIELD].unit = "A/m";

	commands.insert(CMD_SETSTRESS, CommandSpecifier(CMD_SETSTRESS), "setstress");
	commands[CMD_SETSTRESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setstress</b> <i>(meshname) magnitude polar azimuthal</i>";
	commands[CMD_SETSTRESS].limits = { { Any(), Any() }, { DBL3(-MAXSTRESS, -360.0, -360.0), DBL3(MAXSTRESS, 360.0, 360.0) } };
	commands[CMD_SETSTRESS].descr = "[tc0,0.5,0.5,1/tc]Set uniform mechanical stress (Pa) using polar coordinates. If mesh name not specified, this is set for all magnetic meshes - must have MElastic module added.";
	commands[CMD_SETSTRESS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i><Tsig_x, Tsig_y, Tsig_z></i> - applied mechanical stress in Cartesian coordinates for focused mesh, or of meshname if given.";
	commands[CMD_SETSTRESS].unit = "Pa";

	commands.insert(CMD_MODULES, CommandSpecifier(CMD_MODULES), "modules");
	commands[CMD_MODULES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>modules</b> <i>(meshname (modules...))";
	commands[CMD_MODULES].descr = "[tc0,0.5,0.5,1/tc]Show interactive list of available and currently set modules. If meshname specified then return all set modules names, and set given modules if specified.";
	commands[CMD_MODULES].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_MODULES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>modulenames</i>";

	commands.insert(CMD_ADDMODULE, CommandSpecifier(CMD_ADDMODULE), "addmodule");
	commands[CMD_ADDMODULE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addmodule</b> <i>(meshname) handle</i>";
	commands[CMD_ADDMODULE].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_ADDMODULE].descr = "[tc0,0.5,0.5,1/tc]Add module with given handle to named mesh. If mesh name not specified, this is set for all applicable meshes.";

	commands.insert(CMD_DELMODULE, CommandSpecifier(CMD_DELMODULE), "delmodule");
	commands[CMD_DELMODULE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delmodule</b> <i>(meshname) handle</i>";
	commands[CMD_DELMODULE].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_DELMODULE].descr = "[tc0,0.5,0.5,1/tc]Delete module with given handle from named mesh (focused mesh if not specified).";

	commands.insert(CMD_MULTICONV, CommandSpecifier(CMD_MULTICONV), "multiconvolution");
	commands[CMD_MULTICONV].usage = "[tc0,0.5,0,1/tc]USAGE : <b>multiconvolution</b> <i>status</i>";
	commands[CMD_MULTICONV].descr = "[tc0,0.5,0.5,1/tc]Switch between multi-layered convolution (true) and supermesh convolution (false).";

	commands.insert(CMD_2DMULTICONV, CommandSpecifier(CMD_2DMULTICONV), "2dmulticonvolution");
	commands[CMD_2DMULTICONV].usage = "[tc0,0.5,0,1/tc]USAGE : <b>2dmulticonvolution</b> <i>status</i>";
	commands[CMD_2DMULTICONV].descr = "[tc0,0.5,0.5,1/tc]Switch to multi-layered convolution and force it to 2D layering in each mesh (2), or 2D convolution for each mesh (1), or allow 3D (0).";
	commands[CMD_2DMULTICONV].limits = { { int(0), int(2) } };

	commands.insert(CMD_NCOMMONSTATUS, CommandSpecifier(CMD_NCOMMONSTATUS), "ncommonstatus");
	commands[CMD_NCOMMONSTATUS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ncommonstatus</b> <i>status</i>";
	commands[CMD_NCOMMONSTATUS].descr = "[tc0,0.5,0.5,1/tc]Switch to multi-layered convolution and force it to user-defined discretisation (status = true), or default discretisation (status = false).";
	commands[CMD_NCOMMONSTATUS].limits = { { int(0), int(1) } };

	commands.insert(CMD_NCOMMON, CommandSpecifier(CMD_NCOMMON), "ncommon");
	commands[CMD_NCOMMON].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ncommon</b> <i>sizes</i>";
	commands[CMD_NCOMMON].descr = "[tc0,0.5,0.5,1/tc]Switch to multi-layered convolution and force it to user-defined discretisation, specifying sizes as nx ny nz.";
	commands[CMD_NCOMMON].limits = { { INT3(1), Any() } };

	commands.insert(CMD_EXCLUDEMULTICONVDEMAG, CommandSpecifier(CMD_EXCLUDEMULTICONVDEMAG), "excludemulticonvdemag");
	commands[CMD_EXCLUDEMULTICONVDEMAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>excludemulticonvdemag</b> <i>(meshname) status</i>";
	commands[CMD_EXCLUDEMULTICONVDEMAG].descr = "[tc0,0.5,0.5,1/tc]Set exclusion status (0 or 1) of named mesh from multi-layered demag convolution (focused mesh if not specified).";
	commands[CMD_EXCLUDEMULTICONVDEMAG].limits = { { Any(), Any() }, { int(0), int(1) } };

	commands.insert(CMD_GPUKERNELS, CommandSpecifier(CMD_GPUKERNELS), "gpukernels");
	commands[CMD_GPUKERNELS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>gpukernels</b> <i>status</i>";
	commands[CMD_GPUKERNELS].descr = "[tc0,0.5,0.5,1/tc]When in CUDA mode calculate demagnetization kernels initialization on the GPU (1) or on the CPU (0).";
	commands[CMD_GPUKERNELS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_ODE, CommandSpecifier(CMD_ODE), "ode");
	commands[CMD_ODE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ode</b>";
	commands[CMD_ODE].descr = "[tc0,0.5,0.5,1/tc]Show interactive list of available and currently set ODEs and evaluation methods.";

	commands.insert(CMD_EVALSPEEDUP, CommandSpecifier(CMD_EVALSPEEDUP), "evalspeedup");
	commands[CMD_EVALSPEEDUP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>evalspeedup</b> <i>level</i>";
	commands[CMD_EVALSPEEDUP].limits = { { int(EVALSPEEDUP_NONE), int(EVALSPEEDUP_NUMENTRIES) - 1 } };
	commands[CMD_EVALSPEEDUP].descr = "[tc0,0.5,0.5,1/tc]Enable/disable evaluation speedup by extrapolating demag field at evaluation substeps from previous field updates. Status levels: 0 (no speedup), 1 (step), 2 (linear), 3 (quadratic), 4 (cubic), 5 (quartic), 6 (quintic). If enabling speedup strongly recommended to always use quadratic type.";
	commands[CMD_EVALSPEEDUP].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>level</i>";

	commands.insert(CMD_SETODE, CommandSpecifier(CMD_SETODE), "setode");
	commands[CMD_SETODE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setode</b> <i>equation evaluation</i>";
	commands[CMD_SETODE].descr = "[tc0,0.5,0.5,1/tc]Set differential equation to solve in both micromagnetic and atomistic meshes, and method used to solve it (same method is applied to micromagnetic and atomistic meshes).";

	commands.insert(CMD_SETODEEVAL, CommandSpecifier(CMD_SETODEEVAL), "setodeeval");
	commands[CMD_SETODEEVAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setodeeval</b> <i>evaluation</i>";
	commands[CMD_SETODEEVAL].descr = "[tc0,0.5,0.5,1/tc]Set differential equation method used to solve it (same method is applied to micromagnetic and atomistic meshes).";

	commands.insert(CMD_SETATOMODE, CommandSpecifier(CMD_SETATOMODE), "setatomode");
	commands[CMD_SETATOMODE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setatomode</b> <i>equation evaluation</i>";
	commands[CMD_SETATOMODE].descr = "[tc0,0.5,0.5,1/tc]Set differential equation to solve in atomistic meshes, and method used to solve it (same method is applied to micromagnetic and atomistic meshes).";
	
	commands.insert(CMD_SETDT, CommandSpecifier(CMD_SETDT), "setdt");
	commands[CMD_SETDT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdt</b> <i>value</i>";
	commands[CMD_SETDT].limits = { { double(MINTIMESTEP), double(MAXTIMESTEP) } };
	commands[CMD_SETDT].descr = "[tc0,0.5,0.5,1/tc]Set differential equation time-step (only applicable to fixed time-step methods).";
	commands[CMD_SETDT].unit = "s";
	commands[CMD_SETDT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>dT</i>";

	commands.insert(CMD_ASTEPCTRL, CommandSpecifier(CMD_ASTEPCTRL), "astepctrl");
	commands[CMD_ASTEPCTRL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>astepctrl</b> <i>err_fail dT_incr dT_min dT_max</i>";
	commands[CMD_ASTEPCTRL].limits = { 
		{double(MINODERELERROR), double(MAXODERELERROR)},
		{double(1.0), double(100.0)},
		{double(MINTIMESTEP), double(MAXTIMESTEP)},
		{double(MINTIMESTEP), double(MAXTIMESTEP)} };
	commands[CMD_ASTEPCTRL].descr = "[tc0,0.5,0.5,1/tc]Set parameters for adaptive time step control: err_fail - repeat step above this (the tolerance), dT_incr - limit multiplicative increase in dT using this, dT_min, dT_max - dT bounds.";
	commands[CMD_ASTEPCTRL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>err_fail dT_incr dT_min dT_max</i>";

	commands.insert(CMD_SHOWDATA, CommandSpecifier(CMD_SHOWDATA), "showdata");
	commands[CMD_SHOWDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showdata</b> <i>(meshname) dataname (rectangle)</i>";
	commands[CMD_SHOWDATA].descr = "[tc0,0.5,0.5,1/tc]Show value(s) for dataname. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, focused mesh is used with entire mesh rectangle.";
	commands[CMD_SHOWDATA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>varies</i>";
	commands[CMD_SHOWDATA].unit = "m";
	commands[CMD_SHOWDATA].limits = { { Any(), Any() }, { Any(), Any() }, { Rect(), Any() } };

	commands.insert(CMD_DATA, CommandSpecifier(CMD_DATA), "data");
	commands[CMD_DATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>data</b>";
	commands[CMD_DATA].descr = "[tc0,0.5,0.5,1/tc]Shows list of currently set output data and available data.";
	commands[CMD_DATA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>number of set output data fields</i>";

	commands.insert(CMD_ADDDATA, CommandSpecifier(CMD_ADDDATA), "adddata");
	commands[CMD_ADDDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>adddata</b> <i>(meshname) dataname (rectangle)</i>";
	commands[CMD_ADDDATA].limits = { { Any(), Any() }, { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_ADDDATA].descr = "[tc0,0.5,0.5,1/tc]Add dataname to list of output data. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, focused mesh is used with entire mesh rectangle.";
	commands[CMD_ADDDATA].unit = "m";

	commands.insert(CMD_SETDATA, CommandSpecifier(CMD_SETDATA), "setdata");
	commands[CMD_SETDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdata</b> <i>(meshname) dataname (rectangle)</i>";
	commands[CMD_SETDATA].limits = { { Any(), Any() }, { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_SETDATA].descr = "[tc0,0.5,0.5,1/tc]Delete all currently set output data and set dataname to list of output data. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, focused mesh is used with entire mesh rectangle.";
	commands[CMD_SETDATA].unit = "m";

	commands.insert(CMD_DELDATA, CommandSpecifier(CMD_DELDATA), "deldata");
	commands[CMD_DELDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>deldata</b> <i>index</i>";
	commands[CMD_DELDATA].limits = { { int(-1), Any() } };
	commands[CMD_DELDATA].descr = "[tc0,0.5,0.5,1/tc]Delete data from list of output data at index number. If index number is -1 then delete all data fields, leaving just a default time data field - there must always be at least 1 output data field.";

	commands.insert(CMD_EDITDATA, CommandSpecifier(CMD_EDITDATA), "editdata");
	commands[CMD_EDITDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editdata</b> <i>index (meshname) dataname (rectangle)</i>";
	commands[CMD_EDITDATA].limits = { { Any(), Any() }, { int(0), Any() }, { Any(), Any() }, { Rect(), Any() } };		//this order is correct as meshname will always be set first (default focused mesh inserted if not)
	commands[CMD_EDITDATA].descr = "[tc0,0.5,0.5,1/tc]Edit entry in list of output data at given index in list. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, focused mesh is used with entire mesh rectangle.";
	commands[CMD_EDITDATA].unit = "m";

	commands.insert(CMD_SAVEDATA, CommandSpecifier(CMD_SAVEDATA), "savedata");
	commands[CMD_SAVEDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savedata</b> <i>(append_option)</i>";
	commands[CMD_SAVEDATA].limits = { { int(-1), Any() } };
	commands[CMD_SAVEDATA].descr = "[tc0,0.5,0.5,1/tc]Save data (as currently configured in output data - see setdata/adddata commands) to output file (see savedatafile command). append_option = -1 : use default append behaviour (i.e. new file with header created for first save, then append).  append_option = 0 : make new file then save.  append_option = 1 : append to file.";

	commands.insert(CMD_ADDPINNEDDATA, CommandSpecifier(CMD_ADDPINNEDDATA), "addpinneddata");
	commands[CMD_ADDPINNEDDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addpinneddata</b> <i>(meshname) dataname (rectangle))</i>";
	commands[CMD_ADDPINNEDDATA].descr = "[tc0,0.5,0.5,1/tc]Add new entry in data box (at the end) with given dataname and meshname if applicable. A rectangle may also be specified if applicable, however this will not be shown in the data box.";

	commands.insert(CMD_DELPINNEDDATA, CommandSpecifier(CMD_DELPINNEDDATA), "delpinneddata");
	commands[CMD_DELPINNEDDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delpinneddata</b> <i>index</i>";
	commands[CMD_DELPINNEDDATA].limits = { { int(0), Any() } };
	commands[CMD_DELPINNEDDATA].descr = "[tc0,0.5,0.5,1/tc]Delete entry in data box at given index (index in order of appearance in data box from 0 up).";

	commands.insert(CMD_CHDIR, CommandSpecifier(CMD_CHDIR), "chdir");
	commands[CMD_CHDIR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>chdir</b> <i>directory</i>";
	commands[CMD_CHDIR].descr = "[tc0,0.5,0.5,1/tc]Change working directory.";
	commands[CMD_CHDIR].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>directory</i>";

	commands.insert(CMD_SAVEDATAFILE, CommandSpecifier(CMD_SAVEDATAFILE), "savedatafile");
	commands[CMD_SAVEDATAFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savedatafile</b> <i>(directory/)filename</i>";
	commands[CMD_SAVEDATAFILE].descr = "[tc0,0.5,0.5,1/tc]Change output data file (and working directory if specified).";
	commands[CMD_SAVEDATAFILE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>filename</i>";

	commands.insert(CMD_SAVECOMMENT, CommandSpecifier(CMD_SAVECOMMENT), "savecomment");
	commands[CMD_SAVECOMMENT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savecomment</b> <i>(directory/)filename comment</i>";
	commands[CMD_SAVECOMMENT].descr = "[tc0,0.5,0.5,1/tc]Save comment in given file by appending to it.";

	commands.insert(CMD_SAVEIMAGEFILE, CommandSpecifier(CMD_SAVEIMAGEFILE), "saveimagefile");
	commands[CMD_SAVEIMAGEFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveimagefile</b> <i>(directory/)filename</i>";
	commands[CMD_SAVEIMAGEFILE].descr = "[tc0,0.5,0.5,1/tc]Change image file base (and working directory if specified).";
	commands[CMD_SAVEIMAGEFILE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>filename</i>";

	commands.insert(CMD_DATASAVEFLAG, CommandSpecifier(CMD_DATASAVEFLAG), "savedataflag");
	commands[CMD_DATASAVEFLAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savedataflag</b> <i>status</i>";
	commands[CMD_DATASAVEFLAG].descr = "[tc0,0.5,0.5,1/tc]Set data saving flag status.";
	commands[CMD_DATASAVEFLAG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_IMAGESAVEFLAG, CommandSpecifier(CMD_IMAGESAVEFLAG), "saveimageflag");
	commands[CMD_IMAGESAVEFLAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveimageflag</b> <i>status</i>";
	commands[CMD_IMAGESAVEFLAG].descr = "[tc0,0.5,0.5,1/tc]Set image saving flag status.";
	commands[CMD_IMAGESAVEFLAG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_DISKBUFFERLINES, CommandSpecifier(CMD_DISKBUFFERLINES), "diskbufferlines");
	commands[CMD_DISKBUFFERLINES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>diskbufferlines</b> <i>lines</i>";
	commands[CMD_DISKBUFFERLINES].limits = { { int(1), Any() } };
	commands[CMD_DISKBUFFERLINES].descr = "[tc0,0.5,0.5,1/tc]Set output disk buffer number of lines - default is 100. You shouldn't need to adjust this.";
	commands[CMD_DISKBUFFERLINES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>bufferlines</i>";

	commands.insert(CMD_DATAPRECISION, CommandSpecifier(CMD_DATAPRECISION), "dataprecision");
	commands[CMD_DATAPRECISION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dataprecision</b> <i>precision</i>";
	commands[CMD_DATAPRECISION].limits = { { int(SAVEDATAPRECISION), Any() } };
	commands[CMD_DATAPRECISION].descr = "[tc0,0.5,0.5,1/tc]Set number of significant figures used for saving data to file as well as data display. The default (and minimum) is 6.";
	commands[CMD_DATAPRECISION].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>precision</i>";

	commands.insert(CMD_STAGES, CommandSpecifier(CMD_STAGES), "stages");
	commands[CMD_STAGES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>stages</b>";
	commands[CMD_STAGES].descr = "[tc0,0.5,0.5,1/tc]Shows list of currently set simulation stages and available stage types.";
	commands[CMD_STAGES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>number of set stages</i>";

	commands.insert(CMD_ADDSTAGE, CommandSpecifier(CMD_ADDSTAGE), "addstage");
	commands[CMD_ADDSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addstage</b> <i>(meshname) stagetype</i>";
	commands[CMD_ADDSTAGE].descr = "[tc0,0.5,0.5,1/tc]Add a generic stage type to the simulation schedule with name stagetype, specifying a meshname if needed (if not specified and required, focused mesh is used).";
	commands[CMD_ADDSTAGE].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETSTAGE, CommandSpecifier(CMD_SETSTAGE), "setstage");
	commands[CMD_SETSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setstage</b> <i>(meshname) stagetype</i>";
	commands[CMD_SETSTAGE].descr = "[tc0,0.5,0.5,1/tc]Delete all currently set stages, and set a new generic stage type to the simulation schedule with name stagetype, specifying a meshname if needed (if not specified and required, focused mesh is used).";
	commands[CMD_SETSTAGE].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_DELSTAGE, CommandSpecifier(CMD_DELSTAGE), "delstage");
	commands[CMD_DELSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delstage</b> <i>index</i>";
	commands[CMD_DELSTAGE].limits = { { int(-1), Any() } };
	commands[CMD_DELSTAGE].descr = "[tc0,0.5,0.5,1/tc]Delete stage from simulation schedule at index number. If index number is -1 then delete all stages, leaving just a default Relax stage - there must always be at least 1 stage set.";

	commands.insert(CMD_EDITSTAGE, CommandSpecifier(CMD_EDITSTAGE), "editstage");
	commands[CMD_EDITSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editstage</b> <i>index (meshname) stagetype</i>";
	commands[CMD_EDITSTAGE].limits = { { Any(), Any() }, { int(0), Any() }, { Any(), Any() } };		//correct see CMD_EDITDATA
	commands[CMD_EDITSTAGE].descr = "[tc0,0.5,0.5,1/tc]Edit stage type from simulation schedule at index number.";

	commands.insert(CMD_EDITSTAGEVALUE, CommandSpecifier(CMD_EDITSTAGEVALUE), "editstagevalue");
	commands[CMD_EDITSTAGEVALUE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editstagevalue</b> <i>index value</i>";
	commands[CMD_EDITSTAGEVALUE].limits = { { int(0), Any() }, { Any(), Any() } };
	commands[CMD_EDITSTAGEVALUE].descr = "[tc0,0.5,0.5,1/tc]Edit stage setting value in simulation schedule. The value type depends on the stage type.";

	commands.insert(CMD_SETSTAGEVALUE, CommandSpecifier(CMD_SETSTAGEVALUE), "setstagevalue");
	commands[CMD_SETSTAGEVALUE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setstagevalue</b> <i>index</i>";
	commands[CMD_SETSTAGEVALUE].limits = { { int(0), Any() } };
	commands[CMD_SETSTAGEVALUE].descr = "[tc0,0.5,0.5,1/tc]Set configured value for given stage index.";

	commands.insert(CMD_EDITSTAGESTOP, CommandSpecifier(CMD_EDITSTAGESTOP), "editstagestop");
	commands[CMD_EDITSTAGESTOP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editstagestop</b> <i>index stoptype (stopvalue)</i>";
	commands[CMD_EDITSTAGESTOP].descr = "[tc0,0.5,0.5,1/tc]Edit stage/step stopping condition in simulation schedule. Use index < 0 to set condition for all stages.";

	commands.insert(CMD_EDITDATASAVE, CommandSpecifier(CMD_EDITDATASAVE), "editdatasave");
	commands[CMD_EDITDATASAVE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editdatasave</b> <i>index savetype (savevalue)</i>";
	commands[CMD_EDITDATASAVE].descr = "[tc0,0.5,0.5,1/tc]Edit data saving condition in simulation schedule. Use index < 0 to set condition for all stages.";

	commands.insert(CMD_PARAMS, CommandSpecifier(CMD_PARAMS), "params");
	commands[CMD_PARAMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>params</b> <i>(meshname)</i>";
	commands[CMD_PARAMS].descr = "[tc0,0.5,0.5,1/tc]List all material parameters. If meshname not given use the focused mesh.";
	commands[CMD_PARAMS].limits = { { Any(), Any() } };

	commands.insert(CMD_SETPARAM, CommandSpecifier(CMD_SETPARAM), "setparam");
	commands[CMD_SETPARAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparam</b> <i>(meshname) paramname (value)</i>";
	commands[CMD_SETPARAM].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter to given value. If meshname not given use the focused mesh.";
	commands[CMD_SETPARAM].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - return value of named parameter in named mesh (focused mesh if not specified).";
	commands[CMD_SETPARAM].limits = { { Any(), Any() }, { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETKTENS, CommandSpecifier(CMD_SETKTENS), "setktens");
	commands[CMD_SETKTENS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setktens</b> <i>term1 (term2 ...)</i>";
	commands[CMD_SETKTENS].descr = "[tc0,0.5,0.5,1/tc]Set the anisotropy energy terms for tensorial anisotropy. Each term must be of the form dxn1yn2zn3, where d is a number, n1, n2, n3 are integers >= 0. x, y, z are string literals.";

	commands.insert(CMD_PARAMSTEMP, CommandSpecifier(CMD_PARAMSTEMP), "paramstemp");
	commands[CMD_PARAMSTEMP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>paramstemp</b> <i>(meshname)</i>";
	commands[CMD_PARAMSTEMP].descr = "[tc0,0.5,0.5,1/tc]List all material parameters temperature dependence. If meshname not given use the focused mesh.";
	commands[CMD_PARAMSTEMP].limits = { { Any(), Any() } };

	commands.insert(CMD_CLEARPARAMSTEMP, CommandSpecifier(CMD_CLEARPARAMSTEMP), "clearparamstemp");
	commands[CMD_CLEARPARAMSTEMP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearparamstemp</b> <i>(meshname) (paramname)</i>";
	commands[CMD_CLEARPARAMSTEMP].descr = "[tc0,0.5,0.5,1/tc]Clear material parameter temperature dependence in given mesh. If meshname not given clear temperature dependences in all meshes. If paramname not given clear all parameters temperature dependences.";
	commands[CMD_CLEARPARAMSTEMP].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETPARAMTEMPEQUATION, CommandSpecifier(CMD_SETPARAMTEMPEQUATION), "setparamtempequation");
	commands[CMD_SETPARAMTEMPEQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparamtempequation</b> <i>(meshname) paramname text_equation</i>";
	commands[CMD_SETPARAMTEMPEQUATION].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter temperature dependence equation for the named mesh (all applicable meshes if not given).";
	commands[CMD_SETPARAMTEMPEQUATION].limits = { { Any(), Any() } , { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETPARAMTEMPARRAY, CommandSpecifier(CMD_SETPARAMTEMPARRAY), "setparamtemparray");
	commands[CMD_SETPARAMTEMPARRAY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparamtemparray</b> <i>(meshname) paramname filename</i>";
	commands[CMD_SETPARAMTEMPARRAY].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter temperature dependence using an array in the given mesh (all applicable meshes if not given). This must contain temperature values and scaling coefficients. Load directly from a file (tab spaced).";
	commands[CMD_SETPARAMTEMPARRAY].limits = { { Any(), Any() } , { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_COPYPARAMS, CommandSpecifier(CMD_COPYPARAMS), "copyparams");
	commands[CMD_COPYPARAMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>copyparams</b> <i>meshname_from meshname_to (...)</i>";
	commands[CMD_COPYPARAMS].descr = "[tc0,0.5,0.5,1/tc]Copy all mesh parameters from first mesh to all other meshes given - all meshes must be of same type.";

	commands.insert(CMD_COPYMESHDATA, CommandSpecifier(CMD_COPYMESHDATA), "copymeshdata");
	commands[CMD_COPYMESHDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>copymeshdata</b> <i>meshname_from meshname_to (...)</i>";
	commands[CMD_COPYMESHDATA].descr = "[tc0,0.5,0.5,1/tc]Copy all primary mesh data (e.g. magnetization values and shape) from first mesh to all other meshes given - all meshes must be of same type.";

	commands.insert(CMD_PARAMSVAR, CommandSpecifier(CMD_PARAMSVAR), "paramsvar");
	commands[CMD_PARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>paramsvar</b> <i>(meshname)</i>";
	commands[CMD_PARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]List all material parameters spatial variation. If meshname not given use the focused mesh.";
	commands[CMD_PARAMSVAR].limits = { { Any(), Any() } };

	commands.insert(CMD_SETDISPLAYEDPARAMSVAR, CommandSpecifier(CMD_SETDISPLAYEDPARAMSVAR), "setdisplayedparamsvar");
	commands[CMD_SETDISPLAYEDPARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdisplayedparamsvar</b> <i>(meshname) paramname</i>";
	commands[CMD_SETDISPLAYEDPARAMSVAR].limits = { { Any(), Any() }, { int(0), Any() } };
	commands[CMD_SETDISPLAYEDPARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]Set param to display for given mesh (focused mesh if not specified) when ParamVar display is enabled (to show spatial variation if any).";

	commands.insert(CMD_CLEARPARAMSVAR, CommandSpecifier(CMD_CLEARPARAMSVAR), "clearparamsvar");
	commands[CMD_CLEARPARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearparamsvar</b> <i>(meshname) (paramname)</i>";
	commands[CMD_CLEARPARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]Clear parameter spatial dependence in given mesh. If meshname not given clear spatial dependence in all meshes. If paramname not given clear all parameters spatial dependence.";
	commands[CMD_CLEARPARAMSVAR].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETPARAMVAR, CommandSpecifier(CMD_SETPARAMVAR), "setparamvar");
	commands[CMD_SETPARAMVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparamvar</b> <i>(meshname) paramname generatorname (arguments...)</i>";
	commands[CMD_SETPARAMVAR].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter spatial dependence for the named mesh (all applicable meshes if not specified) using the given generator (including any required arguments for the generator - if not given, default values are used).";
	commands[CMD_SETPARAMVAR].limits = { { Any(), Any() }, { Any(), Any() }, { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_GETMESHTYPE, CommandSpecifier(CMD_GETMESHTYPE), "getmeshtype");
	commands[CMD_GETMESHTYPE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>getmeshtype</b> <i>meshname</i>";
	commands[CMD_GETMESHTYPE].descr = "[tc0,0.5,0.5,1/tc]Show mesh type of named mesh.";
	commands[CMD_GETMESHTYPE].limits = { { Any(), Any() } };
	commands[CMD_GETMESHTYPE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshtype</i>.";

	commands.insert(CMD_SAVESIM, CommandSpecifier(CMD_SAVESIM), "savesim");
	commands[CMD_SAVESIM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savesim</b> <i>(directory/)filename</i>";
	commands[CMD_SAVESIM].descr = "[tc0,0.5,0.5,1/tc]Save simulation with given name. If no name given, the last saved/loaded file name will be used.";
	commands[CMD_SAVESIM].limits = { { Any(), Any() } };

	commands.insert(CMD_LOADSIM, CommandSpecifier(CMD_LOADSIM), "loadsim");
	commands[CMD_LOADSIM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadsim</b> <i>(directory/)filename</i>";
	commands[CMD_LOADSIM].descr = "[tc0,0.5,0.5,1/tc]Load simulation with given name.";
	commands[CMD_LOADSIM].limits = { { Any(), Any() } };

	commands.insert(CMD_DEFAULT, CommandSpecifier(CMD_DEFAULT), "default");
	commands[CMD_DEFAULT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>default</b>";
	commands[CMD_DEFAULT].descr = "[tc0,0.5,0.5,1/tc]Reset program to default state.";

	commands.insert(CMD_DISPLAY, CommandSpecifier(CMD_DISPLAY), "display");
	commands[CMD_DISPLAY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>display</b> <i>(meshname) name</i>";
	commands[CMD_DISPLAY].descr = "[tc0,0.5,0.5,1/tc]Change quantity to display for given mesh (focused mesh if not specified).";
	commands[CMD_DISPLAY].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_DISPLAYMODULE, CommandSpecifier(CMD_DISPLAYMODULE), "displaymodule");
	commands[CMD_DISPLAYMODULE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaymodule</b> <i>(meshname) modulename</i>";
	commands[CMD_DISPLAYMODULE].descr = "[tc0,0.5,0.5,1/tc]Select module for effective field and energy density display (focused mesh if not specified). If modulename is none then display total effective field.";
	commands[CMD_DISPLAYMODULE].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_ROTCAMABOUTORIGIN, CommandSpecifier(CMD_ROTCAMABOUTORIGIN), "rotcamaboutorigin");
	commands[CMD_ROTCAMABOUTORIGIN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>rotcamaboutorigin</b> <i>dAzim dPolar</i>";
	commands[CMD_ROTCAMABOUTORIGIN].descr = "[tc0,0.5,0.5,1/tc]Rotate camera about origin using a delta azimuthal and delta polar angle change. This is the same as middle mouse hold and drag on mesh viewer.";
	commands[CMD_ROTCAMABOUTORIGIN].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_ROTCAMABOUTAXIS, CommandSpecifier(CMD_ROTCAMABOUTAXIS), "rotcamaboutaxis");
	commands[CMD_ROTCAMABOUTAXIS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>rotcamaboutaxis</b> <i>dAngle</i>";
	commands[CMD_ROTCAMABOUTAXIS].descr = "[tc0,0.5,0.5,1/tc]Rotate camera about its axis using an angle change. This is the same as right mouse hold and left-right drag on mesh viewer.";
	commands[CMD_ROTCAMABOUTAXIS].limits = { { Any(), Any() } };

	commands.insert(CMD_ADJUSTCAMDISTANCE, CommandSpecifier(CMD_ADJUSTCAMDISTANCE), "adjustcamdistance");
	commands[CMD_ADJUSTCAMDISTANCE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>adjustcamdistance</b> <i>dZ</i>";
	commands[CMD_ADJUSTCAMDISTANCE].descr = "[tc0,0.5,0.5,1/tc]Adjust camera distance from origin using dZ. This is the same as right mouse hold and up-down drag on mesh viewer.";
	commands[CMD_ADJUSTCAMDISTANCE].limits = { { Any(), Any() } };

	commands.insert(CMD_SHIFTCAMORIGIN, CommandSpecifier(CMD_SHIFTCAMORIGIN), "shiftcamorigin");
	commands[CMD_SHIFTCAMORIGIN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shiftcamorigin</b> <i>dX dY</i>";
	commands[CMD_SHIFTCAMORIGIN].descr = "[tc0,0.5,0.5,1/tc]Shift camera origin using dX dY. This is the same as left mouse hold and drag on mesh viewer.";
	commands[CMD_SHIFTCAMORIGIN].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_DISPLAYDETAILLEVEL, CommandSpecifier(CMD_DISPLAYDETAILLEVEL), "displaydetail");
	commands[CMD_DISPLAYDETAILLEVEL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaydetail</b> <i>size</i>";
	commands[CMD_DISPLAYDETAILLEVEL].descr = "[tc0,0.5,0.5,1/tc]Change displayed detail level to given displayed cell size (m).";
	commands[CMD_DISPLAYDETAILLEVEL].limits = { { double(MINMESHSPACE / 2), Any() } };
	commands[CMD_DISPLAYDETAILLEVEL].unit = "m";
	commands[CMD_DISPLAYDETAILLEVEL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>size</i>";

	commands.insert(CMD_DISPLAYRENDERTHRESH, CommandSpecifier(CMD_DISPLAYRENDERTHRESH), "displayrenderthresholds");
	commands[CMD_DISPLAYRENDERTHRESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displayrenderthresholds</b> <i>thresh1 thresh2 thresh3</i>";
	commands[CMD_DISPLAYRENDERTHRESH].descr = "[tc0,0.5,0.5,1/tc]Set display render thresholds of number of displayed cells for faster rendering as: 1) switch to using simpler elements (vectors only), 2) don't display surrounded cells, 3) display on checkerboard pattern even if not surrounded (vectors only). Set to zero to disable.";
	commands[CMD_DISPLAYRENDERTHRESH].limits = { { INT3(), Any() } };
	commands[CMD_DISPLAYRENDERTHRESH].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>thresh1 thresh2 thresh3</i>";

	commands.insert(CMD_DISPLAYBACKGROUND, CommandSpecifier(CMD_DISPLAYBACKGROUND), "displaybackground");
	commands[CMD_DISPLAYBACKGROUND].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaybackground</b> <i>(meshname) name</i>";
	commands[CMD_DISPLAYBACKGROUND].descr = "[tc0,0.5,0.5,1/tc]Change background quantity to display for given mesh (focused mesh if not specified).";
	commands[CMD_DISPLAYBACKGROUND].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_VECREP, CommandSpecifier(CMD_VECREP), "vecrep");
	commands[CMD_VECREP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>vecrep</b> <i>(meshname) vecreptype</i>";
	commands[CMD_VECREP].limits = { { Any(), Any() }, { int(0), int(VEC3REP_NUMOPTIONS) } };
	commands[CMD_VECREP].descr = "[tc0,0.5,0.5,1/tc]Set representation type for vectorial quantities in named mesh, or supermesh (focused mesh if not specified). vecreptype = 0 (full), vecreptype = 1 (x component), vecreptype = 2 (y component), vecreptype = 3 (z component), vecreptype = 4 (direction only), vecreptype = 5 (magnitude only).";

	commands.insert(CMD_SAVEMESHIMAGE, CommandSpecifier(CMD_SAVEMESHIMAGE), "savemeshimage");
	commands[CMD_SAVEMESHIMAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savemeshimage</b> <i>((directory/)filename)</i>";
	commands[CMD_SAVEMESHIMAGE].descr = "[tc0,0.5,0.5,1/tc]Save currently displayed mesh image to given file (as .png). If directory not specified then default directory is used. If filename not specified then default image save file name is used.";

	commands.insert(CMD_SAVEIMAGE, CommandSpecifier(CMD_SAVEIMAGE), "saveimage");
	commands[CMD_SAVEIMAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveimage</b> <i>((directory/)filename)</i>";
	commands[CMD_SAVEIMAGE].descr = "[tc0,0.5,0.5,1/tc]Save currently displayed mesh image, with a centered view of displayed quantity only, to given file (as .png). If directory not specified then default directory is used. If filename not specified then default image save file name is used.";

	commands.insert(CMD_MAKEVIDEO, CommandSpecifier(CMD_MAKEVIDEO), "makevideo");
	commands[CMD_MAKEVIDEO].usage = "[tc0,0.5,0,1/tc]USAGE : <b>makevideo</b> <i>(directory/)filebase fps quality</i>";
	commands[CMD_MAKEVIDEO].limits = { { Any(), Any() }, { double(1), double(120) }, { int(0), int(5) } };
	commands[CMD_MAKEVIDEO].descr = "[tc0,0.5,0.5,1/tc]Make a video from .png files sharing the common filebase name. Make video at given fps and quality (0 to 5 worst to best).";
	
	commands.insert(CMD_IMAGECROPPING, CommandSpecifier(CMD_IMAGECROPPING), "imagecropping");
	commands[CMD_IMAGECROPPING].usage = "[tc0,0.5,0,1/tc]USAGE : <b>imagecropping</b> <i>left bottom right top</i>";
	commands[CMD_IMAGECROPPING].limits = { { double(0), double(1) }, { double(0), double(1) }, { double(0), double(1) }, { double(0), double(1) } };
	commands[CMD_IMAGECROPPING].descr = "[tc0,0.5,0.5,1/tc]Set cropping of saved mesh images using normalized left, bottom, right, top values: 0, 0 point is left, bottom of mesh window and 1, 1 is right, top of mesh window.";

	commands.insert(CMD_DISPLAYTRANSPARENCY, CommandSpecifier(CMD_DISPLAYTRANSPARENCY), "displaytransparency");
	commands[CMD_DISPLAYTRANSPARENCY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaytransparency</b> <i>foreground background</i>";
	commands[CMD_DISPLAYTRANSPARENCY].limits = { { double(0), double(1) }, { double(0), double(1) } };
	commands[CMD_DISPLAYTRANSPARENCY].descr = "[tc0,0.5,0.5,1/tc]Set alpha transparency for display. Values range from 0 (fully transparent) to 1 (opaque). This is applicable in dual display mode when we have a background and foreground for the same mesh.";
	
	commands.insert(CMD_DISPLAYTHRESHOLDS, CommandSpecifier(CMD_DISPLAYTHRESHOLDS), "displaythresholds");
	commands[CMD_DISPLAYTHRESHOLDS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaythresholds</b> <i>minimum maximum</i>";
	commands[CMD_DISPLAYTHRESHOLDS].descr = "[tc0,0.5,0.5,1/tc]Set thresholds for foreground mesh display : magnitude values outside this range are not rendered. If both set to 0 then thresholds are ignored.";

	commands.insert(CMD_DISPLAYTHRESHOLDTRIGGER, CommandSpecifier(CMD_DISPLAYTHRESHOLDTRIGGER), "displaythresholdtrigger");
	commands[CMD_DISPLAYTHRESHOLDTRIGGER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaythresholdtrigger</b> <i>trigtype</i>";
	commands[CMD_DISPLAYTHRESHOLDTRIGGER].limits = { { int(0), int(VEC3REP_NUMOPTIONS) } };
	commands[CMD_DISPLAYTHRESHOLDTRIGGER].descr = "[tc0,0.5,0.5,1/tc]For vector quantities, set component to trigger thresholds on. trigtype = 1 (x component), trigtype = 2 (y component), trigtype = 3 (z component), trigtype = 5 (magnitude only)";

	commands.insert(CMD_MOVINGMESH, CommandSpecifier(CMD_MOVINGMESH), "movingmesh");
	commands[CMD_MOVINGMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>movingmesh</b> <i>status_or_meshname</i>";
	commands[CMD_MOVINGMESH].descr = "[tc0,0.5,0.5,1/tc]Set/unset trigger for movingmesh algorithm. If status_or_meshname = 0 then turn off, if status_or_meshname = 1 then turn on with trigger set on first ferromagnetic mesh, else status_or_meshname should specify the mesh name to use as trigger.";

	commands.insert(CMD_CLEARMOVINGMESH, CommandSpecifier(CMD_CLEARMOVINGMESH), "clearmovingmesh");
	commands[CMD_CLEARMOVINGMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearmovingmesh</b>";
	commands[CMD_CLEARMOVINGMESH].descr = "[tc0,0.5,0.5,1/tc]Clear moving mesh settings made by a prepare command.";

	commands.insert(CMD_MOVINGMESHASYM, CommandSpecifier(CMD_MOVINGMESHASYM), "movingmeshasym");
	commands[CMD_MOVINGMESHASYM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>movingmeshasym</b> <i>status</i>";
	commands[CMD_MOVINGMESHASYM].descr = "[tc0,0.5,0.5,1/tc]Change symmetry type for moving mesh algorithm: 1 for antisymmetric (domain walls), 0 for symmetric (skyrmions).";
	commands[CMD_MOVINGMESHASYM].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_MOVINGMESHTHRESH, CommandSpecifier(CMD_MOVINGMESHTHRESH), "movingmeshthresh");
	commands[CMD_MOVINGMESHTHRESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>movingmeshthresh</b> <i>value</i>";
	commands[CMD_MOVINGMESHTHRESH].limits = { { double(0), double(1) } };
	commands[CMD_MOVINGMESHTHRESH].descr = "[tc0,0.5,0.5,1/tc]Set threshold used to trigger a mesh shift for moving mesh algorithm - normalised value between 0 and 1.";
	commands[CMD_MOVINGMESHTHRESH].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>threshold</i>";

	commands.insert(CMD_PREPAREMOVINGMESH, CommandSpecifier(CMD_PREPAREMOVINGMESH), "preparemovingmesh");
	commands[CMD_PREPAREMOVINGMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>preparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or focused mesh) for moving transverse (or vortex) domain wall simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";
	commands[CMD_PREPAREMOVINGMESH].limits = { { Any(), Any() } };

	commands.insert(CMD_PREPAREMOVINGBLOCHMESH, CommandSpecifier(CMD_PREPAREMOVINGBLOCHMESH), "blochpreparemovingmesh");
	commands[CMD_PREPAREMOVINGBLOCHMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>blochpreparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGBLOCHMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or focused mesh) for moving Bloch domain wall simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";
	commands[CMD_PREPAREMOVINGBLOCHMESH].limits = { { Any(), Any() } };

	commands.insert(CMD_PREPAREMOVINGNEELMESH, CommandSpecifier(CMD_PREPAREMOVINGNEELMESH), "neelpreparemovingmesh");
	commands[CMD_PREPAREMOVINGNEELMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>neelpreparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGNEELMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or focused mesh) for moving Neel domain wall simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";
	commands[CMD_PREPAREMOVINGNEELMESH].limits = { { Any(), Any() } };

	commands.insert(CMD_PREPAREMOVINGSKYRMIONMESH, CommandSpecifier(CMD_PREPAREMOVINGSKYRMIONMESH), "skyrmionpreparemovingmesh");
	commands[CMD_PREPAREMOVINGSKYRMIONMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyrmionpreparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGSKYRMIONMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or focused mesh) for moving skyrmion simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";
	commands[CMD_PREPAREMOVINGSKYRMIONMESH].limits = { { Any(), Any() } };

	commands.insert(CMD_COUPLETODIPOLES, CommandSpecifier(CMD_COUPLETODIPOLES), "coupletodipoles");
	commands[CMD_COUPLETODIPOLES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>coupletodipoles</b> <i>status</i>";
	commands[CMD_COUPLETODIPOLES].descr = "[tc0,0.5,0.5,1/tc]Set/unset coupling to dipoles : if magnetic meshes touch a dipole mesh then interface magnetic cells are coupled to the dipole magnetization direction.";

	commands.insert(CMD_EXCHANGECOUPLEDMESHES, CommandSpecifier(CMD_EXCHANGECOUPLEDMESHES), "exchangecoupledmeshes");
	commands[CMD_EXCHANGECOUPLEDMESHES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>exchangecoupledmeshes</b> <i>(meshname) status</i>";
	commands[CMD_EXCHANGECOUPLEDMESHES].descr = "[tc0,0.5,0.5,1/tc]Set/unset direct exchange coupling to neighboring meshes : if neighboring ferromagnetic meshes touch the named mesh (focused mesh if not specified) then interface magnetic cells are direct exchange coupled to them.";
	commands[CMD_EXCHANGECOUPLEDMESHES].limits = { { Any(), Any() }, { int(0), int(1) } };
	commands[CMD_EXCHANGECOUPLEDMESHES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_ADDELECTRODE, CommandSpecifier(CMD_ADDELECTRODE), "addelectrode");
	commands[CMD_ADDELECTRODE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addelectrode</b> <i>electrode_rect</i>";
	commands[CMD_ADDELECTRODE].unit = "m";
	commands[CMD_ADDELECTRODE].descr = "[tc0,0.5,0.5,1/tc]Add an electrode in given rectangle (m).";

	commands.insert(CMD_DELELECTRODE, CommandSpecifier(CMD_DELELECTRODE), "delelectrode");
	commands[CMD_DELELECTRODE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delelectrode</b> <i>index</i>";
	commands[CMD_DELELECTRODE].limits = { { int(0), Any() } };
	commands[CMD_DELELECTRODE].descr = "[tc0,0.5,0.5,1/tc]Delete electrode with given index.";

	commands.insert(CMD_CLEARELECTRODES, CommandSpecifier(CMD_CLEARELECTRODES), "clearelectrodes");
	commands[CMD_CLEARELECTRODES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearelectrodes</b>";
	commands[CMD_CLEARELECTRODES].descr = "[tc0,0.5,0.5,1/tc]Delete all currently set electrodes.";

	commands.insert(CMD_ELECTRODES, CommandSpecifier(CMD_ELECTRODES), "electrodes");
	commands[CMD_ELECTRODES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>electrodes</b>";
	commands[CMD_ELECTRODES].descr = "[tc0,0.5,0.5,1/tc]Show currently configured electrodes.";

	commands.insert(CMD_SETDEFAULTELECTRODES, CommandSpecifier(CMD_SETDEFAULTELECTRODES), "setdefaultelectrodes");
	commands[CMD_SETDEFAULTELECTRODES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdefaultelectrodes</b> <i>(sides)</i>";
	commands[CMD_SETDEFAULTELECTRODES].limits = { { Any(), Any() } };
	commands[CMD_SETDEFAULTELECTRODES].descr = "[tc0,0.5,0.5,1/tc]Set electrodes at the x-axis ends of the given mesh, both set at 0V. If sides specified (literal x, y, or z) then set at respective axis ends instead of default x-axis ends. Set the left-side electrode as the ground. Delete all other electrodes.";

	commands.insert(CMD_SETELECTRODERECT, CommandSpecifier(CMD_SETELECTRODERECT), "setelectroderect");
	commands[CMD_SETELECTRODERECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setelectroderect</b> <i>electrode_index electrode_rect</i>";
	commands[CMD_SETELECTRODERECT].descr = "[tc0,0.5,0.5,1/tc]Edit rectangle (m) for electrode with given index.";
	commands[CMD_SETELECTRODERECT].unit = "m";

	commands.insert(CMD_SETELECTRODEPOTENTIAL, CommandSpecifier(CMD_SETELECTRODEPOTENTIAL), "setelectrodepotential");
	commands[CMD_SETELECTRODEPOTENTIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setelectrodepotential</b> <i>electrode_index potential</i>";
	commands[CMD_SETELECTRODEPOTENTIAL].descr = "[tc0,0.5,0.5,1/tc]Set potential on electrode with given index.";
	commands[CMD_SETELECTRODEPOTENTIAL].unit = "V";
	commands[CMD_SETELECTRODEPOTENTIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>potential</i>";

	commands.insert(CMD_DESIGNATEGROUND, CommandSpecifier(CMD_DESIGNATEGROUND), "designateground");
	commands[CMD_DESIGNATEGROUND].usage = "[tc0,0.5,0,1/tc]USAGE : <b>designateground</b> <i>electrode_index</i>";
	commands[CMD_DESIGNATEGROUND].descr = "[tc0,0.5,0.5,1/tc]Change ground designation for electrode with given index.";

	commands.insert(CMD_SETPOTENTIAL, CommandSpecifier(CMD_SETPOTENTIAL), "setpotential");
	commands[CMD_SETPOTENTIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setpotential</b> <i>potential</i>";
	commands[CMD_SETPOTENTIAL].descr = "[tc0,0.5,0.5,1/tc]Set a symmetric potential drop : -potential/2 for ground electrode, +potential/2 on all other electrodes.";
	commands[CMD_SETPOTENTIAL].unit = "V";
	commands[CMD_SETPOTENTIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>potential</i>";

	commands.insert(CMD_SETCURRENT, CommandSpecifier(CMD_SETCURRENT), "setcurrent");
	commands[CMD_SETCURRENT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setcurrent</b> <i>current</i>";
	commands[CMD_SETCURRENT].descr = "[tc0,0.5,0.5,1/tc]Set a constant current source with given value. The potential will be adjusted to keep this constant current.";
	commands[CMD_SETCURRENT].unit = "A";
	commands[CMD_SETCURRENT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>current</i>";

	commands.insert(CMD_SETCURRENTDENSITY, CommandSpecifier(CMD_SETCURRENTDENSITY), "setcurrentdensity");
	commands[CMD_SETCURRENTDENSITY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setcurrentdensity</b> <i>(meshname) Jx Jy Jz</i>";
	commands[CMD_SETCURRENTDENSITY].descr = "[tc0,0.5,0.5,1/tc]Set a constant current density vector with given components in given mesh (focused mesh if not specified). Must have transport module added, but transport solver iteration will be disabled.";
	commands[CMD_SETCURRENTDENSITY].unit = "A/m2";
	commands[CMD_SETCURRENTDENSITY].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_OPENPOTENTIAL, CommandSpecifier(CMD_OPENPOTENTIAL), "openpotentialresistance");
	commands[CMD_OPENPOTENTIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>openpotentialresistance</b> <i>resistance</i>";
	commands[CMD_OPENPOTENTIAL].descr = "[tc0,0.5,0.5,1/tc]Set open potential mode for electrodes, by setting a resistance value - set zero to disable open potential mode. In this mode the electrode potential is determined during the simulation, e.g. if thermoelectric effect is enabled, then setting electrodes in open potential mode allows a thermoelectric current to flow through the electrodes, with potential determined as thermoelectric current generated times resistance set.";
	commands[CMD_OPENPOTENTIAL].unit = "Ohm";
	commands[CMD_OPENPOTENTIAL].limits = { { double(0.0), Any() } };
	commands[CMD_OPENPOTENTIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>resistance</i>";

	commands.insert(CMD_TSOLVERCONFIG, CommandSpecifier(CMD_TSOLVERCONFIG), "tsolverconfig");
	commands[CMD_TSOLVERCONFIG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tsolverconfig</b> <i>convergence_error (iters_timeout)</i>";
	commands[CMD_TSOLVERCONFIG].limits = { { double(0.0), double(1.0) }, { int(1), Any() } };
	commands[CMD_TSOLVERCONFIG].descr = "[tc0,0.5,0.5,1/tc]Set transport solver convergence error and iterations for timeout (if given, else use default).";
	commands[CMD_TSOLVERCONFIG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>convergence_error iters_timeout</i>";

	commands.insert(CMD_SSOLVERCONFIG, CommandSpecifier(CMD_SSOLVERCONFIG), "ssolverconfig");
	commands[CMD_SSOLVERCONFIG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ssolverconfig</b> <i>s_convergence_error (s_iters_timeout)</i>";
	commands[CMD_SSOLVERCONFIG].limits = { { double(0.0), double(1.0) }, { int(1), Any() } };
	commands[CMD_SSOLVERCONFIG].descr = "[tc0,0.5,0.5,1/tc]Set spin-transport solver convergence error and iterations for timeout (if given, else use default).";
	commands[CMD_SSOLVERCONFIG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>s_convergence_error s_iters_timeout</i>";
	
	commands.insert(CMD_SETSORDAMPING, CommandSpecifier(CMD_SETSORDAMPING), "setsordamping");
	commands[CMD_SETSORDAMPING].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setsordamping</b> <i>damping_v damping_s</i>";
	commands[CMD_SETSORDAMPING].limits = { { DBL2(MINSORDAMPING), DBL2(MAXSORDAMPING) } };
	commands[CMD_SETSORDAMPING].descr = "[tc0,0.5,0.5,1/tc]Set fixed damping values for SOR algorithm used to solve the Poisson equation for V (electrical potential) and S (spin accumulation) respectively.";
	commands[CMD_SETSORDAMPING].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>damping_v damping_s</i>";

	commands.insert(CMD_STATICTRANSPORTSOLVER, CommandSpecifier(CMD_STATICTRANSPORTSOLVER), "statictransportsolver");
	commands[CMD_STATICTRANSPORTSOLVER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>statictransportsolver</b> <i>status</i>";
	commands[CMD_STATICTRANSPORTSOLVER].limits = { { int(0), int(1) } };
	commands[CMD_STATICTRANSPORTSOLVER].descr = "[tc0,0.5,0.5,1/tc]If static transport solver is set, the transport solver is only iterated at the end of a stage or step. You should set a high iterations timeout if using this mode.";
	commands[CMD_STATICTRANSPORTSOLVER].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_DISABLETRANSPORTSOLVER, CommandSpecifier(CMD_DISABLETRANSPORTSOLVER), "disabletransportsolver");
	commands[CMD_DISABLETRANSPORTSOLVER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>disabletransportsolver</b> <i>status</i>";
	commands[CMD_DISABLETRANSPORTSOLVER].limits = { { int(0), int(1) } };
	commands[CMD_DISABLETRANSPORTSOLVER].descr = "[tc0,0.5,0.5,1/tc]Disable iteration of transport solver so any set current density remains constant.";
	commands[CMD_DISABLETRANSPORTSOLVER].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_TMRTYPE, CommandSpecifier(CMD_TMRTYPE), "tmrtype");
	commands[CMD_TMRTYPE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tmrtype</b> <i>(meshname) setting</i>";
	commands[CMD_TMRTYPE].limits = { {Any(), Any()}, { int(0), int(TMR_NUMOPTIONS) } };
	commands[CMD_TMRTYPE].descr = "[tc0,0.5,0.5,1/tc]Set formula used to calculate TMR angle dependence in named mesh (all meshes if not given). setting = 0: RA = (RAp + (RAap - RAp) * (1 - cos(theta))/2). setting = 1 (Slonczewski form, default): RA = 2*RAp/[(1 + RAp/RAap) + (1 - RAp/RAap)*cos(theta)]";
	commands[CMD_TMRTYPE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>setting</i>";

	commands.insert(CMD_RAPBIAS_EQUATION, CommandSpecifier(CMD_RAPBIAS_EQUATION), "rapbiasequation");
	commands[CMD_RAPBIAS_EQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>rapbiasequation</b> <i>meshname text_equation</i>";
	commands[CMD_RAPBIAS_EQUATION].limits = { {Any(), Any()}, {Any(), Any()} };
	commands[CMD_RAPBIAS_EQUATION].descr = "[tc0,0.5,0.5,1/tc]Set equation for bias dependence of TMR parallel RA value in given mesh (must be insulator with tmr module added) - the text equation variables are RAp (parallel TMR RA), RAap (antiparallel TMR RA), and V (bias across insulator mesh). Call with empty equation string to clear any currently set equation.";

	commands.insert(CMD_RAAPBIAS_EQUATION, CommandSpecifier(CMD_RAAPBIAS_EQUATION), "raapbiasequation");
	commands[CMD_RAAPBIAS_EQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>raapbiasequation</b> <i>meshname text_equation</i>";
	commands[CMD_RAAPBIAS_EQUATION].limits = { {Any(), Any()}, {Any(), Any()} };
	commands[CMD_RAAPBIAS_EQUATION].descr = "[tc0,0.5,0.5,1/tc]Set equation for bias dependence of TMR anti-parallel RA value in given mesh (must be insulator with tmr module added) - the text equation uses V for the bias value. Call with empty equation string to clear any currently set equation.";

	commands.insert(CMD_TAMREQUATION, CommandSpecifier(CMD_TAMREQUATION), "tamrequation");
	commands[CMD_TAMREQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tamrequation</b> <i>meshname text_equation</i>";
	commands[CMD_TAMREQUATION].limits = { {Any(), Any()}, {Any(), Any()} };
	commands[CMD_TAMREQUATION].descr = "[tc0,0.5,0.5,1/tc]Set conductivity equation for TAMR effect (must be a mesh with transport module added) - the text equation uses parameters TAMR, d1, d2, d3. Here TAMR = tamr / 100 is the ratio, where tamr is the material parameter set as a percentage value. di = eai.m, i = 1, 2, 3, where eai is the crystalline axis and m is the magnetization direction. Call with empty equation string to clear any currently set equation and use default conductivity formula: cond / elC = 1 / (1 + TAMR * (1 - d1*d1)), where elC is the base conductivity (material parameter). Note, the formula set always multiplies elC to obtain conductivity due to TAMR.";

	commands.insert(CMD_TEMPERATURE, CommandSpecifier(CMD_TEMPERATURE), "temperature");
	commands[CMD_TEMPERATURE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>temperature</b> <i>(meshname) value</i>";
	commands[CMD_TEMPERATURE].limits = { { Any(), Any() }, { double(0.0), double(MAX_TEMPERATURE) } };
	commands[CMD_TEMPERATURE].descr = "[tc0,0.5,0.5,1/tc]Set mesh base temperature (all meshes if meshname not given) and reset temperature. Also set ambient temperature if Heat module added. If the base temperature setting has a spatial dependence specified through cT, this command will take it into account but only if the Heat module is added. If you want the temperature to remain fixed you can still have the Heat module enabled but disable the heat equation by setting the heat dT to zero (setheatdt 0).";
	commands[CMD_TEMPERATURE].unit = "K";
	commands[CMD_TEMPERATURE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - temperature value for focused mesh.";

	commands.insert(CMD_SETHEATDT, CommandSpecifier(CMD_SETHEATDT), "setheatdt");
	commands[CMD_SETHEATDT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setheatdt</b> <i>value</i>";
	commands[CMD_SETHEATDT].limits = { { double(0), double(MAXTIMESTEP) } };
	commands[CMD_SETHEATDT].descr = "[tc0,0.5,0.5,1/tc]Set heat equation solver time step. Zero value disables solver.";
	commands[CMD_SETHEATDT].unit = "s";
	commands[CMD_SETHEATDT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - heat equation time step.";

	commands.insert(CMD_AMBIENTTEMPERATURE, CommandSpecifier(CMD_AMBIENTTEMPERATURE), "ambient");
	commands[CMD_AMBIENTTEMPERATURE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ambient</b> <i>(meshname) ambient_temperature</i>";
	commands[CMD_AMBIENTTEMPERATURE].limits = { { Any(), Any() }, { double(0.0), double(MAX_TEMPERATURE) } };
	commands[CMD_AMBIENTTEMPERATURE].descr = "[tc0,0.5,0.5,1/tc]Set mesh ambient temperature (all meshes if meshname not given) for Robin boundary conditions : flux normal = alpha * (T_boundary - T_ambient).";
	commands[CMD_AMBIENTTEMPERATURE].unit = "K";
	commands[CMD_AMBIENTTEMPERATURE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>ambient_temperature</i> - ambient temperature for named mesh.";

	commands.insert(CMD_ROBINALPHA, CommandSpecifier(CMD_ROBINALPHA), "robinalpha");
	commands[CMD_ROBINALPHA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>robinalpha</b> <i>(meshname) robin_alpha</i>";
	commands[CMD_ROBINALPHA].limits = { { Any(), Any() }, { double(0.0), double(1e10) } };
	commands[CMD_ROBINALPHA].descr = "[tc0,0.5,0.5,1/tc]Set alpha coefficient (all meshes if meshname not given) for Robin boundary conditions : flux normal = alpha * (T_boundary - T_ambient).";
	commands[CMD_ROBINALPHA].unit = "W/m2K";
	commands[CMD_ROBINALPHA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>robin_alpha</i> - Robin alpha value for named mesh.";

	commands.insert(CMD_INSULATINGSIDES, CommandSpecifier(CMD_INSULATINGSIDES), "insulatingside");
	commands[CMD_INSULATINGSIDES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>insulatingside</b> <i>(meshname) side_literal status</i>";
	commands[CMD_INSULATINGSIDES].limits = { { Any(), Any() }, { Any(), Any() }, { bool(0), bool(1) } };
	commands[CMD_INSULATINGSIDES].descr = "[tc0,0.5,0.5,1/tc]Set temperature insulation (Neumann boundary condition) for named mesh side (focused mesh if not specified). side_literal : x, -x, y, -y, z, -z.";
	commands[CMD_INSULATINGSIDES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status_x status_-x status_y status_-y status_z status_-z</i> - insulating sides status for named mesh.";
	
	commands.insert(CMD_CURIETEMPERATURE, CommandSpecifier(CMD_CURIETEMPERATURE), "curietemperature");
	commands[CMD_CURIETEMPERATURE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>curietemperature</b> <i>(meshname) curie_temperature</i>";
	commands[CMD_CURIETEMPERATURE].limits = { { Any(), Any() }, { double(0.0), double(MAX_TEMPERATURE) } };
	commands[CMD_CURIETEMPERATURE].descr = "[tc0,0.5,0.5,1/tc]Set Curie temperature for ferromagnetic mesh (all ferromagnetic meshes if not specified). This will set default temperature dependencies as: Ms = Ms0*me, A = Ah*me^2, D = D0*me^2, K = K0*me^3 (K1 and K2), damping = damping0*(1-T/3Tc) T < Tc, damping = damping0*2T/3Tc T >= Tc, susrel = dme/d(mu0Hext). Setting the Curie temperature to zero will disable temperature dependence for these parameters.";
	commands[CMD_CURIETEMPERATURE].unit = "K";
	commands[CMD_CURIETEMPERATURE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>curie_temperature</i> - Curie temperature for named mesh.";

	commands.insert(CMD_CURIETEMPERATUREMATERIAL, CommandSpecifier(CMD_CURIETEMPERATUREMATERIAL), "matcurietemperature");
	commands[CMD_CURIETEMPERATUREMATERIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>matcurietemperature</b> <i>(meshname) curie_temperature</i>";
	commands[CMD_CURIETEMPERATUREMATERIAL].limits = { { Any(), Any() }, { double(0.0), double(MAX_TEMPERATURE) } };
	commands[CMD_CURIETEMPERATUREMATERIAL].descr = "[tc0,0.5,0.5,1/tc]Set indicative material Curie temperature for ferromagnetic mesh (focused mesh if not specified). This is not used in calculations, but serves as an indicative value - set the actual Tc value with the <b>curietemperature</b> command.";
	commands[CMD_CURIETEMPERATUREMATERIAL].unit = "K";
	commands[CMD_CURIETEMPERATUREMATERIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>curie_temperature</i> - Indicative material Curie temperature for named mesh.";
	
	commands.insert(CMD_ATOMICMOMENT, CommandSpecifier(CMD_ATOMICMOMENT), "atomicmoment");
	commands[CMD_ATOMICMOMENT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>atomicmoment</b> <i>(meshname) ub_multiple</i>";
	commands[CMD_ATOMICMOMENT].limits = { { Any(), Any() }, { double(0.0), Any() } };
	commands[CMD_ATOMICMOMENT].descr = "[tc0,0.5,0.5,1/tc]Set atomic moment as a multiple of Bohr magnetons for given mesh (all ferromagnetic meshes if not specified). This affects the temperature dependence of 'me' (see curietemperature command). A non-zero value will result in me(T) being dependent on the applied field.";
	commands[CMD_ATOMICMOMENT].unit = "uB";
	commands[CMD_ATOMICMOMENT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>ub_multiple</i> - atomic moment multiple of Bohr magneton for named mesh.";

	commands.insert(CMD_TAU, CommandSpecifier(CMD_TAU), "tau");
	commands[CMD_TAU].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tau</b> <i>(meshname) tau_11 tau_22 (tau_12 tau_21)</i>";
	commands[CMD_TAU].limits = { { Any(), Any() }, { DBL2(),  Any() }, { DBL2(),  Any() } };
	commands[CMD_TAU].descr = "[tc0,0.5,0.5,1/tc]Set ratio of exchange parameters to critical temperature (Neel) for antiferromagnetic mesh (all antiferromagnetic meshes if not specified). tau_11 and tau_22 are the intra-lattice contributions, tau_12 and tau_21 are the inter-lattice contributions.";
	commands[CMD_TAU].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>tau_11 tau_22 tau_12 tau_21</i>";

	commands.insert(CMD_TMODEL, CommandSpecifier(CMD_TMODEL), "tmodel");
	commands[CMD_TMODEL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tmodel</b> <i>(meshname) num_temperatures</i>";
	commands[CMD_TMODEL].limits = { { Any(), Any() }, { int(1), Any() } };
	commands[CMD_TMODEL].descr = "[tc0,0.5,0.5,1/tc]Set temperature model (determined by number of temperatures) in given meshname (focused mesh if not specified). Note insulating meshes only allow a 1-temperature model.";

	commands.insert(CMD_RESETELSOLVER, CommandSpecifier(CMD_RESETELSOLVER), "resetelsolver");
	commands[CMD_RESETELSOLVER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>resetelsolver</b>";
	commands[CMD_RESETELSOLVER].descr = "[tc0,0.5,0.5,1/tc]Reset the elastodynamics solver to unstrained state.";

	commands.insert(CMD_STRAINEQUATION, CommandSpecifier(CMD_STRAINEQUATION), "strainequation");
	commands[CMD_STRAINEQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>strainequation</b> <i>(meshname) text_equation</i>";
	commands[CMD_STRAINEQUATION].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_STRAINEQUATION].descr = "[tc0,0.5,0.5,1/tc]Set equation for diagonal strain components, using a vector text equation (3 equations separated by commas) in given meshname (focused mesh if not specified). This will disable the elastodynamics solver in all meshes.";

	commands.insert(CMD_SHEARSTRAINEQUATION, CommandSpecifier(CMD_SHEARSTRAINEQUATION), "shearstrainequation");
	commands[CMD_SHEARSTRAINEQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shearstrainequation</b> <i>(meshname) text_equation</i>";
	commands[CMD_SHEARSTRAINEQUATION].limits = { { Any(), Any() }, { Any(), Any() } };
	commands[CMD_SHEARSTRAINEQUATION].descr = "[tc0,0.5,0.5,1/tc]Set equation for shear strain components, using a vector text equation (3 equations separated by commas) in given meshname (focused mesh if not specified). This will disable the elastodynamics solver in all meshes.";
	
	commands.insert(CMD_CLEARSTRAINEQUATIONS, CommandSpecifier(CMD_CLEARSTRAINEQUATIONS), "clearstrainequations");
	commands[CMD_CLEARSTRAINEQUATIONS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearstrainequations</b>";
	commands[CMD_CLEARSTRAINEQUATIONS].descr = "[tc0,0.5,0.5,1/tc]Clear strain equations from all meshes.";

	commands.insert(CMD_SETELDT, CommandSpecifier(CMD_SETELDT), "seteldt");
	commands[CMD_SETELDT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>seteldt</b> <i>value</i>";
	commands[CMD_SETELDT].limits = { { double(0), double(MAXTIMESTEP) } };
	commands[CMD_SETELDT].descr = "[tc0,0.5,0.5,1/tc]Set elastodynamics equation solver time step. Zero value disables solver.";
	commands[CMD_SETELDT].unit = "s";
	commands[CMD_SETELDT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - elastodynamics equation time step.";

	commands.insert(CMD_LINKDTELASTIC, CommandSpecifier(CMD_LINKDTELASTIC), "linkdtelastic");
	commands[CMD_LINKDTELASTIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>linkdtelastic</b> <i>flag</i>";
	commands[CMD_LINKDTELASTIC].limits = { { int(0), int(1) } };
	commands[CMD_LINKDTELASTIC].descr = "[tc0,0.5,0.5,1/tc]Links elastodynamics solver time-step to magnetic ODE time-step if set, else elastodynamics time-step is independently controlled.";

	commands.insert(CMD_SURFACEFIX, CommandSpecifier(CMD_SURFACEFIX), "surfacefix");
	commands[CMD_SURFACEFIX].usage = "[tc0,0.5,0,1/tc]USAGE : <b>surfacefix</b> <i>(meshname) rect_or_face</i>";
	commands[CMD_SURFACEFIX].unit = "m";
	commands[CMD_SURFACEFIX].descr = "[tc0,0.5,0.5,1/tc]Add a fixed surface for the elastodynamics solver. A plane rectangle may be specified in absolute coordinates, which should intersect at least one mesh face, which will then be set as fixed. Alternatively an entire mesh face may be specified using a string literal: -x, x, -y, y, -z, or z. In this case a named mesh will also be required (focused mesh if not given).";

	commands.insert(CMD_SURFACESTRESS, CommandSpecifier(CMD_SURFACESTRESS), "surfacestress");
	commands[CMD_SURFACESTRESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>surfacestress</b> <i>(meshname) rect_or_face vector_equation</i>";
	commands[CMD_SURFACESTRESS].unit = "m";
	commands[CMD_SURFACESTRESS].descr = "[tc0,0.5,0.5,1/tc]Add an external stress surface for the elastodynamics solver. A plane rectangle may be specified in absolute coordinates, which should intersect at least one mesh face, which will then be set as fixed. Alternatively an entire mesh face may be specified using a string literal: -x, x, -y, y, -z, or z. In this case a named mesh will also be required (focused mesh if not given). The external stress is specified using a vector equation, with variables x, y, t. Coordinates x, y are relative to specified rectangle, and t is the time.";

	commands.insert(CMD_DELSURFACEFIX, CommandSpecifier(CMD_DELSURFACEFIX), "delsurfacefix");
	commands[CMD_DELSURFACEFIX].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delsurfacefix</b> <i>index</i>";
	commands[CMD_DELSURFACEFIX].limits = { {int(-1), Any()} };
	commands[CMD_DELSURFACEFIX].descr = "[tc0,0.5,0.5,1/tc]Delete a fixed surface with given index (previously added with surfacefix). Use index -1 to delete all fixed surfaces.";

	commands.insert(CMD_DELSURFACESTRESS, CommandSpecifier(CMD_DELSURFACESTRESS), "delsurfacestress");
	commands[CMD_DELSURFACESTRESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delsurfacestress</b> <i>index</i>";
	commands[CMD_DELSURFACESTRESS].limits = { {int(-1), Any()}  };
	commands[CMD_DELSURFACESTRESS].descr = "[tc0,0.5,0.5,1/tc]Delete an external stress surface with given index (previously added with surfacestress). Use index -1 to delete all external stress surfaces.";

	commands.insert(CMD_EDITSURFACEFIX, CommandSpecifier(CMD_EDITSURFACEFIX), "editsurfacefix");
	commands[CMD_EDITSURFACEFIX].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editsurfacefix</b> <i>index rect</i>";
	commands[CMD_EDITSURFACEFIX].limits = { {int(0), Any()},  { Any(), Any() } };
	commands[CMD_EDITSURFACEFIX].unit = "m";
	commands[CMD_EDITSURFACEFIX].descr = "[tc0,0.5,0.5,1/tc]Edit rectangle of fixed surface with given index.";

	commands.insert(CMD_EDITSURFACESTRESS, CommandSpecifier(CMD_EDITSURFACESTRESS), "editsurfacestress");
	commands[CMD_EDITSURFACESTRESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editsurfacestress</b> <i>index rect</i>";
	commands[CMD_EDITSURFACESTRESS].limits = { {int(0), Any()},  { Any(), Any() } };
	commands[CMD_EDITSURFACESTRESS].unit = "m";
	commands[CMD_EDITSURFACESTRESS].descr = "[tc0,0.5,0.5,1/tc]Edit rectangle of stress surface with given index.";

	commands.insert(CMD_EDITSURFACESTRESSEQUATION, CommandSpecifier(CMD_EDITSURFACESTRESSEQUATION), "editsurfacestressequation");
	commands[CMD_EDITSURFACESTRESSEQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editsurfacestressequation</b> <i>index equation</i>";
	commands[CMD_EDITSURFACESTRESSEQUATION].limits = { {int(0), Any()},  { Any(), Any() } };
	commands[CMD_EDITSURFACESTRESSEQUATION].descr = "[tc0,0.5,0.5,1/tc]Edit equation of stress surface with given index.";

	commands.insert(CMD_STOCHASTIC, CommandSpecifier(CMD_STOCHASTIC), "stochastic");
	commands[CMD_STOCHASTIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>stochastic</b>";
	commands[CMD_STOCHASTIC].descr = "[tc0,0.5,0.5,1/tc]Shows stochasticity settings : stochastic cellsize for each mesh and related settings.";

	commands.insert(CMD_LINKSTOCHASTIC, CommandSpecifier(CMD_LINKSTOCHASTIC), "linkstochastic");
	commands[CMD_LINKSTOCHASTIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>linkstochastic</b> <i>(meshname) flag</i>";
	commands[CMD_LINKSTOCHASTIC].limits = { { Any(), Any() }, { int(0), int(1) } };
	commands[CMD_LINKSTOCHASTIC].descr = "[tc0,0.5,0.5,1/tc]Links stochastic cellsize to magnetic cellsize if flag set to 1 for given mesh, else stochastic cellsize is independently controlled. If meshname not given set for all meshes.";

	commands.insert(CMD_SETDTSTOCH, CommandSpecifier(CMD_SETDTSTOCH), "setdtstoch");
	commands[CMD_SETDTSTOCH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdtstoch</b> <i>value</i>";
	commands[CMD_SETDTSTOCH].limits = { { double(MINTIMESTEP), double(MAXTIMESTEP) } };
	commands[CMD_SETDTSTOCH].descr = "[tc0,0.5,0.5,1/tc]Set time step for stochastic field generation.";
	commands[CMD_SETDTSTOCH].unit = "s";
	commands[CMD_SETDTSTOCH].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>dTstoch</i>";

	commands.insert(CMD_LINKDTSTOCHASTIC, CommandSpecifier(CMD_LINKDTSTOCHASTIC), "linkdtstochastic");
	commands[CMD_LINKDTSTOCHASTIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>linkdtstochastic</b> <i>flag</i>";
	commands[CMD_LINKDTSTOCHASTIC].limits = { { int(0), int(1) } };
	commands[CMD_LINKDTSTOCHASTIC].descr = "[tc0,0.5,0.5,1/tc]Links stochastic time-step to ODE time-step if set, else stochastic time-step is independently controlled.";

	commands.insert(CMD_SETDTSPEEDUP, CommandSpecifier(CMD_SETDTSPEEDUP), "setdtspeedup");
	commands[CMD_SETDTSPEEDUP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdtspeedup</b> <i>value</i>";
	commands[CMD_SETDTSPEEDUP].limits = { { double(MINTIMESTEP), double(MAXTIMESTEP) } };
	commands[CMD_SETDTSPEEDUP].descr = "[tc0,0.5,0.5,1/tc]Set time step for evaluation speedup.";
	commands[CMD_SETDTSPEEDUP].unit = "s";
	commands[CMD_SETDTSPEEDUP].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>dTspeedup</i>";

	commands.insert(CMD_LINKDTSPEEDUP, CommandSpecifier(CMD_LINKDTSPEEDUP), "linkdtspeedup");
	commands[CMD_LINKDTSPEEDUP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>linkdtspeedup</b> <i>flag</i>";
	commands[CMD_LINKDTSPEEDUP].limits = { { int(0), int(1) } };
	commands[CMD_LINKDTSPEEDUP].descr = "[tc0,0.5,0.5,1/tc]Links speedup time-step to ODE time-step if set, else speedup time-step is independently controlled.";

	commands.insert(CMD_CUDA, CommandSpecifier(CMD_CUDA), "cuda");
	commands[CMD_CUDA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>cuda</b> <i>status</i>";
	commands[CMD_CUDA].descr = "[tc0,0.5,0.5,1/tc]Switch CUDA GPU computations on/off.";
	commands[CMD_CUDA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_MEMORY, CommandSpecifier(CMD_MEMORY), "memory");
	commands[CMD_MEMORY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>memory</b>";
	commands[CMD_MEMORY].descr = "[tc0,0.5,0.5,1/tc]Show CPU and GPU-addressable memory information (total and free).";

	commands.insert(CMD_SELCUDADEV, CommandSpecifier(CMD_SELCUDADEV), "selectcudadevice");
	commands[CMD_SELCUDADEV].usage = "[tc0,0.5,0,1/tc]USAGE : <b>selectcudadevice</b> <i>number</i>";
	commands[CMD_SELCUDADEV].limits = { { int(0), Any() } };
	commands[CMD_SELCUDADEV].descr = "[tc0,0.5,0.5,1/tc]Select CUDA device to use from available devices. The device CUDA Compute version must match Boris CUDA version. Devices numbered from 0 up, default selection at startup is device 0.";
	commands[CMD_SELCUDADEV].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>number</i>";

	commands.insert(CMD_OPENMANUAL, CommandSpecifier(CMD_OPENMANUAL), "manual");
	commands[CMD_OPENMANUAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>manual</b>";
	commands[CMD_OPENMANUAL].descr = "[tc0,0.5,0.5,1/tc]Opens Boris manual for current version.";

	commands.insert(CMD_REFINEROUGHNESS, CommandSpecifier(CMD_REFINEROUGHNESS), "refineroughness");
	commands[CMD_REFINEROUGHNESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>refineroughness</b> <i>(meshname) value</i>";
	commands[CMD_REFINEROUGHNESS].limits = { { Any(), Any() }, { INT3(1, 1, 1), Any() } };
	commands[CMD_REFINEROUGHNESS].descr = "[tc0,0.5,0.5,1/tc]Set roughness refinement cellsize divider in given mesh (focused mesh if not specified), i.e. cellsize used for roughness initialization is the magnetic cellsize divided by value (3 components, so divide component by component).";
	commands[CMD_REFINEROUGHNESS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - roughness refinement.";

	commands.insert(CMD_CLEARROUGHNESS, CommandSpecifier(CMD_CLEARROUGHNESS), "clearroughness");
	commands[CMD_CLEARROUGHNESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearroughness</b> <i>(meshname)</i>";
	commands[CMD_CLEARROUGHNESS].limits = { { Any(), Any() } };
	commands[CMD_CLEARROUGHNESS].descr = "[tc0,0.5,0.5,1/tc]Clear roughness in given mesh (focused mesh if not specified) by setting the fine shape same as the coarse M shape.";

	commands.insert(CMD_ROUGHENMESH, CommandSpecifier(CMD_ROUGHENMESH), "roughenmesh");
	commands[CMD_ROUGHENMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>roughenmesh</b> <i>(meshname) depth (side (seed))</i>";
	commands[CMD_ROUGHENMESH].limits = { { Any(), Any() }, { double(0), Any() }, { Any(), Any() }, { int(1), Any() } };
	commands[CMD_ROUGHENMESH].unit = "m";
	commands[CMD_ROUGHENMESH].descr = "[tc0,0.5,0.5,1/tc]Roughen given mesh (focused mesh if not specified) to given depth (m) on named side (use side = x, y, z, -x, -y, -z as literal, z by default). The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_SURFROUGHENJAGGED, CommandSpecifier(CMD_SURFROUGHENJAGGED), "surfroughenjagged");
	commands[CMD_SURFROUGHENJAGGED].usage = "[tc0,0.5,0,1/tc]USAGE : <b>surfroughenjagged</b> <i>(meshname) depth spacing (seed (sides))</i>";
	commands[CMD_SURFROUGHENJAGGED].limits = { { Any(), Any() }, { double(0), Any() }, { double(0), Any() }, { int(1), Any() }, { Any(), Any() } };
	commands[CMD_SURFROUGHENJAGGED].unit = "m";
	commands[CMD_SURFROUGHENJAGGED].descr = "[tc0,0.5,0.5,1/tc]Roughen given mesh (focused mesh if not specified) surfaces using a jagged pattern to given depth (m) and peak spacing (m). Roughen both sides by default, unless sides is specified as -z or z (string literal). The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_GENERATE2DGRAINS, CommandSpecifier(CMD_GENERATE2DGRAINS), "generate2dgrains");
	commands[CMD_GENERATE2DGRAINS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>generate2dgrains</b> <i>(meshname) spacing (seed)</i>";
	commands[CMD_GENERATE2DGRAINS].limits = { { Any(), Any() }, { double(0), Any() }, { int(0), Any() } };
	commands[CMD_GENERATE2DGRAINS].unit = "m";
	commands[CMD_GENERATE2DGRAINS].descr = "[tc0,0.5,0.5,1/tc]Generate 2D Voronoi cells in the xy plane at given average spacing in given mesh (focused mesh if not specified). The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_GENERATE3DGRAINS, CommandSpecifier(CMD_GENERATE3DGRAINS), "generate3dgrains");
	commands[CMD_GENERATE3DGRAINS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>generate3dgrains</b> <i>(meshname) spacing (seed)</i>";
	commands[CMD_GENERATE3DGRAINS].limits = { { Any(), Any() }, { double(0), Any() }, { int(0), Any() } };
	commands[CMD_GENERATE3DGRAINS].unit = "m";
	commands[CMD_GENERATE3DGRAINS].descr = "[tc0,0.5,0.5,1/tc]Generate 3D Voronoi cells at given average spacing in given mesh (focused mesh if not specified). The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_BENCHTIME, CommandSpecifier(CMD_BENCHTIME), "benchtime");
	commands[CMD_BENCHTIME].usage = "[tc0,0.5,0,1/tc]USAGE : <b>benchtime</b>";
	commands[CMD_BENCHTIME].descr = "[tc0,0.5,0.5,1/tc]Show the last simulation duration time in ms, between start and stop; used for performance becnhmarking.";
	commands[CMD_BENCHTIME].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_MATERIALSDATABASE, CommandSpecifier(CMD_MATERIALSDATABASE), "materialsdatabase");
	commands[CMD_MATERIALSDATABASE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>materialsdatabase</b> <i>(mdbname)</i>";
	commands[CMD_MATERIALSDATABASE].descr = "[tc0,0.5,0.5,1/tc]Switch materials database in use. This setting is not saved by savesim, so using loadsim doesn't affect this setting; default mdb set on program start.";

	commands.insert(CMD_ADDMATERIAL, CommandSpecifier(CMD_ADDMATERIAL), "addmaterial");
	commands[CMD_ADDMATERIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addmaterial<b> <i>name rectangle</i>";
	commands[CMD_ADDMATERIAL].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDMATERIAL].unit = "m";
	commands[CMD_ADDMATERIAL].descr = "[tc0,0.5,0.5,1/tc]Add a new mesh with material parameters loaded from the materials database. The name is the material name as found in the mdb file (see materialsdatabase command); this also determines the type of mesh to create, as well as the created mesh name. The rectangle (m) can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDMATERIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshname</i> - return name of mesh just added (can differ from the material name).";

	commands.insert(CMD_SETMATERIAL, CommandSpecifier(CMD_SETMATERIAL), "setmaterial");
	commands[CMD_SETMATERIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setmaterial<b> <i>name rectangle</i>";
	commands[CMD_SETMATERIAL].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_SETMATERIAL].unit = "m";
	commands[CMD_SETMATERIAL].descr = "[tc0,0.5,0.5,1/tc]Set a single mesh with material parameters loaded from the materials database (deleting all other meshes). The name is the material name as found in the mdb file (see materialsdatabase command); this also determines the type of mesh to create, as well as the created mesh name. The rectangle (m) can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_SETMATERIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshname</i> - return name of mesh just added (can differ from the material name).";

	commands.insert(CMD_ADDMDBENTRY, CommandSpecifier(CMD_ADDMDBENTRY), "addmdbentry");
	commands[CMD_ADDMDBENTRY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addmdbentry</b> <i>meshname (materialname)</i>";
	commands[CMD_ADDMDBENTRY].descr = "[tc0,0.5,0.5,1/tc]Add new entry in the local materials database from parameters in the given mesh. The name of the new entry is set to materialname if specified, else set to meshname. For a complete entry you should then edit the mdb file manually with all the appropriate fields shown there.";

	commands.insert(CMD_DELMDBENTRY, CommandSpecifier(CMD_DELMDBENTRY), "delmdbentry");
	commands[CMD_DELMDBENTRY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delmdbentry</b> <i>materialname</i>";
	commands[CMD_DELMDBENTRY].descr = "[tc0,0.5,0.5,1/tc]Delete entry in the local materials database (see materialsdatabase for current selection).";

	commands.insert(CMD_REFRESHMDB, CommandSpecifier(CMD_REFRESHMDB), "refreshmdb");
	commands[CMD_REFRESHMDB].usage = "[tc0,0.5,0,1/tc]USAGE : <b>refreshmdb</b>";
	commands[CMD_REFRESHMDB].descr = "[tc0,0.5,0.5,1/tc]Reload the local materials database (see materialsdatabase for current selection). This is useful if you modify the values in the materials database file externally.";

	commands.insert(CMD_REQMDBSYNC, CommandSpecifier(CMD_REQMDBSYNC), "requestmdbsync");
	commands[CMD_REQMDBSYNC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>requestmdbsync</b> <i>materialname (email)</i>";
	commands[CMD_REQMDBSYNC].descr = "[tc0,0.5,0.5,1/tc]Request the given entry in the local materials database is added to the online shared materials database. This must be a completed entry - see manual for instructions. The entry will be checked before being made available to all users through the online materials database. If you want to receive an update about the status of this request include an email address.";

	commands.insert(CMD_UPDATEMDB, CommandSpecifier(CMD_UPDATEMDB), "updatemdb");
	commands[CMD_UPDATEMDB].usage = "[tc0,0.5,0,1/tc]USAGE : <b>updatemdb</b>";
	commands[CMD_UPDATEMDB].descr = "[tc0,0.5,0.5,1/tc]Switch to, and update the local materials database from the online shared materials database.";

	commands.insert(CMD_SHOWLENGHTS, CommandSpecifier(CMD_SHOWLENGHTS), "showlengths");
	commands[CMD_SHOWLENGHTS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showlengths</b> <i>(meshname)</i>";
	commands[CMD_SHOWLENGHTS].descr = "[tc0,0.5,0.5,1/tc]Calculate a number of critical lengths for the focused mesh (must be ferromagnetic; focused mesh if not specified) to inform magnetization cellsize selection. lex = sqrt(2 A / mu0 Ms^2) : exchange length, l_Bloch = sqrt(A / K1) : Bloch wall width, l_sky = PI D / 4 K1 : Neel skyrmion wall width.";
	commands[CMD_SHOWLENGHTS].limits = { { Any(), Any() } };

	commands.insert(CMD_SHOWMCELLS, CommandSpecifier(CMD_SHOWMCELLS), "showmcells");
	commands[CMD_SHOWMCELLS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showmcells</b> <i>(meshname)</i>";
	commands[CMD_SHOWMCELLS].descr = "[tc0,0.5,0.5,1/tc]Show number of discretisation cells for magnetization for focused mesh (must be ferromagnetic; focused mesh if not specified).";
	commands[CMD_SHOWMCELLS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>n</i>";
	commands[CMD_SHOWMCELLS].limits = { { Any(), Any() } };

	commands.insert(CMD_LOADOVF2MAG, CommandSpecifier(CMD_LOADOVF2MAG), "loadovf2mag");
	commands[CMD_LOADOVF2MAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2mag</b> <i>(meshname) (renormalize_value) (directory/)filename</i>";
	commands[CMD_LOADOVF2MAG].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing magnetization data, into the given mesh (which must be magnetic; focused mesh if not specified), mapping the data to the current mesh dimensions. By default the loaded data will not be renormalized: renormalize_value = 0. If a value is specified for renormalize_value, the loaded data will be renormalized to it (e.g. this would be an Ms value).";
	commands[CMD_LOADOVF2MAG].unit = "A/m";
	commands[CMD_LOADOVF2MAG].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_LOADOVF2FIELD, CommandSpecifier(CMD_LOADOVF2FIELD), "loadovf2field");
	commands[CMD_LOADOVF2FIELD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2field</b> <i>(meshname) (directory/)filename</i>";
	commands[CMD_LOADOVF2FIELD].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing applied field data, into the given mesh (which must be magnetic; focused mesh if not specified), mapping the data to the current mesh dimensions. Alternatively if meshname is supermesh, then load a global field with absolute rectangle coordinates and its own discretization.";
	commands[CMD_LOADOVF2FIELD].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_CLEARGLOBALFIELD, CommandSpecifier(CMD_CLEARGLOBALFIELD), "clearglobalfield");
	commands[CMD_CLEARGLOBALFIELD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearglobalfield</b>";
	commands[CMD_CLEARGLOBALFIELD].descr = "[tc0,0.5,0.5,1/tc]Clear a set global field (previously loaded with loadovf2field on supermesh).";

	commands.insert(CMD_SHIFTGLOBALFIELD, CommandSpecifier(CMD_SHIFTGLOBALFIELD), "shiftglobalfield");
	commands[CMD_SHIFTGLOBALFIELD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shiftglobalfield</b> <i>shift</i>";
	commands[CMD_SHIFTGLOBALFIELD].unit = "m";
	commands[CMD_SHIFTGLOBALFIELD].limits = { { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) } };
	commands[CMD_SHIFTGLOBALFIELD].descr = "[tc0,0.5,0.5,1/tc]Shift rectangle of global field (previously loaded with loadovf2field on supermesh).";
	commands[CMD_SHIFTGLOBALFIELD].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>rect</i> (global field rectangle)";

	commands.insert(CMD_LOADOVF2TEMP, CommandSpecifier(CMD_LOADOVF2TEMP), "loadovf2temp");
	commands[CMD_LOADOVF2TEMP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2temp</b> <i>(meshname) (directory/)filename</i>";
	commands[CMD_LOADOVF2TEMP].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing temperature data, into the given mesh (which must have heat module enabled; focused mesh if not specified), mapping the data to the current mesh dimensions.";
	commands[CMD_LOADOVF2TEMP].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_LOADOVF2CURR, CommandSpecifier(CMD_LOADOVF2CURR), "loadovf2curr");
	commands[CMD_LOADOVF2CURR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2curr</b> <i>(meshname) (directory/)filename</i>";
	commands[CMD_LOADOVF2CURR].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing current density data, into the given mesh  (which must have transport module enabled; focused mesh if not specified), mapping the data to the current mesh dimensions.";
	commands[CMD_LOADOVF2CURR].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SAVEOVF2MAG, CommandSpecifier(CMD_SAVEOVF2MAG), "saveovf2mag");
	commands[CMD_SAVEOVF2MAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveovf2mag</b> <i>(meshname) (n) (data_type) (directory/)filename</i>";
	commands[CMD_SAVEOVF2MAG].descr = "[tc0,0.5,0.5,1/tc]Save an OOMMF-style OVF 2.0 file containing magnetization data from the given mesh (which must be magnetic; focused mesh if not specified). You can normalize the data to Ms0 value by specifying the n flag (e.g. saveovf2mag n filename) - by default the data is not normalized. You can specify the data type as data_type = bin4 (single precision 4 bytes per float), data_type = bin8 (double precision 8 bytes per float), or data_type = text. By default bin8 is used.";
	commands[CMD_SAVEOVF2MAG].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SAVEOVF2PARAMVAR, CommandSpecifier(CMD_SAVEOVF2PARAMVAR), "saveovf2param");
	commands[CMD_SAVEOVF2PARAMVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveovf2param</b> <i>(meshname) (data_type) paramname (directory/)filename</i>";
	commands[CMD_SAVEOVF2PARAMVAR].descr = "[tc0,0.5,0.5,1/tc]Save an OOMMF-style OVF 2.0 file containing the named parameter spatial variation data from the given mesh (focused mesh if not specified). You can specify the data type as data_type = bin4 (single precision 4 bytes per float), data_type = bin8 (double precision 8 bytes per float), or data_type = text. By default bin8 is used.";
	commands[CMD_SAVEOVF2PARAMVAR].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SAVEOVF2, CommandSpecifier(CMD_SAVEOVF2), "saveovf2");
	commands[CMD_SAVEOVF2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveovf2</b> <i>(meshname (quantity)) (data_type) (directory/)filename</i>";
	commands[CMD_SAVEOVF2].descr = "[tc0,0.5,0.5,1/tc]Save an OOMMF-style OVF 2.0 file containing data from the given mesh (focused mesh if not specified), depending on currently displayed quantites, or the named quantity if given (see output of display command for possible quantities). You can specify the data type as data_type = bin4 (single precision 4 bytes per float), data_type = bin8 (double precision 8 bytes per float), or data_type = text. By default bin8 is used.";
	commands[CMD_SAVEOVF2].limits = { { Any(), Any() }, { Any(), Any() }, { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_LOADOVF2DISP, CommandSpecifier(CMD_LOADOVF2DISP), "loadovf2disp");
	commands[CMD_LOADOVF2DISP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2disp</b> <i>(meshname) (directory/)filename</i>";
	commands[CMD_LOADOVF2DISP].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing mechanical displacement data, into the given mesh (which must be ferromagnetic and have the melastic module enabled; focused mesh if not specified), mapping the data to the current mesh dimensions. From the mechanical displacement the strain tensor is calculated.";
	commands[CMD_LOADOVF2DISP].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_LOADOVF2STRAIN, CommandSpecifier(CMD_LOADOVF2STRAIN), "loadovf2strain");
	commands[CMD_LOADOVF2STRAIN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2strain</b> <i>(meshname) (directory/)filename_diag filename_odiag</i>";
	commands[CMD_LOADOVF2STRAIN].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing strain tensor data, into the given mesh (which must be ferromagnetic and have the melastic module enabled; focused mesh if not specified), mapping the data to the current mesh dimensions. The symmetric strain tensor is applicable for a cubic crystal, and has 3 diagonal component (specified in filename_diag with vector data as xx, yy, zz), and 3 off-diagonal components (specified in filename_odiag with vector data as yz, xz, xy).";
	commands[CMD_LOADOVF2DISP].limits = { { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SCRIPTSERVER, CommandSpecifier(CMD_SCRIPTSERVER), "scriptserver");
	commands[CMD_SCRIPTSERVER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>scriptserver</b> <i>status</i>";
	commands[CMD_SCRIPTSERVER].descr = "[tc0,0.5,0.5,1/tc]Enable or disable the script communication server. When enabled the program will listen for commands received using network sockets on port 1542.";

	commands.insert(CMD_CHECKUPDATES, CommandSpecifier(CMD_CHECKUPDATES), "checkupdates");
	commands[CMD_CHECKUPDATES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>checkupdates</b>";
	commands[CMD_CHECKUPDATES].descr = "[tc0,0.5,0.5,1/tc]Connect to boris-spintronics.uk to check if updates to program or materials database are available.";

	commands.insert(CMD_EQUATIONCONSTANTS, CommandSpecifier(CMD_EQUATIONCONSTANTS), "equationconstants");
	commands[CMD_EQUATIONCONSTANTS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>equationconstants</b> <i>name value</i>";
	commands[CMD_EQUATIONCONSTANTS].descr = "[tc0,0.5,0.5,1/tc]Create or edit user constant to be used in text equations.";
	
	commands.insert(CMD_DELEQUATIONCONSTANT, CommandSpecifier(CMD_DELEQUATIONCONSTANT), "delequationconstant");
	commands[CMD_DELEQUATIONCONSTANT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delequationconstant</b> <i>name</i>";
	commands[CMD_DELEQUATIONCONSTANT].descr = "[tc0,0.5,0.5,1/tc]Delete named user constant used in text equations.";

	commands.insert(CMD_CLEAREQUATIONCONSTANTS, CommandSpecifier(CMD_CLEAREQUATIONCONSTANTS), "clearequationconstants");
	commands[CMD_CLEAREQUATIONCONSTANTS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearequationconstants</b>";
	commands[CMD_CLEAREQUATIONCONSTANTS].descr = "[tc0,0.5,0.5,1/tc]Clear all user-defined constants for text equations.";

	commands.insert(CMD_FLUSHERRORLOG, CommandSpecifier(CMD_FLUSHERRORLOG), "flusherrorlog");
	commands[CMD_FLUSHERRORLOG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>flusherrorlog</b>";
	commands[CMD_FLUSHERRORLOG].descr = "[tc0,0.5,0.5,1/tc]Clear error log.";

	commands.insert(CMD_ERRORLOG, CommandSpecifier(CMD_ERRORLOG), "errorlog");
	commands[CMD_ERRORLOG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>errorlog</b> <i>status</i>";
	commands[CMD_ERRORLOG].descr = "[tc0,0.5,0.5,1/tc]Set error log status.";

	commands.insert(CMD_STARTUPUPDATECHECK, CommandSpecifier(CMD_STARTUPUPDATECHECK), "startupupdatecheck");
	commands[CMD_STARTUPUPDATECHECK].usage = "[tc0,0.5,0,1/tc]USAGE : <b>startupupdatecheck</b> <i>status</i>";
	commands[CMD_STARTUPUPDATECHECK].descr = "[tc0,0.5,0.5,1/tc]Set startup update check flag.";

	commands.insert(CMD_STARTUPSCRIPTSERVER, CommandSpecifier(CMD_STARTUPSCRIPTSERVER), "startupscriptserver");
	commands[CMD_STARTUPSCRIPTSERVER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>startupscriptserver</b> <i>status</i>";
	commands[CMD_STARTUPSCRIPTSERVER].descr = "[tc0,0.5,0.5,1/tc]Set startup script server flag.";

	commands.insert(CMD_THREADS, CommandSpecifier(CMD_THREADS), "threads");
	commands[CMD_THREADS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>threads</b> <i>number</i>";
	commands[CMD_THREADS].descr = "[tc0,0.5,0.5,1/tc]Set number of threads to use for cuda 0 computations. Value of zero means set maximum available.";
	commands[CMD_THREADS].limits = { { int(0), Any(omp_get_num_procs()) } };
	commands[CMD_THREADS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>num_threads</i>";

	commands.insert(CMD_SERVERPORT, CommandSpecifier(CMD_SERVERPORT), "serverport");
	commands[CMD_SERVERPORT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>serverport</b> <i>port</i>";
	commands[CMD_SERVERPORT].descr = "[tc0,0.5,0.5,1/tc]Set script server port.";
	commands[CMD_SERVERPORT].limits = { { int(0), Any() } };
	commands[CMD_SERVERPORT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>port</i>";

	commands.insert(CMD_SERVERPWD, CommandSpecifier(CMD_SERVERPWD), "serverpassword");
	commands[CMD_SERVERPWD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>serverpassword</b> <i>password</i>";
	commands[CMD_SERVERPWD].descr = "[tc0,0.5,0.5,1/tc]Set/change script server password - this is used to authenticate remote client messages. By default no password is set (blank).";
	commands[CMD_SERVERPWD].limits = { { Any(), Any() } };

	commands.insert(CMD_SERVERSLEEPMS, CommandSpecifier(CMD_SERVERSLEEPMS), "serversleepms");
	commands[CMD_SERVERSLEEPMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>serversleepms</b> <i>time_ms</i>";
	commands[CMD_SERVERSLEEPMS].descr = "[tc0,0.5,0.5,1/tc]Set script server thread sleep time in ms: lower value makes server more responsive, but increases CPU load.";
	commands[CMD_SERVERSLEEPMS].limits = { { int(1), Any() } };
	commands[CMD_SERVERSLEEPMS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>time_ms</i>";

	commands.insert(CMD_NEWINSTANCE, CommandSpecifier(CMD_NEWINSTANCE), "newinstance");
	commands[CMD_NEWINSTANCE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>newinstance</b> <i>port (cudaDevice (password))</i>";
	commands[CMD_NEWINSTANCE].descr = "[tc0,0.5,0.5,1/tc]Start a new local Boris instance with given server port, and optionally cuda device number (0/1/2/3/...); a value of -1 means automatically determine cuda device. If password not blank this should be specified otherwise server will not start.";
	commands[CMD_NEWINSTANCE].limits = { { int(0), Any() }, { int(-1), Any() }, { Any(), Any() } };

	commands.insert(CMD_VERSIONUPDATE, CommandSpecifier(CMD_VERSIONUPDATE), "versionupdate");
	commands[CMD_VERSIONUPDATE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>versionupdate</b> <i>action (target_version)</i>";
	commands[CMD_VERSIONUPDATE].descr = "[tc0,0.5,0.5,1/tc]Incremental version update command. action = check : find if incremental update is available for current program version. action = download : download available incremental updates if any (or up to target version if specified). action = install : install incremental updates if any (or up to target version if specified), downloading if needed. action = rollback : roll back to specified target version if possible.";
	commands[CMD_VERSIONUPDATE].limits = { { Any(), Any() }, { int(0), Any() } };

	commands.insert(CMD_SHOWTC, CommandSpecifier(CMD_SHOWTC), "showtc");
	commands[CMD_SHOWTC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showtc</b> <i>(meshname)</i>";
	commands[CMD_SHOWTC].descr = "[tc0,0.5,0.5,1/tc]Show predicted Tc value (K) for current focused mesh (must be atomistic; focused mesh if not specified), using formula Tc = J*e*z/3kB, where e is the spin-wave correction factor, and z is the coordination number.";
	commands[CMD_SHOWTC].limits = { { Any(), Any() } };
	commands[CMD_SHOWTC].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Tc</i>";

	commands.insert(CMD_SHOWMS, CommandSpecifier(CMD_SHOWMS), "showms");
	commands[CMD_SHOWMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showms</b> <i>(meshname)</i>";
	commands[CMD_SHOWMS].descr = "[tc0,0.5,0.5,1/tc]Show predicted saturation magnetization (A/m) value for current focused mesh (must be atomistic; focused mesh if not specified), using formula Ms = mu_s*n/a^3, where n is the number of atomic moments per unit cell, and a is the atomic cell size.";
	commands[CMD_SHOWMS].limits = { { Any(), Any() } };
	commands[CMD_SHOWMS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Ms</i>";

	commands.insert(CMD_SHOWA, CommandSpecifier(CMD_SHOWA), "showa");
	commands[CMD_SHOWA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showa</b> <i>(meshname)</i>";
	commands[CMD_SHOWA].descr = "[tc0,0.5,0.5,1/tc]Show predicted exchange stiffness (J/m) value for current focused mesh (must be atomistic; focused mesh if not specified), using formula A = J*n/2a, where n is the number of atomic moments per unit cell, and a is the atomic cell size.";
	commands[CMD_SHOWA].limits = { { Any(), Any() } };
	commands[CMD_SHOWA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>A</i>";

	commands.insert(CMD_SHOWK, CommandSpecifier(CMD_SHOWK), "showk");
	commands[CMD_SHOWK].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showk</b> <i>(meshname)</i>";
	commands[CMD_SHOWK].descr = "[tc0,0.5,0.5,1/tc]Show predicted uniaxial anisotropy (J/m^3) constant value for current focused mesh (must be atomistic; focused mesh if not specified), using formula K = k*n/a^3, where n is the number of atomic moments per unit cell, and a is the atomic cell size.";
	commands[CMD_SHOWK].limits = { { Any(), Any() } };
	commands[CMD_SHOWK].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>A</i>";

	commands.insert(CMD_SKYPOSDMUL, CommandSpecifier(CMD_SKYPOSDMUL), "skyposdmul");
	commands[CMD_SKYPOSDMUL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyposdmul</b> <i>(meshname) multiplier</i>";
	commands[CMD_SKYPOSDMUL].descr = "[tc0,0.5,0.5,1/tc]Set skyrmion diameter multiplier for given meshname (focused mesh if not specified) which determines skypos tracking rectangle size. Default is 2.0, i.e. tracking rectangle side is twice the diameter along each axis. Reduce if skyrmion density is too large, but at the risk of losing skyrmion tracking (recommended to keep above 1.2).";
	commands[CMD_SKYPOSDMUL].limits = { {Any(), Any()}, { double(1.0), double(10.0) } };
	commands[CMD_SKYPOSDMUL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>multiplier</i>";

	commands.insert(CMD_DWPOS_COMPONENT, CommandSpecifier(CMD_DWPOS_COMPONENT), "dwposcomponent");
	commands[CMD_DWPOS_COMPONENT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dwposcomponent</b> <i>value</i>";
	commands[CMD_DWPOS_COMPONENT].descr = "[tc0,0.5,0.5,1/tc]Set vector component used to extract dwpos_x, dwpos_y, dwpos_z data: this should be the component with a tanh profile. -1: automatic (detect tanh component), 0: x, 1: y, 2: z. Default is automatic, but you might want to specify component as the automatic detection can fail when very noisy (e.g. close to Tc for atomistic simulations).";
	commands[CMD_DWPOS_COMPONENT].limits = { {int(-1), int(2)} };
	commands[CMD_DWPOS_COMPONENT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_MCSERIAL, CommandSpecifier(CMD_MCSERIAL), "mcserial");
	commands[CMD_MCSERIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcserial</b> <i>(meshname) value</i>";
	commands[CMD_MCSERIAL].descr = "[tc0,0.5,0.5,1/tc]Change Monte-Carlo algorithm type. 0: parallel (default); red-black ordering parallelization 1: serial (testing only); spins picked in random order. If meshname not specified setting is applied to all atomistic meshes.";
	commands[CMD_MCSERIAL].limits = { {Any(), Any()}, { int(0), int(1) } };

	commands.insert(CMD_MCCONSTRAIN, CommandSpecifier(CMD_MCCONSTRAIN), "mcconstrain");
	commands[CMD_MCCONSTRAIN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcconstrain</b> <i>(meshname) value</i>";
	commands[CMD_MCCONSTRAIN].descr = "[tc0,0.5,0.5,1/tc]Set value 0 to revert to classic Monte-Carlo Metropolis for ASD. Set a unit vector direction value (x y z) to switch to constrained Monte Carlo as described in PRB 82, 054415 (2010). If meshname not specified setting is applied to all atomistic meshes.";
	commands[CMD_MCCONSTRAIN].limits = { {Any(), Any()}, { DBL3(), Any() } };

	commands.insert(CMD_MCDISABLE, CommandSpecifier(CMD_MCDISABLE), "mcdisable");
	commands[CMD_MCDISABLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcdisable</b> <i>(meshname) status</i>";
	commands[CMD_MCDISABLE].descr = "[tc0,0.5,0.5,1/tc]Disable or enable Monte Carlo algorithm for given mesh (focused mesh if not specified).";
	commands[CMD_MCDISABLE].limits = { {Any(), Any()}, { int(0), int(1) } };

	commands.insert(CMD_MCCOMPUTEFIELDS, CommandSpecifier(CMD_MCCOMPUTEFIELDS), "mccomputefields");
	commands[CMD_MCCOMPUTEFIELDS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mccomputefields</b> <i>status</i>";
	commands[CMD_MCCOMPUTEFIELDS].descr = "[tc0,0.5,0.5,1/tc]Enable/disable modules update during Monte Carlo stages. Disabled by default, but if you want to save energy density terms you need to enabled this - equivalent to running ComputeFields after every Monte Carlo step.";
	commands[CMD_MCCOMPUTEFIELDS].limits = { { int(0), Any() } };
	commands[CMD_MCCOMPUTEFIELDS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_MCCONEANGLELIMITS, CommandSpecifier(CMD_MCCONEANGLELIMITS), "mcconeangle");
	commands[CMD_MCCONEANGLELIMITS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcconeangle</b> <i>min_angle max_angle</i>";
	commands[CMD_MCCONEANGLELIMITS].descr = "[tc0,0.5,0.5,1/tc]Set Monte Carlo cone angle limits (1 to 180 degrees). Setting min and max values the same results in a fixed cone angle Monte Carlo algorithm.";
	commands[CMD_MCCONEANGLELIMITS].limits = { {DBL2(), DBL2(180, 180)} };
	commands[CMD_MCCONEANGLELIMITS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>min_angle max_angle</i>";

	commands.insert(CMD_PRNGSEED, CommandSpecifier(CMD_PRNGSEED), "prngseed");
	commands[CMD_PRNGSEED].usage = "[tc0,0.5,0,1/tc]USAGE : <b>prngseed</b> <i>(meshname) seed</i>";
	commands[CMD_PRNGSEED].descr = "[tc0,0.5,0.5,1/tc]Set PRNG seed, for computations with stochasticity, in given mesh (all meshes if name not given). If seed = 0 then system tick count is used as seed (default behavior), otherwise the fixed set value is used. NOTE: seed values > 0 only take effect for computations with cuda 1.";
	commands[CMD_PRNGSEED].limits = { {Any(), Any()}, {int(0), Any()} };
	commands[CMD_PRNGSEED].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>seed</i>";

	commands.insert(CMD_SHAPEMOD_ROT, CommandSpecifier(CMD_SHAPEMOD_ROT), "shape_rotation");
	commands[CMD_SHAPEMOD_ROT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_rotation</b> <i>psi theta phi</i>";
	commands[CMD_SHAPEMOD_ROT].descr = "[tc0,0.5,0.5,1/tc]Set modifier for shape generator commands: rotation using psi (around y), theta (around x) and phi (around z) in degrees.";
	commands[CMD_SHAPEMOD_ROT].limits = { { DBL3(-360), DBL3(360) } };

	commands.insert(CMD_SHAPEMOD_REP, CommandSpecifier(CMD_SHAPEMOD_REP), "shape_repetitions");
	commands[CMD_SHAPEMOD_REP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_repetitions</b> <i>x y z</i>";
	commands[CMD_SHAPEMOD_REP].descr = "[tc0,0.5,0.5,1/tc]Set modifier for shape generator commands: number of shape repetitions along x, y, z.";
	commands[CMD_SHAPEMOD_REP].limits = { { INT3(1), Any() } };
	commands[CMD_SHAPEMOD_REP].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_SHAPEMOD_DISP, CommandSpecifier(CMD_SHAPEMOD_DISP), "shape_displacement");
	commands[CMD_SHAPEMOD_DISP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_displacement</b> <i>x y z</i>";
	commands[CMD_SHAPEMOD_DISP].descr = "[tc0,0.5,0.5,1/tc]Set modifier for shape generator commands: displacements along x, y, z to use with shape repetitions (shape_repetitions).";
	commands[CMD_SHAPEMOD_DISP].unit = "m";
	commands[CMD_SHAPEMOD_DISP].limits = { { DBL3(), Any() } };
	commands[CMD_SHAPEMOD_DISP].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_SHAPEMOD_METHOD, CommandSpecifier(CMD_SHAPEMOD_METHOD), "shape_method");
	commands[CMD_SHAPEMOD_METHOD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_method</b> <i>method</i>";
	commands[CMD_SHAPEMOD_METHOD].descr = "[tc0,0.5,0.5,1/tc]Set modifier for shape generator commands: method to use when applying shape as string literal: add or sub.";
	commands[CMD_SHAPEMOD_METHOD].limits = { { Any(), Any() } };
	commands[CMD_SHAPEMOD_METHOD].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>name</i>";

	commands.insert(CMD_SHAPE_DISK, CommandSpecifier(CMD_SHAPE_DISK), "shape_disk");
	commands[CMD_SHAPE_DISK].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_disk</b> <i>(meshname) dia_x dia_y cpos_x cpos_y (z_start z_end)</i>";
	commands[CMD_SHAPE_DISK].descr = "[tc0,0.5,0.5,1/tc]Set a disk shape in given mesh (focused mesh if not specified) with given x and y diameters at given centre position with x and y coordinates. Optionally specify z start and end coordinates, otherwise entire mesh thickness used.";
	commands[CMD_SHAPE_DISK].unit = "m";
	commands[CMD_SHAPE_DISK].limits = { { Any(), Any() }, { DBL2(), Any() }, { DBL2(-MAXSIMSPACE),  DBL2(+MAXSIMSPACE) }, { DBL2(-MAXSIMSPACE),  DBL2(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_DISK].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z</i>";

	commands.insert(CMD_SHAPE_RECT, CommandSpecifier(CMD_SHAPE_RECT), "shape_rect");
	commands[CMD_SHAPE_RECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_rect</b> <i>(meshname) len_x len_y cpos_x cpos_y (z_start z_end)</i>";
	commands[CMD_SHAPE_RECT].descr = "[tc0,0.5,0.5,1/tc]Set a rectangle shape in given mesh (focused mesh if not specified) with given x and y lengths at given centre position with x and y coordinates. Optionally specify z start and end coordinates, otherwise entire mesh thickness used.";
	commands[CMD_SHAPE_RECT].unit = "m";
	commands[CMD_SHAPE_RECT].limits = { { Any(), Any() }, { DBL2(), Any() }, { DBL2(-MAXSIMSPACE),  DBL2(+MAXSIMSPACE) }, { DBL2(-MAXSIMSPACE),  DBL2(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_RECT].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z</i>";

	commands.insert(CMD_SHAPE_TRIANGLE, CommandSpecifier(CMD_SHAPE_TRIANGLE), "shape_triangle");
	commands[CMD_SHAPE_TRIANGLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_triangle</b> <i>(meshname) len_x len_y cpos_x cpos_y (z_start z_end)</i>";
	commands[CMD_SHAPE_TRIANGLE].descr = "[tc0,0.5,0.5,1/tc]Set a triangle shape in given mesh (focused mesh if not specified) with given x and y lengths at given centre position with x and y coordinates. Optionally specify z start and end coordinates, otherwise entire mesh thickness used.";
	commands[CMD_SHAPE_TRIANGLE].unit = "m";
	commands[CMD_SHAPE_TRIANGLE].limits = { { Any(), Any() }, { DBL2(), Any() }, { DBL2(-MAXSIMSPACE),  DBL2(+MAXSIMSPACE) }, { DBL2(-MAXSIMSPACE),  DBL2(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_TRIANGLE].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method</i>";

	commands.insert(CMD_SHAPE_ELLIPSOID, CommandSpecifier(CMD_SHAPE_ELLIPSOID), "shape_ellipsoid");
	commands[CMD_SHAPE_ELLIPSOID].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_ellipsoid</b> <i>(meshname) dia_x dia_y dia_z cpos_x cpos_y cpos_z</i>";
	commands[CMD_SHAPE_ELLIPSOID].descr = "[tc0,0.5,0.5,1/tc]Set a prolate ellipsoid shape in given mesh (focused mesh if not specified) with given diameters (x, y, z) at given centre position coordinates (x, y, z).";
	commands[CMD_SHAPE_ELLIPSOID].unit = "m";
	commands[CMD_SHAPE_ELLIPSOID].limits = { { Any(), Any() }, { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_ELLIPSOID].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method</i>";

	commands.insert(CMD_SHAPE_PYRAMID, CommandSpecifier(CMD_SHAPE_PYRAMID), "shape_pyramid");
	commands[CMD_SHAPE_PYRAMID].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_pyramid</b> <i>(meshname) len_x len_y len_z cpos_x cpos_y cpos_z</i>";
	commands[CMD_SHAPE_PYRAMID].descr = "[tc0,0.5,0.5,1/tc]Set a pyramid shape in given mesh (focused mesh if not specified) with given lengths (x, y, z) at given centre position coordinates (x, y, z).";
	commands[CMD_SHAPE_PYRAMID].unit = "m";
	commands[CMD_SHAPE_PYRAMID].limits = { { Any(), Any() }, { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_PYRAMID].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method</i>";

	commands.insert(CMD_SHAPE_TETRAHEDRON, CommandSpecifier(CMD_SHAPE_TETRAHEDRON), "shape_tetrahedron");
	commands[CMD_SHAPE_TETRAHEDRON].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_tetrahedron</b> <i>(meshname) len_x len_y len_z cpos_x cpos_y cpos_z</i>";
	commands[CMD_SHAPE_TETRAHEDRON].descr = "[tc0,0.5,0.5,1/tc]Set a tetrahedron shape in given mesh (focused mesh if not specified) with given lengths (x, y, z) at given centre position coordinates (x, y, z).";
	commands[CMD_SHAPE_TETRAHEDRON].unit = "m";
	commands[CMD_SHAPE_TETRAHEDRON].limits = { { Any(), Any() }, { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_TETRAHEDRON].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method</i>";

	commands.insert(CMD_SHAPE_CONE, CommandSpecifier(CMD_SHAPE_CONE), "shape_cone");
	commands[CMD_SHAPE_CONE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_cone</b> <i>(meshname) len_x len_y len_z cpos_x cpos_y cpos_z</i>";
	commands[CMD_SHAPE_CONE].descr = "[tc0,0.5,0.5,1/tc]Set a conical shape in given mesh (focused mesh if not specified) with given lengths (x, y, z) at given centre position coordinates (x, y, z).";
	commands[CMD_SHAPE_CONE].unit = "m";
	commands[CMD_SHAPE_CONE].limits = { { Any(), Any() }, { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_CONE].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method</i>";

	commands.insert(CMD_SHAPE_TORUS, CommandSpecifier(CMD_SHAPE_TORUS), "shape_torus");
	commands[CMD_SHAPE_TORUS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_torus</b> <i>(meshname) len_x len_y len_z cpos_x cpos_y cpos_z</i>";
	commands[CMD_SHAPE_TORUS].descr = "[tc0,0.5,0.5,1/tc]Set a torus shape in given mesh (focused mesh if not specified) with given lengths (x, y, z) at given centre position coordinates (x, y, z).";
	commands[CMD_SHAPE_TORUS].unit = "m";
	commands[CMD_SHAPE_TORUS].limits = { { Any(), Any() }, { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) } };
	commands[CMD_SHAPE_TORUS].return_descr = "[tc0,0.5,0,1/tc]Script return values: shape object definition as <i>name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method</i>";

	commands.insert(CMD_SHAPE_SET, CommandSpecifier(CMD_SHAPE_SET), "shape_set");
	commands[CMD_SHAPE_SET].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_set</b> <i>(meshname) name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_psi rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method ...</i>";
	commands[CMD_SHAPE_SET].descr = "[tc0,0.5,0.5,1/tc]General command for setting a shape: meant to be used with scripted simulations, not directly from the console. Can be used with the return data from shape_... commands, and can take multiple elementary shapes to form a composite object. Set a named shape in given mesh (focused mesh if not specified) with given dimensions (x, y, z) at given centre position coordinates (x, y, z).";
	commands[CMD_SHAPE_SET].unit = "m";
	commands[CMD_SHAPE_SET].limits = { { Any(), Any() },  { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) }, { DBL3(-360), DBL3(360) }, { INT3(1), Any() }, { DBL3(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SHAPE_GET, CommandSpecifier(CMD_SHAPE_GET), "shape_get");
	commands[CMD_SHAPE_GET].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_get</b> <i>(meshname (quantity)) name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_psi rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method ...</i>";
	commands[CMD_SHAPE_GET].descr = "[tc0,0.5,0.5,1/tc]Get average value from a shape in given mesh (focused mesh if not specified), depending on currently displayed quantites, or the named quantity if given (see output of display command for possible quantities), where shape definition is shown in shape_set command help. NOTE: this command does not work yet with CUDA enabled, will be finished in a future version.";
	commands[CMD_SHAPE_GET].unit = "m";
	commands[CMD_SHAPE_GET].limits = { { Any(), Any() }, { Any(), Any() }, { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) }, { DBL3(-360), DBL3(360) }, { INT3(1), Any() }, { DBL3(), Any() }, { Any(), Any() } };

	commands.insert(CMD_SETSHAPEANGLE, CommandSpecifier(CMD_SETSHAPEANGLE), "shape_setangle");
	commands[CMD_SETSHAPEANGLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_setangle</b> <i>(meshname) name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_psi rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method ... theta polar</i>";
	commands[CMD_SETSHAPEANGLE].descr = "[tc0,0.5,0.5,1/tc]Set magnetization angle (theta - polar, phi - azimuthal) in degrees for given shape identifier (or composite shape) in given mesh (focused mesh if not specified). The shape identifier is returned by one of the shape_... commands, also same parameters taken by the shape_set command. As with shape_set this command is intended for scripted simulations only.";
	commands[CMD_SETSHAPEANGLE].unit = "m";
	commands[CMD_SETSHAPEANGLE].limits = { { Any(), Any() },  { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) }, { DBL3(-360), DBL3(360) }, { INT3(1), Any() }, { DBL3(), Any() }, { Any(), Any() }, {DBL2(-360), DBL2(360)} };

	commands.insert(CMD_SHAPE_SETPARAM, CommandSpecifier(CMD_SHAPE_SETPARAM), "shape_setparam");
	commands[CMD_SHAPE_SETPARAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>shape_setparam</b> <i>(meshname) paramname name dim_x dim_y dim_z cpos_x cpos_y cpos_z rot_psi rot_theta rot_phi repeat_x repeat_y repeat_z disp_x disp_y disp_z method ... scaling_value</i>";
	commands[CMD_SHAPE_SETPARAM].descr = "[tc0,0.5,0.5,1/tc]Set material parameter scaling value in given shape only, in given mesh (focused mesh if not specified). If you want to set effective parameter value directly, not as a scaling value, then set base parameter value (setparam) to 1 and the real parameter value through this command.";
	commands[CMD_SHAPE_SETPARAM].unit = "m";
	commands[CMD_SHAPE_SETPARAM].limits = { { Any(), Any() }, { Any(), Any() },  { DBL3(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(+MAXSIMSPACE) }, { DBL3(-360), DBL3(360) }, { INT3(1), Any() }, { DBL3(), Any() }, { Any(), Any() }, { Any(), Any() } };

	commands.insert(CMD_RUNCOMMBUFFER, CommandSpecifier(CMD_RUNCOMMBUFFER), "runcommbuffer");
	commands[CMD_RUNCOMMBUFFER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>runcommbuffer</b>";
	commands[CMD_RUNCOMMBUFFER].descr = "[tc0,0.5,0.5,1/tc]Execute commands stored in command buffer.";

	commands.insert(CMD_CLEARCOMMBUFFER, CommandSpecifier(CMD_CLEARCOMMBUFFER), "clearcommbuffer");
	commands[CMD_CLEARCOMMBUFFER].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearcommbuffer</b>";
	commands[CMD_CLEARCOMMBUFFER].descr = "[tc0,0.5,0.5,1/tc]Clear commands stored in command buffer.";

	commands.insert(CMD_BUFFERCOMMAND, CommandSpecifier(CMD_BUFFERCOMMAND), "buffercommand");
	commands[CMD_BUFFERCOMMAND].usage = "[tc0,0.5,0,1/tc]USAGE : <b>buffercommand</b> <i>command (params...)</i>";
	commands[CMD_BUFFERCOMMAND].descr = "[tc0,0.5,0.5,1/tc]Append command with parameters to command buffer.";
	commands[CMD_BUFFERCOMMAND].limits = { { Any(), Any() } };

	commands.insert(CMD_DP_CLEARALL, CommandSpecifier(CMD_DP_CLEARALL), "dp_clearall");
	commands[CMD_DP_CLEARALL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_clearall</b>";
	commands[CMD_DP_CLEARALL].descr = "[tc0,0.5,0.5,1/tc]Clear all dp arrays.";

	commands.insert(CMD_DP_CLEAR, CommandSpecifier(CMD_DP_CLEAR), "dp_clear");
	commands[CMD_DP_CLEAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_clear</b> <i>indexes...</i>";
	commands[CMD_DP_CLEAR].limits = { { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_CLEAR].descr = "[tc0,0.5,0.5,1/tc]Clear dp arrays with specified indexes.";
	
	commands.insert(CMD_DP_SHOWSIZES, CommandSpecifier(CMD_DP_SHOWSIZES), "dp_showsizes");
	commands[CMD_DP_SHOWSIZES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_showsizes</b> <i>(dp_arr)</i>";
	commands[CMD_DP_SHOWSIZES].descr = "[tc0,0.5,0.5,1/tc]List sizes of all non-empty dp arrays, unless a specific dp_arr index is specified, in which case only show the size of dp_arr.";
	commands[CMD_DP_SHOWSIZES].limits = { { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SHOWSIZES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>dp_arr size if specified</i>";
	
	commands.insert(CMD_DP_GET, CommandSpecifier(CMD_DP_GET), "dp_get");
	commands[CMD_DP_GET].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_get</b> <i>dp_arr index</i>";
	commands[CMD_DP_GET].descr = "[tc0,0.5,0.5,1/tc]Show value in dp_arr at given index - the index must be within the dp_arr size.";
	commands[CMD_DP_GET].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), Any() } };
	commands[CMD_DP_GET].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_DP_SET, CommandSpecifier(CMD_DP_SET), "dp_set");
	commands[CMD_DP_SET].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_set</b> <i>dp_arr index value</i>";
	commands[CMD_DP_SET].descr = "[tc0,0.5,0.5,1/tc]Set value in dp_arr at given index - the index must be within the dp_arr size.";
	commands[CMD_DP_SET].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), Any() },{ Any(), Any() } };

	commands.insert(CMD_DP_LOAD, CommandSpecifier(CMD_DP_LOAD), "dp_load");
	commands[CMD_DP_LOAD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_load</b> <i>(directory/)filename file_indexes... dp_indexes...</i>";
	commands[CMD_DP_LOAD].limits = { { Any(), Any() },  { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_LOAD].descr = "[tc0,0.5,0.5,1/tc]Load data columns from filename into dp arrays. file_indexes are the column indexes in filename (.txt termination by default), dp_indexes are used for the dp arrays; count from 0. If directory not specified, the default one is used.";

	commands.insert(CMD_DP_SAVE, CommandSpecifier(CMD_DP_SAVE), "dp_save");
	commands[CMD_DP_SAVE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_save</b> <i>(directory/)filename dp_indexes...</i>";
	commands[CMD_DP_SAVE].limits = { { Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SAVE].descr = "[tc0,0.5,0.5,1/tc]Save specified dp arrays in filename (.txt termination by default). If directory not specified, the default one is used. dp_indexes are used for the dp arrays; count from 0.";

	commands.insert(CMD_DP_SAVEAPPEND, CommandSpecifier(CMD_DP_SAVEAPPEND), "dp_saveappend");
	commands[CMD_DP_SAVEAPPEND].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_saveappend</b> <i>(directory/)filename dp_indexes...</i>";
	commands[CMD_DP_SAVEAPPEND].limits = { { Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SAVEAPPEND].descr = "[tc0,0.5,0.5,1/tc]Save specified dp arrays in filename (.txt termination by default) by appending at the end. If directory not specified, the default one is used. dp_indexes are used for the dp arrays; count from 0.";

	commands.insert(CMD_DP_SAVEASROW, CommandSpecifier(CMD_DP_SAVEASROW), "dp_saveasrow");
	commands[CMD_DP_SAVEASROW].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_saveasrow</b> <i>(directory/)filename dp_index</i>";
	commands[CMD_DP_SAVEASROW].limits = { { Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SAVEASROW].descr = "[tc0,0.5,0.5,1/tc]Save specified dp array in filename (.txt termination by default) as a single row with tab-spaced values. If directory not specified, the default one is used.";

	commands.insert(CMD_DP_SAVEAPPENDASROW, CommandSpecifier(CMD_DP_SAVEAPPENDASROW), "dp_saveappendasrow");
	commands[CMD_DP_SAVEAPPENDASROW].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_saveappendasrow</b> <i>(directory/)filename dp_index</i>";
	commands[CMD_DP_SAVEAPPENDASROW].limits = { { Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SAVEAPPENDASROW].descr = "[tc0,0.5,0.5,1/tc]Save specified dp array in filename (.txt termination by default) as a single row with tab-spaced values, appending to end of file. If directory not specified, the default one is used.";

	commands.insert(CMD_DP_NEWFILE, CommandSpecifier(CMD_DP_NEWFILE), "dp_newfile");
	commands[CMD_DP_NEWFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_newfile</b> <i>(directory/)filename</i>";
	commands[CMD_DP_NEWFILE].limits = { { Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_NEWFILE].descr = "[tc0,0.5,0.5,1/tc]Make new file, erasing any existing file with given name. If directory not specified, the default one is used.";

	commands.insert(CMD_DP_GETPROFILE, CommandSpecifier(CMD_DP_GETPROFILE), "dp_getprofile");
	commands[CMD_DP_GETPROFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_getprofile</b> <i>start end dp_index</i>";
	commands[CMD_DP_GETPROFILE].limits = { { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_GETPROFILE].descr = "[tc0,0.5,0.5,1/tc]Extract profile of physical quantities displayed on screen, at the current display resolution, along the line specified with given start and end cartesian absolute coordinates (m). Place profile in given dp arrays: up to 4 consecutive dp arrays are used, first for distance along line, the next 3 for physical quantity components (e.g. Mx, My, Mz) so allow space for these starting at dp_index. NOTE : if you only want to extract the profile from a single mesh, then use dp_getexactprofile instead; this command should only be used to extract profiles across multiple meshes.";
	commands[CMD_DP_GETPROFILE].unit = "m";

	commands.insert(CMD_DP_GETEXACTPROFILE, CommandSpecifier(CMD_DP_GETEXACTPROFILE), "dp_getexactprofile");
	commands[CMD_DP_GETEXACTPROFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_getexactprofile</b> <i>(meshname (quantity)) start end step dp_index (stencil)</i>";
	commands[CMD_DP_GETEXACTPROFILE].limits = { { Any(), Any() }, { Any(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { MINMESHSPACE, Any() }, { int(-1), int(MAX_ARRAYS - 1) }, { DBL3(MINMESHSPACE), Any() } };
	commands[CMD_DP_GETEXACTPROFILE].descr = "[tc0,0.5,0.5,1/tc]Extract profile of physical quantity displayed on screen, or the named quantity if given (see output of display command for possible quantities), for named mesh (focused mesh if not specified), directly from the mesh (so using the exact mesh resolution not the displayed resolution), along the line specified with given start and end relative cartesian coordinates (m) and with the given step size (m). If stencil specified - as x y z (m) - then obtain profile values using weighted averaging with stencil centered on profile point. If dp_index >= 0 then place profile in given dp arrays: up to 4 consecutive dp arrays are used, first for distance along line, the next 3 for physical quantity components (e.g. Mx, My, Mz) so allow space for these starting at dp_index. If dp_index = -1 then just average profile in internal memory, to be read out with dp_getaveragedprofile. Profile will wrap-around if exceeding mesh size.";
	commands[CMD_DP_GETEXACTPROFILE].unit = "m";

	commands.insert(CMD_DP_GETAVERAGEDPROFILE, CommandSpecifier(CMD_DP_GETAVERAGEDPROFILE), "dp_getaveragedprofile");
	commands[CMD_DP_GETAVERAGEDPROFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_getaveragedprofile</b> <i>(meshname) start end step dp_index</i>";
	commands[CMD_DP_GETAVERAGEDPROFILE].limits = { { Any(), Any() }, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { MINMESHSPACE, Any() }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_GETAVERAGEDPROFILE].descr = "[tc0,0.5,0.5,1/tc]This is used in conjuction with dp_getexactprofile. Get in dp arays starting at dp_index the averaged profile built with dp_getexactprofile. start, end and step are the same parameters used by dp_getexactprofile. Calling this command causes the average counter to be reset to zero. Thus a sequence consists of 1) dp_getexactprofile executed as many times as needed (e.g. to average stochasticity), called with dp_index = -1, 2) dp_getaveragedprofile to read out the profile and reset averaging counter, ready for next batch of dp_getexactprofile calls.";
	commands[CMD_DP_GETAVERAGEDPROFILE].unit = "m";
	commands[CMD_DP_GETAVERAGEDPROFILE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>number of averages</i>";

	commands.insert(CMD_AVERAGEMESHRECT, CommandSpecifier(CMD_AVERAGEMESHRECT), "averagemeshrect");
	commands[CMD_AVERAGEMESHRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>averagemeshrect</b> <i>(meshname (quantity)) (rectangle (dp_index))</i>";
	commands[CMD_AVERAGEMESHRECT].descr = "[tc0,0.5,0.5,1/tc]Calculate the average value depending on currently displayed quantities, or the named quantity if given (see output of display command for possible quantities). The rectangle is specified in relative coordinates to the named mesh (focused mesh if not specified); if not specified average the entire focused mesh. If dp_index is specified then also append value to dp array.";
	commands[CMD_AVERAGEMESHRECT].limits = { { Any(), Any() }, { Any(), Any() }, { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_AVERAGEMESHRECT].unit = "m";
	commands[CMD_AVERAGEMESHRECT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_GETVALUE, CommandSpecifier(CMD_GETVALUE), "getvalue");
	commands[CMD_GETVALUE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>getvalue</b> <i>(meshname (quantity)) position</i>";
	commands[CMD_GETVALUE].descr = "[tc0,0.5,0.5,1/tc]Read value at position depending on currently displayed quantites, or the named quantity if given (see output of display command for possible quantities), in named mesh (focused mesh if not specified).";
	commands[CMD_GETVALUE].limits = { { Any(), Any() },  { Any(), Any() }, { DBL3(), DBL3(MAXSIMSPACE) } };
	commands[CMD_GETVALUE].unit = "m";
	commands[CMD_GETVALUE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";
	
	commands.insert(CMD_DP_TOPOCHARGE, CommandSpecifier(CMD_DP_TOPOCHARGE), "dp_topocharge");
	commands[CMD_DP_TOPOCHARGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_topocharge</b> <i>(meshname) (x y radius)</i>";
	commands[CMD_DP_TOPOCHARGE].descr = "[tc0,0.5,0.5,1/tc]Calculate the topological charge for the given mesh (must be magnetic; focused mesh if not specified), optionally in the given circle with radius and centered at x y (relative values). Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI.";
	commands[CMD_DP_TOPOCHARGE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Q</i> - the calculated topological charge.";
	commands[CMD_DP_TOPOCHARGE].limits = { { Any(), Any() }, {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_TOPOCHARGE].unit = "m";

	commands.insert(CMD_DP_COUNTSKYRMIONS, CommandSpecifier(CMD_DP_COUNTSKYRMIONS), "dp_countskyrmions");
	commands[CMD_DP_COUNTSKYRMIONS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_countskyrmions</b> <i>(meshname) (x y radius)</i>";
	commands[CMD_DP_COUNTSKYRMIONS].descr = "[tc0,0.5,0.5,1/tc]Calculate the number of skyrmions for the given mesh (must be magnetic; focused mesh if not specified), optionally in the given circle with radius and centered at x y (relative values). Use Qmag = Integral(|m.(dm/dx x dm/dy)| dxdy) / 4PI.";
	commands[CMD_DP_COUNTSKYRMIONS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Q</i> - the calculated topological charge.";
	commands[CMD_DP_COUNTSKYRMIONS].limits = { { Any(), Any() }, {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_COUNTSKYRMIONS].unit = "m";

	commands.insert(CMD_DP_HISTOGRAM, CommandSpecifier(CMD_DP_HISTOGRAM), "dp_histogram");
	commands[CMD_DP_HISTOGRAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_histogram</b> <i>(meshname) dp_x dp_y (cx cy cz (numbins min max))</i>";
	commands[CMD_DP_HISTOGRAM].descr = "[tc0,0.5,0.5,1/tc]Calculate a histogram with given number of bins, minimum and maximum bin values, from the magnetization magnitude of the given mesh (must be magnetic; focused mesh if not specified). Save histogram in dp arrays at dp_x, dp_y.If cx cy cz values given, then first average in macrocells containing (cx, cy, cz) individual cells. If histogram parameters not given use 100 bins with minimum and maximum magnetization magnitude values.";
	commands[CMD_DP_HISTOGRAM].limits = { { Any(), Any() }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {INT3(1), Any()}, {int(2), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_HISTOGRAM].unit = "A/m";
	
	commands.insert(CMD_DP_THAVHISTOGRAM, CommandSpecifier(CMD_DP_THAVHISTOGRAM), "dp_thavhistogram");
	commands[CMD_DP_THAVHISTOGRAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_thavhistogram</b> <i>(meshname) dp_x dp_y cx cy cz (numbins min max)</i>";
	commands[CMD_DP_THAVHISTOGRAM].descr = "[tc0,0.5,0.5,1/tc]Calculate a histogram with given number of bins, using thermal averaging in each macrocell containing (cx, cy, cz) individual cells, minimum and maximum bin values, from the magnetization magnitude of the focused mesh (must be magnetic). Save histogram in dp arrays at dp_x, dp_y. If histogram parameters not given use 100 bins with minimum and maximum magnetization magnitude values.";
	commands[CMD_DP_THAVHISTOGRAM].limits = { { Any(), Any() }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {INT3(1), Any()}, {int(2), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_THAVHISTOGRAM].unit = "A/m";

	commands.insert(CMD_DP_ANGHISTOGRAM, CommandSpecifier(CMD_DP_ANGHISTOGRAM), "dp_anghistogram");
	commands[CMD_DP_ANGHISTOGRAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_anghistogram</b> <i>(meshname) dp_x dp_y (cx cy cz (nx ny nz (numbins min max)))</i>";
	commands[CMD_DP_ANGHISTOGRAM].descr = "[tc0,0.5,0.5,1/tc]Calculate an angular deviation histogram with given number of bins, minimum and maximum bin values, from the magnetization of the given mesh (must be magnetic; focused mesh if not specified). Save histogram in dp arrays at dp_x, dp_y. If cx cy cz values given, then first average in macrocells containing (cx, cy, cz) individual cells. If unit vector direction nx ny nz not given, then angular deviation is calculated from the average magnetization direction. If histogram parameters not given use 100 bins with minimum and maximum magnetization magnitude values.";
	commands[CMD_DP_ANGHISTOGRAM].limits = { { Any(), Any() }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {INT3(1), Any()}, {DBL3(-1, -1, -1), DBL3(1, 1, 1)}, {int(2), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_ANGHISTOGRAM].unit = "A/m";

	commands.insert(CMD_DP_THAVANGHISTOGRAM, CommandSpecifier(CMD_DP_THAVANGHISTOGRAM), "dp_thavanghistogram");
	commands[CMD_DP_THAVANGHISTOGRAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_thavanghistogram</b> <i>(meshname) dp_x dp_y cx cy cz (nx ny nz (numbins min max))</i>";
	commands[CMD_DP_THAVANGHISTOGRAM].descr = "[tc0,0.5,0.5,1/tc]Calculate an angular deviation histogram with given number of bins, using thermal averaging in each macrocell containing (cx, cy, cz) individual cells, minimum and maximum bin values, from the magnetization of the given mesh (must be magnetic; focused mesh if not specified). Save histogram in dp arrays at dp_x, dp_y. If unit vector direction nx ny nz not given, then angular deviation is calculated from the average magnetization direction. If histogram parameters not given use 100 bins with minimum and maximum magnetization magnitude values.";
	commands[CMD_DP_THAVANGHISTOGRAM].limits = { { Any(), Any() }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {INT3(1), Any()}, {DBL3(-1, -1, -1), DBL3(1, 1, 1)}, {int(2), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_THAVANGHISTOGRAM].unit = "A/m";

	commands.insert(CMD_DP_HISTOGRAM2, CommandSpecifier(CMD_DP_HISTOGRAM2), "dp_histogram2");
	commands[CMD_DP_HISTOGRAM2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_histogram2</b> <i>(meshname) dp_x dp_y (numbins min max M2 deltaM2)</i>";
	commands[CMD_DP_HISTOGRAM2].descr = "[tc0,0.5,0.5,1/tc]Calculate a histogram for the given 2-sublattice mesh (focused mesh if not specified) with given number of bins, minimum and maximum bin values for sub-lattice A, if the corresponding magnetization magnitude in sub-lattice B equals M2 within the given deltaM2. Save histogram in dp arrays at dp_x, dp_y. If histogram parameters not given use 100 bins with minimum and maximum magnetization magnitude values, with M2 set to MeB and deltaM2 set 0.01*MeB respectively.";
	commands[CMD_DP_HISTOGRAM2].limits = { { Any(), Any() }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {int(2), Any()}, {double(0), Any()}, {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_HISTOGRAM2].unit = "A/m";

	commands.insert(CMD_DP_APPEND, CommandSpecifier(CMD_DP_APPEND), "dp_append");
	commands[CMD_DP_APPEND].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_append</b> <i>dp_original dp_new</i>";
	commands[CMD_DP_APPEND].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_APPEND].descr = "[tc0,0.5,0.5,1/tc]Append data from dp_new to the end of dp_original.";
	
	commands.insert(CMD_DP_SEQUENCE, CommandSpecifier(CMD_DP_SEQUENCE), "dp_sequence");
	commands[CMD_DP_SEQUENCE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_sequence</b> <i>dp_index start_value increment points</i>";
	commands[CMD_DP_SEQUENCE].limits = { { int(0), int(MAX_ARRAYS - 1) },{ Any(), Any() },{ Any(), Any() }, { int(0), Any() } };
	commands[CMD_DP_SEQUENCE].descr = "[tc0,0.5,0.5,1/tc]Generate a sequence of data points in dp_index from start_value using increment.";

	commands.insert(CMD_DP_RAREFY, CommandSpecifier(CMD_DP_RAREFY), "dp_rarefy");
	commands[CMD_DP_RAREFY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_rarefy</b> <i>dp_in dp_out (skip)</i>";
	commands[CMD_DP_RAREFY].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) }, { int(1), Any() } };
	commands[CMD_DP_RAREFY].descr = "[tc0,0.5,0.5,1/tc]Pick elements from dp_in using the skip value (1 by default) and set them in dp_out; e.g. with skip = 2 every 3rd data point is picked. The default skip = 1 picks every other point.";
	
	commands.insert(CMD_DP_EXTRACT, CommandSpecifier(CMD_DP_EXTRACT), "dp_extract");
	commands[CMD_DP_EXTRACT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_extract</b> <i>dp_in dp_out start_index (length)</i>";
	commands[CMD_DP_EXTRACT].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), Any() },{ int(1), Any() } };
	commands[CMD_DP_EXTRACT].descr = "[tc0,0.5,0.5,1/tc]From dp_in array extract a number of points - length - starting at start_index, and place them in dp_out.";

	commands.insert(CMD_DP_ERASE, CommandSpecifier(CMD_DP_ERASE), "dp_erase");
	commands[CMD_DP_ERASE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_erase</b> <i>dp_index start_index length</i>";
	commands[CMD_DP_ERASE].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), Any() },{ int(0), Any() },{ int(1), Any() } };
	commands[CMD_DP_ERASE].descr = "[tc0,0.5,0.5,1/tc]From dp_index array erase a number of points - length - starting at start_index.";

	commands.insert(CMD_DP_ADD, CommandSpecifier(CMD_DP_ADD), "dp_add");
	commands[CMD_DP_ADD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_add</b> <i>dp_source value (dp_dest)</i>";
	commands[CMD_DP_ADD].limits = { { int(0), int(MAX_ARRAYS - 1) }, { Any(), Any() }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_ADD].descr = "[tc0,0.5,0.5,1/tc]Add value to dp array and place it in destination (or at same position if destination not specified).";

	commands.insert(CMD_DP_SUB, CommandSpecifier(CMD_DP_SUB), "dp_sub");
	commands[CMD_DP_SUB].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_sub</b> <i>dp_source value (dp_dest)</i>";
	commands[CMD_DP_SUB].limits = { { int(0), int(MAX_ARRAYS - 1) },{ Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SUB].descr = "[tc0,0.5,0.5,1/tc]Subtract value from dp array and place it in destination (or at same position if destination not specified).";

	commands.insert(CMD_DP_MUL, CommandSpecifier(CMD_DP_MUL), "dp_mul");
	commands[CMD_DP_MUL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_mul</b> <i>dp_source value (dp_dest)</i>";
	commands[CMD_DP_MUL].limits = { { int(0), int(MAX_ARRAYS - 1) },{ Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_MUL].descr = "[tc0,0.5,0.5,1/tc]Multiply dp array with value and place it in destination (or at same position if destination not specified).";

	commands.insert(CMD_DP_DIV, CommandSpecifier(CMD_DP_DIV), "dp_div");
	commands[CMD_DP_DIV].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_div</b> <i>dp_source value (dp_dest)</i>";
	commands[CMD_DP_DIV].limits = { { int(0), int(MAX_ARRAYS - 1) },{ Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_DIV].descr = "[tc0,0.5,0.5,1/tc]Divide dp array by value and place it in destination (or at same position if destination not specified).";

	commands.insert(CMD_DP_POW, CommandSpecifier(CMD_DP_POW), "dp_pow");
	commands[CMD_DP_POW].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_pow</b> <i>dp_source exponent (dp_dest)</i>";
	commands[CMD_DP_POW].limits = { { int(0), int(MAX_ARRAYS - 1) },{ Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_POW].descr = "[tc0,0.5,0.5,1/tc]Exponentiate source aray values and place it in destination (or at same position if destination not specified).";

	commands.insert(CMD_DP_DOTPROD, CommandSpecifier(CMD_DP_DOTPROD), "dp_dotprod");
	commands[CMD_DP_DOTPROD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_dotprod</b> <i>dp_vector ux uy uz dp_out</i>";
	commands[CMD_DP_DOTPROD].limits = { { int(0), int(MAX_ARRAYS - 3) },{ Any(), Any() },{ Any(), Any() },{ Any(), Any() },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_DOTPROD].descr = "[tc0,0.5,0.5,1/tc]Take dot product of (ux, uy, uz) with vectors in dp arrays dp_vector, dp_vector + 1, dp_vector + 2 and place result in dp_out.";

	commands.insert(CMD_DP_ADDDP, CommandSpecifier(CMD_DP_ADDDP), "dp_adddp");
	commands[CMD_DP_ADDDP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_adddp</b> <i>dp_x1 dp_x2 dp_dest</i>";
	commands[CMD_DP_ADDDP].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_ADDDP].descr = "[tc0,0.5,0.5,1/tc]Add dp arrays : dp_dest = dp_x1 + dp_x2";

	commands.insert(CMD_DP_SUBDP, CommandSpecifier(CMD_DP_SUBDP), "dp_subdp");
	commands[CMD_DP_SUBDP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_subdp</b> <i>dp_x1 dp_x2 dp_dest</i>";
	commands[CMD_DP_SUBDP].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SUBDP].descr = "[tc0,0.5,0.5,1/tc]Subtract dp arrays : dp_dest = dp_x1 - dp_x2";

	commands.insert(CMD_DP_MULDP, CommandSpecifier(CMD_DP_MULDP), "dp_muldp");
	commands[CMD_DP_MULDP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_muldp</b> <i>dp_x1 dp_x2 dp_dest</i>";
	commands[CMD_DP_MULDP].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_MULDP].descr = "[tc0,0.5,0.5,1/tc]Multiply dp arrays : dp_dest = dp_x1 * dp_x2";

	commands.insert(CMD_DP_DIVDP, CommandSpecifier(CMD_DP_DIVDP), "dp_divdp");
	commands[CMD_DP_DIVDP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_divdp</b> <i>dp_x1 dp_x2 dp_dest</i>";
	commands[CMD_DP_DIVDP].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_DIVDP].descr = "[tc0,0.5,0.5,1/tc]Divide dp arrays : dp_dest = dp_x1 / dp_x2";

	commands.insert(CMD_DP_DOTPRODDP, CommandSpecifier(CMD_DP_DOTPRODDP), "dp_dotproddp");
	commands[CMD_DP_DOTPRODDP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_dotproddp</b> <i>dp_x1 dp_x2</i>";
	commands[CMD_DP_DOTPRODDP].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_DOTPRODDP].descr = "[tc0,0.5,0.5,1/tc]Take dot product of dp arrays : value = dp_x1.dp_x2";

	commands.insert(CMD_DP_MINMAX, CommandSpecifier(CMD_DP_MINMAX), "dp_minmax");
	commands[CMD_DP_MINMAX].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_minmax</b> <i>dp_index</i>";
	commands[CMD_DP_MINMAX].limits = { { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_MINMAX].descr = "[tc0,0.5,0.5,1/tc]Obtain absolute minimum and maximum values, together with their index position.";
	commands[CMD_DP_MINMAX].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>min_value min_index max_value max_index</i>.";

	commands.insert(CMD_DP_MEAN, CommandSpecifier(CMD_DP_MEAN), "dp_mean");
	commands[CMD_DP_MEAN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_mean</b> <i>dp_index (exclusion_ratio)</i>";
	commands[CMD_DP_MEAN].limits = { { int(0), int(MAX_ARRAYS - 1) }, { double(0), Any() } };
	commands[CMD_DP_MEAN].descr = "[tc0,0.5,0.5,1/tc]Obtain mean value with standard deviation. If exclusion_ratio (>0) included, then exclude any points which are greater than exclusion_ratio normalised distance away from the mean.";
	commands[CMD_DP_MEAN].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>mean stdev</i>.";

	commands.insert(CMD_DP_SUM, CommandSpecifier(CMD_DP_SUM), "dp_sum");
	commands[CMD_DP_SUM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_sum</b> <i>dp_index</i>";
	commands[CMD_DP_SUM].limits = { { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_SUM].descr = "[tc0,0.5,0.5,1/tc]Obtain sum of all elements in dp array.";
	commands[CMD_DP_SUM].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>sum</i>.";

	commands.insert(CMD_DP_CHUNKEDSTD, CommandSpecifier(CMD_DP_CHUNKEDSTD), "dp_chunkedstd");
	commands[CMD_DP_CHUNKEDSTD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_chunkedstd</b> <i>dp_index chunk</i>";
	commands[CMD_DP_CHUNKEDSTD].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), Any() } };
	commands[CMD_DP_CHUNKEDSTD].descr = "[tc0,0.5,0.5,1/tc]Obtain stadard deviation from every chunk number of elements in dp array, then average all the chunk std value to obtain final std. In the special case where chunk is the number of elements in the dp array (set chunk = 0 for this case), the chunked std is the same as the usual std on the mean.";
	commands[CMD_DP_CHUNKEDSTD].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>std</i>.";

	commands.insert(CMD_DP_GETAMPLI, CommandSpecifier(CMD_DP_GETAMPLI), "dp_getampli");
	commands[CMD_DP_GETAMPLI].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_getampli</b> <i>dp_source pointsPeriod</i>";
	commands[CMD_DP_GETAMPLI].limits = { { int(0), int(MAX_ARRAYS - 1) }, {int(1), Any()} };
	commands[CMD_DP_GETAMPLI].descr = "[tc0,0.5,0.5,1/tc]Obtain maximum amplitude obtained every pointsPeriod points.";
	commands[CMD_DP_GETAMPLI].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>amplitude</i>.";

	commands.insert(CMD_DP_LINREG, CommandSpecifier(CMD_DP_LINREG), "dp_linreg");
	commands[CMD_DP_LINREG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_linreg</b> <i>dp_index_x dp_index_y (dp_index_z dp_index_out)</i>";
	commands[CMD_DP_LINREG].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_LINREG].descr = "[tc0,0.5,0.5,1/tc]Fit using linear regression to obtain gradient and intercept with their uncertainties. If dp_index_z is specified multiple linear regressions are performed on adjacent data points with same z value; output in 5 dp arrays starting at dp_index_out as: z g g_err c c_err.";
	commands[CMD_DP_LINREG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>g g_err c c_err</i>.";

	commands.insert(CMD_DP_COERCIVITY, CommandSpecifier(CMD_DP_COERCIVITY), "dp_coercivity");
	commands[CMD_DP_COERCIVITY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_coercivity</b> <i>dp_index_x dp_index_y</i>";
	commands[CMD_DP_COERCIVITY].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_COERCIVITY].descr = "[tc0,0.5,0.5,1/tc]Obtain coercivity from x-y data: find first crossings of x axis in the two possible directions, with uncertainty obtained from step size.";
	commands[CMD_DP_COERCIVITY].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Hc_up Hc_up_err- Hc_up_err+ Hc_dn Hc_dn_err- Hc_dn_err+</i>.";

	commands.insert(CMD_DP_REMANENCE, CommandSpecifier(CMD_DP_REMANENCE), "dp_remanence");
	commands[CMD_DP_REMANENCE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_remanence</b> <i>dp_index_x dp_index_y</i>";
	commands[CMD_DP_REMANENCE].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_REMANENCE].descr = "[tc0,0.5,0.5,1/tc]Obtain remanence from x-y data: find values at zero field in the two possible directions.";
	commands[CMD_DP_REMANENCE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Mr_up Mr_dn</i>.";
	
	commands.insert(CMD_DP_COMPLETEHYSTERLOOP, CommandSpecifier(CMD_DP_COMPLETEHYSTERLOOP), "dp_completehysteresis");
	commands[CMD_DP_COMPLETEHYSTERLOOP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_completehysteresis</b> <i>dp_index_x dp_index_y</i>";
	commands[CMD_DP_COMPLETEHYSTERLOOP].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_COMPLETEHYSTERLOOP].descr = "[tc0,0.5,0.5,1/tc]For a hysteresis loop with only one branch continue it by constructing the other direction branch (invert both x and y data and add it in continuation) - use only with hysteresis loops which are expected to be symmetric.";

	commands.insert(CMD_DP_DUMPTDEP, CommandSpecifier(CMD_DP_DUMPTDEP), "dp_dumptdep");
	commands[CMD_DP_DUMPTDEP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_dumptdep</b> <i>meshname paramname max_temperature dp_index</i>";
	commands[CMD_DP_DUMPTDEP].limits = { { Any(), Any() }, { Any(), Any() }, { double(1.0), double(MAX_TEMPERATURE) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_DUMPTDEP].descr = "[tc0,0.5,0.5,1/tc]Get temperature dependence of named parameter from named mesh up to max_temperature, at dp_index - temperature scaling values obtained.";

	commands.insert(CMD_DP_FITLORENTZ, CommandSpecifier(CMD_DP_FITLORENTZ), "dp_fitlorentz");
	commands[CMD_DP_FITLORENTZ].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitlorentz</b> <i>dp_x dp_y</i>";
	commands[CMD_DP_FITLORENTZ].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_FITLORENTZ].descr = "[tc0,0.5,0.5,1/tc]Fit Lorentz peak function to x y data : f(x) = y0 + S dH / (4(x-H0)^2 + dH^2).";
	commands[CMD_DP_FITLORENTZ].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>S, H0, dH, y0, std_S, std_H0, std_dH, std_y0</i>.";

	commands.insert(CMD_DP_FITLORENTZ2, CommandSpecifier(CMD_DP_FITLORENTZ2), "dp_fitlorentz2");
	commands[CMD_DP_FITLORENTZ2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitlorentz2</b> <i>dp_x dp_y</i>";
	commands[CMD_DP_FITLORENTZ2].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_FITLORENTZ2].descr = "[tc0,0.5,0.5,1/tc]Fit Lorentz peak function with both symmetric and asymmetric parts to x y data : f(x) = y0 + S (dH + A * (x - H0)) / (4(x-H0)^2 + dH^2).";
	commands[CMD_DP_FITLORENTZ2].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>S, A, H0, dH, y0, std_S, std_A, std_H0, std_dH, std_y0</i>.";

	commands.insert(CMD_DP_FITSKYRMION, CommandSpecifier(CMD_DP_FITSKYRMION), "dp_fitskyrmion");
	commands[CMD_DP_FITSKYRMION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitskyrmion</b> <i>dp_x dp_y</i>";
	commands[CMD_DP_FITSKYRMION].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_FITSKYRMION].descr = "[tc0,0.5,0.5,1/tc]Fit skyrmion z component to obtain radius and center position : Mz(x) = Ms * cos(2*arctan(sinh(R/w)/sinh((x-x0)/w))).";
	commands[CMD_DP_FITSKYRMION].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>R, x0, Ms, w, std_R, std_x0, std_Ms, std_w</i>.";

	commands.insert(CMD_DP_FITDW, CommandSpecifier(CMD_DP_FITDW), "dp_fitdw");
	commands[CMD_DP_FITDW].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitdw</b> <i>dp_x dp_y</i>";
	commands[CMD_DP_FITDW].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_FITDW].descr = "[tc0,0.5,0.5,1/tc]Fit domain wall tanh component to obtain width and center position : M(x) = A * tanh(PI * (x - x0) / D).";
	commands[CMD_DP_FITDW].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>D, x0, A, std_D, std_x0, std_A</i>.";

	commands.insert(CMD_DP_FITSTT, CommandSpecifier(CMD_DP_FITSTT), "dp_fitstt");
	commands[CMD_DP_FITSTT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitstt</b> <i>(meshname) (rectangle)</i>";
	commands[CMD_DP_FITSTT].limits = { { Any(), Any() }, {Rect(), Any()} };
	commands[CMD_DP_FITSTT].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity). The fitting is done inside the specified rectangle for the given mesh (focused mesh if not specified), with the rectangle specified using relative coordinates as sx sy sz ex ey ez (entire mesh if not specified). The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not.";
	commands[CMD_DP_FITSTT].unit = "m";
	commands[CMD_DP_FITSTT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>P, beta, std_P, std_beta, Rsq</i>.";

	commands.insert(CMD_DP_FITSOT, CommandSpecifier(CMD_DP_FITSOT), "dp_fitsot");
	commands[CMD_DP_FITSOT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitsot</b> <i>(meshname) (hm_mesh) (rectangle)</i>";
	commands[CMD_DP_FITSOT].limits = { { Any(), Any() }, { Any(), Any() }, {Rect(), Any()} };
	commands[CMD_DP_FITSOT].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent interfacial spin torque using SOT with fitting parameters SHAeff and flST (field-like torque coefficient). hm_mesh specifies the heavy metal mesh from which to obtain the current density - if not specified, then the current density in meshname is used instead. The fitting is done inside the specified rectangle for the given mesh (focused mesh if not specified), with the rectangle specified using relative coordinates as sx sy sz ex ey ez (entire mesh if not specified). This mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and Tsi to have been computed.";
	commands[CMD_DP_FITSOT].unit = "m";
	commands[CMD_DP_FITSOT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>SHAeff, flST, std_SHAeff, std_flST, Rsq</i>.";

	commands.insert(CMD_DP_FITSOTSTT, CommandSpecifier(CMD_DP_FITSOTSTT), "dp_fitsotstt");
	commands[CMD_DP_FITSOTSTT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitsotstt</b> <i>(meshname) (hm_mesh) (rectangle)</i>";
	commands[CMD_DP_FITSOTSTT].limits = { {Any(), Any()}, {Any(), Any()}, {Rect(), Any()} };
	commands[CMD_DP_FITSOTSTT].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity), and simultaneously also using SOT with fitting parameters SHAeff and flST (field-like torque coefficient). hm_mesh specifies the heavy metal mesh from which to obtain the current density for SOT. The fitting is done inside the specified rectangle for the given mesh (focused mesh if not specified), with the rectangle specified using relative coordinates as sx sy sz ex ey ez (entire mesh if not specified). The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not.";
	commands[CMD_DP_FITSOTSTT].unit = "m";
	commands[CMD_DP_FITSOTSTT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>SHAeff, flST, P, beta, std_SHAeff, std_flST, std_P, std_beta, Rsq</i>.";

	commands.insert(CMD_DP_CALCSOT, CommandSpecifier(CMD_DP_CALCSOT), "dp_calcsot");
	commands[CMD_DP_CALCSOT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_calcsot</b> <i>hm_mesh fm_mesh</i>";
	commands[CMD_DP_CALCSOT].limits = { {Any(), Any()}, {Any(), Any()} };
	commands[CMD_DP_CALCSOT].descr = "[tc0,0.5,0.5,1/tc]For the given heavy metal and ferromagnetic meshes calculate the expected effective spin Hall angle and field-like torque coefficient according to analytical equations (see manual).";
	commands[CMD_DP_CALCSOT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>SHAeff, flST</i>.";

	commands.insert(CMD_DP_FITNONADIABATIC, CommandSpecifier(CMD_DP_FITNONADIABATIC), "dp_fitnonadiabatic");
	commands[CMD_DP_FITNONADIABATIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitnonadiabatic</b> <i>(meshname) (abs_err Rsq T_ratio (stencil))</i>";
	commands[CMD_DP_FITNONADIABATIC].limits = { {Any(), Any()}, {double(0.0), Any()}, {double(0.0), Any()}, {double(0.0), Any()}, {int(3), Any()} };
	commands[CMD_DP_FITNONADIABATIC].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity) using a given square in-plane stencil (default size 3) in order to extract the spatial variation of beta. Cut-off values for absolute fitting error (default 0.1), Rsq measure (default 0.9), and normalized torque magnitude (default 0.1) can be set - value of zero disables cutoff. The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not. Output available in Cust_S.";

	commands.insert(CMD_DP_FITADIABATIC, CommandSpecifier(CMD_DP_FITADIABATIC), "dp_fitadiabatic");
	commands[CMD_DP_FITADIABATIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitadiabatic</b> <i>(meshname) (abs_err Rsq T_ratio (stencil))</i>";
	commands[CMD_DP_FITADIABATIC].limits = { {Any(), Any()}, {double(0.0), Any()}, {double(0.0), Any()}, {double(0.0), Any()}, {int(3), Any()} };
	commands[CMD_DP_FITADIABATIC].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity) using a given square in-plane stencil (default size 3) in order to extract the spatial variation of P. Cut-off values for absolute fitting error (default 0.1), Rsq measure (default 0.9), and normalized torque magnitude (default 0.1) can be set - value of zero disables cutoff. The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not. Output available in Cust_S.";

	commands.insert(CMD_DP_CALCTOPOCHARGEDENSITY, CommandSpecifier(CMD_DP_CALCTOPOCHARGEDENSITY), "dp_calctopochargedensity");
	commands[CMD_DP_CALCTOPOCHARGEDENSITY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_calctopochargedensity</b> <i>(meshname)</i>";
	commands[CMD_DP_CALCTOPOCHARGEDENSITY].limits = { {Any(), Any()} };
	commands[CMD_DP_CALCTOPOCHARGEDENSITY].descr = "[tc0,0.5,0.5,1/tc]Calculate topological charge density spatial dependence for the given mesh (must be magnetic; focused mesh if not specified). Output available in Cust_S.";

	commands.insert(CMD_DP_REPLACEREPEATS, CommandSpecifier(CMD_DP_REPLACEREPEATS), "dp_replacerepeats");
	commands[CMD_DP_REPLACEREPEATS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_replacerepeats</b> <i>dp_index (dp_index_out)</i>";
	commands[CMD_DP_REPLACEREPEATS].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_REPLACEREPEATS].descr = "[tc0,0.5,0.5,1/tc]Replace repeated points from data in dp_index using linear interpolation: if two adjacent sets of repeated points found, replace repeats between the mid-points of the sets. If dp_index_out not specified then processed data overwrites dp_index.";

	commands.insert(CMD_DP_REMOVEOFFSET, CommandSpecifier(CMD_DP_REMOVEOFFSET), "dp_removeoffset");
	commands[CMD_DP_REMOVEOFFSET].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_removeoffset</b> <i>dp_index (dp_index_out)</i>";
	commands[CMD_DP_REMOVEOFFSET].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_REMOVEOFFSET].descr = "[tc0,0.5,0.5,1/tc]Subtract the first point (the offset) from all the points in dp_index. If dp_index_out not specified then processed data overwrites dp_index.";

	commands.insert(CMD_CARTESIANTOPOLAR, CommandSpecifier(CMD_CARTESIANTOPOLAR), "dp_cartesiantopolar");
	commands[CMD_CARTESIANTOPOLAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_cartesiantopolar</b> <i>dp_in_x dp_in_y (dp_out_r dp_out_theta)</i>";
	commands[CMD_CARTESIANTOPOLAR].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_CARTESIANTOPOLAR].descr = "[tc0,0.5,0.5,1/tc]Convert from Cartesian coordinates (x,y) to polar (r, theta).";

	commands.insert(CMD_DP_SMOOTH, CommandSpecifier(CMD_DP_SMOOTH), "dp_smooth");
	commands[CMD_DP_SMOOTH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_smooth</b> <i>dp_in dp_out window_size</i>";
	commands[CMD_DP_SMOOTH].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) }, { int(2), Any() } };
	commands[CMD_DP_SMOOTH].descr = "[tc0,0.5,0.5,1/tc]Smooth data in dp_in using nearest-neighbor averaging with given window size, and place result in dp_out (must be different).";

	commands.insert(CMD_DP_MONOTONIC, CommandSpecifier(CMD_DP_MONOTONIC), "dp_monotonic");
	commands[CMD_DP_MONOTONIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_monotonic</b> <i>dp_in_x dp_in_y dp_out_x dp_out_y</i>";
	commands[CMD_DP_MONOTONIC].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_MONOTONIC].descr = "[tc0,0.5,0.5,1/tc]From input x-y data extract monotonic sequence and place it in output x y arrays.";

	commands.insert(CMD_DP_CROSSINGSHISTOGRAM, CommandSpecifier(CMD_DP_CROSSINGSHISTOGRAM), "dp_crossingshistogram");
	commands[CMD_DP_CROSSINGSHISTOGRAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_crossingshistogram</b> <i>dp_in_x dp_in_y dp_level dp_counts (steps)</i>";
	commands[CMD_DP_CROSSINGSHISTOGRAM].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(1), Any() } };
	commands[CMD_DP_CROSSINGSHISTOGRAM].descr = "[tc0,0.5,0.5,1/tc]From input x-y data build a histogram of number of times x-y data crosses a given line (up or down). The line varies between minimum and maximum of y data in given number of steps (100 by default). Output the line values in dp_level with corresponding number of crossings in dp_counts.";

	commands.insert(CMD_DP_CROSSINGSFREQUENCY, CommandSpecifier(CMD_DP_CROSSINGSFREQUENCY), "dp_crossingsfrequency");
	commands[CMD_DP_CROSSINGSFREQUENCY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_crossingsfrequency</b> <i>dp_in_x dp_in_y dp_level dp_freq_up dp_freq_dn (steps)</i>";
	commands[CMD_DP_CROSSINGSFREQUENCY].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(1), Any() } };
	commands[CMD_DP_CROSSINGSFREQUENCY].descr = "[tc0,0.5,0.5,1/tc]From input x-y data build a histogram of average frequency the x-y data crosses a given line (up and down, separated). The line varies between minimum and maximum of y data in given number of steps (100 by default). Output the line values in dp_level with corresponding crossings frequencies in dp_freq_up and dp_freq_dn.";

	commands.insert(CMD_DP_PEAKSFREQUENCY, CommandSpecifier(CMD_DP_PEAKSFREQUENCY), "dp_peaksfrequency");
	commands[CMD_DP_PEAKSFREQUENCY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_peaksfrequency</b> <i>dp_in_x dp_in_y dp_level dp_freq (steps)</i>";
	commands[CMD_DP_PEAKSFREQUENCY].limits = { { int(0), int(MAX_ARRAYS - 1) },{ int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, { int(1), Any() } };
	commands[CMD_DP_PEAKSFREQUENCY].descr = "[tc0,0.5,0.5,1/tc]From input x-y data build a histogram of average frequency of peaks in the x-y data in bands given by the number of steps. The bands vary between minimum and maximum of y data in given number of steps (100 by default). Output the line values in dp_level with corresponding peak frequencies in dp_freq.";

	commands.insert(CMD_TEST, CommandSpecifier(CMD_TEST), "test");
	commands[CMD_TEST].usage = "[tc0,0.5,0,1/tc]USAGE : ";
	commands[CMD_TEST].descr = "[tc0,0.5,0.5,1/tc].";

	//---------------------------------------------------------------- DATA OUTPUTS

	//Build data descriptors, holding labels, handles and type specifier
	dataDescriptor.push_back("sstep", DatumSpecifier("Stage, Step : ", 2), DATA_STAGESTEP);
	dataDescriptor.push_back("time", DatumSpecifier("Time : ", 1, "s"), DATA_TIME);
	dataDescriptor.push_back("stime", DatumSpecifier("Stage time : ", 1, "s"), DATA_STAGETIME);
	dataDescriptor.push_back("iter", DatumSpecifier("Iterations : ", 1), DATA_ITERATIONS);
	dataDescriptor.push_back("siter", DatumSpecifier("Stage Iterations : ", 1), DATA_SITERATIONS);
	dataDescriptor.push_back("dt", DatumSpecifier("ODE dT : ", 1, "s"), DATA_DT);
	dataDescriptor.push_back("heat_dT", DatumSpecifier("heat dT : ", 1, "s"), DATA_HEATDT);
	dataDescriptor.push_back("mxh", DatumSpecifier("|mxh| : ", 1), DATA_MXH);
	dataDescriptor.push_back("dmdt", DatumSpecifier("|dm/dt| : ", 1), DATA_DMDT);
	dataDescriptor.push_back("Ha", DatumSpecifier("Applied Field : ", 3, "A/m", false), DATA_HA);
	dataDescriptor.push_back("<M>", DatumSpecifier("<M> : ", 3, "A/m", false, false), DATA_AVM);
	dataDescriptor.push_back("<M>th", DatumSpecifier("<M>th : ", 3, "A/m", false, false), DATA_THAVM);
	dataDescriptor.push_back("<M2>", DatumSpecifier("<M2> : ", 3, "A/m", false, false), DATA_AVM2);
	dataDescriptor.push_back("|M|mm", DatumSpecifier("|M|mm : ", 2, "A/m", false, false), DATA_M_MINMAX);
	dataDescriptor.push_back("Mx_mm", DatumSpecifier("Mx_mm : ", 2, "A/m", false, false), DATA_MX_MINMAX);
	dataDescriptor.push_back("My_mm", DatumSpecifier("My_mm : ", 2, "A/m", false, false), DATA_MY_MINMAX);
	dataDescriptor.push_back("Mz_mm", DatumSpecifier("Mz_mm : ", 2, "A/m", false, false), DATA_MZ_MINMAX);
	dataDescriptor.push_back("MCparams", DatumSpecifier("MCparams : ", 2, "", false), DATA_MONTECARLOPARAMS);
	dataDescriptor.push_back("<Jc>", DatumSpecifier("<Jc> : ", 3, "A/m^2", false, false), DATA_JC);
	dataDescriptor.push_back("<Jsx>", DatumSpecifier("<Jsx> : ", 3, "A/s", false, false), DATA_JSX);
	dataDescriptor.push_back("<Jsy>", DatumSpecifier("<Jsy> : ", 3, "A/s", false, false), DATA_JSY);
	dataDescriptor.push_back("<Jsz>", DatumSpecifier("<Jsz> : ", 3, "A/s", false, false), DATA_JSZ);
	dataDescriptor.push_back("<mxdmdt>", DatumSpecifier("<mxdmdt> : ", 3, "1/s", false, false), DATA_RESPUMP);
	dataDescriptor.push_back("<dmdt>", DatumSpecifier("<dmdt> : ", 3, "1/s", false, false), DATA_IMSPUMP);
	dataDescriptor.push_back("<mxdmdt2>", DatumSpecifier("<mxdmdt2> : ", 3, "1/s", false, false), DATA_RESPUMP2);
	dataDescriptor.push_back("<dmdt2>", DatumSpecifier("<dmdt2> : ", 3, "1/s", false, false), DATA_IMSPUMP2);
	dataDescriptor.push_back("<mxdm2dt>", DatumSpecifier("<mxdm2dt> : ", 3, "1/s", false, false), DATA_RESPUMP12);
	dataDescriptor.push_back("<m2xdmdt>", DatumSpecifier("<m2xdmdt> : ", 3, "1/s", false, false), DATA_RESPUMP21);
	dataDescriptor.push_back("<V>", DatumSpecifier("<V> : ", 1, "V", false, false), DATA_V);
	dataDescriptor.push_back("<S>", DatumSpecifier("<S> : ", 3, "A/m", false, false), DATA_S);
	dataDescriptor.push_back("<elC>", DatumSpecifier("<elC> : ", 1, "S/m", false, false), DATA_ELC);
	dataDescriptor.push_back("<T>", DatumSpecifier("<Temperature> : ", 1, "K", false, false), DATA_TEMP);
	dataDescriptor.push_back("<T_l>", DatumSpecifier("<Temperature_l> : ", 1, "K", false, false), DATA_TEMP_L);
	dataDescriptor.push_back("V", DatumSpecifier("Potential : ", 1, "V"), DATA_POTENTIAL);
	dataDescriptor.push_back("I", DatumSpecifier("I, Net I : ", 2, "A"), DATA_CURRENT);
	dataDescriptor.push_back("R", DatumSpecifier("Resistance : ", 1, "Ohm"), DATA_RESISTANCE);
	dataDescriptor.push_back("<u>", DatumSpecifier("<u> : ", 3, "m", false, false), DATA_AVU);
	dataDescriptor.push_back("<Sd>", DatumSpecifier("<Sd> : ", 3, "", false, false), DATA_AVSTRAINDIAG);
	dataDescriptor.push_back("<Sod>", DatumSpecifier("<Sod> : ", 3, "", false, false), DATA_AVSTRAINODIAG);
	dataDescriptor.push_back("e_demag", DatumSpecifier("Demag e : ", 1, "J/m3", false, false), DATA_E_DEMAG);
	dataDescriptor.push_back("e_exch", DatumSpecifier("Exchange e : ", 1, "J/m3", false, false), DATA_E_EXCH);
	dataDescriptor.push_back("t_exch", DatumSpecifier("Exchange torque : ", 3, "au", false, false), DATA_T_EXCH);
	dataDescriptor.push_back("e_surfexch", DatumSpecifier("Surface exchange e : ", 1, "J/m3", false, false), DATA_E_SURFEXCH);
	dataDescriptor.push_back("t_surfexch", DatumSpecifier("Surface exchange torque : ", 3, "au", false, false), DATA_T_SURFEXCH);
	dataDescriptor.push_back("e_zee", DatumSpecifier("Zeeman e : ", 1, "J/m3", false, false), DATA_E_ZEE);
	dataDescriptor.push_back("t_zee", DatumSpecifier("Zeeman torque : ", 3, "au", false, false), DATA_T_ZEE);
	dataDescriptor.push_back("e_stray", DatumSpecifier("Strayfield e : ", 1, "J/m3", false, false), DATA_E_STRAY);
	dataDescriptor.push_back("t_stray", DatumSpecifier("Strayfield torque : ", 1, "au", false, false), DATA_T_STRAY);
	dataDescriptor.push_back("e_mo", DatumSpecifier("Magneto-optical e : ", 1, "J/m3", false, false), DATA_E_MOPTICAL);
	dataDescriptor.push_back("e_mel", DatumSpecifier("Magnetoelastic e : ", 1, "J/m3", false, false), DATA_E_MELASTIC);
	dataDescriptor.push_back("e_anis", DatumSpecifier("Anisotropy e : ", 1, "J/m3", false, false), DATA_E_ANIS);
	dataDescriptor.push_back("t_anis", DatumSpecifier("Anisotropy torque : ", 3, "au", false, false), DATA_T_ANIS);
	dataDescriptor.push_back("e_rough", DatumSpecifier("Roughness e : ", 1, "J/m3", false, false), DATA_E_ROUGH);
	dataDescriptor.push_back("e_total", DatumSpecifier("Total e : ", 1, "J/m3", true), DATA_E_TOTAL);
	dataDescriptor.push_back("dwshift", DatumSpecifier("DW shift : ", 1, "m"), DATA_DWSHIFT);
	dataDescriptor.push_back("dwpos_x", DatumSpecifier("DW (pos, width) : ", 2, "m", false, false), DATA_DWPOS_X);
	dataDescriptor.push_back("dwpos_y", DatumSpecifier("DW (pos, width) : ", 2, "m", false, false), DATA_DWPOS_Y);
	dataDescriptor.push_back("dwpos_z", DatumSpecifier("DW (pos, width) : ", 2, "m", false, false), DATA_DWPOS_Z);
	dataDescriptor.push_back("skyshift", DatumSpecifier("Skyrmion shift : ", 2, "m", false, false), DATA_SKYSHIFT);
	dataDescriptor.push_back("skypos", DatumSpecifier("Skyrmion (pos, dia) : ", 4, "m", false, false), DATA_SKYPOS);
	dataDescriptor.push_back("Q_topo", DatumSpecifier("Topological Charge : ", 1, "", false, false), DATA_Q_TOPO);
	dataDescriptor.push_back("v_iter", DatumSpecifier("V Solver Iterations : ", 1), DATA_TRANSPORT_ITERSTOCONV);
	dataDescriptor.push_back("s_iter", DatumSpecifier("S Solver Iterations : ", 1), DATA_TRANSPORT_SITERSTOCONV);
	dataDescriptor.push_back("ts_err", DatumSpecifier("Transport Solver Error : ", 1), DATA_TRANSPORT_CONVERROR);
	dataDescriptor.push_back("TMR", DatumSpecifier("TMR : ", 1, "Ohm", false, false), DATA_TMR);
	dataDescriptor.push_back("commbuf", DatumSpecifier("Command Buffer : ", 1), DATA_COMMBUFFER);

	//---------------------------------------------------------------- MESHES

	meshtypeHandles.push_back("Ferromagnet", MESH_FERROMAGNETIC);
	meshtypeHandles.push_back("AntiFerromagnet", MESH_ANTIFERROMAGNETIC);
	meshtypeHandles.push_back("Dipole", MESH_DIPOLE);
	meshtypeHandles.push_back("Conductor", MESH_METAL);
	meshtypeHandles.push_back("Insulator", MESH_INSULATOR);
	meshtypeHandles.push_back("Atomistic", MESH_ATOM_CUBIC);

	//---------------------------------------------------------------- MODULES

	//Modules
	moduleHandles.push_back("demag_N", MOD_DEMAG_N);
	moduleHandles.push_back("demag", MOD_DEMAG);
	moduleHandles.push_back("exchange", MOD_EXCHANGE);
	moduleHandles.push_back("DMexchange", MOD_DMEXCHANGE);
	moduleHandles.push_back("iDMexchange", MOD_IDMEXCHANGE);
	moduleHandles.push_back("viDMexchange", MOD_VIDMEXCHANGE);
	moduleHandles.push_back("surfexchange", MOD_SURFEXCHANGE);
	moduleHandles.push_back("Zeeman", MOD_ZEEMAN);
	moduleHandles.push_back("moptical", MOD_MOPTICAL);
	moduleHandles.push_back("aniuni", MOD_ANIUNI);
	moduleHandles.push_back("anicubi", MOD_ANICUBI);
	moduleHandles.push_back("anibi", MOD_ANIBI);
	moduleHandles.push_back("anitens", MOD_ANITENS);
	moduleHandles.push_back("melastic", MOD_MELASTIC);
	moduleHandles.push_back("transport", MOD_TRANSPORT);
	moduleHandles.push_back("tmr", MOD_TMR);
	moduleHandles.push_back("heat", MOD_HEAT);
	moduleHandles.push_back("SOTfield", MOD_SOTFIELD);
	moduleHandles.push_back("STfield", MOD_STFIELD);
	moduleHandles.push_back("roughness", MOD_ROUGHNESS);
	moduleHandles.push_back("dipoledipole", MOD_ATOM_DIPOLEDIPOLE);
	moduleHandles.push_back("mstrayfield", MOD_STRAYFIELD_MESH);
	
	//super-mesh modules
	moduleHandles.push_back("sdemag", MODS_SDEMAG);
	moduleHandles.push_back("strayfield", MODS_STRAYFIELD);
	moduleHandles.push_back("stransport", MODS_STRANSPORT);
	moduleHandles.push_back("sheat", MODS_SHEAT);
	moduleHandles.push_back("smelastic", MODS_SMELASTIC);
	moduleHandles.push_back("Oersted", MODS_OERSTED);

	//---------------------------------------------------------------- ODEs

	//ODEs : all (micromagnetic meshes)
	odeHandles.push_back("LLG", ODE_LLG);
	odeHandles.push_back("LLGStatic", ODE_LLGSTATIC);
	odeHandles.push_back("LLGStatic-SA", ODE_LLGSTATICSA);
	odeHandles.push_back("LLG-STT", ODE_LLGSTT);
	odeHandles.push_back("LLB", ODE_LLB);
	odeHandles.push_back("LLB-STT", ODE_LLBSTT);
	odeHandles.push_back("sLLG", ODE_SLLG);
	odeHandles.push_back("sLLG-STT", ODE_SLLGSTT);
	odeHandles.push_back("sLLB", ODE_SLLB);
	odeHandles.push_back("sLLB-STT", ODE_SLLBSTT);
	odeHandles.push_back("LLG-SA", ODE_LLGSA);
	odeHandles.push_back("sLLG-SA", ODE_SLLGSA);
	odeHandles.push_back("LLB-SA", ODE_LLBSA);
	odeHandles.push_back("sLLB-SA", ODE_SLLBSA);

	//ODEs : atomistic meshes only
	atom_odeHandles.push_back("LLG", ODE_LLG);
	atom_odeHandles.push_back("LLGStatic", ODE_LLGSTATIC);
	atom_odeHandles.push_back("LLGStatic-SA", ODE_LLGSTATICSA);
	atom_odeHandles.push_back("LLG-STT", ODE_LLGSTT);
	atom_odeHandles.push_back("sLLG", ODE_SLLG);
	atom_odeHandles.push_back("sLLG-STT", ODE_SLLGSTT);
	atom_odeHandles.push_back("LLG-SA", ODE_LLGSA);
	atom_odeHandles.push_back("sLLG-SA", ODE_SLLGSA);

	//Evaluation methods
	odeEvalHandles.push_back("Euler", EVAL_EULER);
	odeEvalHandles.push_back("TEuler", EVAL_TEULER);
	odeEvalHandles.push_back("AHeun", EVAL_AHEUN);
	odeEvalHandles.push_back("RK4", EVAL_RK4);
	odeEvalHandles.push_back("ABM", EVAL_ABM);
	odeEvalHandles.push_back("RK23", EVAL_RK23);
	odeEvalHandles.push_back("RKF45", EVAL_RKF45);
	odeEvalHandles.push_back("RKF56", EVAL_RKF56);
	odeEvalHandles.push_back("RKCK45", EVAL_RKCK45);
	odeEvalHandles.push_back("RKDP54", EVAL_RKDP54);
	odeEvalHandles.push_back("SDesc", EVAL_SD);

	//Allowed evaluation methods for given ODE
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54), ODE_LLG);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54, EVAL_SD), ODE_LLGSTATIC);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54, EVAL_SD), ODE_LLGSTATICSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54), ODE_LLGSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54), ODE_LLB);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54), ODE_LLBSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54), ODE_LLGSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF45, EVAL_RKF56, EVAL_RKCK45, EVAL_RKDP54), ODE_LLBSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4), ODE_SLLG);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4), ODE_SLLGSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4), ODE_SLLB);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4), ODE_SLLBSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4), ODE_SLLGSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4), ODE_SLLBSA);

	odeDefaultEval.push_back(EVAL_RKF45, ODE_LLG);
	odeDefaultEval.push_back(EVAL_SD, ODE_LLGSTATIC);
	odeDefaultEval.push_back(EVAL_SD, ODE_LLGSTATICSA);
	odeDefaultEval.push_back(EVAL_RKF45, ODE_LLGSTT);
	odeDefaultEval.push_back(EVAL_RKF45, ODE_LLB);
	odeDefaultEval.push_back(EVAL_RKF45, ODE_LLBSTT);
	odeDefaultEval.push_back(EVAL_RKF45, ODE_LLGSA);
	odeDefaultEval.push_back(EVAL_RKF45, ODE_LLBSA);
	odeDefaultEval.push_back(EVAL_AHEUN, ODE_SLLG);
	odeDefaultEval.push_back(EVAL_AHEUN, ODE_SLLGSTT);
	odeDefaultEval.push_back(EVAL_AHEUN, ODE_SLLB);
	odeDefaultEval.push_back(EVAL_AHEUN, ODE_SLLBSTT);
	odeDefaultEval.push_back(EVAL_AHEUN, ODE_SLLGSA);
	odeDefaultEval.push_back(EVAL_AHEUN, ODE_SLLBSA);

	//---------------------------------------------------------------- SIMULATION SCHEDULE

	stageDescriptors.push_back("Relax", StageDescriptor(SS_RELAX), SS_RELAX);
	stageDescriptors.push_back("Hxyz", StageDescriptor(SS_HFIELDXYZ, "A/m", false), SS_HFIELDXYZ);
	stageDescriptors.push_back("Hxyz_seq", StageDescriptor(SS_HFIELDXYZSEQ, "A/m", false), SS_HFIELDXYZSEQ);
	stageDescriptors.push_back("Hpolar_seq", StageDescriptor(SS_HPOLARSEQ, "A/m", false), SS_HPOLARSEQ);
	stageDescriptors.push_back("Hfmr", StageDescriptor(SS_HFMR, "A/m", false), SS_HFMR);
	stageDescriptors.push_back("Hequation", StageDescriptor(SS_HFIELDEQUATION, "", false), SS_HFIELDEQUATION);
	//Don't define equation sequences to reduce crowding of options in the interface : we can set an equation sequence using an equation stage by editing the text itself, e.g. a format of the type "number: ..." means an equation sequence
	stageDescriptors.push_back("Hequation_seq", StageDescriptor(SS_HFIELDEQUATIONSEQ, "", false, true), SS_HFIELDEQUATIONSEQ);
	stageDescriptors.push_back("Hfile", StageDescriptor(SS_HFIELDFILE, "", false), SS_HFIELDFILE);
	stageDescriptors.push_back("V", StageDescriptor(SS_V, "V"), SS_V);
	stageDescriptors.push_back("V_seq", StageDescriptor(SS_VSEQ, "V"), SS_VSEQ);
	//deprecated (use equation instead, but keep SS_ values to keep older simulation files compatible):
	//stageDescriptors.push_back("Vsin", StageDescriptor(SS_VSIN, "V"), SS_VSIN);
	//stageDescriptors.push_back("Vcos", StageDescriptor(SS_VCOS, "V"), SS_VCOS);
	stageDescriptors.push_back("Vequation", StageDescriptor(SS_VEQUATION, ""), SS_VEQUATION);
	//Don't define equation sequences to reduce crowding of options in the interface : we can set an equation sequence using an equation stage by editing the text itself, e.g. a format of the type "number: ..." means an equation sequence
	stageDescriptors.push_back("Vequation_seq", StageDescriptor(SS_VEQUATIONSEQ, "", true, true), SS_VEQUATIONSEQ);
	stageDescriptors.push_back("Vfile", StageDescriptor(SS_VFILE, ""), SS_VFILE);
	stageDescriptors.push_back("I", StageDescriptor(SS_I, "A"), SS_I);
	stageDescriptors.push_back("I_seq", StageDescriptor(SS_ISEQ, "A"), SS_ISEQ);
	//deprecated (use equation instead, but keep SS_ values to keep older simulation files compatible):
	//stageDescriptors.push_back("Isin", StageDescriptor(SS_ISIN, "A"), SS_ISIN);
	//stageDescriptors.push_back("Icos", StageDescriptor(SS_ICOS, "A"), SS_ICOS);
	stageDescriptors.push_back("Iequation", StageDescriptor(SS_IEQUATION, ""), SS_IEQUATION);
	//Don't define equation sequences to reduce crowding of options in the interface : we can set an equation sequence using an equation stage by editing the text itself, e.g. a format of the type "number: ..." means an equation sequence
	stageDescriptors.push_back("Iequation_seq", StageDescriptor(SS_IEQUATIONSEQ, "", true, true), SS_IEQUATIONSEQ);
	stageDescriptors.push_back("Ifile", StageDescriptor(SS_IFILE, ""), SS_IFILE);
	stageDescriptors.push_back("T", StageDescriptor(SS_T, "K", false), SS_T);
	stageDescriptors.push_back("T_seq", StageDescriptor(SS_TSEQ, "K", false), SS_TSEQ);
	stageDescriptors.push_back("Tequation", StageDescriptor(SS_TEQUATION, "", false), SS_TEQUATION);
	//Don't define equation sequences to reduce crowding of options in the interface : we can set an equation sequence using an equation stage by editing the text itself, e.g. a format of the type "number: ..." means an equation sequence
	stageDescriptors.push_back("Tequation_seq", StageDescriptor(SS_TEQUATIONSEQ, "", false, true), SS_TEQUATIONSEQ);
	stageDescriptors.push_back("Tfile", StageDescriptor(SS_TFILE, "", false), SS_TFILE);
	stageDescriptors.push_back("Q", StageDescriptor(SS_Q, "W/m3", false), SS_Q);
	stageDescriptors.push_back("Q_seq", StageDescriptor(SS_QSEQ, "W/m3", false), SS_QSEQ);
	stageDescriptors.push_back("Qequation", StageDescriptor(SS_QEQUATION, "", false), SS_QEQUATION);
	//Don't define equation sequences to reduce crowding of options in the interface : we can set an equation sequence using an equation stage by editing the text itself, e.g. a format of the type "number: ..." means an equation sequence
	stageDescriptors.push_back("Qequation_seq", StageDescriptor(SS_QEQUATIONSEQ, "", false, true), SS_QEQUATIONSEQ);
	stageDescriptors.push_back("Qfile", StageDescriptor(SS_QFILE, "", false), SS_QFILE);
	stageDescriptors.push_back("Sunif", StageDescriptor(SS_TSIGPOLAR, "Pa", false), SS_TSIGPOLAR);
	stageDescriptors.push_back("MonteCarlo", StageDescriptor(SS_MONTECARLO), SS_MONTECARLO);

	stageStopDescriptors.push_back("nostop", StageStopDescriptor(STOP_NOSTOP), STOP_NOSTOP);
	stageStopDescriptors.push_back("iter", StageStopDescriptor(STOP_ITERATIONS), STOP_ITERATIONS);
	stageStopDescriptors.push_back("mxh", StageStopDescriptor(STOP_MXH), STOP_MXH);
	stageStopDescriptors.push_back("dmdt", StageStopDescriptor(STOP_DMDT), STOP_DMDT);
	stageStopDescriptors.push_back("time", StageStopDescriptor(STOP_TIME, "s"), STOP_TIME);
	stageStopDescriptors.push_back("mxh_iter", StageStopDescriptor(STOP_MXH_ITER), STOP_MXH_ITER);
	stageStopDescriptors.push_back("dmdt_iter", StageStopDescriptor(STOP_DMDT_ITER), STOP_DMDT_ITER);

	dataSaveDescriptors.push_back("none", DataSaveDescriptor(DSAVE_NONE), DSAVE_NONE);
	dataSaveDescriptors.push_back("stage", DataSaveDescriptor(DSAVE_STAGE), DSAVE_STAGE);
	dataSaveDescriptors.push_back("step", DataSaveDescriptor(DSAVE_STEP), DSAVE_STEP);
	dataSaveDescriptors.push_back("iter", DataSaveDescriptor(DSAVE_ITER), DSAVE_ITER);
	dataSaveDescriptors.push_back("time", DataSaveDescriptor(DSAVE_TIME, "s"), DSAVE_TIME);

	//---------------------------------------------------------------- DEFAULT STARTING STATE

	//Default simulation stages - this must be added first so simStages is not an empty vector
	//Default stage: nothing set, not stopping, nothing saved
	AddGenericStage(SS_RELAX);
	EditStageStopCondition(0, STOP_NOSTOP);

	//default save data settings - note, saveDataList should not be empty
	NewSaveDataEntry(DATA_STAGESTEP);
	NewSaveDataEntry(DATA_ITERATIONS);
	NewSaveDataEntry(DATA_TIME);
	NewSaveDataEntry(DATA_HA, "permalloy");
	NewSaveDataEntry(DATA_AVM, "permalloy");

	//Default data box entries
	NewDataBoxField(DatumConfig(DATA_STAGESTEP));
	NewDataBoxField(DatumConfig(DATA_ITERATIONS));
	NewDataBoxField(DatumConfig(DATA_TIME));
	NewDataBoxField(DatumConfig(DATA_MXH));

	//---------------------------------------------------------------- STARTUP OPTIONS
	
	//set error log file with path
	errorlog_fileName = GetUserDocumentsPath() + boris_data_directory + "errorlog.txt";
	//set startup options file with path
	startup_options_file = GetUserDocumentsPath() + boris_data_directory + "startup.txt";

	//Load options for startup first
	Load_Startup_Flags();

	//---------------------------------------------------------------- SERVER START

	server_port = server_port_;
	server_pwd = server_pwd_;

	//start network sockets thread to listen for incoming messages
	commSocket.Change_Password(server_pwd);
	if (server_port.length()) commSocket.Change_Port(server_port);
	commSocket.Change_RecvSleep(server_recv_sleep_ms);
	if (start_scriptserver && !server_port_.length()) {

		Script_Server_Control(true);
	}
	else if (server_port_.length()) {

		//if a server_port_ was passed in, then over-ride start_scriptserver flag
		//in this case the provided password must match the local password
		if (server_pwd == server_pwd_ || !server_pwd.length()) {

			server_port = server_port_;
			commSocket.Change_Port(server_port);

			Script_Server_Control(true);
		}
	}
	
	//---------------------------------------------------------------- FINAL SETTINGS

	//Set number of OpenMP threads (default is use all available but user can configure this)
	if (!OmpThreads) OmpThreads = omp_get_num_procs();

	//Update display - do not animate starting view
	UpdateScreen_AutoSet_Sudden();

	//make sure the disk buffer has the correct size
	if (savedata_diskbuffer.size() != savedata_diskbuffer_size) savedata_diskbuffer.resize(savedata_diskbuffer_size);
	savedata_diskbuffer_position = 0;

	//---------------------------------------------------------------- EMBEDDED PYTHON SCRIPT

	python_script = python_script_;
	python_script_mGPU = python_script_mGPU_;
	if (python_script_parallel_.size() >= 2) python_script_parallel = python_script_parallel_;
	python_script_terminate = python_script_terminate_;
	python_script_deletesource = python_script_deletesource_;

	///////////////////////////////////////////////////////////// FROM THIS POINT NO SETTINGS ESSENTIAL FOR CORRECT INITIALIZATION ALLOWED

	//---------------------------------------------------------------- UPDATES

	//Check with "www.boris-spintronics.uk" if program version is up to date, and get latest update time for materials database
	if (start_check_updates && !python_script.length()) single_call_launch(&Simulation::CheckUpdate, THREAD_HANDLEMESSAGE2);
		
	//---------------------------------------------------------------- MESSAGES

	//if this directory doesn't exist create it
	if (!MakeDirectory(directory)) BD.DisplayConsoleError("ERROR : Couldn't create user directory.");

	BD.DisplayConsoleMessage("Console activated...");
	BD.DisplayFormattedConsoleMessage("[tc0,0.5,0,1/tc]To open manual use the <b>manual</b> command.");

	//---------------------------------------------------------------- SAVE DEFAULT and CTOR FINISH

	//save the default state in program directory for loading with "default" command (update this automatically in case the program version has changed)
	if(!python_script.length()) SaveSimulation(directory + "default");
}

Simulation::~Simulation()
{
	BD.DisplayConsoleMessage("Shutting down...");

#if PYTHON_EMBEDDING == 1
	if (is_thread_running(THREAD_PYTHONSCRIPT)) python_script_running = false;
#endif

	while (is_thread_running(THREAD_LOOP)) {
		StopSimulation();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}

	Stop_All_Threads();

	BD.DisplayConsoleMessage("All threads stopped. Clean-up...");
}
