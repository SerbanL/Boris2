#include "stdafx.h"
#include "Simulation.h"

#if GRAPHICS == 1
Simulation::Simulation(HWND hWnd, int Program_Version, std::string server_port_, std::string server_pwd_, int cudaDevice) :
	err_hndl(this),
	BD(hWnd, new SimTOFunct(this, &Simulation::ConsoleActionHandler, &Simulation::ConsoleInteractiveObjectState)),
	SimulationSharedData(true),
	mdb(params_for_meshtype(MESH_SUPERMESH)),
	ProgramStateNames(this,
		{
			VINFO(BD),
			VINFO(directory), VINFO(savedataFile), VINFO(imageSaveFileBase), VINFO(currentSimulationFile), VINFO(appendToDataFile), VINFO(saveDataFlag), VINFO(saveImageFlag),
			VINFO(saveDataList), VINFO(dataBoxList),
			VINFO(stage_step),
			VINFO(simStages), VINFO(iterUpdate), VINFO(autocomplete),
			VINFO(SMesh),
			VINFO(cudaEnabled), VINFO(cudaDeviceSelect),
			VINFO(shape_change_individual),
			VINFO(static_transport_solver),
			VINFO(image_cropping), VINFO(displayTransparency), VINFO(displayThresholds), VINFO(displayThresholdTrigger),
			VINFO(userConstants)
		}, {})
#else
Simulation::Simulation(int Program_Version, std::string server_port_, std::string server_pwd_, int cudaDevice) :
	err_hndl(this),
	BD(),
	SimulationSharedData(true),
	mdb(params_for_meshtype(MESH_SUPERMESH)),
	ProgramStateNames(this,
		{
			VINFO(BD),
			VINFO(directory), VINFO(savedataFile), VINFO(imageSaveFileBase), VINFO(currentSimulationFile), VINFO(appendToDataFile), VINFO(saveDataFlag), VINFO(saveImageFlag),
			VINFO(saveDataList), VINFO(dataBoxList),
			VINFO(stage_step),
			VINFO(simStages), VINFO(iterUpdate), VINFO(autocomplete),
			VINFO(SMesh),
			VINFO(cudaEnabled), VINFO(cudaDeviceSelect),
			VINFO(shape_change_individual),
			VINFO(static_transport_solver),
			VINFO(image_cropping), VINFO(displayTransparency), VINFO(displayThresholds), VINFO(displayThresholdTrigger),
			VINFO(userConstants)
		}, {})
#endif
{
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

	if (cudaAvailable) cudaDeviceReset();

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

	commands.insert(CMD_COMPUTEFIELDS, CommandSpecifier(CMD_COMPUTEFIELDS), "computefields");
	commands[CMD_COMPUTEFIELDS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>computefields</b>";
	commands[CMD_COMPUTEFIELDS].descr = "[tc0,0.5,0.5,1/tc]Run simulation from current state for a single iteration without advancing the simulation time.";

	commands.insert(CMD_ISRUNNING, CommandSpecifier(CMD_ISRUNNING), "isrunning");
	commands[CMD_ISRUNNING].usage = "[tc0,0.5,0,1/tc]USAGE : <b>isrunning</b>";
	commands[CMD_ISRUNNING].descr = "[tc0,0.5,0.5,1/tc]Checks if the simulation is running and sends state value to the calling script.";
	commands[CMD_ISRUNNING].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>state</i> - return simulation running state.";

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

	commands.insert(CMD_ADDDIPOLEMESH, CommandSpecifier(CMD_ADDDIPOLEMESH), "adddipole");
	commands[CMD_ADDDIPOLEMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDDIPOLEMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>adddipole</b> <i>name rectangle</i>";
	commands[CMD_ADDDIPOLEMESH].descr = "[tc0,0.5,0.5,1/tc]Add a rectangular dipole with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDDIPOLEMESH].unit = "m";

	commands.insert(CMD_ADDINSULATORMESH, CommandSpecifier(CMD_ADDINSULATORMESH), "addinsulator");
	commands[CMD_ADDINSULATORMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addinsulator</b> <i>name rectangle</i>";
	commands[CMD_ADDINSULATORMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDINSULATORMESH].descr = "[tc0,0.5,0.5,1/tc]Add an insulator mesh with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDINSULATORMESH].unit = "m";

	commands.insert(CMD_ADDDIAMAGNETICMESH, CommandSpecifier(CMD_ADDDIAMAGNETICMESH), "adddiamagnet");
	commands[CMD_ADDDIAMAGNETICMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>adddiamagnet</b> <i>name rectangle</i>";
	commands[CMD_ADDDIAMAGNETICMESH].limits = { { Any(), Any() },{ Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_ADDDIAMAGNETICMESH].descr = "[tc0,0.5,0.5,1/tc]Add a diamagnetic mesh with given name and rectangle (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_ADDDIAMAGNETICMESH].unit = "m";

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
	commands[CMD_DELMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delmesh</b> <i>name</i>";
	commands[CMD_DELMESH].descr = "[tc0,0.5,0.5,1/tc]Delete mesh with given name.";

	commands.insert(CMD_RENAMEMESH, CommandSpecifier(CMD_RENAMEMESH), "renamemesh");
	commands[CMD_RENAMEMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>renamemesh</b> <i>(old_name) new_name</i>";
	commands[CMD_RENAMEMESH].descr = "[tc0,0.5,0.5,1/tc]Rename mesh. If old_name not specified then the mesh in focus is renamed.";

	commands.insert(CMD_MESHFOCUS, CommandSpecifier(CMD_MESHFOCUS), "meshfocus");
	commands[CMD_MESHFOCUS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshfocus</b> <i>meshname</i>";
	commands[CMD_MESHFOCUS].descr = "[tc0,0.5,0.5,1/tc]Change mesh focus to given mesh name.";
	commands[CMD_MESHFOCUS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshname</i> - return name of mesh in focus.";

	commands.insert(CMD_MESHFOCUS2, CommandSpecifier(CMD_MESHFOCUS2), "meshfocus2");
	commands[CMD_MESHFOCUS2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshfocus2</b> <i>meshname</i>";
	commands[CMD_MESHFOCUS2].descr = "[tc0,0.5,0.5,1/tc]Change mesh focus to given mesh name but do not change camera orientation.";
	commands[CMD_MESHFOCUS2].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>meshname</i> - return name of mesh in focus.";

	commands.insert(CMD_MESHRECT, CommandSpecifier(CMD_MESHRECT), "meshrect");
	commands[CMD_MESHRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>meshrect</b> <i>rectangle</i>";
	commands[CMD_MESHRECT].limits = { { Rect(DBL3(-MAXSIMSPACE/2), DBL3(-MAXSIMSPACE/2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_MESHRECT].descr = "[tc0,0.5,0.5,1/tc]Change rectangle of mesh in focus (m). The rectangle can be specified as: <i>sx sy sz ex ey ez</i> for the start and end points in Cartesian coordinates, or as: <i>ex ey ez</i> with the start point as the origin.";
	commands[CMD_MESHRECT].unit = "m";
	commands[CMD_MESHRECT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>rectangle</i> - return rectangle of mesh in focus.";

	commands.insert(CMD_SCALEMESHRECTS, CommandSpecifier(CMD_SCALEMESHRECTS), "scalemeshrects");
	commands[CMD_SCALEMESHRECTS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>scalemeshrects</b> <i>status</i>";
	commands[CMD_SCALEMESHRECTS].descr = "[tc0,0.5,0.5,1/tc]When changing a mesh rectangle scale and shift all other mesh rectangles in proportion if status set.";
	commands[CMD_SCALEMESHRECTS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_MESH, CommandSpecifier(CMD_MESH), "mesh");
	commands[CMD_MESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mesh</b>";
	commands[CMD_MESH].descr = "[tc0,0.5,0.5,1/tc]Display information for all meshes.";

	commands.insert(CMD_CELLSIZE, CommandSpecifier(CMD_CELLSIZE), "cellsize");
	commands[CMD_CELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>cellsize</b> <i>value</i>";
	commands[CMD_CELLSIZE].limits = { { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_CELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of mesh in focus (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_CELLSIZE].unit = "m";
	commands[CMD_CELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return cellsize of mesh in focus.";

	commands.insert(CMD_ECELLSIZE, CommandSpecifier(CMD_ECELLSIZE), "ecellsize");
	commands[CMD_ECELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ecellsize</b> <i>value</i>";
	commands[CMD_ECELLSIZE].limits = { { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_ECELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of mesh in focus for electrical conduction (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_ECELLSIZE].unit = "m";
	commands[CMD_ECELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return electrical conduction cellsize of mesh in focus.";

	commands.insert(CMD_TCELLSIZE, CommandSpecifier(CMD_TCELLSIZE), "tcellsize");
	commands[CMD_TCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tcellsize</b> <i>value</i>";
	commands[CMD_TCELLSIZE].limits = { { DBL3(MINMESHSPACE/2), Any() } };
	commands[CMD_TCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of mesh in focus for thermal conduction (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_TCELLSIZE].unit = "m";
	commands[CMD_TCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return thermal conduction cellsize of mesh in focus.";

	commands.insert(CMD_SCELLSIZE, CommandSpecifier(CMD_SCELLSIZE), "scellsize");
	commands[CMD_SCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>scellsize</b> <i>value</i>";
	commands[CMD_SCELLSIZE].limits = { { DBL3(MINMESHSPACE / 2), Any() } };
	commands[CMD_SCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of mesh in focus for stochastic properties (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_SCELLSIZE].unit = "m";
	commands[CMD_SCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return stochastic properties cellsize of mesh in focus.";

	commands.insert(CMD_MCELLSIZE, CommandSpecifier(CMD_MCELLSIZE), "mcellsize");
	commands[CMD_MCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcellsize</b> <i>value</i>";
	commands[CMD_MCELLSIZE].limits = { { DBL3(MINMESHSPACE / 2), Any() } };
	commands[CMD_MCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change cellsize of mesh in focus for mechanical solver (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_MCELLSIZE].unit = "m";
	commands[CMD_MCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return mechanical cellsize of mesh in focus.";

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
	commands[CMD_ATOMDMCELLSIZE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dmcellsize</b> <i>value</i>";
	commands[CMD_ATOMDMCELLSIZE].limits = { { DBL3(MINMESHSPACE / 2), Any() } };
	commands[CMD_ATOMDMCELLSIZE].descr = "[tc0,0.5,0.5,1/tc]Change demagnetizing field macrocell size of mesh in focus, for atomistic meshes (m). The cellsize can be specified as: <i>hx hy hz</i>, or as: <i>hxyz</i>";
	commands[CMD_ATOMDMCELLSIZE].unit = "m";
	commands[CMD_ATOMDMCELLSIZE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>cellsize</i> - return demagnetizing field macrocell size.";

	commands.insert(CMD_DELRECT, CommandSpecifier(CMD_DELRECT), "delrect");
	commands[CMD_DELRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delrect</b> <i>rectangle (meshname)</i>";
	commands[CMD_DELRECT].limits = { { Rect(), Any() }, { Any(), Any() } };
	commands[CMD_DELRECT].descr = "[tc0,0.5,0.5,1/tc]Void rectangle (m) within given mesh (active mesh if name not given). The rectangle coordinates are relative to specified mesh.";
	commands[CMD_DELRECT].unit = "m";

	commands.insert(CMD_ADDRECT, CommandSpecifier(CMD_ADDRECT), "addrect");
	commands[CMD_ADDRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addrect</b> <i>rectangle (meshname)</i>";
	commands[CMD_ADDRECT].limits = { { Rect(), Any() }, { Any(), Any() } };
	commands[CMD_ADDRECT].descr = "[tc0,0.5,0.5,1/tc]Fill rectangle (m) within given mesh (active mesh if name not given). The rectangle coordinates are relative to specified mesh.";
	commands[CMD_ADDRECT].unit = "m";

	commands.insert(CMD_RESETMESH, CommandSpecifier(CMD_RESETMESH), "resetmesh");
	commands[CMD_RESETMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>resetmesh</b> <i>(meshname)</i>";
	commands[CMD_RESETMESH].descr = "[tc0,0.5,0.5,1/tc]Reset to constant magnetization in given mesh (active mesh if name not given).";

	commands.insert(CMD_LOADMASKFILE, CommandSpecifier(CMD_LOADMASKFILE), "loadmaskfile");
	commands[CMD_LOADMASKFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadmaskfile</b> <i>(z_depth) (directory/)filename</i>";
	commands[CMD_LOADMASKFILE].descr = "[tc0,0.5,0.5,1/tc]Apply .png mask file to magnetization in active mesh (i.e. transfer shape from .png file to mesh - white means empty cells). If image is in grayscale then void cells up to given depth top down (z_depth > 0) or down up (z_depth < 0). If z-depth = 0 then void top down up to all z cells.";
	commands[CMD_LOADMASKFILE].unit = "m";

	commands.insert(CMD_INDIVIDUALMASKSHAPE, CommandSpecifier(CMD_INDIVIDUALMASKSHAPE), "individualshape");
	commands[CMD_INDIVIDUALMASKSHAPE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>individualmaskshape</b> <i>status</i>";
	commands[CMD_INDIVIDUALMASKSHAPE].limits = { { int(0), int(1) } };
	commands[CMD_INDIVIDUALMASKSHAPE].descr = "[tc0,0.5,0.5,1/tc]When changing the shape inside a mesh, e.g. through a mask file, set this flag to true so the shape is applied only to the primary displayed physical quantity. If set to false then all relevant physical quantities are shaped.";
	commands[CMD_INDIVIDUALMASKSHAPE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_SETANGLE, CommandSpecifier(CMD_SETANGLE), "setangle");
	commands[CMD_SETANGLE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setangle</b> <i>polar azimuthal (meshname)</i>";
	commands[CMD_SETANGLE].limits = { { double(-360.0), double(360.0) },{ double(-360.0), double(360.0) }, { Any(), Any() } };
	commands[CMD_SETANGLE].descr = "[tc0,0.5,0.5,1/tc]Set magnetization angle in mesh uniformly using polar coordinates. If mesh name not specified, this is set for all ferromagnetic meshes.";

	commands.insert(CMD_RANDOM, CommandSpecifier(CMD_RANDOM), "random");
	commands[CMD_RANDOM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>random</b> <i>(meshname, (seed))</i>";
	commands[CMD_RANDOM].limits = { { Any(), Any() }, { int(1), Any() } };
	commands[CMD_RANDOM].descr = "[tc0,0.5,0.5,1/tc]Set random magnetization distribution in mesh, with pseud-random number generator seed (default 1 if not specified). If mesh name not specified, set for focused mesh.";

	commands.insert(CMD_SETRECT, CommandSpecifier(CMD_SETRECT), "setrect");
	commands[CMD_SETRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setrect</b> <i>polar azimuthal rectangle (meshname)</i>";
	commands[CMD_SETRECT].limits = { { double(-360.0), double(360.0) },{ double(-360.0), double(360.0) }, { Rect(), Any() }, { Any(), Any() } };
	commands[CMD_SETRECT].descr = "[tc0,0.5,0.5,1/tc]Set magnetization angle in given rectangle of mesh (relative coordinates) uniformly using polar coordinates. If mesh name not specified, the active mesh is used.";
	commands[CMD_SETRECT].unit = "m";

	commands.insert(CMD_INVERTMAG, CommandSpecifier(CMD_INVERTMAG), "invertmag");
	commands[CMD_INVERTMAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>invertmag</b> <i>(components) (meshname)</i>";
	commands[CMD_INVERTMAG].limits = { { Any(), Any() } };
	commands[CMD_INVERTMAG].descr = "[tc0,0.5,0.5,1/tc]Invert magnetization direction. If mesh name not specified, the active mesh is used. You can choose to invert just one or two components instead of the entire vector: specify components as x y z, e.g. invertmag x";

	commands.insert(CMD_MIRRORMAG, CommandSpecifier(CMD_MIRRORMAG), "mirrormag");
	commands[CMD_MIRRORMAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mirrormag</b> <i>axis (meshname)</i>";
	commands[CMD_MIRRORMAG].limits = { { Any(), Any() } };
	commands[CMD_MIRRORMAG].descr = "[tc0,0.5,0.5,1/tc]Mirror magnetization in a given axis, specified as x, y, z, e.g. mirrormag x. If mesh name not specified, the active mesh is used";

	commands.insert(CMD_DWALL, CommandSpecifier(CMD_DWALL), "dwall");
	commands[CMD_DWALL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dwall</b> <i>longitudinal transverse width position (meshname)</i>";
	commands[CMD_DWALL].limits = { { Any(), Any() } , { Any(), Any() }, { double(0.0), Any() }, { double(0.0), Any() }, { Any(), Any() } };
	commands[CMD_DWALL].unit = "m";
	commands[CMD_DWALL].descr = "[tc0,0.5,0.5,1/tc]Create an idealised domain wall (tanh profile for longitudinal component, 1/cosh profile for transverse component) along the x-axis direction in the given mesh (active mesh if name not specified). For longitudinal and transverse specify the components of magnetization as x, -x, y, -y, z, -z, i.e. specify using these string literals. For width and position use metric units.";

	commands.insert(CMD_VORTEX, CommandSpecifier(CMD_VORTEX), "vortex");
	commands[CMD_VORTEX].usage = "[tc0,0.5,0,1/tc]USAGE : <b>vortex</b> <i>longitudinal rotation core (rectangle) (meshname)</i>";
	commands[CMD_VORTEX].limits = { { int(-1), int(+1) } , { int(-1), int(+1) }, { int(-1), int(+1) }, { Rect(), Any() }, { Any(), Any() } };
	commands[CMD_VORTEX].unit = "m";
	commands[CMD_VORTEX].descr = "[tc0,0.5,0.5,1/tc]Create a vortex domain wall with settings: longitudinal (-1: tail-to-tail, 1: head-to-head), rotation (-1: clockwise, 1: counter-clockwise), core (-1: down, 1: up). The vortex may be set in the given rectangle (entire mesh if not given), in the given mesh (focused mesh if not given).";

	commands.insert(CMD_SKYRMION, CommandSpecifier(CMD_SKYRMION), "skyrmion");
	commands[CMD_SKYRMION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyrmion</b> <i>core chirality diameter position (meshname)</i>";
	commands[CMD_SKYRMION].limits = { { int(-1), int(1) } , { int(-1), int(1) }, { double(1e-9), Any() },{ Any(), Any() }, { Any(), Any() } };
	commands[CMD_SKYRMION].unit = "m";
	commands[CMD_SKYRMION].descr = "[tc0,0.5,0.5,1/tc]Create an idealised Neel-type skyrmion with given diameter and centre position in the x-y plane (2 relative coordinates needed only) of the given mesh (active mesh if name not specified). Core specifies the skyrmion core direction: -1 for down, 1 for up. Chirality specifies the radial direction rotation: 1 for towards core, -1 away from core. For diameter and position use metric units.";

	commands.insert(CMD_SKYRMIONBLOCH, CommandSpecifier(CMD_SKYRMIONBLOCH), "skyrmionbloch");
	commands[CMD_SKYRMIONBLOCH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyrmionbloch</b> <i>core chirality diameter position (meshname)</i>";
	commands[CMD_SKYRMIONBLOCH].limits = { { int(-1), int(1) } ,{ int(-1), int(1) },{ double(1e-9), Any() },{ Any(), Any() },{ Any(), Any() } };
	commands[CMD_SKYRMIONBLOCH].unit = "m";
	commands[CMD_SKYRMIONBLOCH].descr = "[tc0,0.5,0.5,1/tc]Create an idealised Bloch-type skyrmion with given diameter and centre position in the x-y plane (2 relative coordinates needed only) of the given mesh (active mesh if name not specified). Core specifies the skyrmion core direction: -1 for down, 1 for up. Chirality specifies the radial direction rotation: 1 for clockwise, -1 for anti-clockwise. For diameter and position use metric units.";

	commands.insert(CMD_PBC, CommandSpecifier(CMD_PBC), "pbc");
	commands[CMD_PBC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>pbc</b> <i>meshname flag images</i>";
	commands[CMD_PBC].descr = "[tc0,0.5,0.5,1/tc]Set periodic boundary conditions for magnetization in given mesh (must be ferromagnetic). Flags specify types of perodic boundary conditions: x, y, or z; images specify the number of mesh images to use either side for the given direction when calculating the demagnetising kernel - a value of zero disables pbc. e.g. ""pbc x 10"" sets x periodic boundary conditions with 10 images either side for the focused mesh; ""pbc x 0"" clears pbc for the x axis.";
	commands[CMD_PBC].limits = { { Any(), Any() },{ Any(), Any() }, {int(0), int(1000)} };

	commands.insert(CMD_SETFIELD, CommandSpecifier(CMD_SETFIELD), "setfield");
	commands[CMD_SETFIELD].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setfield</b> <i>magnitude polar azimuthal (meshname)</i>";
	commands[CMD_SETFIELD].limits = { { DBL3(-MAXFIELD, -360.0, -360.0), DBL3(MAXFIELD, 360.0, 360.0) }, { Any(), Any() } };
	commands[CMD_SETFIELD].descr = "[tc0,0.5,0.5,1/tc]Set uniform magnetic field (A/m) using polar coordinates. If mesh name not specified, this is set for all magnetic meshes - must have Zeeman module added.";
	commands[CMD_SETFIELD].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i><Ha_x, Ha_y, Ha_z></i> - applied field in Cartesian coordinates for mesh in focus.";
	commands[CMD_SETFIELD].unit = "A/m";

	commands.insert(CMD_SETSTRESS, CommandSpecifier(CMD_SETSTRESS), "setstress");
	commands[CMD_SETSTRESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setstress</b> <i>magnitude polar azimuthal (meshname)</i>";
	commands[CMD_SETSTRESS].limits = { { DBL3(-MAXSTRESS, -360.0, -360.0), DBL3(MAXSTRESS, 360.0, 360.0) }, { Any(), Any() } };
	commands[CMD_SETSTRESS].descr = "[tc0,0.5,0.5,1/tc]Set uniform mechanical stress (Pa) using polar coordinates. If mesh name not specified, this is set for all magnetic meshes - must have MElastic module added.";
	commands[CMD_SETSTRESS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i><Tsig_x, Tsig_y, Tsig_z></i> - applied mechanical stress in Cartesian coordinates for mesh in focus.";
	commands[CMD_SETSTRESS].unit = "Pa";

	commands.insert(CMD_MODULES, CommandSpecifier(CMD_MODULES), "modules");
	commands[CMD_MODULES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>modules</b>";
	commands[CMD_MODULES].descr = "[tc0,0.5,0.5,1/tc]Show interactive list of available and currently set modules.";

	commands.insert(CMD_ADDMODULE, CommandSpecifier(CMD_ADDMODULE), "addmodule");
	commands[CMD_ADDMODULE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addmodule</b> <i>meshname handle</i>";
	commands[CMD_ADDMODULE].descr = "[tc0,0.5,0.5,1/tc]Add module with given handle to named mesh.";

	commands.insert(CMD_DELMODULE, CommandSpecifier(CMD_DELMODULE), "delmodule");
	commands[CMD_DELMODULE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delmodule</b> <i>meshname handle</i>";
	commands[CMD_DELMODULE].descr = "[tc0,0.5,0.5,1/tc]Delete module with given handle from named mesh.";

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
	commands[CMD_EXCLUDEMULTICONVDEMAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>excludemulticonvdemag</b> <i>status meshname</i>";
	commands[CMD_EXCLUDEMULTICONVDEMAG].descr = "[tc0,0.5,0.5,1/tc]Set exclusion status (0 or 1) of named mesh from multi-layered demag convolution.";
	commands[CMD_EXCLUDEMULTICONVDEMAG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status</i>";

	commands.insert(CMD_ODE, CommandSpecifier(CMD_ODE), "ode");
	commands[CMD_ODE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ode</b>";
	commands[CMD_ODE].descr = "[tc0,0.5,0.5,1/tc]Show interactive list of available and currently set ODEs and evaluation methods.";

	commands.insert(CMD_EVALSPEEDUP, CommandSpecifier(CMD_EVALSPEEDUP), "evalspeedup");
	commands[CMD_EVALSPEEDUP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>evalspeedup</b> <i>status</i>";
	commands[CMD_EVALSPEEDUP].limits = { { int(EVALSPEEDUP_NONE), int(EVALSPEEDUP_NUMENTRIES) - 1 } };
	commands[CMD_EVALSPEEDUP].descr = "[tc0,0.5,0.5,1/tc]<b>!!!Experimental!!!</b>Status levels: 0 (no speedup), 1 (accurate), 2 (aggressive), 3 (extreme).";

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
	commands[CMD_ASTEPCTRL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>astepctrl</b> <i>err_fail err_high err_low dT_incr dT_min dT_max</i>";
	commands[CMD_ASTEPCTRL].limits = { 
		{double(MINODERELERROR), double(MAXODERELERROR)},
		{double(MINODERELERROR), double(MAXODERELERROR)}, 
		{double(MINODERELERROR), double(MAXODERELERROR)},
		{double(1.0), double(2.0)},
		{double(MINTIMESTEP), double(MAXTIMESTEP)},
		{double(MINTIMESTEP), double(MAXTIMESTEP)} };
	commands[CMD_ASTEPCTRL].descr = "[tc0,0.5,0.5,1/tc]Set parameters for adaptive time step control: err_fail - repeat step above this, err_high - decrease dT above this, err_low - increase dT below this, dT_incr - increase dT using fixed multiplier, dT_min, dT_max - dT bounds.";

	commands.insert(CMD_SHOWDATA, CommandSpecifier(CMD_SHOWDATA), "showdata");
	commands[CMD_SHOWDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showdata</b> <i>dataname (meshname, (rectangle))</i>";
	commands[CMD_SHOWDATA].descr = "[tc0,0.5,0.5,1/tc]Show value(s) for dataname. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, active mesh is used with entire mesh rectangle.";
	commands[CMD_SHOWDATA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>varies</i>";
	commands[CMD_SHOWDATA].unit = "m";

	commands.insert(CMD_DATA, CommandSpecifier(CMD_DATA), "data");
	commands[CMD_DATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>data</b>";
	commands[CMD_DATA].descr = "[tc0,0.5,0.5,1/tc]Shows list of currently set output data and available data.";
	commands[CMD_DATA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>number of set output data fields</i>";

	commands.insert(CMD_ADDDATA, CommandSpecifier(CMD_ADDDATA), "adddata");
	commands[CMD_ADDDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>adddata</b> <i>dataname (meshname, (rectangle))</i>";
	commands[CMD_ADDDATA].limits = { { Any(), Any() }, { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_ADDDATA].descr = "[tc0,0.5,0.5,1/tc]Add dataname to list of output data. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, active mesh is used with entire mesh rectangle.";
	commands[CMD_ADDDATA].unit = "m";

	commands.insert(CMD_SETDATA, CommandSpecifier(CMD_SETDATA), "setdata");
	commands[CMD_SETDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdata</b> <i>dataname (meshname, (rectangle))</i>";
	commands[CMD_SETDATA].limits = { { Any(), Any() }, { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_SETDATA].descr = "[tc0,0.5,0.5,1/tc]Delete all currently set output data and set dataname to list of output data. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, active mesh is used with entire mesh rectangle.";
	commands[CMD_SETDATA].unit = "m";

	commands.insert(CMD_DELDATA, CommandSpecifier(CMD_DELDATA), "deldata");
	commands[CMD_DELDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>deldata</b> <i>index</i>";
	commands[CMD_DELDATA].limits = { { int(-1), Any() } };
	commands[CMD_DELDATA].descr = "[tc0,0.5,0.5,1/tc]Delete data from list of output data at index number. If index number is -1 then delete all data fields, leaving just a default time data field - there must always be at least 1 output data field.";

	commands.insert(CMD_EDITDATA, CommandSpecifier(CMD_EDITDATA), "editdata");
	commands[CMD_EDITDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editdata</b> <i>index dataname (meshname, (rectangle))</i>";
	commands[CMD_EDITDATA].limits = { { int(0), Any() }, { Any(), Any() }, { Any(), Any() }, { Rect(), Any() } };
	commands[CMD_EDITDATA].descr = "[tc0,0.5,0.5,1/tc]Edit entry in list of output data at given index in list. If applicable specify meshname and rectangle (m) in mesh. If not specified and required, active mesh is used with entire mesh rectangle.";
	commands[CMD_EDITDATA].unit = "m";

	commands.insert(CMD_ADDPINNEDDATA, CommandSpecifier(CMD_ADDPINNEDDATA), "addpinneddata");
	commands[CMD_ADDPINNEDDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addpinneddata</b> <i>dataname (meshname, (rectangle))</i>";
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

	commands.insert(CMD_STAGES, CommandSpecifier(CMD_STAGES), "stages");
	commands[CMD_STAGES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>stages</b>";
	commands[CMD_STAGES].descr = "[tc0,0.5,0.5,1/tc]Shows list of currently set simulation stages and available stage types.";
	commands[CMD_STAGES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>number of set stages</i>";

	commands.insert(CMD_ADDSTAGE, CommandSpecifier(CMD_ADDSTAGE), "addstage");
	commands[CMD_ADDSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>addstage</b> <i>stagetype (meshname)</i>";
	commands[CMD_ADDSTAGE].descr = "[tc0,0.5,0.5,1/tc]Add a generic stage type to the simulation schedule with name stagetype, specifying a meshname if needed (if not specified and required, active mesh is used).";

	commands.insert(CMD_SETSTAGE, CommandSpecifier(CMD_SETSTAGE), "setstage");
	commands[CMD_SETSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setstage</b> <i>stagetype (meshname)</i>";
	commands[CMD_SETSTAGE].descr = "[tc0,0.5,0.5,1/tc]Delete all currently set stages, and set a new generic stage type to the simulation schedule with name stagetype, specifying a meshname if needed (if not specified and required, active mesh is used).";

	commands.insert(CMD_DELSTAGE, CommandSpecifier(CMD_DELSTAGE), "delstage");
	commands[CMD_DELSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>delstage</b> <i>index</i>";
	commands[CMD_DELSTAGE].limits = { { int(-1), Any() } };
	commands[CMD_DELSTAGE].descr = "[tc0,0.5,0.5,1/tc]Delete stage from simulation schedule at index number. If index number is -1 then delete all stages, leaving just a default Relax stage - there must always be at least 1 stage set.";

	commands.insert(CMD_EDITSTAGE, CommandSpecifier(CMD_EDITSTAGE), "editstage");
	commands[CMD_EDITSTAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editstage</b> <i>index stagetype (meshname)</i>";
	commands[CMD_EDITSTAGE].limits = { { int(0), Any() }, { Any(), Any() }, { Any(), Any() } };
	commands[CMD_EDITSTAGE].descr = "[tc0,0.5,0.5,1/tc]Edit stage type from simulation schedule at index number.";

	commands.insert(CMD_EDITSTAGEVALUE, CommandSpecifier(CMD_EDITSTAGEVALUE), "editstagevalue");
	commands[CMD_EDITSTAGEVALUE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editstagevalue</b> <i>index value</i>";
	commands[CMD_EDITSTAGEVALUE].limits = { { int(0), Any() }, { Any(), Any() } };
	commands[CMD_EDITSTAGEVALUE].descr = "[tc0,0.5,0.5,1/tc]Edit stage setting value in simulation schedule. The value type depends on the stage type.";

	commands.insert(CMD_EDITSTAGESTOP, CommandSpecifier(CMD_EDITSTAGESTOP), "editstagestop");
	commands[CMD_EDITSTAGESTOP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editstagestop</b> <i>index stoptype (stopvalue)</i>";
	commands[CMD_EDITSTAGESTOP].descr = "[tc0,0.5,0.5,1/tc]Edit stage/step stopping condition in simulation schedule. Use index < 0 to set condition for all stages.";

	commands.insert(CMD_EDITDATASAVE, CommandSpecifier(CMD_EDITDATASAVE), "editdatasave");
	commands[CMD_EDITDATASAVE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>editdatasave</b> <i>index savetype (savevalue)</i>";
	commands[CMD_EDITDATASAVE].descr = "[tc0,0.5,0.5,1/tc]Edit data saving condition in simulation schedule. Use index < 0 to set condition for all stages.";

	commands.insert(CMD_PARAMS, CommandSpecifier(CMD_PARAMS), "params");
	commands[CMD_PARAMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>params</b> <i>(meshname)</i>";
	commands[CMD_PARAMS].descr = "[tc0,0.5,0.5,1/tc]List all material parameters. If meshname not given use the active mesh.";

	commands.insert(CMD_SETPARAM, CommandSpecifier(CMD_SETPARAM), "setparam");
	commands[CMD_SETPARAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparam</b> <i>meshname paramname (value)</i>";
	commands[CMD_SETPARAM].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter to given value.";
	commands[CMD_SETPARAM].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - return value of named parameter in named mesh.";

	commands.insert(CMD_PARAMSTEMP, CommandSpecifier(CMD_PARAMSTEMP), "paramstemp");
	commands[CMD_PARAMSTEMP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>paramstemp</b> <i>(meshname)</i>";
	commands[CMD_PARAMSTEMP].descr = "[tc0,0.5,0.5,1/tc]List all material parameters temperature dependence. If meshname not given use the active mesh.";

	commands.insert(CMD_CLEARPARAMSTEMP, CommandSpecifier(CMD_CLEARPARAMSTEMP), "clearparamstemp");
	commands[CMD_CLEARPARAMSTEMP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearparamstemp</b> <i>(meshname, (paramname))</i>";
	commands[CMD_CLEARPARAMSTEMP].descr = "[tc0,0.5,0.5,1/tc]Clear material parameter temperature dependence in given mesh. If meshname not given clear temperature dependences in all meshes for all parameters. If paramname not given clear all parameters temperature dependences in named mesh.";
	
	commands.insert(CMD_SETPARAMTEMPEQUATION, CommandSpecifier(CMD_SETPARAMTEMPEQUATION), "setparamtempequation");
	commands[CMD_SETPARAMTEMPEQUATION].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparamtempequation</b> <i>meshname paramname text_equation</i>";
	commands[CMD_SETPARAMTEMPEQUATION].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter temperature dependence equation for the named mesh.";

	commands.insert(CMD_SETPARAMTEMPARRAY, CommandSpecifier(CMD_SETPARAMTEMPARRAY), "setparamtemparray");
	commands[CMD_SETPARAMTEMPARRAY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparamtemparray</b> <i>meshname paramname filename</i>";
	commands[CMD_SETPARAMTEMPARRAY].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter temperature dependence using an array in the given mesh. This must contain temperature values and scaling coefficients. Load directly from a file (tab spaced).";

	commands.insert(CMD_COPYPARAMS, CommandSpecifier(CMD_COPYPARAMS), "copyparams");
	commands[CMD_COPYPARAMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>copyparams</b> <i>meshname_from meshname_to (...)</i>";
	commands[CMD_COPYPARAMS].descr = "[tc0,0.5,0.5,1/tc]Copy all mesh parameters from first mesh to all other meshes given - all meshes must be of same type.";

	commands.insert(CMD_COPYMESHDATA, CommandSpecifier(CMD_COPYMESHDATA), "copymeshdata");
	commands[CMD_COPYMESHDATA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>copymeshdata</b> <i>meshname_from meshname_to (...)</i>";
	commands[CMD_COPYMESHDATA].descr = "[tc0,0.5,0.5,1/tc]Copy all primary mesh data (e.g. magnetization values and shape) from first mesh to all other meshes given - all meshes must be of same type.";

	commands.insert(CMD_PARAMSVAR, CommandSpecifier(CMD_PARAMSVAR), "paramsvar");
	commands[CMD_PARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>paramsvar</b> <i>(meshname)</i>";
	commands[CMD_PARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]List all material parameters spatial variation. If meshname not given use the active mesh.";

	commands.insert(CMD_SETDISPLAYEDPARAMSVAR, CommandSpecifier(CMD_SETDISPLAYEDPARAMSVAR), "setdisplayedparamsvar");
	commands[CMD_SETDISPLAYEDPARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdisplayedparamsvar</b> <i>meshname paramname</i>";
	commands[CMD_SETDISPLAYEDPARAMSVAR].limits = { { int(0), Any() } };
	commands[CMD_SETDISPLAYEDPARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]Set param to display for given mesh when ParamVar display is enabled (to show spatial variation if any).";

	commands.insert(CMD_CLEARPARAMSVAR, CommandSpecifier(CMD_CLEARPARAMSVAR), "clearparamsvar");
	commands[CMD_CLEARPARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearparamsvar</b> <i>(meshname)</i>";
	commands[CMD_CLEARPARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]Clear all material parameters spatial dependence in given mesh. If meshname not given clear spatial dependence in all meshes.";

	commands.insert(CMD_CLEARPARAMVAR, CommandSpecifier(CMD_CLEARPARAMVAR), "clearparamvar");
	commands[CMD_CLEARPARAMVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearparamvar</b> <i>meshname paramname</i>";
	commands[CMD_CLEARPARAMVAR].descr = "[tc0,0.5,0.5,1/tc]Clear parameter spatial dependence in given mesh.";
	
	commands.insert(CMD_PARAMSVAR, CommandSpecifier(CMD_PARAMSVAR), "paramsvar");
	commands[CMD_PARAMSVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>paramsvar</b> <i>(meshname)</i>";
	commands[CMD_PARAMSVAR].descr = "[tc0,0.5,0.5,1/tc]List all material parameters spatial variation. If meshname not given use the active mesh.";

	commands.insert(CMD_SETPARAMVAR, CommandSpecifier(CMD_SETPARAMVAR), "setparamvar");
	commands[CMD_SETPARAMVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setparamvar</b> <i>meshname paramname generatorname (arguments...)</i>";
	commands[CMD_SETPARAMVAR].descr = "[tc0,0.5,0.5,1/tc]Set the named parameter spatial dependence for the named mesh using the given generator (including any required arguments for the generator - if not given, default values are used).";

	commands.insert(CMD_SAVESIM, CommandSpecifier(CMD_SAVESIM), "savesim");
	commands[CMD_SAVESIM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savesim</b> <i>(directory/)filename</i>";
	commands[CMD_SAVESIM].descr = "[tc0,0.5,0.5,1/tc]Save simulation with given name. If no name given, the last saved/loaded file name will be used.";

	commands.insert(CMD_LOADSIM, CommandSpecifier(CMD_LOADSIM), "loadsim");
	commands[CMD_LOADSIM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadsim</b> <i>(directory/)filename</i>";
	commands[CMD_LOADSIM].descr = "[tc0,0.5,0.5,1/tc]Load simulation with given name.";

	commands.insert(CMD_DEFAULT, CommandSpecifier(CMD_DEFAULT), "default");
	commands[CMD_DEFAULT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>default</b>";
	commands[CMD_DEFAULT].descr = "[tc0,0.5,0.5,1/tc]Reset program to default state.";

	commands.insert(CMD_DISPLAY, CommandSpecifier(CMD_DISPLAY), "display");
	commands[CMD_DISPLAY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>display</b> <i>name (meshname)</i>";
	commands[CMD_DISPLAY].descr = "[tc0,0.5,0.5,1/tc]Change quantity to display for given mesh (active mesh if name not given).";

	commands.insert(CMD_DISPLAYDETAILLEVEL, CommandSpecifier(CMD_DISPLAYDETAILLEVEL), "displaydetail");
	commands[CMD_DISPLAYDETAILLEVEL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaydetail</b> <i>size</i>";
	commands[CMD_DISPLAYDETAILLEVEL].descr = "[tc0,0.5,0.5,1/tc]Change displayed detail level to given displayed cell size (m).";
	commands[CMD_DISPLAYDETAILLEVEL].unit = "m";

	commands.insert(CMD_DISPLAYBACKGROUND, CommandSpecifier(CMD_DISPLAYBACKGROUND), "displaybackground");
	commands[CMD_DISPLAYBACKGROUND].usage = "[tc0,0.5,0,1/tc]USAGE : <b>displaybackground</b> <i>name (meshname)</i>";
	commands[CMD_DISPLAYBACKGROUND].descr = "[tc0,0.5,0.5,1/tc]Change background quantity to display for given mesh (active mesh if name not given).";

	commands.insert(CMD_VECREP, CommandSpecifier(CMD_VECREP), "vecrep");
	commands[CMD_VECREP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>vecrep</b> <i>meshname vecreptype</i>";
	commands[CMD_VECREP].limits = { { Any(), Any() }, { int(0), int(VEC3REP_NUMOPTIONS) } };
	commands[CMD_VECREP].descr = "[tc0,0.5,0.5,1/tc]Set representation type for vectorial quantities in named mesh (or supermesh). vecreptype = 0 (full), vecreptype = 1 (x component), vecreptype = 2 (y component), vecreptype = 3 (z component), vecreptype = 4 (direction only), vecreptype = 5 (magnitude only).";

	commands.insert(CMD_SAVEMESHIMAGE, CommandSpecifier(CMD_SAVEMESHIMAGE), "savemeshimage");
	commands[CMD_SAVEMESHIMAGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>savemeshimage</b> <i>((directory/)filename)</i>";
	commands[CMD_SAVEMESHIMAGE].descr = "[tc0,0.5,0.5,1/tc]Save currently displayed mesh image to given file (as .png). If directory not specified then default directory is used. If filename not specified then default image save file name is used.";

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
	commands[CMD_PREPAREMOVINGMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or active mesh) for moving transverse (or vortex) domain wall simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";

	commands.insert(CMD_PREPAREMOVINGBLOCHMESH, CommandSpecifier(CMD_PREPAREMOVINGBLOCHMESH), "blochpreparemovingmesh");
	commands[CMD_PREPAREMOVINGBLOCHMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>blochpreparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGBLOCHMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or active mesh) for moving Bloch domain wall simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";

	commands.insert(CMD_PREPAREMOVINGNEELMESH, CommandSpecifier(CMD_PREPAREMOVINGNEELMESH), "neelpreparemovingmesh");
	commands[CMD_PREPAREMOVINGNEELMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>neelpreparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGNEELMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or active mesh) for moving Neel domain wall simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";

	commands.insert(CMD_PREPAREMOVINGSKYRMIONMESH, CommandSpecifier(CMD_PREPAREMOVINGSKYRMIONMESH), "skyrmionpreparemovingmesh");
	commands[CMD_PREPAREMOVINGSKYRMIONMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyrmionpreparemovingmesh</b> <i>(meshname)</i>";
	commands[CMD_PREPAREMOVINGSKYRMIONMESH].descr = "[tc0,0.5,0.5,1/tc]Setup the named mesh (or active mesh) for moving skyrmion simulations: 1) set movingmesh trigger, 2) set domain wall structure, 3) set dipoles left and right to remove end magnetic charges, 4) enable strayfield module.";

	commands.insert(CMD_COUPLETODIPOLES, CommandSpecifier(CMD_COUPLETODIPOLES), "coupletodipoles");
	commands[CMD_COUPLETODIPOLES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>coupletodipoles</b> <i>status</i>";
	commands[CMD_COUPLETODIPOLES].descr = "[tc0,0.5,0.5,1/tc]Set/unset coupling to dipoles : if ferromagnetic meshes touch a dipole mesh then interface magnetic cells are exchange coupled to the dipole magnetization direction.";

	commands.insert(CMD_EXCHANGECOUPLEDMESHES, CommandSpecifier(CMD_EXCHANGECOUPLEDMESHES), "exchangecoupledmeshes");
	commands[CMD_EXCHANGECOUPLEDMESHES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>exchangecoupledmeshes</b> <i>status (meshname)</i>";
	commands[CMD_EXCHANGECOUPLEDMESHES].descr = "[tc0,0.5,0.5,1/tc]Set/unset direct exchange coupling to neighboring meshes : if neighboring ferromagnetic meshes touch the named mesh (set for focused mesh if meshname not given) then interface magnetic cells are direct exchange coupled to them.";
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
	commands[CMD_SETDEFAULTELECTRODES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setdefaultelectrodes</b>";
	commands[CMD_SETDEFAULTELECTRODES].descr = "[tc0,0.5,0.5,1/tc]Set electrodes at the x-axis ends of the given mesh, both set at 0V. Set the left-side electrode as the ground. Delete all other electrodes.";

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

	commands.insert(CMD_TSOLVERCONFIG, CommandSpecifier(CMD_TSOLVERCONFIG), "tsolverconfig");
	commands[CMD_TSOLVERCONFIG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tsolverconfig</b> <i>convergence_error (iters_timeout)</i>";
	commands[CMD_TSOLVERCONFIG].limits = { { double(1e-10), double(1e-1) }, { int(1), int(500000) } };
	commands[CMD_TSOLVERCONFIG].descr = "[tc0,0.5,0.5,1/tc]Set transport solver convergence error and iterations for timeout (if given, else use default).";
	commands[CMD_TSOLVERCONFIG].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>convergence_error iters_timeout</i>";

	commands.insert(CMD_SSOLVERCONFIG, CommandSpecifier(CMD_SSOLVERCONFIG), "ssolverconfig");
	commands[CMD_SSOLVERCONFIG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ssolverconfig</b> <i>s_convergence_error (s_iters_timeout)</i>";
	commands[CMD_SSOLVERCONFIG].limits = { { double(1e-10), double(1e-1) },{ int(1), int(50000) } };
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

	commands.insert(CMD_TEMPERATURE, CommandSpecifier(CMD_TEMPERATURE), "temperature");
	commands[CMD_TEMPERATURE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>temperature</b> <i>value (meshname)</i>";
	commands[CMD_TEMPERATURE].limits = { { double(0.0), double(MAX_TEMPERATURE) }, { Any(), Any() } };
	commands[CMD_TEMPERATURE].descr = "[tc0,0.5,0.5,1/tc]Set mesh base temperature (all meshes if meshname not given) and reset temperature. Also set ambient temperature if Heat module added. If the base temperature setting has a spatial dependence specified through cT, this command will take it into account but only if the Heat module is added. If you want the temperature to remain fixed you can still have the Heat module enabled but disable the heat equation by setting the heat dT to zero (setheatdt 0).";
	commands[CMD_TEMPERATURE].unit = "K";
	commands[CMD_TEMPERATURE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - temperature value for mesh in focus.";

	commands.insert(CMD_SETHEATDT, CommandSpecifier(CMD_SETHEATDT), "setheatdt");
	commands[CMD_SETHEATDT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>setheatdt</b> <i>value</i>";
	commands[CMD_SETHEATDT].limits = { { double(0), double(MAXTIMESTEP) } };
	commands[CMD_SETHEATDT].descr = "[tc0,0.5,0.5,1/tc]Set heat equation solver time step.";
	commands[CMD_SETHEATDT].unit = "s";
	commands[CMD_SETHEATDT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - heat equation time step.";

	commands.insert(CMD_AMBIENTTEMPERATURE, CommandSpecifier(CMD_AMBIENTTEMPERATURE), "ambient");
	commands[CMD_AMBIENTTEMPERATURE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>ambient</b> <i>ambient_temperature (meshname)</i>";
	commands[CMD_AMBIENTTEMPERATURE].limits = { { double(0.0), double(MAX_TEMPERATURE) }, { Any(), Any() } };
	commands[CMD_AMBIENTTEMPERATURE].descr = "[tc0,0.5,0.5,1/tc]Set mesh ambient temperature (all meshes if meshname not given) for Robin boundary conditions : flux normal = alpha * (T_boundary - T_ambient).";
	commands[CMD_AMBIENTTEMPERATURE].unit = "K";
	commands[CMD_AMBIENTTEMPERATURE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>ambient_temperature</i> - ambient temperature for mesh in focus.";

	commands.insert(CMD_ROBINALPHA, CommandSpecifier(CMD_ROBINALPHA), "robinalpha");
	commands[CMD_ROBINALPHA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>robinalpha</b> <i>robin_alpha (meshname)</i>";
	commands[CMD_ROBINALPHA].limits = { { double(0.0), double(1e10) }, { Any(), Any() } };
	commands[CMD_ROBINALPHA].descr = "[tc0,0.5,0.5,1/tc]Set alpha coefficient (all meshes if meshname not given) for Robin boundary conditions : flux normal = alpha * (T_boundary - T_ambient).";
	commands[CMD_ROBINALPHA].unit = "W/m2K";
	commands[CMD_ROBINALPHA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>robin_alpha</i> - Robin alpha value for mesh in focus.";

	commands.insert(CMD_INSULATINGSIDES, CommandSpecifier(CMD_INSULATINGSIDES), "insulatingside");
	commands[CMD_INSULATINGSIDES].usage = "[tc0,0.5,0,1/tc]USAGE : <b>insulatingside</b> <i>side_literal status (meshname)</i>";
	commands[CMD_INSULATINGSIDES].limits = { { Any(), Any() }, { bool(0), bool(1) }, { Any(), Any() } };
	commands[CMD_INSULATINGSIDES].descr = "[tc0,0.5,0.5,1/tc]Set temperature insulation (Neumann boundary condition) for named mesh side (active mesh if not given). side_literal : x, -x, y, -y, z, -z.";
	commands[CMD_INSULATINGSIDES].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>status_x status_-x status_y status_-y status_z status_-z</i> - insulating sides status for mesh in focus.";
	
	commands.insert(CMD_CURIETEMPERATURE, CommandSpecifier(CMD_CURIETEMPERATURE), "curietemperature");
	commands[CMD_CURIETEMPERATURE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>curietemperature</b> <i>curie_temperature (meshname)</i>";
	commands[CMD_CURIETEMPERATURE].limits = { { double(0.0), double(MAX_TEMPERATURE) }, { Any(), Any() } };
	commands[CMD_CURIETEMPERATURE].descr = "[tc0,0.5,0.5,1/tc]Set Curie temperature (all ferromagnetic meshes if meshname not given) for ferromagnetic mesh. This will set default temperature dependencies as: Ms = Ms0*me, A = Ah*me^2, D = D0*me^2, K = K0*me^3 (K1 and K2), damping = damping0*(1-T/3Tc) T < Tc, damping = damping0*2T/3Tc T >= Tc, susrel = dme/d(mu0Hext). Setting the Curie temperature to zero will disable temperature dependence for these parameters.";
	commands[CMD_CURIETEMPERATURE].unit = "K";
	commands[CMD_CURIETEMPERATURE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>curie_temperature</i> - Curie temperature for mesh in focus.";

	commands.insert(CMD_CURIETEMPERATUREMATERIAL, CommandSpecifier(CMD_CURIETEMPERATUREMATERIAL), "matcurietemperature");
	commands[CMD_CURIETEMPERATUREMATERIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>matcurietemperature</b> <i>curie_temperature (meshname)</i>";
	commands[CMD_CURIETEMPERATUREMATERIAL].limits = { { double(0.0), double(MAX_TEMPERATURE) },{ Any(), Any() } };
	commands[CMD_CURIETEMPERATUREMATERIAL].descr = "[tc0,0.5,0.5,1/tc]Set indicative material Curie temperature for ferromagnetic mesh (focused ferromagnetic mesh if meshname not given). This is not used in calculations, but serves as an indicative value - set the actual Tc value with the <b>curietemperature</b> command.";
	commands[CMD_CURIETEMPERATUREMATERIAL].unit = "K";
	commands[CMD_CURIETEMPERATUREMATERIAL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>curie_temperature</i> - Indicative material Curie temperature for mesh in focus.";
	
	commands.insert(CMD_ATOMICMOMENT, CommandSpecifier(CMD_ATOMICMOMENT), "atomicmoment");
	commands[CMD_ATOMICMOMENT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>atomicmoment</b> <i>ub_multiple (meshname)</i>";
	commands[CMD_ATOMICMOMENT].limits = { { double(0.0), Any() },{ Any(), Any() } };
	commands[CMD_ATOMICMOMENT].descr = "[tc0,0.5,0.5,1/tc]Set atomic moment as a multiple of Bohr magnetons (all applicable meshes if meshname not given) for given mesh. This affects the temperature dependence of 'me' (see curietemperature command). A non-zero value will result in me(T) being dependent on the applied field.";
	commands[CMD_ATOMICMOMENT].unit = "uB";
	commands[CMD_ATOMICMOMENT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>ub_multiple</i> - atomic moment multiple of Bohr magneton for mesh in focus.";

	commands.insert(CMD_TAU, CommandSpecifier(CMD_TAU), "tau");
	commands[CMD_TAU].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tau</b> <i>tau_11 tau_22 (tau_12 tau_21) (meshname)</i>";
	commands[CMD_TAU].limits = { { double(0.0),  Any() }, { double(0.0),  Any() }, { double(0.0),  Any() }, { double(0.0),  Any() } };
	commands[CMD_TAU].descr = "[tc0,0.5,0.5,1/tc]Set ratio of exchange parameters to critical temperature (Neel) (all antiferromagnetic meshes if meshname not given) for antiferromagnetic mesh. tau_11 and tau_22 are the intra-lattice contributions, tau_12 and tau_21 are the inter-lattice contributions.";
	commands[CMD_TAU].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>tau_11 tau_22 tau_12 tau_21</i>";

	commands.insert(CMD_TMODEL, CommandSpecifier(CMD_TMODEL), "tmodel");
	commands[CMD_TMODEL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>tmodel</b> <i>num_temperatures (meshname)</i>";
	commands[CMD_TMODEL].limits = { { int(1), Any() } };
	commands[CMD_TMODEL].descr = "[tc0,0.5,0.5,1/tc]Set temperature model (determined by number of temperatures) in given meshname (focused mesh if meshname not given). Note insulating meshes only allow a 1-temperature model.";

	commands.insert(CMD_STOCHASTIC, CommandSpecifier(CMD_STOCHASTIC), "stochastic");
	commands[CMD_STOCHASTIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>stochastic</b>";
	commands[CMD_STOCHASTIC].descr = "[tc0,0.5,0.5,1/tc]Shows stochasticity settings : stochastic cellsize for each mesh and related settings.";

	commands.insert(CMD_LINKSTOCHASTIC, CommandSpecifier(CMD_LINKSTOCHASTIC), "linkstochastic");
	commands[CMD_LINKSTOCHASTIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>linkstochastic</b> <i>flag (meshname)</i>";
	commands[CMD_LINKSTOCHASTIC].limits = { { int(0), int(1) },{ Any(), Any() } };
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
	commands[CMD_SETDTSPEEDUP].descr = "[tc0,0.5,0.5,1/tc]Set time step for evaluation speedup, to be used in when in extreme mode.";
	commands[CMD_SETDTSPEEDUP].unit = "s";
	commands[CMD_SETDTSPEEDUP].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>dTspeedup</i>";

	commands.insert(CMD_LINKDTSPEEDUP, CommandSpecifier(CMD_LINKDTSPEEDUP), "linkdtspeedup");
	commands[CMD_LINKDTSPEEDUP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>linkdtspeedup</b> <i>flag</i>";
	commands[CMD_LINKDTSPEEDUP].limits = { { int(0), int(1) } };
	commands[CMD_LINKDTSPEEDUP].descr = "[tc0,0.5,0.5,1/tc]Links speedup time-step to ODE time-step if set, else speedup time-step is independently controlled. Applicable in extreme mode only.";

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
	commands[CMD_REFINEROUGHNESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>refineroughness</b> <i>value (meshname)</i>";
	commands[CMD_REFINEROUGHNESS].limits = { { INT3(1, 1, 1), Any() } };
	commands[CMD_REFINEROUGHNESS].descr = "[tc0,0.5,0.5,1/tc]Set roughness refinement cellsize divider in given mesh, i.e. cellsize used for roughness initialization is the ferromagnetic cellsize divided by value (3 components, so divide component by component).";
	commands[CMD_REFINEROUGHNESS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i> - roughness refinement.";

	commands.insert(CMD_CLEARROUGHNESS, CommandSpecifier(CMD_CLEARROUGHNESS), "clearroughness");
	commands[CMD_CLEARROUGHNESS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>clearroughness</b> <i>(meshname)</i>";
	commands[CMD_CLEARROUGHNESS].limits = { { Any(), Any() } };
	commands[CMD_CLEARROUGHNESS].descr = "[tc0,0.5,0.5,1/tc]Clear roughness by setting the fine shape same as the coarse M shape.";

	commands.insert(CMD_ROUGHENMESH, CommandSpecifier(CMD_ROUGHENMESH), "roughenmesh");
	commands[CMD_ROUGHENMESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>roughenmesh</b> <i>depth (side, (seed))</i>";
	commands[CMD_ROUGHENMESH].limits = { { double(0), Any() }, { Any(), Any() }, { int(1), Any() } };
	commands[CMD_ROUGHENMESH].unit = "m";
	commands[CMD_ROUGHENMESH].descr = "[tc0,0.5,0.5,1/tc]Roughen active mesh to given depth (m) on named side (use side = x, y, z, -x, -y, -z as literal, z by default). The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_SURFROUGHENJAGGED, CommandSpecifier(CMD_SURFROUGHENJAGGED), "surfroughenjagged");
	commands[CMD_SURFROUGHENJAGGED].usage = "[tc0,0.5,0,1/tc]USAGE : <b>surfroughenjagged</b> <i>depth spacing (seed, (sides))</i>";
	commands[CMD_SURFROUGHENJAGGED].limits = { { double(0), Any() },{ double(0), Any() },{ int(1), Any() },{ Any(), Any() } };
	commands[CMD_SURFROUGHENJAGGED].unit = "m";
	commands[CMD_SURFROUGHENJAGGED].descr = "[tc0,0.5,0.5,1/tc]Roughen active mesh surfaces using a jagged pattern to given depth (m) and peak spacing (m). Roughen both sides by default, unless sides is specified as -z or z (string literal). The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_GENERATE2DGRAINS, CommandSpecifier(CMD_GENERATE2DGRAINS), "generate2dgrains");
	commands[CMD_GENERATE2DGRAINS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>generate2dgrains</b> <i>spacing (seed)</i>";
	commands[CMD_GENERATE2DGRAINS].limits = { { double(0), Any() } };
	commands[CMD_GENERATE2DGRAINS].unit = "m";
	commands[CMD_GENERATE2DGRAINS].descr = "[tc0,0.5,0.5,1/tc]Generate 2D Voronoi cells in the xy plane at given average spacing. The seed is used for the pseudo-random number generator, 1 by default.";

	commands.insert(CMD_GENERATE3DGRAINS, CommandSpecifier(CMD_GENERATE3DGRAINS), "generate3dgrains");
	commands[CMD_GENERATE3DGRAINS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>generate3dgrains</b> <i>spacing (seed)</i>";
	commands[CMD_GENERATE3DGRAINS].limits = { { double(0), Any() } };
	commands[CMD_GENERATE3DGRAINS].unit = "m";
	commands[CMD_GENERATE3DGRAINS].descr = "[tc0,0.5,0.5,1/tc]Generate 3D Voronoi cells at given average spacing. The seed is used for the pseudo-random number generator, 1 by default.";

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
	commands[CMD_SHOWLENGHTS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showlengths</b>";
	commands[CMD_SHOWLENGHTS].descr = "[tc0,0.5,0.5,1/tc]Calculate a number of critical lengths for the focused mesh (must be ferromagnetic) to inform magnetization cellsize selection. lex = sqrt(2 A / mu0 Ms^2) : exchange length, l_Bloch = sqrt(A / K1) : Bloch wall width, l_sky = PI D / 4 K1 : Neel skyrmion wall width.";

	commands.insert(CMD_SHOWMCELLS, CommandSpecifier(CMD_SHOWMCELLS), "showmcells");
	commands[CMD_SHOWMCELLS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showmcells</b>";
	commands[CMD_SHOWMCELLS].descr = "[tc0,0.5,0.5,1/tc]Show number of discretisation cells for magnetization for focused mesh (must be ferromagnetic).";
	commands[CMD_SHOWMCELLS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>n</i>";

	commands.insert(CMD_LOADOVF2MESH, CommandSpecifier(CMD_LOADOVF2MESH), "loadovf2mesh");
	commands[CMD_LOADOVF2MESH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2mesh</b> <i>(renormalize_value) (directory/)filename</i>";
	commands[CMD_LOADOVF2MESH].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing 3-component vector data. This will create a new permalloy ferromagnetic mesh with dimensions and magnetization data obtained from the OVF 2.0 file. By default the loaded data will not be renormalized: renormalize_value = 0. If a value is specified for renormalize_value, the loaded data will be renormalized to it (e.g. this would be an Ms value).";

	commands.insert(CMD_LOADOVF2MAG, CommandSpecifier(CMD_LOADOVF2MAG), "loadovf2mag");
	commands[CMD_LOADOVF2MAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2mag</b> <i>(renormalize_value) (directory/)filename</i>";
	commands[CMD_LOADOVF2MAG].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing magnetization data, into the currently focused mesh (which must be ferromagnetic), mapping the data to the current mesh dimensions. By default the loaded data will not be renormalized: renormalize_value = 0. If a value is specified for renormalize_value, the loaded data will be renormalized to it (e.g. this would be an Ms value).";
	commands[CMD_LOADOVF2MAG].unit = "A/m";

	commands.insert(CMD_SAVEOVF2MAG, CommandSpecifier(CMD_SAVEOVF2MAG), "saveovf2mag");
	commands[CMD_SAVEOVF2MAG].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveovf2mag</b> <i>(n) (data_type) (directory/)filename</i>";
	commands[CMD_SAVEOVF2MAG].descr = "[tc0,0.5,0.5,1/tc]Save an OOMMF-style OVF 2.0 file containing magnetization data from the currently focused mesh (which must be ferromagnetic). You can normalize the data to Ms0 value by specifying the n flag (e.g. saveovf2mag n filename) - by default the data is not normalized. You can specify the data type as data_type = bin4 (single precision 4 bytes per float), data_type = bin8 (double precision 8 bytes per float), or data_type = text. By default bin8 is used.";

	commands.insert(CMD_SAVEOVF2PARAMVAR, CommandSpecifier(CMD_SAVEOVF2PARAMVAR), "saveovf2param");
	commands[CMD_SAVEOVF2PARAMVAR].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveovf2param</b> <i>(data_type) (meshname) paramname (directory/)filename</i>";
	commands[CMD_SAVEOVF2PARAMVAR].descr = "[tc0,0.5,0.5,1/tc]Save an OOMMF-style OVF 2.0 file containing the named parameter spatial variation data from the named mesh (currently focused mesh if not specified). You can specify the data type as data_type = bin4 (single precision 4 bytes per float), data_type = bin8 (double precision 8 bytes per float), or data_type = text. By default bin8 is used.";
	
	commands.insert(CMD_SAVEOVF2, CommandSpecifier(CMD_SAVEOVF2), "saveovf2");
	commands[CMD_SAVEOVF2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>saveovf2</b> <i>(data_type) (directory/)filename</i>";
	commands[CMD_SAVEOVF2].descr = "[tc0,0.5,0.5,1/tc]Save an OOMMF-style OVF 2.0 file containing data from the currently focused mesh. You can specify the data type as data_type = bin4 (single precision 4 bytes per float), data_type = bin8 (double precision 8 bytes per float), or data_type = text. By default bin8 is used.";

	commands.insert(CMD_LOADOVF2DISP, CommandSpecifier(CMD_LOADOVF2DISP), "loadovf2disp");
	commands[CMD_LOADOVF2DISP].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2disp</b> <i>(directory/)filename</i>";
	commands[CMD_LOADOVF2DISP].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing mechanical displacement data, into the currently focused mesh (which must be ferromagnetic and have the melastic module enabled), mapping the data to the current mesh dimensions. From the mechanical displacement the strain tensor is calculated.";

	commands.insert(CMD_LOADOVF2STRAIN, CommandSpecifier(CMD_LOADOVF2STRAIN), "loadovf2strain");
	commands[CMD_LOADOVF2STRAIN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>loadovf2strain</b> <i>(directory/)filename_diag filename_odiag</i>";
	commands[CMD_LOADOVF2STRAIN].descr = "[tc0,0.5,0.5,1/tc]Load an OOMMF-style OVF 2.0 file containing strain tensor data, into the currently focused mesh (which must be ferromagnetic and have the melastic module enabled), mapping the data to the current mesh dimensions. The symmetric strain tensor is applicable for a cubic crystal, and has 3 diagonal component (specified in filename_diag with vector data as xx, yy, zz), and 3 off-diagonal components (specified in filename_odiag with vector data as yz, xz, xy).";

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

	commands.insert(CMD_SHOWTC, CommandSpecifier(CMD_SHOWTC), "showtc");
	commands[CMD_SHOWTC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showtc</b>";
	commands[CMD_SHOWTC].descr = "[tc0,0.5,0.5,1/tc]Show predicted Tc value (K) for current mesh in focus (must be atomistic), using formula Tc = J*e*z/3kB, where e is the spin-wave correction factor, and z is the coordination number.";
	commands[CMD_SHOWTC].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Tc</i>";

	commands.insert(CMD_SHOWMS, CommandSpecifier(CMD_SHOWMS), "showms");
	commands[CMD_SHOWMS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showms</b>";
	commands[CMD_SHOWMS].descr = "[tc0,0.5,0.5,1/tc]Show predicted saturation magnetization (A/m) value for current mesh in focus (must be atomistic), using formula Ms = mu_s*n/a^3, where n is the number of atomic moments per unit cell, and a is the atomic cell size.";
	commands[CMD_SHOWMS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Ms</i>";

	commands.insert(CMD_SHOWA, CommandSpecifier(CMD_SHOWA), "showa");
	commands[CMD_SHOWA].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showa</b>";
	commands[CMD_SHOWA].descr = "[tc0,0.5,0.5,1/tc]Show predicted exchange stiffness (J/m) value for current mesh in focus (must be atomistic), using formula A = J*n/2a, where n is the number of atomic moments per unit cell, and a is the atomic cell size.";
	commands[CMD_SHOWA].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>A</i>";

	commands.insert(CMD_SHOWK, CommandSpecifier(CMD_SHOWK), "showk");
	commands[CMD_SHOWK].usage = "[tc0,0.5,0,1/tc]USAGE : <b>showk</b>";
	commands[CMD_SHOWK].descr = "[tc0,0.5,0.5,1/tc]Show predicted uniaxial anisotropy (J/m^3) constant value for current mesh in focus (must be atomistic), using formula K = k*n/a^3, where n is the number of atomic moments per unit cell, and a is the atomic cell size.";
	commands[CMD_SHOWK].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>A</i>";

	commands.insert(CMD_SKYPOSDMUL, CommandSpecifier(CMD_SKYPOSDMUL), "skyposdmul");
	commands[CMD_SKYPOSDMUL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>skyposdmul</b> <i>multiplier (meshname)</i>";
	commands[CMD_SKYPOSDMUL].descr = "[tc0,0.5,0.5,1/tc]Set skyrmion diameter multiplier for given meshname (focused mesh if not given) which determines skypos tracking rectangle size. Default is 2.0, i.e. tracking rectangle side is twice the diameter along each axis. Reduce if skyrmion density is too large, but at the risk of losing skyrmion tracking (recommended to keep above 1.2).";
	commands[CMD_SKYPOSDMUL].limits = { { double(1.0), double(10.0) }, {Any(), Any()} };
	commands[CMD_SKYPOSDMUL].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>multiplier</i>";

	commands.insert(CMD_MCSERIAL, CommandSpecifier(CMD_MCSERIAL), "mcserial");
	commands[CMD_MCSERIAL].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcserial</b> <i>value (meshname)</i>";
	commands[CMD_MCSERIAL].descr = "[tc0,0.5,0.5,1/tc]Change Monte-Carlo algorithm type. 0: parallel (default); red-black ordering parallelization 1: serial (testing only); spins picked in random order. If meshname not specified setting is applied to all atomistic meshes.";
	commands[CMD_MCSERIAL].limits = { { int(0), int(1) }, {Any(), Any()} };

	commands.insert(CMD_MCCONSTRAIN, CommandSpecifier(CMD_MCCONSTRAIN), "mcconstrain");
	commands[CMD_MCCONSTRAIN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>mcconstrain</b> <i>value (meshname)</i>";
	commands[CMD_MCCONSTRAIN].descr = "[tc0,0.5,0.5,1/tc]Set value 0 to revert to classic Monte-Carlo Metropolis for ASD. Set a unit vector direction value (x y z) to switch to constrained Monte Carlo as described in PRB 82, 054415 (2010). If meshname not specified setting is applied to all atomistic meshes.";
	commands[CMD_MCCONSTRAIN].limits = { { DBL3(), Any() }, {Any(), Any()} };

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
	commands[CMD_DP_GETPROFILE].descr = "[tc0,0.5,0.5,1/tc]Extract profile of physical quantity displayed on screen, at the current display resolution, along the line specified with given start and end cartesian absolute coordinates (m). Place profile in given dp arrays: up to 4 consecutive dp arrays are used, first for distance along line, the next 3 for physical quantity components (e.g. Mx, My, Mz) so allow space for these starting at dp_index.";
	commands[CMD_DP_GETPROFILE].unit = "m";

	commands.insert(CMD_DP_GETEXACTPROFILE, CommandSpecifier(CMD_DP_GETEXACTPROFILE), "dp_getexactprofile");
	commands[CMD_DP_GETEXACTPROFILE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_getexactprofile</b> <i>start end step dp_index (stencil)</i>";
	commands[CMD_DP_GETEXACTPROFILE].limits = { { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) }, { MINMESHSPACE, Any() }, { int(0), int(MAX_ARRAYS - 1) }, { DBL3(MINMESHSPACE), Any() } };
	commands[CMD_DP_GETEXACTPROFILE].descr = "[tc0,0.5,0.5,1/tc]Extract profile of physical quantity displayed on screen, directly from the mesh so using the exact mesh resolution not the displayed resolution, along the line specified with given start and end cartesian absolute coordinates (m), and with the given step size (m). If stencil specified - as x y z (m) - then obtain profile values using weighted averaging with stencil centered on profile point. Place profile in given dp arrays: up to 3 consecutive dp arrays are used for physical quantity components (e.g. Mx, My, Mz) so allow space for these starting at dp_index.";
	commands[CMD_DP_GETEXACTPROFILE].unit = "m";

	commands.insert(CMD_DP_GETPATH, CommandSpecifier(CMD_DP_GETPATH), "dp_getpath");
	commands[CMD_DP_GETPATH].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_getpath</b> <i>dp_index_in dp_index_out</i>";
	commands[CMD_DP_GETPATH].limits = { { int(0), int(MAX_ARRAYS - 3) }, { int(0), int(MAX_ARRAYS - 3) } };
	commands[CMD_DP_GETPATH].descr = "[tc0,0.5,0.5,1/tc]Extract profile of physical quantity displayed on screen, directly from stored mesh data thus independent of display resolution, along the path specified in Cartesian absolute coordinates (m) through dp arrays at dp_index_in, dp_index_in + 1, dp_index_in + 2 (x, y, z coordinates resp.). Place extracted profile in given dp arrays dp_index_out, dp_index_out + 1, dp_index_out + 2 (x, y, z components for vector data).";

	commands.insert(CMD_GETVALUE, CommandSpecifier(CMD_GETVALUE), "getvalue");
	commands[CMD_GETVALUE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>getvalue</b> <i>abspos</i>";
	commands[CMD_GETVALUE].limits = { { DBL3(-MAXSIMSPACE), DBL3(MAXSIMSPACE) } };
	commands[CMD_GETVALUE].descr = "[tc0,0.5,0.5,1/tc]Get data value at abspos (absolute position in Cartesian coordinates) depending on currently displayed quantities.";
	commands[CMD_GETVALUE].unit = "m";
	commands[CMD_GETVALUE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_AVERAGEMESHRECT, CommandSpecifier(CMD_AVERAGEMESHRECT), "averagemeshrect");
	commands[CMD_AVERAGEMESHRECT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>averagemeshrect</b> <i>(rectangle)</i>";
	commands[CMD_AVERAGEMESHRECT].descr = "[tc0,0.5,0.5,1/tc]Calculate the average value depending on currently displayed quantities. The rectangle is specified in relative coordinates to the currently focused mesh; if not specified average the entire focused mesh.";
	commands[CMD_AVERAGEMESHRECT].limits = { { Rect(DBL3(-MAXSIMSPACE / 2), DBL3(-MAXSIMSPACE / 2) + DBL3(MINMESHSPACE)), Rect(DBL3(-MAXSIMSPACE / 2), DBL3(MAXSIMSPACE / 2)) } };
	commands[CMD_AVERAGEMESHRECT].unit = "m";
	commands[CMD_AVERAGEMESHRECT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>value</i>";

	commands.insert(CMD_DP_TOPOCHARGE, CommandSpecifier(CMD_DP_TOPOCHARGE), "dp_topocharge");
	commands[CMD_DP_TOPOCHARGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_topocharge</b> <i>(x y radius)</i>";
	commands[CMD_DP_TOPOCHARGE].descr = "[tc0,0.5,0.5,1/tc]Calculate the topological charge for focused mesh (must be magnetic), optionally in the given circle with radius and centered at x y (relative values). Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI.";
	commands[CMD_DP_TOPOCHARGE].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Q</i> - the calculated topological charge.";
	commands[CMD_DP_TOPOCHARGE].limits = { {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_TOPOCHARGE].unit = "m";

	commands.insert(CMD_DP_COUNTSKYRMIONS, CommandSpecifier(CMD_DP_COUNTSKYRMIONS), "dp_countskyrmions");
	commands[CMD_DP_COUNTSKYRMIONS].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_countskyrmions</b> <i>(x y radius)</i>";
	commands[CMD_DP_COUNTSKYRMIONS].descr = "[tc0,0.5,0.5,1/tc]Calculate the number of skyrmions for focused mesh (must be magnetic), optionally in the given circle with radius and centered at x y (relative values). Use Qmag = Integral(|m.(dm/dx x dm/dy)| dxdy) / 4PI.";
	commands[CMD_DP_COUNTSKYRMIONS].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>Q</i> - the calculated topological charge.";
	commands[CMD_DP_COUNTSKYRMIONS].limits = { {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_COUNTSKYRMIONS].unit = "m";

	commands.insert(CMD_DP_HISTOGRAM, CommandSpecifier(CMD_DP_HISTOGRAM), "dp_histogram");
	commands[CMD_DP_HISTOGRAM].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_histogram</b> <i>dp_x dp_y (bin min max)</i>";
	commands[CMD_DP_HISTOGRAM].descr = "[tc0,0.5,0.5,1/tc]Calculate a histogram with given bin, minimum and maximum values, from the magnetization magnitude of the focused mesh (must be magnetic). Save histogram in dp arrays at dp_x, dp_y. If histogram parameters not given use a bin with 100 steps between minimum and maximum magnetization magnitude.";
	commands[CMD_DP_HISTOGRAM].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
	commands[CMD_DP_HISTOGRAM].unit = "A/m";

	commands.insert(CMD_DP_HISTOGRAM2, CommandSpecifier(CMD_DP_HISTOGRAM2), "dp_histogram2");
	commands[CMD_DP_HISTOGRAM2].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_histogram2</b> <i>dp_x dp_y (bin min max M2 deltaM2)</i>";
	commands[CMD_DP_HISTOGRAM2].descr = "[tc0,0.5,0.5,1/tc]Calculate a histogram for a 2-sublattice mesh with given bin, minimum and maximum values for sub-lattice A, if the corresponding magnetization magnitude in sub-lattice B equals M2 within the given deltaM2. Save histogram in dp arrays at dp_x, dp_y. If histogram parameters not given use a bin with 100 steps between minimum and maximum magnetization magnitude, with M2 set to MeB and deltaM2 set 0.01*MeB respectively.";
	commands[CMD_DP_HISTOGRAM2].limits = { { int(0), int(MAX_ARRAYS - 1) }, { int(0), int(MAX_ARRAYS - 1) }, {double(0), Any()}, {double(0), Any()}, {double(0), Any()}, {double(0), Any()}, {double(0), Any()} };
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
	commands[CMD_DP_MEAN].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_mean</b> <i>dp_index</i>";
	commands[CMD_DP_MEAN].limits = { { int(0), int(MAX_ARRAYS - 1) } };
	commands[CMD_DP_MEAN].descr = "[tc0,0.5,0.5,1/tc]Obtain mean value with standard deviation.";
	commands[CMD_DP_MEAN].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>mean stdev</i>.";

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

	commands.insert(CMD_DP_FITSTT, CommandSpecifier(CMD_DP_FITSTT), "dp_fitstt");
	commands[CMD_DP_FITSTT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitstt</b> <i>(rectangle)</i>";
	commands[CMD_DP_FITSTT].limits = { {Rect(), Any()} };
	commands[CMD_DP_FITSTT].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity). The fitting is done inside the specified rectangle for the focused mesh, with the rectangle specified using relative coordinates as sx sy sz ex ey ez (entire mesh if not specified). The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not.";
	commands[CMD_DP_FITSTT].unit = "m";
	commands[CMD_DP_FITSTT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>P, beta, std_P, std_beta, Rsq</i>.";

	commands.insert(CMD_DP_FITSOT, CommandSpecifier(CMD_DP_FITSOT), "dp_fitsot");
	commands[CMD_DP_FITSOT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitsot</b> <i>hm_mesh (rectangle)</i>";
	commands[CMD_DP_FITSOT].limits = { {Any(), Any()}, {Rect(), Any()} };
	commands[CMD_DP_FITSOT].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent interfacial spin torque using SOT with fitting parameters SHAeff and flST (field-like torque coefficient). hm_mesh specifies the heavy metal mesh from which to obtain the current density. The fitting is done inside the specified rectangle for the focused mesh, with the rectangle specified using relative coordinates as sx sy sz ex ey ez (entire mesh if not specified). The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and Tsi to have been computed.";
	commands[CMD_DP_FITSOT].unit = "m";
	commands[CMD_DP_FITSOT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>SHAeff, flST, std_SHAeff, std_flST, Rsq</i>.";

	commands.insert(CMD_DP_FITSOTSTT, CommandSpecifier(CMD_DP_FITSOTSTT), "dp_fitsotstt");
	commands[CMD_DP_FITSOTSTT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitsotstt</b> <i>hm_mesh (rectangle)</i>";
	commands[CMD_DP_FITSOTSTT].limits = { {Any(), Any()}, {Rect(), Any()} };
	commands[CMD_DP_FITSOTSTT].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity), and simultaneously also using SOT with fitting parameters SHAeff and flST (field-like torque coefficient). hm_mesh specifies the heavy metal mesh from which to obtain the current density for SOT. The fitting is done inside the specified rectangle for the focused mesh, with the rectangle specified using relative coordinates as sx sy sz ex ey ez (entire mesh if not specified). The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not.";
	commands[CMD_DP_FITSOTSTT].unit = "m";
	commands[CMD_DP_FITSOTSTT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>SHAeff, flST, P, beta, std_SHAeff, std_flST, std_P, std_beta, Rsq</i>.";

	commands.insert(CMD_DP_CALCSOT, CommandSpecifier(CMD_DP_CALCSOT), "dp_calcsot");
	commands[CMD_DP_CALCSOT].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_calcsot</b> <i>hm_mesh fm_mesh</i>";
	commands[CMD_DP_CALCSOT].limits = { {Any(), Any()}, {Any(), Any()} };
	commands[CMD_DP_CALCSOT].descr = "[tc0,0.5,0.5,1/tc]For the given heavy metal and ferromagnetic meshes calculate the expected effective spin Hall angle and field-like torque coefficient according to analytical equations (see manual).";
	commands[CMD_DP_CALCSOT].return_descr = "[tc0,0.5,0,1/tc]Script return values: <i>SHAeff, flST</i>.";

	commands.insert(CMD_DP_FITNONADIABATIC, CommandSpecifier(CMD_DP_FITNONADIABATIC), "dp_fitnonadiabatic");
	commands[CMD_DP_FITNONADIABATIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitnonadiabatic</b> <i>(abs_err Rsq T_ratio (stencil))</i>";
	commands[CMD_DP_FITNONADIABATIC].limits = { {double(0.0), Any()}, {double(0.0), Any()}, {double(0.0), Any()}, {int(3), Any()} };
	commands[CMD_DP_FITNONADIABATIC].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity) using a given square in-plane stencil (default size 3) in order to extract the spatial variation of beta. Cut-off values for absolute fitting error (default 0.1), Rsq measure (default 0.9), and normalized torque magnitude (default 0.1) can be set - value of zero disables cutoff. The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not. Output available in Cust_S.";

	commands.insert(CMD_DP_FITADIABATIC, CommandSpecifier(CMD_DP_FITADIABATIC), "dp_fitadiabatic");
	commands[CMD_DP_FITADIABATIC].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_fitadiabatic</b> <i>(abs_err Rsq T_ratio (stencil))</i>";
	commands[CMD_DP_FITADIABATIC].limits = { {double(0.0), Any()}, {double(0.0), Any()}, {double(0.0), Any()}, {int(3), Any()} };
	commands[CMD_DP_FITADIABATIC].descr = "[tc0,0.5,0.5,1/tc]Fit the computed self-consistent spin torque (see below) using Zhang-Li STT with fitting parameters P and beta (non-adiabaticity) using a given square in-plane stencil (default size 3) in order to extract the spatial variation of P. Cut-off values for absolute fitting error (default 0.1), Rsq measure (default 0.9), and normalized torque magnitude (default 0.1) can be set - value of zero disables cutoff. The focused mesh must be ferromagnetic, have the transport module set with spin solver enabled, and we also require Jc and either Ts or Tsi to have been computed. The fitting is done on Ts, Tsi, or on their sum depending if they’ve been enabled or not. Output available in Cust_S.";

	commands.insert(CMD_DP_CALCEXCHANGE, CommandSpecifier(CMD_DP_CALCEXCHANGE), "dp_calcexchange");
	commands[CMD_DP_CALCEXCHANGE].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_calcexchange</b>";
	commands[CMD_DP_CALCEXCHANGE].descr = "[tc0,0.5,0.5,1/tc]Calculate spatial dependence of exchange energy density for the focused mesh (must be magnetic and have an exchange module enabled). Output available in Cust_S.";

	commands.insert(CMD_DP_CALCTOPOCHARGEDENSITY, CommandSpecifier(CMD_DP_CALCTOPOCHARGEDENSITY), "dp_calctopochargedensity");
	commands[CMD_DP_CALCTOPOCHARGEDENSITY].usage = "[tc0,0.5,0,1/tc]USAGE : <b>dp_calctopochargedensity</b>";
	commands[CMD_DP_CALCTOPOCHARGEDENSITY].descr = "[tc0,0.5,0.5,1/tc]Calculate topological charge density spatial dependence for the focused mesh (must be magnetic). Output available in Cust_S.";

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
	dataDescriptor.push_back("<Mxsq>", DatumSpecifier("<Mxsq> : ", 1, "A/m", false, false), DATA_AVMXSQ);
	dataDescriptor.push_back("<Mysq>", DatumSpecifier("<Mysq> : ", 1, "A/m", false, false), DATA_AVMYSQ);
	dataDescriptor.push_back("<Mzsq>", DatumSpecifier("<Mzsq> : ", 1, "A/m", false, false), DATA_AVMZSQ);
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
	dataDescriptor.push_back("<V>", DatumSpecifier("<V> : ", 1, "V", false, false), DATA_V);
	dataDescriptor.push_back("<S>", DatumSpecifier("<S> : ", 3, "A/m", false, false), DATA_S);
	dataDescriptor.push_back("<elC>", DatumSpecifier("<elC> : ", 1, "S/m", false, false), DATA_ELC);
	dataDescriptor.push_back("<T>", DatumSpecifier("<Temperature> : ", 1, "K", false, false), DATA_TEMP);
	dataDescriptor.push_back("<T_l>", DatumSpecifier("<Temperature_l> : ", 1, "K", false, false), DATA_TEMP_L);
	dataDescriptor.push_back("V", DatumSpecifier("Potential : ", 1, "V"), DATA_POTENTIAL);
	dataDescriptor.push_back("I", DatumSpecifier("I, Net I : ", 2, "A"), DATA_CURRENT);
	dataDescriptor.push_back("R", DatumSpecifier("Resistance : ", 1, "Ohm"), DATA_RESISTANCE);
	dataDescriptor.push_back("e_demag", DatumSpecifier("Demag e : ", 1, "J/m3", false), DATA_E_DEMAG);
	dataDescriptor.push_back("e_exch", DatumSpecifier("Exchange e : ", 1, "J/m3", false, false), DATA_E_EXCH);
	dataDescriptor.push_back("e_exch_max", DatumSpecifier("Exchange Maximum e : ", 1, "J/m3", false, false), DATA_E_EXCH_MAX);
	dataDescriptor.push_back("e_surfexch", DatumSpecifier("Surface exchange e : ", 1, "J/m3", false), DATA_E_SURFEXCH);
	dataDescriptor.push_back("e_zee", DatumSpecifier("Zeeman e : ", 1, "J/m3", false, false), DATA_E_ZEE);
	dataDescriptor.push_back("e_mel", DatumSpecifier("Magnetoelastic e : ", 1, "J/m3", false), DATA_E_MELASTIC);
	dataDescriptor.push_back("e_anis", DatumSpecifier("Anisotropy e : ", 1, "J/m3", false, false), DATA_E_ANIS);
	dataDescriptor.push_back("e_rough", DatumSpecifier("Roughness e : ", 1, "J/m3", false), DATA_E_ROUGH);
	dataDescriptor.push_back("e_total", DatumSpecifier("Total e : ", 1, "J/m3", true), DATA_E_TOTAL);
	dataDescriptor.push_back("dwshift", DatumSpecifier("DW shift : ", 1, "m"), DATA_DWSHIFT);
	dataDescriptor.push_back("skyshift", DatumSpecifier("Skyrmion shift : ", 2, "m", false, false), DATA_SKYSHIFT);
	dataDescriptor.push_back("skypos", DatumSpecifier("Skyrmion (pos, dia) : ", 4, "m", false, false), DATA_SKYPOS);
	dataDescriptor.push_back("Q_topo", DatumSpecifier("Topological Charge : ", 1, "", false, false), DATA_Q_TOPO);
	dataDescriptor.push_back("v_iter", DatumSpecifier("V Solver Iterations : ", 1), DATA_TRANSPORT_ITERSTOCONV);
	dataDescriptor.push_back("s_iter", DatumSpecifier("S Solver Iterations : ", 1), DATA_TRANSPORT_SITERSTOCONV);
	dataDescriptor.push_back("ts_err", DatumSpecifier("Transport Solver Error : ", 1), DATA_TRANSPORT_CONVERROR);
	
	//---------------------------------------------------------------- MODULES

	//Modules
	moduleHandles.push_back("demag_N", MOD_DEMAG_N);
	moduleHandles.push_back("demag", MOD_DEMAG);
	moduleHandles.push_back("exchange", MOD_EXCHANGE);
	moduleHandles.push_back("DMexchange", MOD_DMEXCHANGE);
	moduleHandles.push_back("iDMexchange", MOD_IDMEXCHANGE);
	moduleHandles.push_back("surfexchange", MOD_SURFEXCHANGE);
	moduleHandles.push_back("Zeeman", MOD_ZEEMAN);
	moduleHandles.push_back("moptical", MOD_MOPTICAL);
	moduleHandles.push_back("aniuni", MOD_ANIUNI);
	moduleHandles.push_back("anicubi", MOD_ANICUBI);
	moduleHandles.push_back("melastic", MOD_MELASTIC);
	moduleHandles.push_back("transport", MOD_TRANSPORT);
	moduleHandles.push_back("heat", MOD_HEAT);
	moduleHandles.push_back("SOTfield", MOD_SOTFIELD);
	moduleHandles.push_back("STfield", MOD_STFIELD);
	moduleHandles.push_back("roughness", MOD_ROUGHNESS);
	moduleHandles.push_back("dipoledipole", MOD_ATOM_DIPOLEDIPOLE);
	
	//super-mesh modules
	moduleHandles.push_back("sdemag", MODS_SDEMAG);
	moduleHandles.push_back("strayfield", MODS_STRAYFIELD);
	moduleHandles.push_back("stransport", MODS_STRANSPORT);
	moduleHandles.push_back("sheat", MODS_SHEAT);
	moduleHandles.push_back("Oersted", MODS_OERSTED);

	//---------------------------------------------------------------- ODEs

	//ODEs : all (micromagnetic meshes)
	odeHandles.push_back("LLG", ODE_LLG);
	odeHandles.push_back("LLGStatic", ODE_LLGSTATIC);
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
	odeEvalHandles.push_back("RKF45", EVAL_RKF);
	odeEvalHandles.push_back("RKCK45", EVAL_RKCK);
	odeEvalHandles.push_back("RKDP54", EVAL_RKDP);
	odeEvalHandles.push_back("SDesc", EVAL_SD);

	//Allowed evaluation methods for given ODE
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP), ODE_LLG);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP, EVAL_SD), ODE_LLGSTATIC);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP), ODE_LLGSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP), ODE_LLB);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP), ODE_LLBSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP), ODE_LLGSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN, EVAL_RK4, EVAL_ABM, EVAL_RK23, EVAL_RKF, EVAL_RKCK, EVAL_RKDP), ODE_LLBSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN), ODE_SLLG);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN), ODE_SLLGSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN), ODE_SLLB);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN), ODE_SLLBSTT);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN), ODE_SLLGSA);
	odeAllowedEvals.push_back(make_vector(EVAL_EULER, EVAL_TEULER, EVAL_AHEUN), ODE_SLLBSA);

	odeDefaultEval.push_back(EVAL_RKF, ODE_LLG);
	odeDefaultEval.push_back(EVAL_SD, ODE_LLGSTATIC);
	odeDefaultEval.push_back(EVAL_RKF, ODE_LLGSTT);
	odeDefaultEval.push_back(EVAL_RKF, ODE_LLB);
	odeDefaultEval.push_back(EVAL_RKF, ODE_LLBSTT);
	odeDefaultEval.push_back(EVAL_RKF, ODE_LLGSA);
	odeDefaultEval.push_back(EVAL_RKF, ODE_LLBSA);
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

	dataSaveDescriptors.push_back("none", DataSaveDescriptor(DSAVE_NONE), DSAVE_NONE);
	dataSaveDescriptors.push_back("stage", DataSaveDescriptor(DSAVE_STAGE), DSAVE_STAGE);
	dataSaveDescriptors.push_back("step", DataSaveDescriptor(DSAVE_STEP), DSAVE_STEP);
	dataSaveDescriptors.push_back("iter", DataSaveDescriptor(DSAVE_ITER), DSAVE_ITER);
	dataSaveDescriptors.push_back("time", DataSaveDescriptor(DSAVE_TIME, "s"), DSAVE_TIME);

	//---------------------------------------------------------------- DEFAULT STARTING STATE

	//Default simulation stages - this must be added first so simStages is not an empty vector
	AddGenericStage(SS_RELAX);

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

	//start network sockets thread to listen for incoming messages
	commSocket.Change_Password(server_pwd);
	commSocket.Change_Port(server_port);
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

	//---------------------------------------------------------------- UPDATES

	//Check with "www.boris-spintronics.uk" if program version is up to date, and get latest update time for materials database
	if (start_check_updates) single_call_launch(&Simulation::CheckUpdate, THREAD_HANDLEMESSAGE2);

	//---------------------------------------------------------------- FINAL SETTINGS

	//Set number of OpenMP threads (default is use all available but user can configure this)
	if (!OmpThreads) OmpThreads = omp_get_num_procs();

	//Update display - do not animate starting view
	UpdateScreen_AutoSet_Sudden();

	//---------------------------------------------------------------- MESSAGES

	//if this directory doesn't exist create it
	if (!MakeDirectory(directory)) BD.DisplayConsoleError("ERROR : Couldn't create user directory.");

	BD.DisplayConsoleMessage("Console activated...");
	BD.DisplayFormattedConsoleMessage("[tc0,0.5,0,1/tc]To open manual use the <b>manual</b> command.");

	//---------------------------------------------------------------- SAVE DEFAULT and CTOR FINISH

	Save_Startup_Flags();

	//save the default state in program directory for loading with "default" command (update this automatically in case the program version has changed)
	SaveSimulation(directory + "default");
}

Simulation::~Simulation()
{
	BD.DisplayConsoleMessage("Shutting down...");

	StopSimulation();

	Stop_All_Threads();

	BD.DisplayConsoleMessage("All threads stopped. Clean-up...");
}
