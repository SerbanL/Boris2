#include "stdafx.h"
#include "Simulation.h"

//MAIN SIMULATION LOOP. Runs in SimulationThread
void Simulation::Simulate(void)
{
	//if simulation is being stopped (done from THREAD_LOOP_STOP), then stop doing any more work/saving etc.
	if (is_thread_running(THREAD_LOOP_STOP)) return;

	//stop other parts of the program from changing simulation parameters in the middle of an interation
	//non-blocking std::mutex is needed here so we can stop the simulation from HandleCommand - it also uses the simulationMutex. If Simulation thread gets blocked by this std::mutex they'll wait on each other forever.
	if (simulationMutex.try_lock()) {

		CheckSaveDataConditions();

		if (simStages[stage_step.major].stage_type() == SS_MONTECARLO) {

			//Monte-Carlo stages are special - use Iterate_MonteCarlo to advance simulation instead
#if COMPILECUDA == 1
			if (cudaEnabled) {

				SMesh.Iterate_MonteCarloCUDA(simStages[stage_step.major].get_value<double>(stage_step.minor));
				if (SMesh.Get_MonteCarlo_ComputeFields()) SMesh.ComputeFieldsCUDA();
			}
			else {

				SMesh.Iterate_MonteCarlo(simStages[stage_step.major].get_value<double>(stage_step.minor));
				if (SMesh.Get_MonteCarlo_ComputeFields()) SMesh.ComputeFields();
			}
#else
			SMesh.Iterate_MonteCarlo(simStages[stage_step.major].get_value<double>(stage_step.minor));
			if (SMesh.Get_MonteCarlo_ComputeFields()) SMesh.ComputeFields();
#endif
		}
		else {

			//advance time for this iteration
#if COMPILECUDA == 1
			if (cudaEnabled) SMesh.AdvanceTimeCUDA();
			else SMesh.AdvanceTime();
#else
			SMesh.AdvanceTime();
#endif
		}

		//Display update (asynchronous only if cuda is enabled)
		if (iterUpdate && SMesh.GetIteration() % iterUpdate == 0) UpdateScreen_Quick(cudaEnabled);

		//Check conditions for advancing simulation schedule
		CheckSimulationSchedule();

		//finished this iteration
		simulationMutex.unlock();

		//THREAD_HANDLEMESSAGE is used to run HandleCommand, which also uses simulationMutex to guard access.
		//With Visual Studio 2017 v141 toolset : without the short wait below, when HandleCommand has been called, simulationMutex will block access for a long time as this Simulate method gets called over and over again on its thread.
		//This means the command gets executed very late (ten seconds not unusual) - not good!
		//This wasn't a problem with Visual Studio 2012, v110 or v120 toolset. Maybe with the VS2017 compiler the calls to Simulate on the infinite loop thread are all inlined. 
		//Effectively there is almost no delay between unlocking and locking the std::mutex again on the next iteration - THREAD_HANDLEMESSAGE cannot sneak in to lock simulationMutex easily!
		if (is_thread_running(THREAD_HANDLEMESSAGE)) std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}
}

//Similar to Simulate but only runs for one iteration and does not advance time
void Simulation::ComputeFields(void)
{
	if (is_thread_running(THREAD_LOOP)) {

		StopSimulation();
	}
	else {

		BError error;
		BD.DisplayConsoleMessage("Initializing modules...");

		bool initialization_error;

		if (!cudaEnabled) {

			initialization_error = err_hndl.qcall(error, &SuperMesh::InitializeAllModules, &SMesh);
		}
		else {

#if COMPILECUDA == 1
			initialization_error = err_hndl.qcall(error, &SuperMesh::InitializeAllModulesCUDA, &SMesh);
#endif
		}

		if (initialization_error) {

			BD.DisplayConsoleError("Failed to initialize simulation.");
			err_hndl.show_error(error, true);
			return;
		}
	}

	BD.DisplayConsoleMessage("Initialized. Updating fields.");

#if COMPILECUDA == 1
	if (cudaEnabled) SMesh.ComputeFieldsCUDA();
	else SMesh.ComputeFields();
#else
	SMesh.ComputeFields();
#endif

	//Display update
	UpdateScreen();

	BD.DisplayConsoleMessage("Fields updated.");
}

void Simulation::RunSimulation(void)
{
	//make sure the disk buffer has the correct size
	if (savedata_diskbuffer.size() != savedata_diskbuffer_size) savedata_diskbuffer.resize(savedata_diskbuffer_size);
	if (savedata_diskoverflowbuffer.size() != savedata_diskbuffer_size) savedata_diskoverflowbuffer.resize(savedata_diskbuffer_size);
	savedata_diskbuffer_position = 0;
	savedata_diskoverflowbuffer_position = 0;

	if (is_thread_running(THREAD_LOOP)) {

		BD.DisplayConsoleMessage("Simulation already running.");
		return;
	}

	BD.DisplayConsoleMessage("Initializing modules...");

	BError error;
	bool initialization_error;

	if (!cudaEnabled) {

		initialization_error = err_hndl.qcall(error, &SuperMesh::InitializeAllModules, &SMesh);
	}
	else {
#if COMPILECUDA == 1
		initialization_error = err_hndl.qcall(error, &SuperMesh::InitializeAllModulesCUDA, &SMesh);
#endif
	}

	if (initialization_error) {

		BD.DisplayConsoleError("Failed to initialize simulation.");
		err_hndl.show_error(error, true);
		return;
	}

	//show any warnings, and continue
	if (error.warning_set()) err_hndl.show_error(error, true);

	//set initial stage values if at the beginning (stage = 0, step = 0, and stageiteration = 0)
	if (Check_and_GetStageStep() == INT2()) {

		if (SMesh.GetStageIteration() == 0) {

			SetSimulationStageValue();
			appendToDataFile = false;
		}
	}

	infinite_loop_launch(&Simulation::Simulate, &Simulation::SetupRunSimulation, THREAD_LOOP);
	BD.DisplayConsoleMessage("Initialized. Simulation running. Started at: " + Get_Date_Time());

	sim_start_ms = GetSystemTickCount();
}

//SetupRunSimulation sets cuda device and number of OpenMP threads for the RunSimulation, called on the same thread as RunSimulation
void Simulation::SetupRunSimulation(void)
{
	//Set number of OpenMP threads to use on this thread
	if (OmpThreads && OmpThreads <= omp_get_num_procs()) omp_set_num_threads(OmpThreads);

#if COMPILECUDA == 1
	//Commands are executed on newly spawned threads, so if cuda is on and we are not using device 0 (default device) we must switch to required device, otherwise 0 will be used
	if (cudaEnabled && cudaDeviceSelect != 0) cudaSetDevice(cudaDeviceSelect);
#endif
}

void Simulation::StopSimulation(void)
{
	if (is_thread_running(THREAD_LOOP)) {

		stop_thread(THREAD_LOOP);

		sim_end_ms = GetSystemTickCount();

		//flush disk buffer
		if (savedata_diskoverflowbuffer_position) SaveData_DiskBufferFlush(&savedata_diskoverflowbuffer, &savedata_diskoverflowbuffer_position);
		if (savedata_diskbuffer_position) SaveData_DiskBufferFlush(&savedata_diskbuffer, &savedata_diskbuffer_position);

		BD.DisplayConsoleMessage("Simulation stopped. " + Get_Date_Time());

		if (commSocket.ClientConnected()) {

			//if client connected, signal simulation has finished
			commSocket.SetSendData({ "stopped" });
			commSocket.SendDataParams();
		}

		UpdateScreen();
	}
}

void Simulation::ResetSimulation(void)
{
	StopSimulation();

	stage_step = INT2();
	SMesh.ResetODE();
	
	UpdateScreen();
}

//Execute the commands stored in the command buffer : when you call this must ensure simulationMutex is unlocked, and THREAD_HANDLEMESSAGE is free (or else launch it asynchonously)
void Simulation::RunCommandBuffer(void)
{
	for (int idx = 0; idx < command_buffer.size(); idx++) {

		HandleCommand(command_buffer[idx]);
	}
}