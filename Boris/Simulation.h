//Defines the top level object

//IMPORTANT NOTES

//1.
//Stop DiagTrack and DPS (Diagnostic Policy Service) from Windows services. These sometime insert a hook into the executable which slows down the program very significantly. Also causes unexplained crashes (!!!).

//2.
//If using Visual Studio 2017 with CUDA 9.2:
//CUDA 9.2 only officially suports VS2017 up to _MSC_VER = 1913 (Visual Studio 2017 version 15.6)
//If you work with higher versions of VS2017 (e.g. 15.9 where _MSV_VER = 1916) then you won't be able to compile.
//To keep using 9.2 you can manually modify the check in host_config.h file located at:
//C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.2\include\crt
//replace #if _MSC_VER < 1600 || _MSC_VER > 1913
//with
//#if _MSC_VER < 1600 || _MSC_VER > 1916

//3.
//Cannot import BorisLib.h in .cu files. This is due to nvcc not supporting C++17 standards at the time of writing (only C++14).
//I got around this by keeping all compilation units in .cu files with C++14 code only, which didn't interfere with program design as I didn't really need any C++17 features there (but very useful in many other compilation units!).

//4.
//On Windows 10 it's possible CUDA nvcc will keep coming up with "exit error code 1", including when you simply try to clean the project.
//This is a generic error, and if you enable maximum verbosity in VS (Tools > Options > Projects & Solutions > Build and Run) you'll see this is because of "Access is denied." for a process which writes some temporary data.
//In particular when the process is trying to write temporary data to "C:\Users\<UserName>\AppData\Local" it won't be able to because of denied access.

//5.
//Don't set the number of CUDA threads per block too large. Currently using 128.
//If the amount of code that is included in a CUDA kernel is too large (through inlining of various functions, etc. etc.) the program will start to exhibit very strange bugs.

#pragma once

#include "BorisLib.h"

#include "Boris_Enums_Defs.h"
#include "ErrorHandler.h"
#include "OVF2_Handlers.h"

#include "BorisDisplay.h"
#include "BorisDisplayNonGraphical.h"
#include "BorisInteractiveObjects.h"

#include "Commands.h"

#include "SimSharedData.h"
#include "SimulationData.h"
#include "SimSchedule.h"
#include "DataProcessing.h"
#include "MaterialsDataBase.h"

#include "Mesh.h"
#include "Atom_Mesh.h"
#include "SuperMesh.h"


#if COMPILECUDA == 1
#include "BorisCUDALib.h"
#endif

//Top level object
class Simulation : 
	public SimulationSharedData, 
	public Threads<Simulation>,
	public ProgramState<Simulation, 
	std::tuple<
	BorisDisplay, 
	std::string, std::string, std::string, std::string, bool, bool, bool, int, int,
	vector_lut<DatumConfig>, vector_lut<DatumConfig>, 
	INT2, 
	vector_lut<StageConfig>, int, bool, 
	SuperMesh, 
	bool, int, 
	bool, 
	bool, bool,
	DBL4, DBL2, DBL2, int,
	DBL3, INT3, DBL3, std::string,
	vector_key<double>,
	std::vector<std::string>>,
	std::tuple<> >
{
private:

	////////////////////////////////////////////////////////////////////////////////////
	//IDs

#if GRAPHICS == 1 && OPERATING_SYSTEM == OS_WIN
	//Save window handle if in graphics mode on windows
	HWND hWnd;
#endif

	int Program_Version;

	////////////////////////////////////////////////////////////////////////////////////
	//VARIOUS STANDARD FILES AND DIRECTORIES

	//saved simulation file header and program version number
	std::string simfile_header = "Boris simulation file version ";
	
	std::string domain_name = "www.boris-spintronics.uk";
	std::string download_page = "/download";
	
	//check if current program version is the latest
	std::string version_checker = "version.php";
	
	//check latest avilable incremental update for current program version
	std::string delta_update_checker = "delta_update_check.php";
	std::string installers_dir = "/installers";

	//patcher program utility
	std::string patch_utility = "PatchApply.exe";

	//allow users to send a new entry for the online materials database
	std::string mdb_entry_handler = "mdbentry.php";
	//update the local materials database from the online materials database
	std::string mdb_update_handler = "mdb.php";
	//check when the online materials database was last updated
	std::string mdb_lastupdate = "mdbtime.php";

	//Embedded Python utility files (located in executable directory)
	std::string embeddedPythonWraper_fileName = "EmbeddedPythonWrapper.py";

	//directory name with Boris Data (e.g. error log output, default.bsm, and other important folders
#if OPERATING_SYSTEM == OS_WIN
	std::string boris_data_directory = "Boris Data/";
#elif OPERATING_SYSTEM == OS_LIN
	std::string boris_data_directory = "Boris_Data/";
#endif
	//Default simulations directory folder, located within boris_data_directory
	std::string boris_simulations_directory = "Simulations/";

	//BorisPythonScripts directory, containing pre-configured Python scripts and NetSocks, located within boris_data_directory
	std::string BorisPythonScripts_directory = "BorisPythonScripts/";

	//program name (from argv[0])
	std::string progName;

	////////////////////////////////////////////////////////////////////////////////////
	//STARTUP OPTIONS

	//check for updates on program startup?
	bool start_check_updates = false;

	//start script server on program startup?
	bool start_scriptserver = true;

	//flag used to interrupt embedded Python script : set to false
	bool python_script_running = false;

	//save/load startup flags in this file (e.g. log_errors, start_check_updates, start_scriptserver)
	std::string startup_options_file;

	//value set by version update checker :
	//-1: attempting to connect
	//0: connection failure
	//1: program up to date
	//2: update available
	int version_checking = -1;

	//number of threads to use for OpenMP (0 means determine maximum available and set that)
	int OmpThreads = 0;

	//log errors generated by command input?
	bool log_errors = true;
	std::string errorlog_fileName;

	////////////////////////////////////////////////////////////////////////////////////
	//SERVER OPTIONS

	//Default netsocks server port
	std::string server_port = DEFAULT_PORT;
	//default server receiver sleep time : lower value improves server responsivity but increase CPU load for server thread.
	int server_recv_sleep_ms = RECVSLEEPMS;
	//default server password : no password by default (mostly comms will happen locally behind a firewall, or in any case a firewall would be on by default, but user can set a password anyway)
	std::string server_pwd = "";

	////////////////////////////////////////////////////////////////////////////////////
	//EMBEDDED PYTHON SCRIPT OPTIONS

	//script to execute when passed through command line
	std::string python_script = "";

	//script multi-GPU option for #pragma parallel for directives
	int python_script_mGPU = 1;
	
	//workload splitting for #pragma parallel for directive (first one is the stride, after we have the offsets)
	std::vector<int> python_script_parallel = {};
	
	//terminate program or not after execution of command line script
	int python_script_terminate = 0;

	//delete source file on completion
	int python_script_deletesource = 0;

	////////////////////////////////////////////////////////////////////////////////////
	//VARIOUS AUXILIARY OBJECTS

	ErrorHandler<Simulation> err_hndl;
	//the display master object - all display updates done through this object
	BorisDisplay BD;

	//communication with scripting clients - listen for incoming commands and respond with requested data
	NetSocks commSocket;

	//Materials data base management - local and online
	MaterialsDB mdb;

	////////////////////////////////////////////////////////////////////////////////////
	//SUPERMESH

	//Simulation super-mesh
	SuperMesh SMesh;

	////////////////////////////////////////////////////////////////////////////////////
	//DESCRIPTORS

	//collection of commands, with the CommandSpecifier indexed by the command name (name as typed in console). The CommandSpecifier describes the properties common to all commands (vector_key more convenient than a map)
	vector_key<CommandSpecifier> commands;

	//lut-indexed output data descriptor : can be indexed using a value from DATA_ enum; can also be indexed using a key: this is the data handle as used in the console
	vector_key_lut<DatumSpecifier> dataDescriptor;

	//lut-indexed vector of module handles (handle as entered in the console) : can be indexed using a value from MOD_ enum
	vector_lut<std::string> moduleHandles;

	//lut-indexed vector of mesh type handles can be indexed using a value from MESH_ enum
	vector_lut<std::string> meshtypeHandles;

	//For each ODE_ entry specify allowed EVAL_ entries (i.e. for each differential equation specify allowed evaluation methods)
	vector_lut< std::vector<EVAL_> > odeAllowedEvals;

	//The default evaluation method for each ODE : this is the evaluation method set when changing ODEs
	vector_lut<EVAL_> odeDefaultEval;

	//Link ODE_ entries with text handles (ODE_ is the major id) : micromagnetic meshes
	vector_lut<std::string> odeHandles;

	//Link ODE_ entries with text handles (ODE_ is the major id) : atomistic meshes (smaller subset of ODEs allowed, but same evaluation methods)
	vector_lut<std::string> atom_odeHandles;

	//Link EVAL_ entries with text handles (EVAL_ is the major id)
	vector_lut<std::string> odeEvalHandles;

	//handles and descriptors for simulation stages types, stop conditions and data saving conditions, indexed by SS_ and STOP_ respectively, as well as keys (handles)
	vector_key_lut<StageDescriptor> stageDescriptors;
	vector_key_lut<StageStopDescriptor> stageStopDescriptors;
	vector_key_lut<DataSaveDescriptor> dataSaveDescriptors;

	//interactive console object info : index with majorId (IOI_ entry) and minorId (subtype). If only one entry for majorId then it applies irrespective of minorId.
	vector_lut<std::string> ioInfo;

	////////////////////////////////////////////////////////////////////////////////////
	//OUTPUT DATA SETTINGS

	//files for saving output data
	std::string savedataFile = "out_data.txt";
	//image save file base (during a simulation this is complemented by _iteration.png termination)
	std::string imageSaveFileBase = "mesh_image";

	//currently loaded simulation file, including directory (empty if none)
	std::string currentSimulationFile;
	//append or create new data file?
	bool appendToDataFile = false;
	//saveDataList is in SimulationSharedData

	//flags for enabling data and image saving (they both use the same saving condition in the simulation schedule)
	bool saveDataFlag = true, saveImageFlag = false;

	//precision (number of significant figures) of data in saved files and display (mesh viewer, databox, console: affects ToString precision)
	int dataprecision = SAVEDATAPRECISION;

	//output buffered data saving to disk
	int savedata_diskbuffer_size = DISKBUFFERLINES;
	//savedata_diskbuffer_position maximum value is savedata_diskbuffer_size, and indicates the next position to write data to buffer; when maximum reached then flush buffer
	int savedata_diskbuffer_position = 0;
	std::vector<std::string> savedata_diskbuffer;

	//data to display in data box
	vector_lut<DatumConfig> dataBoxList;

	////////////////////////////////////////////////////////////////////////////////////
	//SHAPE OPTIONS

	//modifiers for shape generator commands (shape_...)
	DBL3 shape_rotation = DBL3(0, 0, 0);
	INT3 shape_repetitions = INT3(1, 1, 1);
	DBL3 shape_displacement = DBL3(0, 0, 0);
	//mesh shape generation method
	std::string shape_method = MeshShape(MSHAPEMETHOD_ADD).method_name();

	////////////////////////////////////////////////////////////////////////////////////
	//OTHER OPTIONS
	
	//display updating on number of iterations during simulation
	int iterUpdate = 100;

	//autocomplete commands in console
	bool autocomplete = true;

	//default no cropping of saved mesh images : left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
	//use a small back-off (0.001) for default settings to avoid capturing lines at the edges
	DBL4 image_cropping = DBL4(0.001, 0.001, 0.999, 0.999);

	////////////////////////////////////////////////////////////////////////////////////
	//COMMANDS AUXILIARY
	
	//vector of command strings (commands with parameters, as would be typed into the console or received by command server)
	//this is used to execute a sequennce of commands, e.g. for more advanced data extraction during a data saving schedule
	std::vector<std::string> command_buffer;

	////////////////////////////////////////////////////////////////////////////////////
	//DP ARRAYS

	//internal data processing arrays : used to load values into columns and process data using console commands
	DPArrays dpArr = DPArrays(MAX_ARRAYS);
	
	////////////////////////////////////////////////////////////////////////////////////
	//CUDA

	//are cuda computations available for this program version and hardware?
#if COMPILECUDA == 0
	bool cudaAvailable = false;
#else
	//still need to test for this in the constructor
	bool cudaAvailable;
#endif

	////////////////////////////////////////////////////////////////////////////////////
	//MUTEXES

	//thread-safe access to simulation state changes - e.g. don't want to change parameters in the middle of an iteration
	std::mutex simulationMutex;

	//thread-safe access for external input: must make sure there is no clash between direct user input and scripted input, otherwise crashes can result
	std::mutex userInputMutex;

	std::mutex diskMutex;

	////////////////////////////////////////////////////////////////////////////////////
	//BENCHMARKING

	//timing values used for benchmarking - use benchtime command to return sim_end_ms - sim_start_ms
	unsigned int sim_start_ms = 0, sim_end_ms = 0;

public:

private:

	//-------------------------------------Startup Options

	//save/load flags: start_check_updates, start_scriptserver, log_errors
	void Load_Startup_Flags(void);
	void Save_Startup_Flags(void);

	//-------------------------------------Program update checker (in BorisUpdate.cpp)

	//Check with "www.boris-spintronics.uk" if program version is up to date
	void CheckUpdate(void);
	void OpenDownloadPage(void);

	//check latest version available
	int CheckLatestVersion(void);
	//for current program version check available incremental update versions
	std::vector<int> CheckDeltaUpdate(void);
	//check and download incremental update files; if target_version specified download only incremental update files up to target_version
	std::vector<std::string> DownloadDeltaUpdate(int target_version = 0);
	//check, download, and install latest update file
	int InstallDeltaUpdate(int target_version = 0);
	//Roll back to a backup version using the PatchApply.exe utility
	void RollbackUpdate(int target_version);

	//-------------------------------------Communication with interactive console objects (TextObject with action handler set)

	//These are not needed in non-graphical mode
#if GRAPHICS == 1

	//Interactive objects in the console can generate messages, which must be handled here. This handler is passed using a functionoid (made in the constructor) to the BorisConsole object	 
	//Iteractive objects are constructed using the passed functionoid. On user interaction this handler is called with their properties and interaction actionCode (thus specifying the action requested).
	InteractiveObjectActionOutcome ConsoleActionHandler(int actionCode, InteractiveObjectProperties& iop, TextObject *pTO);

	//called by TextObject when it draws itself: check if the textobject needs to change. return true if any changes made so the caller knows to recalculate things (e.g. text placements, etc.)
	InteractiveObjectStateChange ConsoleInteractiveObjectState(InteractiveObjectProperties &iop, TextObject *pTO);

#endif

	//-------------------------------------Files loading and saving

	//Save / Load simulation from .bsm file : use ProgramState methods.
	BError SaveSimulation(std::string fileName);
	BError LoadSimulation(std::string fileName);

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) 
	{ 
		Conversion::tostringconversion::precision() = dataprecision;
		dpArr.clear_all(); 
	}

	//-------------------------------------Simulation control

	//Called by RunSimulation to prepare everything, just before Simulate() is launched
	void PrepareRunSimulation(void);

	//Start and stop Simulate method in a separate thread
	void RunSimulation(void);
	//SetupRunSimulation sets cuda device and number of OpenMP threads for the RunSimulation, called on the same thread as RunSimulation
	void SetupRunSimulation(void);
	void StopSimulation(void);
	//stop and reset simulation back to starting point
	void ResetSimulation(void);

	//main simulation method - computes a single complete iteration (a full ode time-step)
	void Simulate(void);
	//a dummy function which does no work, to keep THREAD_LOOP busy when needed
	void Simulate_Dummy(void);

	//Similar to Simulate but only runs for one iteration and does not advance time
	void ComputeFields(void);

	//Command handler
	void HandleCommand(std::string command);
	
	//Execute the commands stored in the command buffer : when you call this must ensure simulationMutex is unlocked, and THREAD_HANDLEMESSAGE is free (or else launch it asynchonously)
	void RunCommandBuffer(void);

	//-------------------------------------Simulation schedule

	//add a new simulation stage with default settings
	void AddGenericStage(SS_ stageType, std::string meshName = "");
	void DeleteStage(int stageIndex);

	//set stop condition for index-th entry in simStages (the index-th stage) with default stop value.
	void SetGenericStopCondition(int index, STOP_ stopType);
	//set data save conditon for index-th entry in simStages (the index-th stage) with default saving value.
	void SetGenericDataSaveCondition(int index, DSAVE_ dsaveType);

	//edit set stage type, stop condition, or data save condition
	void EditStageType(int index, SS_ stageType, std::string meshName = "");
	void EditStageValue(int stageIndex, std::string value_string);
	void EditStageStopCondition(int index, STOP_ stopType, std::string stopValueString = "");
	void EditDataSaveCondition(int index, DSAVE_ dsaveType, std::string dsaveValueString = "");

	//update mesh names in simStages
	void UpdateStageMeshNames(std::string oldMeshName, std::string newMeshName);

	//get stage step, step as an INT2, also checking if they are valid
	INT2 Check_and_GetStageStep();

	//Check if conditions for ending stage, step have been met and set values for next stage if so.
	void CheckSimulationSchedule(void);

	//when stage step stopping conditions hve been met call this. Next stage step will be set. Simulation is stopped if schedule finishes.
	void AdvanceSimulationSchedule(void);

	//set the values for current stage_step value (if stage_index < 0), as configured in simStages
	//If stage_index >= 0 then set value for given stage
	void SetSimulationStageValue(int stage_index = -1);

	//check if conditions for saving data have been met for curent stage
	void CheckSaveDataConditions();

	//-------------------------------------Console messages helper methods

	//show usage for given console command
	void PrintCommandUsage(std::string command_name) { BD.DisplayFormattedConsoleMessage(commands[command_name].usage); BD.DisplayFormattedConsoleMessage(commands[command_name].descr); BD.DisplayFormattedConsoleMessage(commands[command_name].return_descr); }

	//file name can have data specifiers, e.g. %iter% means '%iter%' should be replaced by the actual value of the iter data field, etc.
	std::string ScanFileNameData(std::string& fileName);

	//---------------------------------------------------- MESH LIST

	//NOTE, in general we have two main types of methods:

	//1. void Print_something_List(void)			: this displays in the console required formatted text involving a list of things
	//2. std::string Build_something_ListLine(int index) : this makes the formatted text for a single line in the list for the above method. Also needed in state handler to update the list when something changes.

	//Sometimes methods like this are also use : std::string Build_something_Text(int index) : this makes unformatted text used in interactive objects textId and displayed text.

	//show list of meshes
	void Print_Mesh_List(void);
	//build formatted std::string, describing the given mesh
	std::string Build_Mesh_ListLine(int meshIndex);

	//---------------------------------------------------- MESH DISPLAY

	void Print_MeshDisplay_List(void);
	std::string Build_MeshDisplay_ListLine(int meshIndex);

	//---------------------------------------------------- MODULES LIST

	//show list of modules (active and not active) for all meshes
	void Print_Modules_List(void);
	//build formatted std::string, describing the modules available for a given named mesh
	std::string Build_Modules_ListLine(int meshIndex);

	//---------------------------------------------------- MODULES EFFCTIVE FIELD DISPLAY LIST

	//show list of modules (active and not active) for all meshes
	void Print_DisplayModules_List(void);
	//build formatted std::string, describing the modules available for a given named mesh
	std::string Build_DisplayModules_ListLine(int meshIndex);

	//---------------------------------------------------- ODES LIST

	//show ODE configuration
	void Print_ODEs(void);

	//---------------------------------------------------- SHOW DATA LIST (for console display)

	//show data available for showdata command
	void Print_ShowData(void);

	//---------------------------------------------------- OUTPUT DATA LIST (for saving during a simulation)

	//Print available and currently set output data as interactive lists
	void Print_AvailableOutputData(void);

	//build formatted std::string for interactive object describing an entry in saveDataList at index_in_list (also helper method to build the actual unformatted display text in the object)
	void Print_SetOutputData_List(void);
	std::string Build_SetOutputData_ListLine(int index_in_list);
	std::string Build_SetOutputData_Text(int index_in_list);

	//---------------------------------------------------- SIMULATION STAGES

	//Print available and currently set simulation stages as interactive lists
	void Print_SetStages_List(void);

	//build formatted std::string for interactive object describing an entry in simStages at index_in_list (also helper method to build the actual unformatted display text in the object)
	std::string Build_SetStages_ListLine(int index_in_list);
	std::string Build_SetStages_Text(int index_in_list);

	//helpers for interactive objects for set stages lines
	std::string Build_SetStages_ValueText(int index_in_list);
	std::string Build_SetStages_StopConditionText(int index_in_list);
	std::string Build_SetStages_SaveConditionText(int index_in_list, int dsaveIdx);

	//---------------------------------------------------- MESH PARAMETERS

	//print all mesh parameters for the given mesh name
	void Print_MeshParams(std::string meshName);

	//build an entire line with mesh parameters for a given mesh index
	std::string Build_MeshParams_Line(int meshIndex);

	//build mesh parameter text for a given mesh index and param identifier - helper method
	std::string Build_MeshParams_Text(int meshIdx, PARAM_ paramId);

	//---------------------------------------------------- MESH PARAMETERS TEMPERATURE DEPENDENCE

	//print all mesh parameters temperature dependence for the given mesh name
	void Print_MeshParamsTemperature(std::string meshName);

	//build an entire line with mesh parameters temperature dependence for a given mesh index
	std::string Build_MeshParamsTemp_Text(int meshIndex);

	//---------------------------------------------------- MESH PARAMETERS SPATIAL DEPENDENCE

	//print all mesh parameters spatial dependence for the given mesh name
	void Print_MeshParamsVariation(std::string meshName);

	//build an entire line with mesh parameters spatial dependence for a given mesh index
	std::string Build_MeshParamsVariation_Text(int meshIndex);

	//---------------------------------------------------- MOVING MESH SETTINGS

	void PrintMovingMeshSettings(void);

	//---------------------------------------------------- ELECTRODES and TRANSPORT SETTINGS

	void Print_Electrodes_List(void);
	std::string Build_Electrodes_ListLine(int el_index);

	void PrintTransportSolverConfig(void);

	//---------------------------------------------------- TMR SETTINGS

	void Print_TMRType_List(void);
	std::string Print_TMRType_ListLine(int meshIndex);

	//---------------------------------------------------- TEMPERATURE

	void Print_MeshTemperature_List(void);
	std::string Build_MeshTemperature_ListLine(int meshIndex);

	void Print_HeatBoundaries_List(void);
	std::string Build_HeatBoundaries_ListLine(int meshIndex);

	//---------------------------------------------------- CURIE TEMPERATURE and ATOMIC MOMENT

	void Print_CurieandMoment_List(void);
	std::string Build_CurieandMoment_ListLine(int meshIndex);
	
	//---------------------------------------------------- TEMPERATURE MODEL TYPE

	void Print_TemperatureModel_List(void);
	std::string Build_TemperatureModel_ListLine(int meshIndex);

	//---------------------------------------------------- STOCHASTICITY SETIINGS

	void Print_Stochasticity_List(void);
	std::string Build_Stochasticity_ListLine(int meshIndex);

	//---------------------------------------------------- EVALUATION SPEEDUP SETIINGS

	void Print_Speedup_List(void);
	std::string Build_Speedup_ListLine(int meshIndex);

	//---------------------------------------------------- CUDA and MEMORY INFO

	void Print_CUDAStatus(void);
	void Print_MemoryInfo(void);

	//---------------------------------------------------- SCALE RECTS STATUS

	void Print_Scale_Rects_Status(void);

	//---------------------------------------------------- COUPLED-To-DIPOLES STATUS

	void Print_CoupledToDipoles_Settings(void);

	//---------------------------------------------------- DIPOLE SHIFTING ALGORITHM

	void Print_DipoleShiftingAlgorithm_List(void);
	std::string Build_DipoleShiftingAlgorithm_ListLine(int meshIndex);

	//---------------------------------------------------- ERROR LOG STATUS and other STARTUP OPTIONS

	void Print_ErrorLogStatus(void);

	void Print_StartupUpdateCheckStatus(void);
	void Print_StartupScriptServerStatus(void);

	//print all startup options
	void Print_StartupOptions(void);

	//---------------------------------------------------- NUMBER OF THREADS

	void Print_Threads(void);

	//---------------------------------------------------- SCRIPT SERVER INFO

	void Print_ServerInfo(void);

	//---------------------------------------------------- NEIGHBORING MESHES EXCHANGE COUPLING

	void Print_ExchangeCoupledMeshes_List(void);
	std::string Build_ExchangeCoupledMeshes_ListLine(int meshIndex);

	//---------------------------------------------------- MESH ROUGHNESS REFINEMENT

	void Print_MeshRoughnessRefinement(std::string meshName);

	//---------------------------------------------------- MULTILAYERED CONVOLUTION CONFIGURATION

	void Print_MultiConvolution_Config(void);

	//---------------------------------------------------- GPU KERNELS DEMAG INITIALIZATION CONFIGURATION

	void Print_GPUKernels_Config(void);

	//---------------------------------------------------- MATERIALS DATABASE

	void Print_MaterialsDatabase(void);

	//---------------------------------------------------- ADAPTIVE TIME STEP CONTROL

	void Print_AStepCtrl(void);

	//---------------------------------------------------- PERIODIC BOUNDARY CONDITIONS

	void Print_PBC(void);
	std::string Build_PBC_ListLine(int meshIndex);

	//---------------------------------------------------- INDIVIDUAL SHAPE CONTROL

	void Print_IndividualShapeStatus(void);

	//---------------------------------------------------- USER EQUATION CONSTANTS

	//Print currently set equation constants
	void Print_EquationConstants(void);

	//build formatted std::string for interactive objects describing user constants at index_in_list from userConstants (also helper method to build the actual unformatted display text in the object)
	std::string Build_EquationConstants_ListLine(int index_in_list);
	std::string Build_EquationConstants_Text(int index_in_list);

	//---------------------------------------------------- SKYPOS SETTINGS
	
	//Print skypos multiplier value
	void Print_skypos_dmul(void);
	std::string Build_skypos_dmul_ListLine(int meshIndex);

	//---------------------------------------------------- DWPOS SETTINGS

	//Print dwpos fitting component value
	void Print_DWPos_Component(void);

	//---------------------------------------------------- MONTE-CARLO SETTINGS

	void Print_MCSettings(void);
	std::string Build_MCSettings_ListLine(int meshIndex);

	//---------------------------------------------------- SHAPE MODIFIERS

	void Print_ShapeSettings(void);

	//---------------------------------------------------- SHAPE MODIFIERS
	
	void Print_DisplayRenderSettings(void);

	//---------------------------------------------------- MAKE INTERACTIVE OBJECT : Auxiliary method

	//Generate a formatted std::string depending on the interactive object identifier
	template <typename ... PType>
	std::string MakeIO(IOI_ identifier, PType ... params);

	//construct ioInfo object in constructor
	void MakeIOInfo(void);

	//set entry line in console to given text
	void SetConsoleEntryLineText(std::string text) { BD.SetConsoleEntryLineText(text); }

	//-------------------------------------Script client

	//Scripted communication: listen for commands and dispatch them to HandleCommand; respond with data after command handled
	void Listen_Incoming_Message(void);

	//disable / enable script server
	void Script_Server_Control(bool status);

#if PYTHON_EMBEDDING == 1
	//This method loads and runs the Python script - when called from HandleCommand must be launched on new thread to avoid mutex blocking
	void RunPythonScript(std::string fileName);
#endif

	//-------------------------------------Data helper methods

	//Make new entry in data box
	void NewDataBoxField(DatumConfig dConfig);

	//Update all currently set data box entries - do not call directly, call it through UpdateScreen(); or UpdateDataBox_Refresh(); Doesn't refresh screen.
	void UpdateDataBox(void);

	//Clear and remake all data box fields. Doesn't refresh screen.
	void RebuildDataBox(void);

	//delete any data box fields with given mesh name set
	void DeleteDataBoxFields(std::string meshName);

	//change names of mesh names in data box fields
	void ChangeDataBoxLabels(std::string oldMeshName, std::string newMeshName);

	//update rectangles in dataBoxList (rect_old and rect_new are the old and new rectangles for the given mesh. All dataBoxList rectangles are scaled accordingly.
	void UpdateDataBoxEntries(Rect rect_old, Rect rect_new, std::string meshName);

	//from given data identifier obtain the corresponding simulation value
	Any GetDataValue(DatumConfig dConfig);

	//this is GetDataValue but with std::string conversion
	std::string GetDataValueString(DatumConfig dConfig, bool ignore_unit = false);

	//make a new entry in saveDataList
	void NewSaveDataEntry(DATA_ dataId, std::string meshName = "", Rect dataRect = Rect());
	void EditSaveDataEntry(int index, DATA_ dataId, std::string meshName = "", Rect dataRect = Rect());

	//update mesh names in saveDataList
	void UpdateSaveDataEntries(std::string oldMeshName, std::string newMeshName);

	//update rectangles in saveDataList (rect_old and rect_new are the old and new rectangles for the given mesh. All data rectangles are scaled accordingly.
	void UpdateSaveDataEntries(Rect rect_old, Rect rect_new, std::string meshName);

	//delete save data entries which depend on given meshName
	void DeleteSaveDataEntries(std::string meshName);

	//save currently configured data for saving (in saveDataList) to save data file (savedataFile in directory)
	void SaveData(void);
	//This method does the actual writing to disk : can be launched asynchronously from SaveData method
	void SaveData_DiskBufferFlush(std::vector<std::string>* pdiskbuffer, int* diskbuffer_position);

#if GRAPHICS == 1

	//-------------------------------------Mesh display methods

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. Uses animation to change view
	void AutoSetMeshDisplay(void) { BD.AutoSetMeshDisplaySettings(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. Uses animation to change view but keeps camera view orientation.
	void AutoSetMeshDisplay_KeepOrientation(void) { BD.AutoSetMeshDisplaySettings_KeepOrientation(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. No animation used, just set view
	void AutoSetMeshDisplay_Sudden(void) { BD.AutoSetMeshDisplaySettings_Sudden(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//update mesh viewer with physical quantity - do not call directly, call it through UpdateScreen(); Also refreshes screen, updating interactive objects
	void UpdateMeshDisplay(bool asynchronous = false) { BD.UpdateMeshDisplay(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel()), asynchronous); }
	//update mesh viewer with physical quantity - do not call directly, call it through UpdateScreen_Quick(); Also refreshes screen, but does not update interactive objects, so slightly quicker than UpdateMeshDisplay
	void UpdateMeshDisplay_Quick(bool asynchronous = false) { BD.UpdateMeshDisplay_Quick(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel()), asynchronous); }

	//-------------------------------------Display update methods

	//refresh screen (so similar to RefreshScreen();) but also update displayed quantites. 
	void UpdateScreen(bool asynchronous = false) { UpdateDataBox(); UpdateMeshDisplay(asynchronous); }
	//quicker version where interactive objects are not updated
	void UpdateScreen_Quick(bool asynchronous = false) { UpdateDataBox(); UpdateMeshDisplay_Quick(asynchronous); }

	//update screen and set default view settings - animate view change from current to new settings
	void UpdateScreen_AutoSet(void) { UpdateDataBox(); AutoSetMeshDisplay(); RefreshScreen(); }

	void UpdateScreen_AutoSet_KeepOrientation(void) { UpdateDataBox(); AutoSetMeshDisplay_KeepOrientation(); RefreshScreen(); }

	//update screen and set default view settings - no animation, just set new view
	void UpdateScreen_AutoSet_Sudden(void) { UpdateDataBox(); AutoSetMeshDisplay_Sudden(); RefreshScreen(); }

	//only update data box, with refresh (needed in very special cases, e.g. inter-object interactions in the data box)
	void UpdateDataBox_Refresh(void) { UpdateDataBox(); RefreshScreen(); }

	//clear console and make sure entries in data box match those stored in dataBoxList
	void ClearScreen(void) { BD.ClearScreen(); RebuildDataBox(); RefreshScreen(); }

#else
	
	//DUMMY METHODS TO KEEP THE REST OF THE CODE CLEAN
	//-------------------------------------Mesh display methods

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. Uses animation to change view
	void AutoSetMeshDisplay(void) {}

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. Uses animation to change view but keeps camera view orientation.
	void AutoSetMeshDisplay_KeepOrientation(void) {}

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. No animation used, just set view
	void AutoSetMeshDisplay_Sudden(void) {}

	//update mesh viewer with physical quantity - do not call directly, call it through UpdateScreen(); Also refreshes screen, updating interactive objects
	void UpdateMeshDisplay(bool asynchronous = false) {}
	//update mesh viewer with physical quantity - do not call directly, call it through UpdateScreen_Quick(); Also refreshes screen, but does not update interactive objects, so slightly quicker than UpdateMeshDisplay
	void UpdateMeshDisplay_Quick(bool asynchronous = false) {}

	//-------------------------------------Display update methods

	//refresh screen (so similar to RefreshScreen();) but also update displayed quantites. 
	void UpdateScreen(bool asynchronous = false) {}
	void UpdateScreen_Quick(bool asynchronous = false) {}

	//update screen and set default view settings - animate view change from current to new settings
	void UpdateScreen_AutoSet(void) {}

	void UpdateScreen_AutoSet_KeepOrientation(void) {}

	//update screen and set default view settings - no animation, just set new view
	void UpdateScreen_AutoSet_Sudden(void) {}

	//only update data box, with refresh (needed in very special cases, e.g. inter-object interactions in the data box)
	void UpdateDataBox_Refresh(void) {}

	//clear console and make sure entries in data box match those stored in dataBoxList
	void ClearScreen(void) {}

#endif

public:

	//-------------------------------------Error Handler methods

	//intended to create a restore point to which we can go back if a critical error occurs
	void create_restore(void);

	//restore state after a critical error occurs
	void restore_state(void);

	//show the error
	void show_error(BError error, std::string error_message, bool verbose = true);

	//-------------------------------------Constructor

#if GRAPHICS == 1
	Simulation(
		HWND hWnd, 
		int Program_Version, std::string progName_,
		std::string server_port_ = "", std::string server_pwd_ = "", 
		int cudaDevice = -1,
		std::string python_script_ = "", int python_script_mGPU_ = 1, std::vector<int> python_script_parallel_ = {}, int python_script_terminate_ = 1, int python_script_deletesource = 0);
#else
	Simulation(
		int Program_Version, std::string progName_,
		std::string server_port_ = "", std::string server_pwd_ = "",
		int cudaDevice = -1,
		std::string python_script_ = "", int python_script_mGPU_ = 1, std::vector<int> python_script_parallel_ = {}, int python_script_terminate_ = 1, int python_script_deletesource = 0);
#endif
	~Simulation();

	//-------------------------------------Drawing

#if GRAPHICS == 1
	//Draw screen on request -> delegated to BorisDisplay using thread-safe access 
	//Refresh : full refresh including interactive objects
	void RefreshScreen(void) { BD.Refresh_ThreadSafe(); }
	//Draw : only draw, do not refresh interactive objects
	void DrawScreen(void) { BD.Draw_ThreadSafe(); }
#else
	//DUMMY METHODS TO KEEP THE REST OF THE CODE CLEAN
	void RefreshScreen(void) {}
	void DrawScreen(void) {}
#endif

	//-------------------------------------Handle User Messages

#if GRAPHICS == 1
	//Process message received from WndProc -> delegated to BorisDisplay using thread-safe access
	void NewMessage(AC_ aCode, INT2 mouse, std::string data = "");
#else
	void NewMessage(std::string message);
#endif

	//-------------------------------------Console Messages

	void DisplayConsoleMessage(std::string text) { BD.DisplayConsoleMessage(text); }
	void DisplayConsoleListing(std::string text) { BD.DisplayConsoleListing(text); }
	void DisplayConsoleError(std::string text) { BD.DisplayConsoleError(text); }

	//-------------------------------------Embedded Python Communication

#if PYTHON_EMBEDDING == 1
	//Receive and execute a message from embedded Python script, and return response. Special messages handled separately (see below)
	std::string NewPythonMessage(std::string message);
	
	//Received run command from embedded Python script - special handling required.
	void NewPythonMessage_Run(void);
	void NewPythonMessage_RunStage(int stage);

	//Received run command from embedded Python script, which also required further Python code to be executed every simulation iteration - special handling required.
	//return 1 if simulation is still running, 0 if not
	int NewPythonMessage_RunWithCode(void);
	//before we can use NewPythonMessage_RunWithCode (in a loop), must call NewPythonMessage_PrepareRunWithCode
	void NewPythonMessage_PrepareRunWithCode(void);

	//Embedded Python script will poll this flag to see if it needs to terminate execution (stopping embedded Python script early is surprisingly not straightforward at all if you don't want to crash everything!!!)
	//To avoid user having to worry about these technical details, Python script to be executed will need to be modified by Boris prior to execution (in a manner invisible to end user)
	bool StopPythonScriptFlag(void) { return !python_script_running; }

	std::string Get_Simulations_Directory(void) { return GetUserDocumentsPath() + boris_data_directory + boris_simulations_directory; }
	std::string Get_BorisPythonScripts_Directory(void) { return GetUserDocumentsPath() + boris_data_directory + BorisPythonScripts_directory; }
#endif

	//-------------------------------------Menu Auxiliary Functions

	//Functions called from menu to display various things

	void Show_StartupOptions(void) { Print_StartupOptions(); }

	//Various property checkers required for menu ticks

	//-1 : N/A, 0 : Off, 1 : On
	int GetCUDAStatus(void) { if (!cudaAvailable) return -1; else return cudaEnabled; }
	
	bool GetIndividualShapeStatus(void) { return shape_change_individual; }
	
	bool GetScaleRectsStatus(void) { return SMesh.Get_Scale_Rects(); }
	
	bool GetCoupleToDipolesStatus(void) { return SMesh.Get_Coupled_To_Dipoles(); }
	
	bool GetMovingMeshStatus(void) { return SMesh.IsMovingMeshSet(); }

	bool GetDisabledTransportStatus(void) { return disabled_transport_solver; }

	bool GetScriptServerStatus(void) { return is_thread_running(THREAD_NETWORK); }
};