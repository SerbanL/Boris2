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
//
//Did this and didn't cause any problems, but need to be aware this isn't officially supported and could potentially cause problems.
//I understand CUDA 10 officially supports all VS2017 versions but didn't try it yet (maybe I should!). VS2019 also available but no need to upgrade yet.

//3.
//Cannot import BorisLib.h in .cu files. I believe this may be due to CUDA 9.2 nvcc not supporting C++14 standards (only C++11). Tried CUDA 10.2 - same story.
//I got around this by keeping all compilation units in .cu files with C++11 code only, which didn't interfere with program design as I didn't really need any C++14 features there (but very useful in many other compilation units!).

//4.
//On Windows 10 it's possible CUDA nvcc will keep coming up with "exit error code 1", including when you simply try to clean the project.
//This is a generic error, and if you enable maximum verbosity in VS (Tools > Options > Projects & Solutions > Build and Run) you'll see this is because of "Access is denied." for a process which writes some temporary data.
//In particular when the process is trying to write temporary data to "C:\Users\<UserName>\AppData\Local" it won't be able to because of denied access.
//I got around this by giving full control for this folder in Security settings. This problem appeared out of the blue and I don't really know why, I suspect it's something to do with a Windows Defender update or something similar.

//5.
//Don't set the number of CUDA threads per block too large. Currently using 128.
//If the amount of code that is included in a CUDA kernel is too large (through inlining of various functions, etc. etc.) the program will start to exhibit very strange bugs with no apparent solution, and completely defying logic.
//I've struggled with this for some time, until I realised it's because of too many threads per block used when launching kernels. Transport solver kernels are particularly notorious.
//As soon as I reduced them - hey presto! - logic restored and everything works again. 
//Used to be 512 threads per block, then reduced to 256. When some kernels got too big again, same type of bugs appeared. This time I knew not to waste time and reduced to the current value of 128 - good again!
//Tested this on CUDA 5 and CUDA 6 architectures and didn't find any variation in this behaviour.

//BUGS

//NOT SOLVED:

//LIKELY SOLVED:

//SOLVED:
//1. Drag and drop simulation file sometimes crashes program. Found bad conversion function - I'm certain that was the problem, so consider this solved but keep an eye on this for a while.
//2. Saving simulation file sometimes sets dT to zero (to a floating point error). I've only seen it happen with CUDA enabled. Solution: Checks in ODECommonCUDA destructor fixed. 
//3. Using a python script may result in program hanging if issuing a flood of commands. Solution : reworked WinSocks module.
//4. If heat solver diverges (e.g. due to too high a time step), and at least 1 material parameter has a temperature dependence, when in CUDA mode out of gpu memory errors can result requiring a program restart. 
//I seem to have fixed it using extra checks on Temperature when getting updated parameter values, but I don't understand why this happens without the checks so the solution seems like a hack. Not happy with this!

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

using namespace std;

//Top level object
class Simulation : 
	public SimulationSharedData, 
	public Threads<Simulation>,
	public ProgramState<Simulation, 
	tuple<BorisDisplay, string, string, string, string, bool, bool, bool, vector_lut<DatumConfig>, vector_lut<DatumConfig>, INT2, vector_lut<StageConfig>, int, bool, SuperMesh, bool, bool, bool, DBL4, DBL2, DBL2, int, vector_key<double>>,
	tuple<> >
{
private:

	//saved simulation file header and program version number
	string simfile_header = "Boris simulation file version ";
	
	string domain_name = "www.boris-spintronics.uk";
	string download_page = "/download";
	//check if current program version is the latest
	string version_checker = "version.php";
	//allow users to send a new entry for the online materials database
	string mdb_entry_handler = "mdbentry.php";
	//update the local materials database from the online materials database
	string mdb_update_handler = "mdb.php";
	//check when the online materials database was last updated
	string mdb_lastupdate = "mdbtime.php";

	string boris_data_directory = "Boris Data\\";

	//check for updates on program startup?
	bool start_check_updates = true;

	//start script server on program startup?
	bool start_scriptserver = true;

	//save/load startup flags in this file (e.g. log_errors, start_check_updates, start_scriptserver)
	string startup_options_file;

	int Program_Version;

	//value set by version update checker :
	//-1: attempting to connect
	//0: connection failure
	//1: program up to date
	//2: update available
	int version_checking = -1;

	ErrorHandler<Simulation> err_hndl;

	//log errors generated by command input?
	bool log_errors = true;
	string errorlog_fileName;

	//the display master object - all display updates done through this object
	BorisDisplay BD;

	//communication with scripting clients - listen for incoming commands and respond with requested data
	WinSocks commSocket;

	//Materials data base management - local and online
	MaterialsDB mdb;

	//collection of commands, with the CommandSpecifier indexed by the command name (name as typed in console). The CommandSpecifier describes the properties common to all commands (vector_key more convenient than a map)
	vector_key<CommandSpecifier> commands;

	//lut-indexed output data descriptor : can be indexed using a value from DATA_ enum; can also be indexed using a key: this is the data handle as used in the console
	vector_key_lut<DatumSpecifier> dataDescriptor;

	//files for saving output data
	string savedataFile = "out_data.txt";
	//image save file base (during a simulation this is complemented by _iteration.png termination)
	string imageSaveFileBase = "mesh_image";

	//currently loaded simulation file, including directory (empty if none)
	string currentSimulationFile;
	//append or create new data file?
	bool appendToDataFile = false;
	//saveDataList is in SimulationSharedData

	//flags for enabling data and image saving (they both use the same saving condition in the simulation schedule)
	bool saveDataFlag = true, saveImageFlag = false;

	//data to display in data box
	vector_lut<DatumConfig> dataBoxList;

	//lut-indexed vector of module handles (handle as entered in the console) : can be indexed using a value from MOD_ enum
	vector_lut<string> moduleHandles;

	//For each ODE_ entry specify allowed EVAL_ entries (i.e. for each differential equation specify allowed evaluation methods)
	vector_lut< vector<EVAL_> > odeAllowedEvals;
	//The default evaluation method for each ODE : this is the evaluation method set when changing ODEs
	vector_lut<EVAL_> odeDefaultEval;
	//Link ODE_ entries with text handles (ODE_ is the major id) : micromagnetic meshes
	vector_lut<string> odeHandles;
	//Link ODE_ entries with text handles (ODE_ is the major id) : atomistic meshes (smaller subset of ODEs allowed, but same evaluation methods)
	vector_lut<string> atom_odeHandles;
	//Link EVAL_ entries with text handles (EVAL_ is the major id)
	vector_lut<string> odeEvalHandles;

	//handles and descriptors for simulation stages types, stop conditions and data saving conditions, indexed by SS_ and STOP_ respectively, as well as keys (handles)
	vector_key_lut<StageDescriptor> stageDescriptors;
	vector_key_lut<StageStopDescriptor> stageStopDescriptors;
	vector_key_lut<DataSaveDescriptor> dataSaveDescriptors;

	//interactive console object info : index with majorId (IOI_ entry) and minorId (subtype). If only one entry for majorId then it applies irrespective of minorId.
	vector_lut<string> ioInfo;

	//simulation stages describing the simulation schedule
	vector_lut<StageConfig> simStages;

	//display updating on number of iterations during simulation
	int iterUpdate = 100;

	//autocomplete commands in console
	bool autocomplete = true;

	//default no cropping of saved mesh images : left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
	//use a small back-off (0.001) for default settings to avoid capturing lines at the edges
	DBL4 image_cropping = DBL4(0.001, 0.001, 0.999, 0.999);

	//Simulation super-mesh
	SuperMesh SMesh;

	//internal data processing arrays : used to load values into columns and process data using console commands
	DPArrays dpArr = DPArrays(MAX_ARRAYS);
	
	//are cuda computations available for this program version and hardware?
#if COMPILECUDA == 0
	bool cudaAvailable = false;
#else
	//still need to test for this in the constructor
	bool cudaAvailable;
#endif

	//thread-safe access to simulation state changes - e.g. don't want to change parameters in the middle of an iteration
	mutex simulationMutex;

	//timing values used for benchmarking - use benchtime command to return sim_end_ms - sim_start_ms
	unsigned int sim_start_ms = 0, sim_end_ms = 0;

public:

private:

	//-------------------------------------Startup Options

	//save/load flags: start_check_updates, start_scriptserver, log_errors
	void Load_Startup_Flags(void);
	void Save_Startup_Flags(void);

	//-------------------------------------Program update checker (in BorisIOGenerators.cpp)

	//Check with "www.boris-spintronics.uk" if program version is up to date
	void CheckUpdate(void);
	void OpenDownloadPage(void);

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
	BError SaveSimulation(string fileName);
	BError LoadSimulation(string fileName);

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) { dpArr.clear_all(); }

	//-------------------------------------Simulation control

	//Start and stop Simulate method in a separate thread
	void RunSimulation(void);
	void StopSimulation(void);
	//stop and reset simulation back to starting point
	void ResetSimulation(void);

	//main simulation method - computes a single complete iteration (a full ode time-step)
	void Simulate(void);

	//Similar to Simulate but only runs for one iteration and does not advance time
	void ComputeFields(void);

	//Command handler
	void HandleCommand(string command);

	//-------------------------------------Simulation schedule

	//add a new simulation stage with default settings
	void AddGenericStage(SS_ stageType, string meshName = "");
	void DeleteStage(int stageIndex);

	//set stop condition for index-th entry in simStages (the index-th stage) with default stop value.
	void SetGenericStopCondition(int index, STOP_ stopType);
	//set data save conditon for index-th entry in simStages (the index-th stage) with default saving value.
	void SetGenericDataSaveCondition(int index, DSAVE_ dsaveType);

	//edit set stage type, stop condition, or data save condition
	void EditStageType(int index, SS_ stageType, string meshName = "");
	void EditStageValue(int stageIndex, string value_string);
	void EditStageStopCondition(int index, STOP_ stopType, string stopValueString = "");
	void EditDataSaveCondition(int index, DSAVE_ dsaveType, string dsaveValueString = "");

	//update mesh names in simStages
	void UpdateStageMeshNames(string oldMeshName, string newMeshName);

	//get stage step, step as an INT2, also checking if they are valid
	INT2 Check_and_GetStageStep();

	//Check if conditions for ending stage, step have been met and set values for next stage if so.
	void CheckSimulationSchedule(void);

	//when stage step stopping conditions hve been met call this. Next stage step will be set. Simulation is stopped if schedule finishes.
	void AdvanceSimulationSchedule(void);

	//set the values for current stage_step value, if any, as configured in simStages
	void SetSimulationStageValue(void);

	//check if conditions for saving data have been met for curent stage
	void CheckSaveDataConditions();

	//-------------------------------------Console messages helper methods

	//show usage for given console command
	void PrintCommandUsage(string command_name) { BD.DisplayFormattedConsoleMessage(commands[command_name].usage); BD.DisplayFormattedConsoleMessage(commands[command_name].descr); BD.DisplayFormattedConsoleMessage(commands[command_name].return_descr); }

	//---------------------------------------------------- MESH LIST

	//NOTE, in general we have two main types of methods:

	//1. void Print_something_List(void)			: this displays in the console required formatted text involving a list of things
	//2. string Build_something_ListLine(int index) : this makes the formatted text for a single line in the list for the above method. Also needed in state handler to update the list when something changes.

	//Sometimes methods like this are also use : string Build_something_Text(int index) : this makes unformatted text used in interactive objects textId and displayed text.

	//show list of meshes
	void Print_Mesh_List(void);
	//build formatted string, describing the given mesh
	string Build_Mesh_ListLine(int meshIndex);

	//---------------------------------------------------- MESH DISPLAY

	void Print_MeshDisplay_List(void);
	string Build_MeshDisplay_ListLine(int meshIndex);

	//---------------------------------------------------- MODULES LIST

	//show list of modules (active and not active) for all meshes
	void Print_Modules_List(void);
	//build formatted string, describing the modules available for a given named mesh
	string Build_Modules_ListLine(int meshIndex);

	//---------------------------------------------------- ODES LIST

	//show ODE configuration
	void Print_ODEs(void);

	//---------------------------------------------------- SHOW DATA LIST (for console display)

	//show data available for showdata command
	void Print_ShowData(void);

	//---------------------------------------------------- OUTPUT DATA LIST (for saving during a simulation)

	//Print available and currently set output data as interactive lists
	void Print_AvailableOutputData(void);

	//build formatted string for interactive object describing an entry in saveDataList at index_in_list (also helper method to build the actual unformatted display text in the object)
	void Print_SetOutputData_List(void);
	string Build_SetOutputData_ListLine(int index_in_list);
	string Build_SetOutputData_Text(int index_in_list);

	//---------------------------------------------------- SIMULATION STAGES

	//Print available and currently set simulation stages as interactive lists
	void Print_SetStages_List(void);

	//build formatted string for interactive object describing an entry in simStages at index_in_list (also helper method to build the actual unformatted display text in the object)
	string Build_SetStages_ListLine(int index_in_list);
	string Build_SetStages_Text(int index_in_list);

	//helpers for interactive objects for set stages lines
	string Build_SetStages_ValueText(int index_in_list);
	string Build_SetStages_StopConditionText(int index_in_list);
	string Build_SetStages_SaveConditionText(int index_in_list, int dsaveIdx);

	//---------------------------------------------------- MESH PARAMETERS

	//print all mesh parameters for the given mesh name
	void Print_MeshParams(string meshName);

	//build an entire line with mesh parameters for a given mesh index
	string Build_MeshParams_Line(int meshIndex);

	//build mesh parameter text for a given mesh index and param identifier - helper method
	string Build_MeshParams_Text(int meshIdx, PARAM_ paramId);

	//---------------------------------------------------- MESH PARAMETERS TEMPERATURE DEPENDENCE

	//print all mesh parameters temperature dependence for the given mesh name
	void Print_MeshParamsTemperature(string meshName);

	//build an entire line with mesh parameters temperature dependence for a given mesh index
	string Build_MeshParamsTemp_Text(int meshIndex);

	//---------------------------------------------------- MESH PARAMETERS SPATIAL DEPENDENCE

	//print all mesh parameters spatial dependence for the given mesh name
	void Print_MeshParamsVariation(string meshName);

	//build an entire line with mesh parameters spatial dependence for a given mesh index
	string Build_MeshParamsVariation_Text(int meshIndex);

	//---------------------------------------------------- MOVING MESH SETTINGS

	void PrintMovingMeshSettings(void);

	//---------------------------------------------------- ELECTRODES and TRANSPORT SETTINGS

	void Print_Electrodes_List(void);
	string Build_Electrodes_ListLine(int el_index);

	void PrintTransportSolverConfig(void);

	//---------------------------------------------------- TEMPERATURE

	void Print_MeshTemperature_List(void);
	string Build_MeshTemperature_ListLine(int meshIndex);

	void Print_HeatBoundaries_List(void);
	string Build_HeatBoundaries_ListLine(int meshIndex);

	//---------------------------------------------------- CURIE TEMPERATURE and ATOMIC MOMENT

	void Print_CurieandMoment_List(void);
	string Build_CurieandMoment_ListLine(int meshIndex);
	
	//---------------------------------------------------- TEMPERATURE MODEL TYPE

	void Print_TemperatureModel_List(void);
	string Build_TemperatureModel_ListLine(int meshIndex);

	//---------------------------------------------------- STOCHASTICITY SETIINGS

	void Print_Stochasticity_List(void);
	string Build_Stochasticity_ListLine(int meshIndex);

	//---------------------------------------------------- EVALUATION SPEEDUP SETIINGS

	void Print_Speedup_List(void);
	string Build_Speedup_ListLine(int meshIndex);

	//---------------------------------------------------- CUDA and MEMORY INFO

	void Print_CUDAStatus(void);
	void Print_MemoryInfo(void);

	//---------------------------------------------------- SCALE RECTS STATUS

	void Print_Scale_Rects_Status(void);

	//---------------------------------------------------- COUPLED-To-DIPOLES STATUS

	void Print_CoupledToDipoles_Settings(void);

	//---------------------------------------------------- ERROR LOG STATUS and other STARTUP OPTIONS

	void Print_ErrorLogStatus(void);

	void Print_StartupUpdateCheckStatus(void);
	void Print_StartupScriptServerStatus(void);

	//print all startup options
	void Print_StartupOptions(void);

	//---------------------------------------------------- NEIGHBORING MESHES EXCHANGE COUPLING

	void Print_ExchangeCoupledMeshes_List(void);
	string Build_ExchangeCoupledMeshes_ListLine(int meshIndex);

	//---------------------------------------------------- MESH ROUGHNESS REFINEMENT

	void Print_MeshRoughnessRefinement(string meshName);

	//---------------------------------------------------- MULTILAYERED CONVOLUTION CONFIGURATION

	void Print_MultiConvolution_Config(void);

	//---------------------------------------------------- MATERIALS DATABASE

	void Print_MaterialsDatabase(void);

	//---------------------------------------------------- ADAPTIVE TIME STEP CONTROL

	void Print_AStepCtrl(void);

	//---------------------------------------------------- PERIODIC BOUNDARY CONDITIONS

	void Print_PBC(void);
	string Build_PBC_ListLine(int meshIndex);

	//---------------------------------------------------- INDIVIDUAL SHAPE CONTROL

	void Print_IndividualShapeStatus(void);

	//---------------------------------------------------- USER EQUATION CONSTANTS

	//Print currently set equation constants
	void Print_EquationConstants(void);

	//build formatted string for interactive objects describing user constants at index_in_list from userConstants (also helper method to build the actual unformatted display text in the object)
	string Build_EquationConstants_ListLine(int index_in_list);
	string Build_EquationConstants_Text(int index_in_list);

	//---------------------------------------------------- MAKE INTERACTIVE OBJECT : Auxiliary method

	//Generate a formatted string depending on the interactive object identifier
	template <typename ... PType>
	string MakeIO(IOI_ identifier, PType ... params);

	//construct ioInfo object in constructor
	void MakeIOInfo(void);

	//set entry line in console to given text
	void SetConsoleEntryLineText(string text) { BD.SetConsoleEntryLineText(text); }

	//-------------------------------------Script client

	//Scripted communication: listen for commands and dispatch them to HandleCommand; respond with data after command handled
	void Listen_Incoming_Message(void);

	//disable / enable script server
	void Script_Server_Control(bool status);

	//-------------------------------------Data helper methods

	//Make new entry in data box
	void NewDataBoxField(DatumConfig dConfig);

	//Update all currently set data box entries - do not call directly, call it through UpdateScreen(); or UpdateDataBox_Refresh(); Doesn't refresh screen.
	void UpdateDataBox(void);

	//Clear and remake all data box fields. Doesn't refresh screen.
	void RebuildDataBox(void);

	//delete any data box fields with given mesh name set
	void DeleteDataBoxFields(string meshName);

	//change names of mesh names in data box fields
	void ChangeDataBoxLabels(string oldMeshName, string newMeshName);

	//update rectangles in dataBoxList (rect_old and rect_new are the old and new rectangles for the given mesh. All dataBoxList rectangles are scaled accordingly.
	void UpdateDataBoxEntries(Rect rect_old, Rect rect_new, string meshName);

	//from given data identifier obtain the corresponding simulation value
	Any GetDataValue(DatumConfig dConfig);

	//this is GetDataValue but with string conversion
	string GetDataValueString(DatumConfig dConfig, bool ignore_unit = false);

	//make a new entry in saveDataList
	void NewSaveDataEntry(DATA_ dataId, string meshName = "", Rect dataRect = Rect());
	void EditSaveDataEntry(int index, DATA_ dataId, string meshName = "", Rect dataRect = Rect());

	//update mesh names in saveDataList
	void UpdateSaveDataEntries(string oldMeshName, string newMeshName);

	//update rectangles in saveDataList (rect_old and rect_new are the old and new rectangles for the given mesh. All data rectangles are scaled accordingly.
	void UpdateSaveDataEntries(Rect rect_old, Rect rect_new, string meshName);

	//delete save data entries which depend on given meshName
	void DeleteSaveDataEntries(string meshName);

	//save currently configured data for saving (in saveDataList) to save data file (savedataFile in directory)
	void SaveData(void);

#if GRAPHICS == 1
	//-------------------------------------Mesh display methods

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. Uses animation to change view
	void AutoSetMeshDisplay(void) { BD.AutoSetMeshDisplaySettings(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. Uses animation to change view but keeps camera view orientation.
	void AutoSetMeshDisplay_KeepOrientation(void) { BD.AutoSetMeshDisplaySettings_KeepOrientation(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//set default settings for mesh viewer - do not call directly, call it through UpdateScreen_AutoSet();. Doesn't refresh screen. No animation used, just set view
	void AutoSetMeshDisplay_Sudden(void) { BD.AutoSetMeshDisplaySettings_Sudden(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//update mesh viewer with physical quantity - do not call directly, call it through UpdateScreen();. Doesn't refresh screen.
	void UpdateMeshDisplay(void) { BD.UpdateMeshDisplay(SMesh.FetchOnScreenPhysicalQuantity(BD.Get_MeshDisplay_DetailLevel())); }

	//-------------------------------------Display update methods

	//refresh screen (so similar to RefreshScreen();) but also update displayed quantites. 
	void UpdateScreen(void) { UpdateDataBox(); UpdateMeshDisplay(); RefreshScreen(); }
	//quicker version where interactive objects are not updated
	void UpdateScreen_Quick(void) { UpdateDataBox(); UpdateMeshDisplay(); DrawScreen(); }

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

	//update mesh viewer with physical quantity - do not call directly, call it through UpdateScreen();. Doesn't refresh screen.
	void UpdateMeshDisplay(void) {}

	//-------------------------------------Display update methods

	//refresh screen (so similar to RefreshScreen();) but also update displayed quantites. 
	void UpdateScreen(void) {}
	void DrawScreen(void) {}

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
	void show_error(BError error, string error_message, bool verbose = true);

	//-------------------------------------Other public methods

#if GRAPHICS == 1
	Simulation(HWND hWnd, int Program_Version);
#else
	Simulation(int Program_Version);
#endif
	~Simulation();

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

#if GRAPHICS == 1
	//Process message received from WndProc -> delegated to BorisDisplay using thread-safe access
	void NewMessage(AC_ aCode, INT2 mouse, string data = "");
#else
	void NewMessage(string message);
#endif

	//Functions called from menu to display various things

	void Show_StartupOptions(void) { Print_StartupOptions(); }

	//Various property checkers required for menu ticks

	//-1 : N/A, 0 : Off, 1 : On
	int GetCUDAStatus(void) { if (!cudaAvailable) return -1; else return cudaEnabled; }
	
	bool GetIndividualShapeStatus(void) { return shape_change_individual; }
	
	bool GetScaleRectsStatus(void) { return SMesh.Get_Scale_Rects(); }
	
	bool GetCoupleToDipolesStatus(void) { return SMesh.Get_Coupled_To_Dipoles(); }
	
	bool GetMovingMeshStatus(void) { return SMesh.IsMovingMeshSet(); }

	bool GetScriptServerStatus(void) { return is_thread_running(THREAD_NETWORK); }
};