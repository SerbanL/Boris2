#pragma once

#include <omp.h>

#include "CompileFlags.h"
#include "ErrorHandler.h"

#include "MeshDefs.h"
#include "PhysQRep.h"
#include "SimSharedData.h"

#include "Modules.h"
#include "ParametersDefs.h"
#include "Transport_Defs.h"

#include "MeshParamsBase.h"

#if COMPILECUDA == 1
#include "MeshBaseCUDA.h"
#endif

#include "MeshDefs.h"

class SuperMesh;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Abstract base Mesh class - implement various types of meshes using this base
//
//	This is a base both for micromagnetic and atomistic meshes, so contains only data and methods that are applicable to both
//	Some methods of course need different implementations depending on the context, so these are pure virtual.
//
//	From MeshBase the next derived classes are Mesh and Atom_Mesh, which specialise this interface for a Micromagnetic or Atomistic context.
//	These in turn are also interfaces, and are implemented by the specific type of mesh needed, e.g. Ferromagnetic, Antiferromagnetic, etc.
//	or Simple Cubic, BCC, FCC, HCP etc for atomistic meshes
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MeshBase :
	public SimulationSharedData,
	virtual public MeshParamsBase			//need virtual inheritance : see MeshParamsBase.h comments
{

private:

#if COMPILECUDA == 1
	friend MeshBaseCUDA;
#endif

protected:

	//--------Auxiliary

	int OmpThreads;

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//storage for extracted mesh profiles
	int num_profile_averages = 0;
	std::vector<double> profile_storage_dbl;
	std::vector<DBL3> profile_storage_dbl3;

	//auxiliary VEC for computations
	VEC<DBL3> auxVEC;

	//--------Mesh ID

	//the type of mesh
	int meshType;

	//unique mesh id generated when the object is created (meshIdCounter simply increments every time a new mesh is created and sets a new meshId)
	static int meshIdCounter;

	//a unique mesh id number generated when the object is created, using meshIdCounter.
	int meshId = ++meshIdCounter;

	//--------Display Data

	//current quantity displayed on screen
	int displayedPhysicalQuantity;

	//we can also display, for the same mesh, a secondary physical quantity, at the same time as the primary one.
	//in this case we'll want to use transparency for one or both so we can see them at the same time
	int displayedBackgroundPhysicalQuantity = MESHDISPLAY_NONE;

	//type of representation to use for vectorial quantities (see VEC3REP_ enum in PhysQDefs.h)
	int vec3rep = (int)VEC3REP_FULL;

	//if displaying a parameter spatial variation, this value will say which parameter (a value from PARAM_ enum)
	int displayedParamVar = PARAM_GREL;

	//select module with which we calculate Heff and energy density display (MOD_ERROR means none selected, MOD_ALL means show total effective field)
	int Module_Heff_Display = MOD_ALL;
	int Module_Energy_Display = MOD_ERROR;

	//--------Modules

	//effective field modules currently assigned to this mesh
	vector_lut<Modules*> pMod;

	//--------Configuration

	//if this mesh can participate in multilayered demag convolution, you have the option of excluding it (e.g. antiferromagnetic mesh) - by default all meshes with magnetic computation enabled are included.
	bool exclude_from_multiconvdemag = false;

	//----- MONTE-CARLO Data

	// MONTE-CARLO DATA

	//random number generator - used by Monte Carlo methods
	BorisRand prng;

	//Monte-Carlo current cone angle (vary to reach MONTECARLO_TARGETACCEPTANCE)
	double mc_cone_angledeg = 30.0;

	//last Monte-Carlo step acceptance probability (save it so we can read it out)
	double mc_acceptance_rate = 0.0;

	//use parallel Monte-Carlo algorithms?
	bool mc_parallel = true;

	//can disable mc iteration in given mesh
	bool mc_disabled = false;

	// Constrained MONTE-CARLO DATA

	//use constrained Monte-Carlo?
	bool mc_constrain = false;

	//Constrained Monte-Carlo direction
	DBL3 cmc_n = DBL3(1.0, 0.0, 0.0);

public:

#if COMPILECUDA == 1
	//the CUDA version of this Mesh
	MeshBaseCUDA* pMeshBaseCUDA = nullptr;
#endif

	//--------Auxiliary

	//pointer to the supermesh holding this mesh (some modules in this mesh will need to see other meshes).
	SuperMesh* pSMesh;

	//--------Mesh Dimensions

	//the mesh rectangle in meters : this defines the mesh position in the world (in the super-mesh). cellsizes must divide the mesh rectangle, giving the number of cells in each dimension
	Rect meshRect;

	//-----Magnetic properties

	//number of cells (n.x, n.y, n.z)
	SZ3 n = SZ3(1);

	//cellsizes (h.x, h.y, h.z)
	DBL3 h = DBL3(5e-9);

	//-----Electric conduction properties (Electron charge and spin Transport)

	//number of cells for electrical properties
	SZ3 n_e = SZ3(1);

	//cellsize for electrical properties
	DBL3 h_e = DBL3(5e-9);

	//electrical potential - on n_e, h_e mesh
	VEC_VC<double> V;

	//electrical conductivity - on n_e, h_e mesh
	VEC_VC<double> elC;

	//electric field - on n_e, h_e mesh
	VEC_VC<DBL3> E;

	//spin accumulation - on n_e, h_e mesh
	VEC_VC<DBL3> S;

	//-----Thermal conduction properties

	//number of cells for thermal properties
	SZ3 n_t = SZ3(1);

	//cellsize for thermal properties
	DBL3 h_t = DBL3(5e-9);

	//temperature calculated by Heat module (primary temperature, always used for 1-temperature model; for multi-temperature models in metals this is the itinerant electron temperature)
	VEC_VC<double> Temp;

	//lattice temperature used in many-T models
	VEC_VC<double> Temp_l;

	//-----Mechanical properties

	//number of cells for mechanical properties
	SZ3 n_m = SZ3(1);

	//cellsize for mechanical properties
	DBL3 h_m = DBL3(5e-9);

	//mechanical displacement vectors - on n_m, h_m mesh
	VEC_VC<DBL3> u_disp;

	//strain tensor (symmetric):
	//diagonal and off-diagonal components - on n_m, h_m mesh
	//xx, yy, zz
	VEC_VC<DBL3> strain_diag;
	//yz, xz, xy
	VEC_VC<DBL3> strain_odiag;

	//-----Custom Mesh Display

	VEC<DBL3> displayVEC_VEC;
	VEC<double> displayVEC_SCA;

protected:

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS AUXILIARY

	//average into profile_storage_dbl / profile_storage_dbl3
	void average_mesh_profile(std::vector<double>& line_profile);
	void average_mesh_profile(std::vector<DBL3>& line_profile);

public:

	//------------------------CTOR/DTOR

	MeshBase(MESH_ meshType, SuperMesh *pSMesh_);

	virtual ~MeshBase();

	//obtain error_on_create from this mesh, as well as any set modules - return first error found
	BError Error_On_Create(void);

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) = 0;

	//This is a "softer" version of UpdateConfiguration, which can be used any time and doesn't require the object to be Uninitialized; 
	//this will typically involve changing a value across multiple objects, thus better to call this method rather than try to remember which objects need the value changed.
	virtual void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) = 0;

	//switch CUDA state on/off
	virtual BError SwitchCUDAState(bool cudaState) = 0;

	//at the start of each iteration the mesh might have to be prepared (e.g. state flags set)
	virtual void PrepareNewIteration(void) = 0;

#if COMPILECUDA == 1
	virtual void PrepareNewIterationCUDA(void) = 0;
#endif

	//Take a Monte Carlo step in this mesh (overloaded by mesh implementations) using settings in each mesh
	virtual void Iterate_MonteCarlo(double acceptance_rate) {}

#if COMPILECUDA == 1
	virtual void Iterate_MonteCarloCUDA(double acceptance_rate) {}
#endif

	//----------------------------------- OTHER CONTROL METHODS

	//used by move mesh algorithm : shift mesh quantities (e.g. magnetization) by the given shift (metric units) value within this mesh. The shift is along the x-axis direction only (+ve or -ve).
	virtual void MoveMesh(double x_shift) = 0;

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	virtual double CheckMoveMesh(void) = 0;

	//set PBC for required VECs : should only be called from a demag module
	virtual BError Set_Magnetic_PBC(INT3 pbc_images) = 0;
	virtual INT3 Get_Magnetic_PBC(void) = 0;

	void Set_MonteCarlo_Serial(bool status) { mc_parallel = !status; }
	bool Get_MonteCarlo_Serial(void) { return !mc_parallel; }

	void Set_MonteCarlo_Disabled(bool status) { mc_disabled = status; }
	bool Get_MonteCarlo_Disabled(void) { return mc_disabled; }

	bool Get_MonteCarlo_Constrained(void) { return mc_constrain; }

	void Set_MonteCarlo_Constrained(DBL3 cmc_n_);
	DBL3 Get_MonteCarlo_Constrained_Direction(void) { return cmc_n; }

	//----------------------------------- MODULES indexing

	//index by actual index in pMod
	Modules* operator[](const int &modIdx) { return pMod[modIdx]; }

	//index by module ID
	Modules* operator()(const int &modID) { return pMod(modID); }

	//return entire pMod vector (used to access methods in vector_lut on pMod)
	vector_lut<Modules*>& operator()(void) { return pMod; }

	//----------------------------------- MODULES CONTROL :  MeshBaseModules.cpp

	bool IsModuleSet(MOD_ moduleID) { return pMod.is_ID_set(moduleID); }

	//get number of active modules with given ID
	int GetNumModules(MOD_ moduleID) { return pMod.get_num_IDs(moduleID); }

	Modules* GetModule(MOD_ moduleID, int module_number = 0) { if (pMod.is_id_set(INT2(moduleID, module_number))) return pMod[INT2(moduleID, module_number)]; else return nullptr; }

#if COMPILECUDA == 1
	ModulesCUDA* GetCUDAModule(MOD_ moduleId) { if (IsModuleSet(moduleId)) return pMod(moduleId)->GetCUDAModule(); else return nullptr; }
#endif

	//Add module to list of set modules, also deleting any exclusive modules to this one
	//If you set force_add = true then duplicate modules are allowed : in this case you must keep track of them
	virtual BError AddModule(MOD_ moduleID, bool force_add = false) = 0;

	//delete all modules with given id
	void DelModule(MOD_ moduleID);

	//Delete this particular module, if found.
	//This method is intended to be called from within the module itself, asking for it to be deleted - a call to its dtor will be issued.
	//When using this method in this way, the calling module should immediately return as any further data access there will result in bad memory access.
	//The return addess is still valid as it will be stored on the stack.
	void DelModule(Modules* pModule);

	//get all set modules IDs
	std::vector<MOD_> GetModulesIDs(void);

	//---- INITIALIZE MODULES -> Calls Modules::Initialize

	//initialize all modules prior to computations starting
	BError InitializeAllModules(void);

#if COMPILECUDA == 1
	BError InitializeAllModulesCUDA(void);
#endif

	//---- UPDATE MODULES -> Call Modules::UpdateField

	//update computational state of all modules in this mesh; return total energy density -> each module will have a contribution, so sum it
	double UpdateModules(void);

	//update MOD_TRANSPORT module only if set
	virtual void UpdateTransportSolver(void) = 0;

#if COMPILECUDA == 1
	void UpdateModulesCUDA(void);
	virtual void UpdateTransportSolverCUDA(void) = 0;
#endif

	//----------------------------------- MONTE CARLO AUXILIARY :  MeshBaseMonteCarlo.cpp

	//set cone_angle value depending on required traget acceptance rate (acceptance_rate) and current value of mc_acceptance_rate
	void MonteCarlo_AdaptiveAngle(double& cone_angle, double acceptance_rate);

	//----------------------------------- PARAMETERS CONTROL/INFO : MeshBaseParamsControl.cpp

	//set/get mesh base temperature; by default any text equation dependence will be cleared unless indicated specifically not to (e.g. called when setting base temperature value after evaluating the text equation)
	virtual void SetBaseTemperature(double Temperature, bool clear_equation = true) = 0;
	double GetBaseTemperature(void) { return base_temperature; }

	//set text equation for base temperature : when iterating, the base temperature will be evaluated and set using the text equation
	BError SetBaseTemperatureEquation(std::string equation_string, int step);
	void UpdateTEquationUserConstants(void);

	//parameters spatial variation setters

	//update all mesh parameters spatial variation (if needed)
	bool update_meshparam_var(void);

	//update text equations for mesh parameters with user constants, mesh dimensions, Curie temperature, base temperature
	bool update_all_meshparam_equations(void);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a std::string to be interpreted and converted using ToNum)
	BError set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, std::string generatorArgs, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader);

	//set parameter spatial variation using a shape : set value in given shape only
	BError set_meshparam_shape(PARAM_ paramID, std::vector<MeshShape> shapes, std::string value_text);

	//others	

	//copy all parameters from another Mesh
	virtual BError copy_mesh_parameters(MeshBase& copy_this) = 0;

	//each parameter has a certain type (PARAMTYPE_) - return cellsize in this mesh associated with this parameter type (e.g. ferromagnetic, electric, thermal, elastic)
	DBL3 get_paramtype_cellsize(PARAM_ paramID);

	//set tensorial anisotropy terms
	virtual BError set_tensorial_anisotropy(std::vector<DBL4> Kt) { return BError(); }

	DBL2 Get_MonteCarlo_Params(void) { return DBL2(mc_cone_angledeg, mc_acceptance_rate); }

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	//Return quantity currently set to display on screen.
	//Detail_level : the number of displayed elements in each mesh is obtained by dividing its rectangle by this cubic cell (to the nearest integer in each dimension).
	//detail_level is only needed here when CUDA is enabled so we can fetch data from gpu memory to a coarser array.
	//getBackground : if true then use displayedBackgroundPhysicalQuantity instead of displayedPhysicalQuantity
	virtual PhysQ FetchOnScreenPhysicalQuantity(double detail_level = 0.0, bool getBackground = false) = 0;

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	virtual BError SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType, MESHDISPLAY_ quantity) = 0;

	//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
	//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
	//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
	virtual void GetPhysicalQuantityProfile(
		DBL3 start, DBL3 end, double step, DBL3 stencil, 
		std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, 
		bool do_average, bool read_average, MESHDISPLAY_ quantity) = 0;
	
	//read number of averages performed on profile
	int GetProfileAverages(void) { return num_profile_averages; }

	//return average value for currently displayed mesh quantity in the given relative rectangle
	virtual Any GetAverageDisplayedMeshValue(Rect rel_rect, std::vector<MeshShape> shapes, MESHDISPLAY_ quantity) = 0;

	int GetDisplayedPhysicalQuantity(void) { return displayedPhysicalQuantity; }
	int GetDisplayedBackgroundPhysicalQuantity(void) { return displayedBackgroundPhysicalQuantity; }

	int GetDisplayedParamVar(void) { return displayedParamVar; }

	void SetDisplayedPhysicalQuantity(int displayedPhysicalQuantity_)
	{
		if (displayedPhysicalQuantity_ != displayedBackgroundPhysicalQuantity || displayedPhysicalQuantity_ == MESHDISPLAY_NONE) displayedPhysicalQuantity = displayedPhysicalQuantity_;
	}

	void SetDisplayedBackgroundPhysicalQuantity(int displayedBackgroundPhysicalQuantity_)
	{
		if (displayedBackgroundPhysicalQuantity_ != displayedPhysicalQuantity || displayedBackgroundPhysicalQuantity_ == MESHDISPLAY_NONE) {

			if (displayedBackgroundPhysicalQuantity_ == displayedBackgroundPhysicalQuantity) displayedBackgroundPhysicalQuantity = MESHDISPLAY_NONE;
			else displayedBackgroundPhysicalQuantity = displayedBackgroundPhysicalQuantity_;
		}
	}

	void SetDisplayedParamVar(int displayedParamVar_) { displayedParamVar = displayedParamVar_; }

	void SetVEC3Rep(int vec3rep_) { vec3rep = vec3rep_; }
	int GetVEC3Rep(void) { return vec3rep; }

	bool IsDisplayBackgroundEnabled(void) { return displayedBackgroundPhysicalQuantity != MESHDISPLAY_NONE; }

	//Get settings for module display data 
	//Return module displaying its effective field (MOD_ALL means total Heff)
	int Get_Module_Heff_Display(void) { return Module_Heff_Display; }
	//Return module displaying its energy density spatial variation (MOD_ERROR means none set)
	int Get_Module_Energy_Display(void) { return Module_Energy_Display; }

	//this is similar to Get_Module_Heff_Display, but also takes into account which module is actually set, since e.g. if Module_Heff_Display == MOD_EXCHANGE, but we actually have DM type module set, then need to return the DM module not MOD_EXCHANGE
	//i.e. the meaning of the Module_Heff_Display value changes depending on what modules are actually set. (MOD_EXCHANGE means total exchange value, but MOD_DMEXCHANGE means just the DMI exchange).
	//thus: normally this function returns same as Get_Module_Heff_Display, but exceptions handled here
	int Get_ActualModule_Heff_Display(void)
	{
		if (Module_Heff_Display == MOD_EXCHANGE) {
			if (IsModuleSet(MOD_DMEXCHANGE)) return MOD_DMEXCHANGE;
			else if (IsModuleSet(MOD_IDMEXCHANGE)) return MOD_IDMEXCHANGE;
			else if (IsModuleSet(MOD_VIDMEXCHANGE)) return MOD_VIDMEXCHANGE;
		}
		if (Module_Heff_Display == MOD_DEMAG) {
			if (IsModuleSet(MOD_SDEMAG_DEMAG)) return MOD_SDEMAG_DEMAG;
		}
		return Module_Heff_Display;
	}

	//as above but for energy
	int Get_ActualModule_Energy_Display(void)
	{
		if (Module_Energy_Display == MOD_EXCHANGE) {
			if (IsModuleSet(MOD_DMEXCHANGE)) return MOD_DMEXCHANGE;
			else if (IsModuleSet(MOD_IDMEXCHANGE)) return MOD_IDMEXCHANGE;
			else if (IsModuleSet(MOD_VIDMEXCHANGE)) return MOD_VIDMEXCHANGE;
		}
		if (Module_Energy_Display == MOD_DEMAG) {
			if (IsModuleSet(MOD_SDEMAG_DEMAG)) return MOD_SDEMAG_DEMAG;
		}
		return Module_Energy_Display;
	}

	//Set module display quantity : must call UpdateConfiguration so modules can allocate / free energy as required
	void Set_Module_Heff_Display(int Module_Heff_Display_) { Module_Heff_Display = Module_Heff_Display_; }
	void Set_Module_Energy_Display(int Module_Energy_Display_) { Module_Energy_Display = Module_Energy_Display_; }
	void Set_Module_Display(int Module_Display_) { Module_Heff_Display = Module_Display_; Module_Energy_Display = Module_Display_; }

	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS : MeshBase.cpp

	int get_id(void) { return meshId; }
	void swap_ids(MeshBase* pMesh_swap);

	MESH_ GetMeshType(void) { return (MESH_)meshType; }

	Rect GetMeshRect(void) { return meshRect; }
	DBL3 GetMeshDimensions(void) { return meshRect.size(); }
	DBL3 GetOrigin(void) { return meshRect.s; }

	//set new mesh rectangle and adjust cellsizes to default values if needed - the cellsizes can be manually set after changing the mesh rectangle
	//there's the option of adjusting the number of cells (adjust_num_cells) so set limits are not exceeded (useful when creating a new large mesh which might exceed memory availability with too small a cellsize)
	virtual BError SetMeshRect(Rect meshRect_, bool adjust_num_cells) = 0;

	//ferromagnetic properties
	SZ3 GetMeshSize(void) { return n; }
	//electrical conduction properties
	SZ3 GetMeshESize(void) { return n_e; }
	//thermal conduction properties
	SZ3 GetMeshTSize(void) { return n_t; }

	//ferromagnetic properties
	virtual BError SetMeshCellsize(DBL3 h_) = 0;
	DBL3 GetMeshCellsize(void) { return h; }

	//electrical conduction properties
	virtual BError SetMeshECellsize(DBL3 h_e_) = 0;
	DBL3 GetMeshECellsize(void) { return h_e; }

	//thermal conduction properties
	virtual BError SetMeshTCellsize(DBL3 h_t_) = 0;
	DBL3 GetMeshTCellsize(void) { return h_t; }

	//mechanical properties
	virtual BError SetMeshMCellsize(DBL3 h_m_) = 0;
	DBL3 GetMeshMCellsize(void) { return h_m; }

	//others

	virtual double Get_NonEmpty_Magnetic_Volume(void) = 0;

	bool is_atomistic(void) { return meshBaseType == MESHTYPE_ATOMISTIC; }

	//Set/Get multilayered demag exclusion : will need to call UpdateConfiguration from the supermesh when the flag is changed, so the correct SDemag_Demag modules and related settings are set from the SDemag module.
	//Thus only call Set_Demag_Exclusion directly from SMesh.
	void Set_Demag_Exclusion(bool exclude_from_multiconvdemag_) { exclude_from_multiconvdemag = exclude_from_multiconvdemag_; }
	bool Get_Demag_Exclusion(void) { return exclude_from_multiconvdemag; }

	//search save data list (saveDataList) for given dataID set for this mesh. Return true if found and its rectangle is not Null or is not the entire mesh; else return false.
	bool IsOutputDataSet_withRect(int datumId);
	//return true if data is set (with any rectangle)
	bool IsOutputDataSet(int datumId);

	//check if given stage is set
	bool IsStageSet(int stageType);

	//set computefields_if_MC flag on SuperMesh
	void Set_Force_MonteCarlo_ComputeFields(bool status);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//is magnetic dynamics computation enabled? (check Heff not M - e.g. M is not empty for dipole meshes)
	virtual bool MComputation_Enabled(void) = 0;

	//slighly more general than MComputation_Enabled - this means the mesh can have magnetic properties but doesn't necessarily have an ODE set, e.g. dipole mesh
	virtual bool Magnetism_Enabled(void) = 0;

	//is electrical conduction computation enabled?
	virtual bool EComputation_Enabled(void) = 0;

	//is thermal conduction computation enabled?
	virtual bool TComputation_Enabled(void) = 0;

	//is mechanical computation enabled?
	virtual bool MechComputation_Enabled(void) = 0;

	//check if interface conductance is enabled (for spin transport solver)
	virtual bool GInterface_Enabled(void) = 0;

	//are periodic boundary conditions set?
	virtual int Is_PBC_x(void) = 0;
	virtual int Is_PBC_y(void) = 0;
	virtual int Is_PBC_z(void) = 0;

	//is there a demag-type module set for this mesh? (SDemag not included as this is a SuperMesh module)
	virtual bool Is_Demag_Enabled(void) = 0;

	virtual bool iSHA_nonzero(void) = 0;
	virtual bool SHA_nonzero(void) = 0;

	//----------------------------------- VALUE GETTERS

	//get energy value for given module or one of its exclusive modules (if none active return 0); call it with MOD_ALL to return total energy density in this mesh;
	//If avRect is null then get energy density in entire mesh
	double GetEnergyDensity(MOD_ moduleType, Rect avRect = Rect());

	//as above but for the torque
	DBL3 GetTorque(MOD_ moduleType, Rect avRect = Rect());

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	virtual double GetTopologicalCharge(Rect rectangle = Rect()) = 0;

	virtual double GetAverageElectricalPotential(Rect rectangle = Rect()) = 0;
	virtual DBL3 GetAverageSpinAccumulation(Rect rectangle = Rect()) = 0;

	virtual double GetAverageElectricalConductivity(Rect rectangle = Rect()) = 0;

	//get base temperature or average temperature (if Temp enabled)
	virtual double GetAverageTemperature(Rect rectangle = Rect()) = 0;
	virtual double GetAverageLatticeTemperature(Rect rectangle = Rect()) = 0;

	//----------------------------------- OTHER CALCULATION METHODS

	//compute topological charge density spatial dependence and have it available to display in Cust_S
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	virtual void Compute_TopoChargeDensity(void) = 0;

	//----------------------------------- OTHER MESH SHAPE CONTROL

	virtual BError copy_mesh_data(MeshBase& copy_this) = 0;

	//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
	virtual BError applymask(double zDepth_m, std::string fileName, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader) = 0;

	//set cells to empty in given rectangle (delete by setting entries to zero). The rectangle is relative to this mesh.
	virtual BError delrect(Rect rectangle) = 0;

	//set cells to non-empty in given box
	virtual BError setrect(Rect rectangle) = 0;

	//roughen mesh sides (side = "x", "y", "z", "-x", "-y", or "-z") to given depth (same units as h) with prng instantiated with given seed.
	virtual BError RoughenMeshSides(std::string side, double depth, int seed) = 0;

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	virtual BError RoughenMeshSurfaces_Jagged(double depth, double spacing, int seed, std::string sides) = 0;

	//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	virtual BError GenerateGrains2D(double spacing, int seed) = 0;

	//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	virtual BError GenerateGrains3D(double spacing, int seed) = 0;

	//Advanced mesh shape control methods
	//Method: or (add shape) / not (delete shape) / xor (add and delete overlaps)
	
	//Disk with dimensions (x, y diameters, thickness), centre position relative to mesh, rotation angles, number of repetitions along x, y, z (1, 1, 1 for no repetitions), displacement values if repetitions used
	virtual BError shape_disk(MeshShape shape) = 0;

	//rectangle shape
	virtual BError shape_rect(MeshShape shape) = 0;

	//isosceles triangle shape
	virtual BError shape_triangle(MeshShape shape) = 0;

	//prolate ellipsoid
	virtual BError shape_ellipsoid(MeshShape shape) = 0;

	//pyramid
	virtual BError shape_pyramid(MeshShape shape) = 0;

	//tetrahedron
	virtual BError shape_tetrahedron(MeshShape shape) = 0;

	//cone
	virtual BError shape_cone(MeshShape shape) = 0;

	//torus
	virtual BError shape_torus(MeshShape shape) = 0;

	//general shape setting function, can set composite shape using combination of the above elementary shapes
	virtual BError shape_set(std::vector<MeshShape> shapes) = 0;

	//----------------------------------- METHODS REDEFINED IN SOME IMPLEMENTATIONS (virtual here)

	//use virtual to allow calling using base pointer - if called on "incorrect" mesh then nothing happens (the versions defined here used instead)

	virtual bool GetMoveMeshTrigger(void) { return false; }
	virtual void SetMoveMeshTrigger(bool status) {}

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in output data, not with ShowData or with data box)
	virtual DBL2 Get_skyshift(Rect skyRect) { return DBL2(); }

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in output data, not with ShowData or with data box), as well as diameters along x and y directions.
	virtual DBL4 Get_skypos_diameters(Rect skyRect) { return DBL4(); }

	//set/get skypos tracker rect size diameter multiplier
	virtual double Get_skypos_dmul(void) { return 0.0; }
	virtual void Set_skypos_dmul(double dia_mul_) {}

	//couple this mesh to touching dipoles by setting skip cells as required : used for domain wall moving mesh algorithms
	virtual void CoupleToDipoles(bool status) {}

	//Fit domain wall along the x direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
	virtual DBL2 FitDomainWall_X(Rect rectangle) { return DBL2(); }
	//Fit domain wall along the y direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
	virtual DBL2 FitDomainWall_Y(Rect rectangle) { return DBL2(); }
	//Fit domain wall along the z direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
	virtual DBL2 FitDomainWall_Z(Rect rectangle) { return DBL2(); }

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	virtual DBL3 Average_dmdt(Rect avRect) { return DBL3(); }
	virtual DBL3 Average_dmdt2(Rect avRect) { return DBL3(); }

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	virtual DBL3 Average_mxdmdt(Rect avRect) { return DBL3(); }
	//for sub-lattice B
	virtual DBL3 Average_mxdmdt2(Rect avRect) { return DBL3(); }
	//mixed sub-lattices A and B
	virtual DBL3 Average_mxdm2dt(Rect avRect) { return DBL3(); }
	virtual DBL3 Average_m2xdmdt(Rect avRect) { return DBL3(); }

	//compute magnitude histogram data
	//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first. 
	//output probabilities in histogram_p, corresponding to values set in histogram_x min, min + bin, ..., max, where bin = (max - min) / (num_bins - 1)
	//if macrocell_dims is not INT3(1) then first average in macrocells containing given number of individual mesh cells, then obtain histogram
	virtual bool Get_Histogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims) { return true; }

	//As for Get_Histogram, but use thermal averaging in each macrocell
	virtual bool Get_ThAvHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims) { return true; }

	//angular deviation histogram computed from ndir unit vector direction. If ndir not given (DBL3()), then angular deviation computed from average magnetization direction
	virtual bool Get_AngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir) { return true; }

	//As for Get_AngHistogram, but use thermal averaging in each macrocell
	virtual bool Get_ThAvAngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir) { return true; }

	virtual DBL3 GetThermodynamicAverageMagnetization(Rect rectangle) { return 0.0; }

	//shift a dipole mesh rectangle by given amount (only overload in dipole meshes)
	virtual void Shift_Dipole(DBL3 shift) {}

	//methods for dipole shifting algorithm
	virtual void Set_Dipole_Velocity(DBL3 velocity, DBL3 clipping) {}
	virtual DBL3 Get_Dipole_Velocity(void) { return DBL3(); }
	virtual DBL3 Get_Dipole_Clipping(void) { return DBL3(); }

	virtual void SetTMRType(TMR_ type) {}
	virtual TMR_ GetTMRType(void) { return TMR_NONE; }

	//----------------------------------- MESH QUANTITIES CONTROL

	//this method is also used by the dipole mesh where it does something else - sets the dipole direction
	virtual void SetMagAngle(double polar, double azim, Rect rectangle = Rect()) {}

	virtual void SetMagAngle_Shape(double polar, double azim, std::vector<MeshShape> shapes) {}

	//Set magnetization angle in solid object only containing given relative position uniformly using polar coordinates
	virtual void SetMagAngle_Object(double polar, double azim, DBL3 position) {}

	//Flower state magnetization
	virtual void SetMagFlower(int direction, DBL3 centre, double radius, double thickness) {}

	//Onion state magnetization
	virtual void SetMagOnion(int direction, DBL3 centre, double radius1, double radius2, double thickness) {}

	//Crosstie state magnetization
	virtual void SetMagCrosstie(int direction, DBL3 centre, double radius, double thickness) {}

	//Invert magnetization direction in given mesh (must be magnetic)
	virtual void SetInvertedMag(bool x, bool y, bool z) {}

	//Mirror magnetization in given axis (literal x, y, or z) in given mesh (must be magnetic)
	virtual void SetMirroredMag(std::string axis) {}

	//Set random magentisation distribution in given mesh (must be magnetic)
	virtual void SetRandomMag(int seed) {}
	virtual void SetRandomXYMag(int seed) {}

	//set a domain wall with given width (metric units) at position within mesh (metric units). 
	//Longitudinal and transverse are magnetization componets as: 1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
	virtual void SetMagDomainWall(int longitudinal, int transverse, double width, double position) {}

	//set Neel skyrmion with given orientation (core is up: 1, core is down: -1), chirality (1 for towards centre, -1 away from it) in given rectangle (relative to mesh), calculated in the x-y plane
	virtual void SetSkyrmion(int orientation, int chirality, Rect skyrmion_rect) {}

	//set Bloch skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
	virtual void SetSkyrmionBloch(int orientation, int chirality, Rect skyrmion_rect) {}

	//set M from given data VEC (0 values mean empty points) -> stretch data to M dimensions if needed.
	virtual void SetMagFromData(VEC<DBL3>& data, const Rect& dstRect = Rect()) {}

	//Set Temp from given data VEC -> stretch data to mesh dimensions if needed.
	void SetTempFromData(VEC<double>& data, const Rect& dstRect = Rect());

	//Set E from current density data, by dividing by electrical conductivity (elC). Set V and S to zero. This is meant for computations with fixed Jc and transport solver iteration disabled.
	void SetEFromJcData(VEC<DBL3>& data);
	
	//set electric field VEC from a constant Jc value
	void SetEFromJcValue(DBL3 Jcvalue);

	//Set/Get mesh exchange coupling status to other meshes
	virtual void SetMeshExchangeCoupling(bool status) {}
	virtual bool GetMeshExchangeCoupling(void) { return false; }

	//----------------------------------- MODULE METHODS TEMPLATED CALLERS

	//IMPORTANT NOTE: read note for CallModuleMethod in SuperMesh, repeated here:
	//IMPORTANT NOTE: These templated callers work by trying to cast the Modules* (Module is the abstract base class) to a derived implementation type (i.e. to Owner*) - dynamic_cast results in nullptr if couldn't cast.
	//				  Thus it is important that Owner is actually the implementation and not the base class Modules. If runThisMethod points to a method in Modules, e.g. GetEnergyDensity, don't use the template deduction mechanism!
	//				  e.g. if you use CallModuleMethod(&STransport::GetEnergyDensity), then Owner will not be STransport but Modules, so dynamic_cast will succeed on first attempt which may not be the right module!
	//				  In this case explicitly specify the template parameters as : CallModuleMethod<double, STransport>(&STransport::GetEnergyDensity)
	//FURTHER NOTE:   Some modules have an additional ABC, in addition to Modules, e.g. ZeemanBase ABC used by Zeeman and Atom_Zeeman implementations which result in polymorphic behaviour of micromagnetic and atomistic meshes
	//				  CallModuleMethod also works in this case by calling a method defined in the secondary base class, e.g. &ZeemanBase::SetField:
	//				  In this case dynamic cast will be attempted from Modules* to ZeemanBase*, and this only succeeds if the Modules ABC implementation is also an implementation of the additional ABC (ZeemanBase in this example).

	//Call a method on a mesh module (if module available, which is automatically determined here)
	template <typename RType, typename Owner>
	RType CallModuleMethod(RType(Owner::*runThisMethod)())
	{
		for (int idx = 0; idx < pMod.size(); idx++) {

			Owner* pOwner = dynamic_cast<Owner*>(pMod[idx]);
			if (pOwner) return CALLFP(pOwner, runThisMethod)();
		}

		return RType();
	}

	template <typename RType, typename Owner, typename ... PType>
	RType CallModuleMethod(RType(Owner::*runThisMethod)(PType ...), PType ... params)
	{
		for (int idx = 0; idx < pMod.size(); idx++) {

			Owner* pOwner = dynamic_cast<Owner*>(pMod[idx]);
			if (pOwner) return CALLFP(pOwner, runThisMethod)(params...);
		}

		return RType();
	}

	//similar to CallModuleMethod, but instead of calling just the first module found, call all matching modules in this mesh. No return value possible.
	template <typename Owner>
	void CallAllModulesMethod(void(Owner::*runThisMethod)())
	{
		for (int idx = 0; idx < pMod.size(); idx++) {

			Owner* pOwner = dynamic_cast<Owner*>(pMod[idx]);
			if (pOwner) CALLFP(pOwner, runThisMethod)();
		}
	}

	template <typename Owner, typename ... PType>
	void CallAllModulesMethod(void(Owner::*runThisMethod)(PType ...), PType ... params)
	{
		for (int idx = 0; idx < pMod.size(); idx++) {

			Owner* pOwner = dynamic_cast<Owner*>(pMod[idx]);
			if (pOwner) CALLFP(pOwner, runThisMethod)(params...);
		}
	}
};