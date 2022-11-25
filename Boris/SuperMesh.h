#pragma once

#include "ErrorHandler.h"

#include "SimSharedData.h"

#include "Mesh_Ferromagnetic.h"
#include "Mesh_AntiFerromagnetic.h"
#include "Mesh_Dipole.h"
#include "Mesh_Metal.h"
#include "Mesh_Insulator.h"

#include "Atom_Mesh_Cubic.h"

#include "SDemag.h"
#include "StrayField.h"
#include "STransport.h"
#include "Oersted.h"
#include "SHeat.h"
#include "SMElastic.h"

#if COMPILECUDA == 1
#include "SuperMeshCUDA.h"
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	This is a container of simulation meshes and modules which span two or more meshes (supermesh modules).
//	Coordinates the simulation flow between the different meshes.
//	Performs calculations which cannot be done at individual mesh level (e.g. demag field over multiple ferromagnetic meshes).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SuperMesh : 
	public SimulationSharedData, 
	public ProgramState<SuperMesh,
	std::tuple<
	int, int, 
	SZ3, DBL3, Rect, SZ3, DBL3, Rect,
	VEC<DBL3>,
	ODECommon, Atom_ODECommon,
	vector_key<MeshBase*>, 
	vector_lut<Modules*>, 
	std::string, std::string, 
	bool, bool, int,
	bool,
	bool, DBL2>,
	std::tuple<
	//Micromagnetic Meshes
	FMesh, DipoleMesh, MetalMesh, InsulatorMesh, AFMesh,
	//Atomistic Meshes
	Atom_Mesh_Cubic,
	//supermesh modules
	SDemag, StrayField, STransport, Oersted, SHeat, SMElastic> >
{
	//all supermesh modules are friends
	friend SDemag;
	friend StrayField;
	friend STransport;
	friend SHeat;
	friend SMElastic;
	friend Oersted;

	friend Mesh;
	friend Atom_Mesh;

#if COMPILECUDA == 1
	friend SuperMeshCUDA;
	friend MeshCUDA;
	friend Atom_MeshCUDA;

	friend SDemagCUDA;
	friend StrayFieldCUDA;
	friend STransportCUDA;
	friend SHeatCUDA;
	friend OerstedCUDA;

public:

	//the CUDA version of this Mesh
	SuperMeshCUDA* pSMeshCUDA = nullptr;
#endif

private:

	//-----

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//-----Display data

	//current quantity displayed on screen for super-mesh
	int displayedPhysicalQuantity = MESHDISPLAY_NONE;

	//type of representation to use for vectorial quantities (see VEC3REP_ enum in PhysQDefs.h)
	int vec3rep = (int)VEC3REP_FULL;

	//name of currently active mesh in the pMesh vector (used for quicker data inputting in the console)
	std::string activeMeshName = "permalloy";

	//storage for extracted supermesh profiles
	std::vector<double> profile_storage_dbl;
	std::vector<DBL3> profile_storage_dbl3;

	//-----Supermesh dimensions

	//the entire super-mesh rectangle (rectangle containing all meshes)
	Rect sMeshRect;

	//-----Supermesh dimensions : magnetic

	//ferromagnetic super-mesh dimensions (this is the super-mesh that encompasses all ferromagnetic meshes)
	SZ3 n_fm;

	//ferromagnetic super-mesh discretization cellsize
	DBL3 h_fm = DBL3(5e-9);

	//ferromagnetic super-mesh rectangle
	Rect sMeshRect_fm;

	//-----Supermesh dimensions : electric

	//electric super-mesh dimensions (this is the super-mesh that encompasses all meshes with electrical computations enabled)
	SZ3 n_e;

	//electric super-mesh discretization cellsize
	DBL3 h_e = DBL3(5e-9);

	//electric super-mesh rectangle
	Rect sMeshRect_e;

	//-----Special Supermesh Quantities

	//a set global field with its own discretization and rectangle (independent of supermesh, only held by it)
	//global field is used by Zeeman modules as an additional contribution
	VEC<DBL3> globalField;

	//-----Micromagnetics ODE

	//the micromagnetics equation to solve, together with its evaluation method is common to all meshes.
	//This is the "static" base class (i.e. contains a bunch of properties which are shared by all particular ODE solvers in the different ferromagnetic meshes, together with methods to manage them)
	//There are particular ODE objects in the ferromagnetic meshes, derived from ODECommon, since they differ in dimensions.
	ODECommon odeSolver;

	//-----Atomistic ODE

	//the atomistic equation to solve (in atomistic meshes), together with its evaluation method is common to all atomistic meshes.
	//This is the "static" base class (i.e. contains a bunch of properties which are shared by all particular atomistic ODE solvers in the different atomistic meshes, together with methods to manage them)
	//There are particular ODE objects in the atomistic meshes, derived from Atom_ODECommon, since they differ in dimensions.
	Atom_ODECommon atom_odeSolver;

	//-----Meshes

	//individual simulation meshes (key is the mesh name)
	vector_key<MeshBase*> pMesh;

	//-----Supermesh Modules

	//super-mesh effective field modules (used for long range interactions)
	vector_lut<Modules*> pSMod;

	//-----Options

	//when changing a mesh rectangle, scale all other rectangles in proportion
	bool scale_rects = false;

	//if ferromagnetic meshes touch a dipole mesh then interface magnetic cells are frozen (ode doesn't update them - use skip cell flags) if this flag is true
	//Moreover the interface cells are set with magnetization direction along the dipole magnetization direction. This is an easy way of simulating exchange coupling to the dipoles.
	bool coupled_dipoles = false;
	
	//in CUDA mode initialize kernels on GPU. Can be set to false to initialize on CPU (slightly more accurate, but not enough to make this the default option)
	bool kernel_initialize_on_gpu = true;

	//-----Mesh data settings

	//select which component to use when fitting to obtain domain wall width and position for dwpos_x, dwpos_y, dwpos_z parameters
	//-1: automatic (detect tanh component), 0: x, 1: y, 2: z. You might want to set to defined component if too noisy to reliably auto detect.
	int dwpos_component = -1;

	//-----Monte Carlo settings

	//if running a Monte Carlo algorithm normally we don't want to update fields as well, but set this true if needed (e.g. if you also want to save energy density values for which we need to update fields)
	bool computefields_if_MC = false;

	//if certain modules are added then we need to force effective fields to be computed (e.g. dipole-dipole or demag modules)
	//this flag is similar to computefields_if_MC, but is purely controlled by the respective modules on initialization
	bool force_computefields_if_MC = false;

	//MC cone angle limits. If min and max are the same this turn the adaptive MC algorithms into fixed cone angle.
	DBL2 cone_angle_minmax = DBL2(MONTECARLO_CONEANGLEDEG_MIN, MONTECARLO_CONEANGLEDEG_MAX);

	//-----Auxiliary

	//the total calculated energy density -> return by UpdateModules() or UpdateField() methods.
	double total_energy_density = 0.0;
	//different meshes have different weights when contributing to the total energy density -> ratio of their non-empty volume to total non-empty volume
	//this vector is calculated at initialization and has same size as pMesh vector
	std::vector<double> energy_density_weights;

public:

	//name of super-mesh for use in console (e.g. addmodule supermesh sdemag). It is also a reserved name : no other module can be named with this handle
	std::string superMeshHandle = "supermesh";

private:

public:

	//--------------------------------------------------------- CTOR/DTOR

	SuperMesh(void);
	~SuperMesh();

	//obtain error_on_create from supermesh, as well as any meshes and super-mesh modules added - return first error found
	BError Error_On_Create(void);

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) { error_on_create = UpdateConfiguration(UPDATECONFIG_REPAIROBJECTSTATE); }

	//--------------------------------------------------------- MESH INDEXING and related methods

	//index the super-mesh using mesh name to return a mesh pointer (assumes the referenced mesh is contained in the pMesh vector).
	MeshBase* operator[](const std::string &key) { return pMesh[key]; }
	MeshBase* operator[](const int &meshIdx) { return pMesh[meshIdx]; }

	//obtain mesh which contains the given coordinate (units m)
	MeshBase* operator[](DBL3 coordinate)
	{
		for (int idx = 0; idx < pMesh.size(); idx++)
			if (pMesh[idx]->GetMeshRect().contains(coordinate)) return pMesh[idx];
		return nullptr;
	}

	MeshBase* active_mesh(void) { return pMesh[activeMeshName]; }

	//return entire pMesh vector (used to access methods in vector_key on pMesh)
	vector_key<MeshBase*>& operator()(void) { return pMesh; }

	//check if SMesh contains the given mesh : identify it by its key, or by its unique mesh id - for the latter return the index in pMesh (-1 if not found)
	bool contains(std::string &key) { return pMesh.has_key(key); }
	int contains_id(int meshId) { for (int idx = 0; idx < (int)pMesh.size(); idx++) if (pMesh[idx]->get_id() == meshId) return idx; return -1; }

	//from unique mesh id get key of entry in pMesh vector
	std::string key_from_meshId(int meshId) { for (int idx = 0; idx < (int)pMesh.size(); idx++) if (pMesh[idx]->get_id() == meshId) return pMesh.get_key_from_index(idx); return ""; }
	std::string key_from_meshIdx(int meshIdx) { return pMesh.get_key_from_index(meshIdx); }

	int size(void) { return (int)pMesh.size(); }

	//---------------------------------------------------------IMPORTANT CONTROL METHODS : SuperMeshControl.cpp

	//call whenever changes are made to meshes (e.g. dimensions, cellsize, modules added/deleted)
	//This is the master UpdateConfiguration method, which then cascades to meshes, modules and ODE solvers
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	
	//This is a "softer" version of UpdateConfiguration, which can be used any time and doesn't require any simulation objects to be Uninitialized
	//this will typically involve changing a value across multiple objects, thus better to call this method rather than try to remember which objects need the value changed.
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	//switch CUDA state on/off; when switching on set selected cuda device number (from 1 up)
	BError SwitchCUDAState(bool cudaState, int __cudaDeviceSelect);

	//---------------------------------------------------------MULTI-MESH CONTROL METHODS : SuperMeshControl.cpp

	//couple ferromagnetic meshes to any touching dipole meshes, setting interface cell values and flags
	void CoupleToDipoles(void);

	//--------------------------------------------------------- SIMULATION CONTROL : SuperMeshSimulation.cpp

	//initialize all modules in all meshes, including the ones in super-mesh
	BError InitializeAllModules(void);

#if COMPILECUDA == 1
	BError InitializeAllModulesCUDA(void);
#endif

	//called by Simulation to advance simulation by a time step
	void AdvanceTime(void);

	//Similar to AdvanceTime but only computes effective fields and does not run the ODE solver
	void ComputeFields(void);

	//iterate transport solver only, if available
	void UpdateTransportSolver(void);

#if COMPILECUDA == 1
	void AdvanceTimeCUDA(void);
	void ComputeFieldsCUDA(void);
	void UpdateTransportSolverCUDA(void);
#endif

	//Take a Monte Carlo step over all atomistic meshes using settings in each mesh; increase the iterations counters.
	void Iterate_MonteCarlo(double acceptance_rate);

#if COMPILECUDA == 1
	void Iterate_MonteCarloCUDA(double acceptance_rate);
#endif

	//----------------------------------- ODE SOLVER CONTROL  : SuperMeshODE.cpp

	//Reset all ODE solvers in meshes with on ODE
	void ResetODE(void);

	//set new stage for ODE solvers
	void NewStageODE(void);

	//set the ode and evaluation method. Any new ODE in a magnetic mesh will use these settings. Currently Micromagnetic and Atomistic ODEs use the same evaluation method.
	BError SetODE(ODE_ setOde, EVAL_ evalMethod);
	//same for the atomistic ODE. Currently Micromagnetic and Atomistic ODEs use the same evaluation method.
	BError SetAtomisticODE(ODE_ setOde, EVAL_ evalMethod);

	//set ODE evaluation method, applicable to both micromagnetic and atomistic solvers
	BError SetODEEval(EVAL_ evalMethod);

	//set the time step for the magnetization solver
	void SetTimeStep(double dT);

	//set parameters for adaptive time step control
	void SetAdaptiveTimeStepCtrl(double err_fail, double dT_incr, double dT_min, double dT_max);

	void SetStochTimeStep(double dTstoch);
	double GetStochTimeStep(void);
	void SetLink_dTstoch(bool flag);
	bool GetLink_dTstoch(void);

	void SetSpeedupTimeStep(double dTspeedup);
	double GetSpeedupTimeStep(void);
	void SetLink_dTspeedup(bool flag);
	bool GetLink_dTspeedup(void);
	
	//set evaluation speedup type in ode solver
	void SetEvaluationSpeedup(int status);
	//check evaluation speedup settings in ode solver
	int GetEvaluationSpeedup(void);

	//is the current time step fully finished? - most evaluation schemes need multiple sub-steps
	bool CurrentTimeStepSolved(void);

	//check in ODECommon the type of field update we need to do depending on the ODE evaluation step
	int Check_Step_Update(void);
	//get total time with evaluation step resolution level
	double Get_EvalStep_Time(void);

	//check if ODE solver needs spin accumulation solved
	bool SolveSpinCurrent(void);

	void SetMoveMeshAntisymmetric(bool antisymmetric);
	void SetMoveMeshThreshold(double threshold);

	void SetMoveMeshTrigger(bool status, std::string meshName = "");

	//---Other ODE Getters

	//get set ODE and evaluation method
	void QueryODE(ODE_ &setODE, EVAL_ &evalMethod);
	void QueryODE(ODE_ &setODE);

	//get set atomistic ODE and evaluation method
	void QueryAtomODE(ODE_ &setODE, EVAL_ &evalMethod);
	void QueryAtomODE(ODE_ &setODE);

	int GetIteration(void);
	int GetStageIteration(void);

	double GetTime(void);
	double GetStageTime(void);
	
	double GetTimeStep(void);
	
	double Get_mxh(void);
	double Get_dmdt(void);

	double Get_AStepRelErrCtrl(void);
	DBL3 Get_AStepdTCtrl(void);

	bool IsMovingMeshSet(void);
	int GetId_of_MoveMeshTrigger(void);
	double Get_dwshift(void);

	bool MoveMeshAntisymmetric(void);
	double MoveMeshThreshold(void);

#if COMPILECUDA == 1
	cu_obj<ManagedDiffEq_CommonCUDA>& Get_ManagedDiffEq_CommonCUDA(void) { return odeSolver.Get_pODECUDA()->Get_ManagedDiffEq_CommonCUDA(); }
	cu_obj<ManagedAtom_DiffEq_CommonCUDA>& Get_ManagedAtom_DiffEq_CommonCUDA(void) { return atom_odeSolver.Get_pODECUDA()->Get_ManagedAtom_DiffEq_CommonCUDA(); }
#endif

	//---------------------------------------------------------MONTE-CARLO SOLVER CONTROL : SuperMesh_MonteCarlo.cpp

	//switch to serial (true) or parallel (false) in given mesh - all if meshName is the supermesh handle
	BError Set_MonteCarlo_Serial(bool status, std::string meshName);

	//switch to constrained Monnte-Carlo (true) or classical (false) in given mesh - all if meshName is the supermesh handle; if constrained, then use cmc_n direction.
	BError Set_MonteCarlo_Constrained(bool status, DBL3 cmc_n, std::string meshName);

	//Disable/enable MC iteration in named mesh
	BError Set_MonteCarlo_Disabled(bool status, std::string meshName);

	void Set_MonteCarlo_ComputeFields(bool status) { computefields_if_MC = status; }
	bool Get_MonteCarlo_ComputeFields(void) { return computefields_if_MC || force_computefields_if_MC; }

	void Set_Force_MonteCarlo_ComputeFields(bool status) { force_computefields_if_MC = status; }

	void Set_MonteCarlo_ConeAngleLimits(DBL2 cone_angle_minmax_) { cone_angle_minmax = cone_angle_minmax_; }
	DBL2 Get_MonteCarlo_ConeAngleLimits(void) { return cone_angle_minmax; }

	//--------------------------------------------------------- MESH HANDLING - COMPONENTS : SuperMeshMeshes.cpp

	//Add a new mesh of given type, name and dimensions
	BError AddMesh(std::string meshName, MESH_ meshType, Rect meshRect);

	//delete a ferromagnetic mesh
	BError DelMesh(std::string meshName);

	//rename a mesh
	BError RenameMesh(std::string oldName, std::string newName);

	//change mesh focus
	BError SetMeshFocus(std::string meshName);
	//get name of mesh in focus
	std::string GetMeshFocus(void) { return activeMeshName; }

	//--------------------------------------------------------- MESH HANDLING - SHAPES : SuperMeshMeshes_Shapes.cpp

	//set mesh rect for named mesh (any if applicable any other dependent meshes) and update dependent save data rects by calling the provided function.
	BError SetMeshRect(std::string meshName, Rect meshRect, std::function<void(std::string, Rect)> save_data_updater);

	//copy all primary mesh data (magnetization, elC, Temp, etc.) but do not change dimensions or discretisation
	BError copy_mesh_data(std::string meshName_from, std::string meshName_to);

	//delete or add rectangle for given mesh
	BError delrect(std::string meshName, Rect rectangle);
	BError setrect(std::string meshName, Rect rectangle);
	BError resetrect(std::string meshName);

	//roughen mesh sides (side = "x", "y", "z", "-x", "-y", or "-z") to given depth (same units as h) with prng instantiated with given seed.
	BError RoughenMeshSides(std::string meshName, std::string side, double depth, int seed);

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	BError RoughenMeshSurfaces_Jagged(std::string meshName, double depth, double spacing, int seed, std::string sides);

	//clear roughness: set fine shape to coarse shape
	BError ClearMeshRoughness(std::string meshName);

	//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains2D(std::string meshName, double spacing, int seed);

	//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains3D(std::string meshName, double spacing, int seed);

	//--------------------------------------------------------- MESH HANDLING - SETTINGS : SuperMeshMeshes_Settings.cpp

	//set magnetization angle in given mesh (all if meshName not given)
	BError SetMagAngle(std::string meshName, double polar, double azim);

	//set magnetization angle in given mesh (must be specified), and only in given rectangle of meshName (relative coordinates)
	BError SetMagAngle_Rect(std::string meshName, double polar, double azim, Rect rectangle);

	//set magnetization angle in given mesh using given shape
	BError SetMagAngle_Shape(std::string meshName, double polar, double azim, std::vector<MeshShape> shapes);

	//Set magnetization angle in solid object only containing given relative position uniformly using polar coordinates
	BError SetMagAngle_Object(std::string meshName, double polar, double azim, DBL3 position);

	//Flower state magnetization
	BError SetMagFlower(std::string meshName, int direction, DBL3 centre, double radius, double thickness);

	//Onion state magnetization
	BError SetMagOnion(std::string meshName, int direction, DBL3 centre, double radius1, double radius2, double thickness);

	//Crosstie state magnetization
	BError SetMagCrosstie(std::string meshName, int direction, DBL3 centre, double radius, double thickness);
	
	//Invert magnetization direction in given mesh (must be magnetic)
	BError SetInvertedMag(std::string meshName, bool x = true, bool y = true, bool z = true);

	//Mirror magnetization in given axis (literal x, y, or z) in given mesh (must be magnetic)
	BError SetMirroredMag(std::string meshName, std::string axis);

	//Set random magnetization distribution in given mesh (must be magnetic)
	BError SetRandomMag(std::string meshName, int seed);
	BError SetRandomXYMag(std::string meshName, int seed);

	//longitudinal and transverse are the components specified as std::string literals : "-z", "-y", "-x", "x", "y", "z"
	BError SetMagDomainWall(std::string meshName, std::string longitudinal, std::string transverse, double width, double position);

	//Set Neel skyrmion in given mesh with chirality, diameter and relative x-y plane position
	BError SetSkyrmion(std::string meshName, int orientation, int chirality, double diameter, DBL2 position);
	
	//Set Bloch skyrmion in given mesh with chirality, diameter and relative x-y plane position
	BError SetSkyrmionBloch(std::string meshName, int orientation, int chirality, double diameter, DBL2 position);

	//Set applied magnetic field
	BError SetField(std::string meshName, DBL3 field_cartesian);

	//Set uniform applied mechanical stress
	BError SetUniformStress(std::string meshName, DBL3 stress_cartesian);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - transverse wall
	BError PrepareMovingMesh(std::string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - bloch wall
	BError PrepareMovingMesh_Bloch(std::string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - neel wall
	BError PrepareMovingMesh_Neel(std::string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - skyrmion
	BError PrepareMovingMesh_Skyrmion(std::string meshName);

	//Clear moving mesh settings made by a prepare method
	void ClearMovingMesh(void);

	//set ferromagnetic mesh roughness refinement if Roughness module enabled in given mesh
	BError SetMeshRoughnessRefinement(std::string meshName, INT3 refine);

	//Set periodic boundary conditions for magnetization
	//possible flags: x, y, z
	BError Set_PBC(std::string meshName, std::string flag, int images);

	//set exchange coupling to neighboring meshes - an exchange-type module (i.e. inherit from ExchangeBase) must be enabled in the named mesh
	BError Set_ExchangeCoupledMeshes(bool status, std::string meshName);

	//Set/Get multilayered demag exclusion : will need to call UpdateConfiguration when the flag is changed, so the correct SDemag_Demag modules and related settings are set from the SDemag module.
	BError Set_Demag_Exclusion(bool exclude_from_multiconvdemag, std::string meshName);

	//set link_stochastic flag in named mesh, or all meshes if supermesh handle given
	BError SetLinkStochastic(bool link_stochastic, std::string meshName);

	//set electric field VEC from a constant Jc value in named mesh
	BError SetEFromJcValue(DBL3 Jcvalue, std::string meshName);

	//Set TMR type in named mesh (must be an insulator mesh, or leave blank to apply to all meshes)
	BError SetTMRType(std::string meshName, TMR_ TMR_type);

	//set text equation for RAp and RAap in insulator mesh with tmr module added
	BError SetTMR_BiasEquationParallel(std::string meshName, std::string equation_string);
	BError SetTMR_BiasEquationAntiParallel(std::string meshName, std::string equation_string);

	//set text equation for TAMR conductivity in transport module
	BError Set_TAMR_Conductivity_Equation(std::string meshName, std::string equation_string);

	//load field in supermesh (globalField) or in named mesh
	BError LoadOVF2Field(std::string meshName, std::string fileName);
	BError ClearGlobalField(void);
	//shift globalField rectangle
	void ShiftGlobalField(DBL3 shift);

	//set prng_seed value in given mesh (all meshes if name empty)
	BError Set_PRNG_Seed(std::string meshName, unsigned seed);

	//--------------------------------------------------------- MESH PARAMETERS : SuperMeshParams.cpp

	//these set parameter values and temperature dependence in the indicated mesh - call through these since it's important to call UpdateConfiguration also
	BError set_meshparam_value(std::string meshName, std::string paramHandle, std::string value_text);
	
	//get named parameter value from given mesh. Set value as a std::string in value_text, without units
	BError get_meshparam_value(std::string meshName, std::string paramHandle, std::string& value_text);

	//temperature dependence

	BError set_meshparam_t_equation(std::string meshName, std::string paramHandle, std::string equationText);
	
	BError set_meshparam_tscaling_array(std::string meshName, std::string paramHandle, std::vector<double>& temp, std::vector<double>& scaling_x, std::vector<double>& scaling_y, std::vector<double>& scaling_z, std::string fileName_info = "");

	//clear parameters temperature dependence in given mesh (all meshes if empty std::string)
	BError clear_meshparam_temp(std::string meshName, std::string paramHandle);
	
	//spatial dependence

	BError set_meshparam_s_equation(std::string meshName, std::string paramHandle, std::string equationText);

	//clear parameter spatial dependence (variation) in given mesh for named parameter only (if meshName empty then in all meshes; if paramHandle empty then for all parameters)
	BError clear_meshparam_variation(std::string meshName, std::string paramHandle);

	//set parameter to display in given mesh when ParamVar spatial variation display is enabled
	BError set_meshparamvar_display(std::string meshName, std::string paramHandle);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a std::string to be interpreted and converted using ToNum)
	BError set_meshparam_var(std::string meshName, std::string paramHandle, std::string generatorHandle, std::string generatorArgs, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader);

	//set parameter spatial variation using a shape : set value in given shape only
	BError set_meshparam_shape(std::string meshName, std::string paramHandle, std::vector<MeshShape> shapes, std::string value_text);

	//copy all parameters from another Mesh
	BError copy_mesh_parameters(std::string meshName_from, std::string meshName_to);

	//--------------------------------------------------------- TEMPERATURE / HEAT SOLVER CONTROL : SuperMeshTemperature.cpp

	//set mesh base temperature. If spatial variation set and Heat module enabled then non-uniform base temperature will be set
	BError SetBaseTemperature(std::string meshName, double Temperature);

	//ambient and alpha boundary coefficient for Robin boundary conditions - set in Heat module if active
	BError SetAmbientTemperature(std::string meshName, double T_ambient);
	BError SetAlphaHeatBoundary(std::string meshName, double alpha_boundary);
	
	//insulating mesh sides for heat equation (Neumann boundary conditions). literal can be "x", "-x", "y", "-y", "z", "-z"
	BError SetInsulatingSides(std::string meshName, std::string literal, bool status);

	//set Curie temperature/atomic moment as Bohr magneton multiple for named mesh or all meshes (if meshName is the supermesh handle)
	//this is for the actually set Tc value
	//applicable for micromagnetic meshes only
	BError SetCurieTemperature(std::string meshName, double T_Curie);
	
	//this is for the indicative material Tc value
	//applicable for micromagnetic meshes only
	BError SetCurieTemperatureMaterial(std::string meshName, double T_Curie_material);
	BError SetAtomicMagneticMoment(std::string meshName, DBL2 atomic_moment);

	//set Tc (critical temperature) coupling terms for 2-sublattice model
	//applicable for micromagnetic meshes only
	BError SetTcCoupling(std::string meshName, DBL2 tau_ii, DBL2 tau_ij);
	BError SetTcCoupling_Intra(std::string meshName, DBL2 tau_ii);
	BError SetTcCoupling_Inter(std::string meshName, DBL2 tau_ij);

	//Set temperature model
	BError SetTemperatureModel(std::string meshName, int tmtype);

	//--------------------------------------------------------- ELASTODYNAMICS SOLVER CONTROL : SuperMeshElastodynamics.cpp

	//reset elastodynamics solver state
	BError Reset_ElSolver(void);

	//set text equation in given mesh for diagonal strain
	BError Set_Sd_Equation(std::string meshName, std::string text_equation);

	//set text equation in given mesh for shear strain
	BError Set_Sod_Equation(std::string meshName, std::string text_equation);

	//clear text equations in given mesh (or all meshes if meshName is empty)
	BError Clear_Sd_Sod_Equations(std::string meshName);

	//add a fixed surface for elastodynamics solver (rectangle in absolute coordinates)
	//can also specify a given face in a given mesh (face -x, x, -y, y, -z, or z)
	BError Add_Fixed_Surface(std::string meshName, std::string face, Rect surface_rect);

	//add a stress surface for elastodynamics solver (rectangle in absolute coordinates) with set vector equation
	//can also specify a given face in a given mesh (face -x, x, -y, y, -z, or z)
	BError Add_Stress_Surface(std::string meshName, std::string face, Rect surface_rect, std::string equation);

	//Delete fixed or stressed surfaces with given index (-1 to delete all)
	BError Del_Fixed_Surface(int index);
	BError Del_Stress_Surface(int index);

	//edit values of existing fixed or stress surface
	BError Edit_Fixed_Surface(int index, Rect rect);
	BError Edit_Stress_Surface_Rectangle(int index, Rect rect);
	BError Edit_Stress_Surface_Equation(int index, std::string equation);

	//----------------------------------- MODULES CONTROL : SuperMeshModules.cpp

	//Add module to given mesh (delegate implementation to AddModule in the referenced mesh), checking for super-mesh module clashes
	BError AddModule(std::string meshName, MOD_ moduleId);

	//Delete module from given mesh (delegate implementation to DelModule in the referenced mesh), checking for super-mesh module clashes
	BError DelModule(std::string meshName, MOD_ moduleId);

	bool IsSuperMeshModuleSet(MOD_ moduleId) { return pSMod.is_ID_set(moduleId); }

	Modules* GetSuperMeshModule(MOD_ moduleId) { if (IsSuperMeshModuleSet(moduleId)) return pSMod(moduleId); else return nullptr; }

	//--------------------------------------------------------- GET/SET PROPERTIES / VALUES at SuperMesh level : SuperMesh_Settings.cpp

	//--------Getters

	Rect GetSMeshRect(void) { return sMeshRect; }

	DBL3 GetFMSMeshCellsize(void) { return h_fm; }
	SZ3 GetFMSMeshsize(void) { return n_fm; }
	Rect GetFMSMeshRect(void) { return sMeshRect_fm; }

	DBL3 GetESMeshCellsize(void) { return h_e; }
	SZ3 GetESMeshsize(void) { return n_e; }
	Rect GetESMeshRect(void) { return sMeshRect_e; }

	bool Get_Scale_Rects(void) { return scale_rects; }

	bool Get_Coupled_To_Dipoles(void) { return coupled_dipoles; }

	bool Get_Kernel_Initialize_on_GPU(void) { return kernel_initialize_on_gpu; }

	int Get_DWPos_Component(void) { return dwpos_component; }

	//get total volume energy density
	double GetTotalEnergyDensity(void);

	//search save data list (saveDataList) for given dataID set for given mesh. Return true if found and its rectangle is not Null; else return false.
	bool IsOutputDataSet_withRect(int datumId, MeshBase* pmesh);

	//return true if data is set (with any rectangle)
	bool IsOutputDataSet(int datumId, MeshBase* pmesh);

	//check if given stage is set
	bool IsStageSet(int stageType);

	VEC<DBL3>& GetGlobalField(void) { return globalField; }
#if COMPILECUDA == 1
	//before calling this, make sure cuda is switched on (so pSMeshCUDA not nullptr)
	cu_obj<cuVEC<cuReal3>>& GetGlobalFieldCUDA(void) { return pSMeshCUDA->GetGlobalField(); }
#endif

	//--------Getters for supermesh modules specific properties

	//status for multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
	int Get_Multilayered_Convolution_Status(void);

	//status for force 2d multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
	int Get_2D_Multilayered_Convolution_Status(void);

	//status for use default n_common for multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
	int Use_Default_n_Status(void);

	//get n_common for multilayerd convolution (if return = SZ3() : N/A)
	SZ3 Get_n_common(void);

	//--------Setters

	BError SetFMSMeshCellsize(DBL3 h_fm_);
	BError SetESMeshCellsize(DBL3 h_e_);

	void Set_Scale_Rects(bool status) { scale_rects = status; }

	void Set_Coupled_To_Dipoles(bool status) { coupled_dipoles = status; CoupleToDipoles(); }

	BError Set_Kernel_Initialize_on_GPU(bool status) { kernel_initialize_on_gpu = status; return UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE); }

	void Set_DWPos_Component(int component) { dwpos_component = component; }

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS : SuperMeshDisplay.cpp

	std::vector<PhysQ> FetchOnScreenPhysicalQuantity(double detail_level = 0.0);
	
	//save the quantity currently displayed on screen for named mesh in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(std::string meshName, std::string fileName, std::string ovf2_dataType, MESHDISPLAY_ quantity = MESHDISPLAY_NONE);

	//extract profile from named mesh, from currently display mesh quantity, but reading directly from the quantity
	//Displayed mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
	//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
	void GetPhysicalQuantityProfile(
		DBL3 start, DBL3 end, double step, DBL3 stencil, 
		std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, 
		std::string meshName, bool do_average, bool read_average, MESHDISPLAY_ quantity = MESHDISPLAY_NONE);

	//return average value for currently displayed mesh quantity for named mesh in the given relative rectangle
	Any GetAverageDisplayedMeshValue(std::string meshName, Rect rel_rect, std::vector<MeshShape> shapes = {}, MESHDISPLAY_ quantity = MESHDISPLAY_NONE);

	int GetDisplayedPhysicalQuantity(void) { return displayedPhysicalQuantity; }

	BError SetDisplayedPhysicalQuantity(std::string meshName, int displayedPhysicalQuantity_);
	BError SetDisplayedBackgroundPhysicalQuantity(std::string meshName, int displayedBackgroundPhysicalQuantity_);

	//Get/Set vectorial quantity representation options in named mesh (which could be the supermesh)
	BError SetVEC3Rep(std::string meshName, int vec3rep_);
	int GetVEC3Rep(std::string meshName);

	//--------------------------------------------------------- MODULE METHODS TEMPLATED CALLERS

	//IMPORTANT NOTE: These templated callers work by trying to cast the Modules* (Module is the abstract base class) to a derived implementation type (i.e. to Owner*) - dynamic_cast results in nullptr if couldn't cast.
	//				  Thus it is important that Owner is actually the implementation and not the base class Modules. If runThisMethod points to a method in Modules, e.g. GetEnergyDensity, don't use the template deduction mechanism!
	//				  e.g. if you use CallModuleMethod(&STransport::GetEnergyDensity), then Owner will not be STransport but Modules, so dynamic_cast will succeed on first attempt which may not be the right module!
	//				  In this case explicitly specify the template parameters as : CallModuleMethod<double, STransport>(&STransport::GetEnergyDensity)

	//Call a method on a supermesh module (if module available)
	template <typename RType, typename Owner>
	RType CallModuleMethod(RType(Owner::*runThisMethod)())
	{
		for (int idx = 0; idx < pSMod.size(); idx++) {

			Owner* pOwner = dynamic_cast<Owner*>(pSMod[idx]);
			if (pOwner) return CALLFP(pOwner, runThisMethod)();
		}

		return RType();
	}

	template <typename RType, typename Owner, typename ... PType>
	RType CallModuleMethod(RType(Owner::*runThisMethod)(PType ...), PType ... params)
	{ 
		for (int idx = 0; idx < pSMod.size(); idx++) {

			Owner* pOwner = dynamic_cast<Owner*>(pSMod[idx]);
			if (pOwner) return CALLFP(pOwner, runThisMethod)(params...);
		}

		return RType();
	}
};