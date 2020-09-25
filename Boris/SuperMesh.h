#pragma once

#include "ErrorHandler.h"

#include "SimSharedData.h"

#include "Mesh_Ferromagnetic.h"
#include "Mesh_AntiFerromagnetic.h"
#include "Mesh_Diamagnetic.h"
#include "Mesh_Dipole.h"
#include "Mesh_Metal.h"
#include "Mesh_Insulator.h"

#include "Atom_Mesh_Cubic.h"

#include "SDemag.h"
#include "StrayField.h"
#include "STransport.h"
#include "Oersted.h"
#include "SHeat.h"

#if COMPILECUDA == 1
#include "SuperMeshCUDA.h"
#endif

using namespace std;

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
	tuple<
	int, int, 
	SZ3, DBL3, Rect, SZ3, DBL3, Rect, 
	ODECommon, Atom_ODECommon,
	vector_key<MeshBase*>, 
	vector_lut<Modules*>, 
	string, string, 
	bool, bool>,
	tuple<
	//Micromagnetic Meshes
	FMesh, DipoleMesh, MetalMesh, InsulatorMesh, AFMesh, DiaMesh,
	//Atomistic Meshes
	Atom_Mesh_Cubic,
	//supermesh modules
	SDemag, StrayField, STransport, Oersted, SHeat> >
{
	//all supermesh modules are friends
	friend SDemag;
	friend StrayField;
	friend STransport;
	friend SHeat;
	friend Oersted;

	friend Mesh;

#if COMPILECUDA == 1
	friend SuperMeshCUDA;
	friend MeshCUDA;

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
	string activeMeshName = "permalloy";

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
	//Moreover the interface cells are set with magnetisation direction along the dipole magnetisation direction. This is an easy way of simulating exchange coupling to the dipoles.
	bool coupled_dipoles = false;

	//the total calculated energy density -> return by UpdateModules() or UpdateField() methods.
	double total_energy_density = 0.0;
	//different meshes have different weights when contributing to the total energy density -> ratio of their non-empty volume to total non-empty volume
	//this vector is calculated at initialization and has same size as pMesh vector
	vector<double> energy_density_weights;

public:

	//name of super-mesh for use in console (e.g. addmodule supermesh sdemag). It is also a reserved name : no other module can be named with this handle
	string superMeshHandle = "supermesh";

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
	MeshBase* operator[](const string &key) { return pMesh[key]; }
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
	bool contains(string &key) { return pMesh.has_key(key); }
	int contains_id(int meshId) { for (int idx = 0; idx < (int)pMesh.size(); idx++) if (pMesh[idx]->get_id() == meshId) return idx; return -1; }

	//from unique mesh id get key of entry in pMesh vector
	string key_from_meshId(int meshId) { for (int idx = 0; idx < (int)pMesh.size(); idx++) if (pMesh[idx]->get_id() == meshId) return pMesh.get_key_from_index(idx); return ""; }
	string key_from_meshIdx(int meshIdx) { return pMesh.get_key_from_index(meshIdx); }

	int size(void) { return (int)pMesh.size(); }

	//---------------------------------------------------------IMPORTANT CONTROL METHODS : SuperMeshControl.cpp

	//call whenever changes are made to meshes (e.g. dimensions, cellsize, modules added/deleted)
	//This is the master UpdateConfiguration method, which then cascades to meshes, modules and ODE solvers
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	
	//This is a "softer" version of UpdateConfiguration, which can be used any time and doesn't require any simulation objects to be Uninitialized
	//this will typically involve changing a value across multiple objects, thus better to call this method rather than try to remember which objects need the value changed.
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

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

	//set the time step for the magnetisation solver
	void SetTimeStep(double dT);

	//set parameters for adaptive time step control
	void SetAdaptiveTimeStepCtrl(double err_fail, double err_high, double err_low, double dT_incr, double dT_min, double dT_max);

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

	//check if ODE solver needs spin accumulation solved
	bool SolveSpinCurrent(void);

	void SetMoveMeshAntisymmetric(bool antisymmetric);
	void SetMoveMeshThreshold(double threshold);

	void SetMoveMeshTrigger(bool status, string meshName = "");

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

	DBL3 Get_AStepRelErrCtrl(void);
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
	BError Set_MonteCarlo_Serial(bool status, string meshName);

	//switch to constrained Monnte-Carlo (true) or classical (false) in given mesh - all if meshName is the supermesh handle; if constrained, then use cmc_n direction.
	BError Set_MonteCarlo_Constrained(bool status, DBL3 cmc_n, string meshName);

	//--------------------------------------------------------- MESH HANDLING - COMPONENTS : SuperMeshMeshes.cpp

	//Add a new mesh of given type, name and dimensions
	BError AddMesh(string meshName, MESH_ meshType, Rect meshRect);

	//delete a ferromagnetic mesh
	BError DelMesh(string meshName);

	//rename a mesh
	BError RenameMesh(string oldName, string newName);

	//change mesh focus
	BError SetMeshFocus(string meshName);
	//get name of mesh in focus
	string GetMeshFocus(void) { return activeMeshName; }

	//--------------------------------------------------------- MESH HANDLING - SHAPES : SuperMeshMeshes_Shapes.cpp

	//set mesh rect for named mesh (any if applicable any other dependent meshes) and update dependent save data rects by calling the provided function.
	BError SetMeshRect(string meshName, Rect meshRect, std::function<void(string, Rect)> save_data_updater);

	//copy all primary mesh data (magnetisation, elC, Temp, etc.) but do not change dimensions or discretisation
	BError copy_mesh_data(string meshName_from, string meshName_to);

	//delete or add rectangle for given mesh
	BError delrect(string meshName, Rect rectangle);
	BError setrect(string meshName, Rect rectangle);
	BError resetrect(string meshName);

	//roughen mesh sides perpendicular to a named axis (axis = "x", "y", "z") to given depth (same units as h) with prng instantiated with given seed.
	BError RoughenMeshSides(string meshName, string axis, double depth, int seed);

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	BError RoughenMeshSurfaces_Jagged(string meshName, double depth, double spacing, int seed, string sides);

	//clear roughness: set fine shape to coarse shape
	BError ClearMeshRoughness(string meshName);

	//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains2D(string meshName, double spacing, int seed);

	//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains3D(string meshName, double spacing, int seed);

	//--------------------------------------------------------- MESH HANDLING - SETTINGS : SuperMeshMeshes_Settings.cpp

	//set magnetisation angle in given mesh (all if meshName not given)
	BError SetMagAngle(string meshName, double polar, double azim);

	//set magnetisation angle in given mesh (must be specified), and only in given rectangle of meshName (relative coordinates)
	BError SetMagAngle_Rect(string meshName, double polar, double azim, Rect rectangle);

	//Invert magnetisation direction in given mesh (must be magnetic)
	BError SetInvertedMag(string meshName, bool x = true, bool y = true, bool z = true);

	//Mirror magnetisation in given axis (literal x, y, or z) in given mesh (must be magnetic)
	BError SetMirroredMag(string meshName, string axis);

	//Set random magentisation distribution in given mesh (must be magnetic)
	BError SetRandomMag(string meshName);

	//longitudinal and transverse are the components specified as string literals : "-z", "-y", "-x", "x", "y", "z"
	BError SetMagDomainWall(string meshName, string longitudinal, string transverse, double width, double position);

	//Set Neel skyrmion in given mesh with chirality, diameter and relative x-y plane position
	BError SetSkyrmion(string meshName, int orientation, int chirality, double diameter, DBL2 position);
	
	//Set Bloch skyrmion in given mesh with chirality, diameter and relative x-y plane position
	BError SetSkyrmionBloch(string meshName, int orientation, int chirality, double diameter, DBL2 position);

	//Set applied magnetic field
	BError SetField(string meshName, DBL3 field_cartesian);

	//Set uniform applied mechanical stress
	BError SetUniformStress(string meshName, DBL3 stress_cartesian);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - transverse wall
	BError PrepareMovingMesh(string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - bloch wall
	BError PrepareMovingMesh_Bloch(string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - neel wall
	BError PrepareMovingMesh_Neel(string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - skyrmion
	BError PrepareMovingMesh_Skyrmion(string meshName);

	//Clear moving mesh settings made by a prepare method
	void ClearMovingMesh(void);

	//set ferromagnetic mesh roughness refinement if Roughness module enabled in given mesh
	BError SetMeshRoughnessRefinement(string meshName, INT3 refine);

	//Set periodic boundary conditions for magnetization
	//possible flags: x, y, z
	BError Set_PBC(string meshName, string flag, int images);

	//set exchange coupling to neighboring meshes - an exchange-type module (i.e. inherit from ExchangeBase) must be enabled in the named mesh
	BError Set_ExchangeCoupledMeshes(bool status, string meshName);

	//Set/Get multilayered demag exclusion : will need to call UpdateConfiguration when the flag is changed, so the correct SDemag_Demag modules and related settings are set from the SDemag module.
	BError Set_Demag_Exclusion(bool exclude_from_multiconvdemag, string meshName);

	//set link_stochastic flag in named mesh, or all meshes if supermesh handle given
	BError SetLinkStochastic(bool link_stochastic, string meshName);

	//--------------------------------------------------------- MESH PARAMETERS : SuperMeshParams.cpp

	//these set parameter values and temperature dependence in the indicated mesh - call through these since it's important to call UpdateConfiguration also
	BError set_meshparam_value(string meshName, string paramHandle, string value_text);
	
	//get named parameter value from given mesh. Set value as a string in value_text, without units
	BError get_meshparam_value(string meshName, string paramHandle, string& value_text);

	//temperature dependence

	BError set_meshparam_t_equation(string meshName, string paramHandle, string equationText);
	
	BError set_meshparam_tscaling_array(string meshName, string paramHandle, vector<double>& temp, vector<double>& scaling_x, vector<double>& scaling_y, vector<double>& scaling_z, string fileName_info = "");

	//clear parameters temperature dependence in given mesh (all meshes if empty string)
	BError clear_meshparam_temp(string meshName, string paramHandle);
	
	//spatial dependence

	BError set_meshparam_s_equation(string meshName, string paramHandle, string equationText);

	//clear parameters spatial dependence (variation) in given mesh (all meshes if empty string)
	BError clear_meshparam_variation(string meshName);

	//clear parameter spatial dependence (variation) in given mesh for named parameter only
	BError clear_meshparam_variation(string meshName, string paramHandle);

	//set parameter to display in given mesh when ParamVar spatial variation display is enabled
	BError set_meshparamvar_display(string meshName, string paramHandle);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
	BError set_meshparam_var(string meshName, string paramHandle, string generatorHandle, string generatorArgs, function<vector<unsigned char>(string, INT2)>& bitmap_loader);

	//copy all parameters from another Mesh
	BError copy_mesh_parameters(string meshName_from, string meshName_to);

	//--------------------------------------------------------- TEMPERATURE / HEAT SOLVER CONTROL : SuperMeshTemperature.cpp

	//set mesh base temperature. If spatial variation set and Heat module enabled then non-uniform base temperature will be set
	BError SetBaseTemperature(string meshName, double Temperature);

	//ambient and alpha boundary coefficient for Robin boundary conditions - set in Heat module if active
	BError SetAmbientTemperature(string meshName, double T_ambient);
	BError SetAlphaHeatBoundary(string meshName, double alpha_boundary);
	
	//insulating mesh sides for heat equation (Neumann boundary conditions). literal can be "x", "-x", "y", "-y", "z", "-z"
	BError SetInsulatingSides(string meshName, string literal, bool status);

	//set Curie temperature/atomic moment as Bohr magneton multiple for named mesh or all meshes (if meshName is the supermesh handle)
	//this is for the actually set Tc value
	//applicable for micromagnetic meshes only
	BError SetCurieTemperature(string meshName, double T_Curie);
	
	//this is for the indicative material Tc value
	//applicable for micromagnetic meshes only
	BError SetCurieTemperatureMaterial(string meshName, double T_Curie_material);
	BError SetAtomicMagneticMoment(string meshName, DBL2 atomic_moment);

	//set Tc (critical temperature) coupling terms for 2-sublattice model
	//applicable for micromagnetic meshes only
	BError SetTcCoupling(string meshName, DBL2 tau_ii, DBL2 tau_ij);
	BError SetTcCoupling_Intra(string meshName, DBL2 tau_ii);
	BError SetTcCoupling_Inter(string meshName, DBL2 tau_ij);

	//Set temperature model
	BError SetTemperatureModel(string meshName, int tmtype);

	//----------------------------------- MODULES CONTROL : SuperMeshModules.cpp

	//Add module to given mesh (delegate implementation to AddModule in the referenced mesh), checking for super-mesh module clashes
	BError AddModule(string meshName, MOD_ moduleId);

	//Delete module from given mesh (delegate implementation to DelModule in the referenced mesh), checking for super-mesh module clashes
	BError DelModule(string meshName, MOD_ moduleId);

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

	//get total volume energy density
	double GetTotalEnergyDensity(void);

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

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS : SuperMeshDisplay.cpp

	vector<PhysQ> FetchOnScreenPhysicalQuantity(double detail_level = 0.0);
	
	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType);

	//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
	void PrepareDisplayedMeshValue(void);

	//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
	//Return an Any as the displayed quantity could be either a scalar or a vector.
	Any GetDisplayedMeshValue(DBL3 abs_pos);

	//return average value for currently displayed mesh quantity in the given relative rectangle
	Any GetAverageDisplayedMeshValue(Rect rel_rect);

	int GetDisplayedPhysicalQuantity(void) { return displayedPhysicalQuantity; }

	BError SetDisplayedPhysicalQuantity(string meshName, int displayedPhysicalQuantity_);
	BError SetDisplayedBackgroundPhysicalQuantity(string meshName, int displayedBackgroundPhysicalQuantity_);

	//Get/Set vectorial quantity representation options in named mesh (which could be the supermesh)
	BError SetVEC3Rep(string meshName, int vec3rep_);
	int GetVEC3Rep(string meshName);

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