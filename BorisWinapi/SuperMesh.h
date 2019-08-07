#pragma once

#include "ErrorHandler.h"

#include "SimSharedData.h"

#include "Mesh_Ferromagnetic.h"
#include "Mesh_Dipole.h"
#include "Mesh_Metal.h"
#include "Mesh_Insulator.h"

#include "SDemag.h"
#include "StrayField.h"
#include "STransport.h"
#include "Oersted.h"
#include "SHeat.h"

#if COMPILECUDA == 1
#include "SuperMeshCUDA.h"
#endif

using namespace std;

//This is a container of simulation meshes.
//Coordinates the simulation flow between the different meshes.
//Performs calculations which cannot be done at individual mesh level (e.g. demag field over multiple ferromagnetic meshes).

class SuperMesh : 
	public SimulationSharedData, 
	public ProgramState<SuperMesh,
	tuple<int, SZ3, DBL3, Rect, SZ3, DBL3, Rect, ODECommon, vector_key<Mesh*>, vector_lut<Modules*>, string, string, bool, bool>,
	tuple<FMesh, DipoleMesh, MetalMesh, InsulatorMesh,
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

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//current quantity displayed on screen for super-mesh
	int displayedPhysicalQuantity = MESHDISPLAY_NONE;

	//the entire super-mesh rectangle (rectangle containing all meshes)
	Rect sMeshRect;

	//ferromagnetic super-mesh dimensions (this is the super-mesh that encompasses all ferromagnetic meshes)
	SZ3 n_fm;

	//ferromagnetic super-mesh discretization cellsize
	DBL3 h_fm = DBL3(5e-9);

	//ferromagnetic super-mesh rectangle
	Rect sMeshRect_fm;

	//electric super-mesh dimensions (this is the super-mesh that encompasses all meshes with electrical computations enabled)
	SZ3 n_e;

	//electric super-mesh discretization cellsize
	DBL3 h_e = DBL3(5e-9);

	//electric super-mesh rectangle
	Rect sMeshRect_e;

	//the micromagnetics equation to solve, together with its evaluation method is common to all meshes.
	//This is the "static" base class (i.e. contains a bunch of properties which are shared by all particular ODE solvers in the different ferromagnetic meshes, together with methods to manage them)
	//There are particular ODE objects in the ferromagnetic meshes, derived from ODECommon, since they differ in dimensions.
	ODECommon odeSolver;

	//individual simulation meshes (key is the mesh name)
	vector_key<Mesh*> pMesh;

	//super-mesh effective field modules (used for long range interactions)
	vector_lut<Modules*> pSMod;

	//name of currently active mesh in the pMesh vector (used for quicker data inputting in the console)
	string activeMeshName = "permalloy";

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

	SuperMesh(void);
	~SuperMesh();

	//obtain error_on_create from supermesh, as well as any meshes and super-mesh modules added - return first error found
	BError Error_On_Create(void);

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) { error_on_create = UpdateConfiguration(); }

	//--------------------------------------------------------- MESH INDEXING and related methods

	//index the super-mesh using mesh name to return a mesh pointer (assumes the referenced mesh is contained in the pMesh vector).
	Mesh* operator[](const string &key) { return pMesh[key]; }
	Mesh* operator[](const int &meshIdx) { return pMesh[meshIdx]; }

	//obtain mesh which contains the given coordinate (units m)
	Mesh* operator[](DBL3 coordinate)
	{
		for (int idx = 0; idx < pMesh.size(); idx++)
			if (pMesh[idx]->GetMeshRect().contains(coordinate)) return pMesh[idx];
		return nullptr;
	}

	Mesh* active_mesh(void) { return pMesh[activeMeshName]; }

	//return entire pMesh vector (used to access methods in vector_key on pMesh)
	vector_key<Mesh*>& operator()(void) { return pMesh; }

	//check if SMesh contains the given mesh : identify it by its key, or by its unique mesh id - for the latter return the index in pMesh (-1 if not found)
	bool contains(string &key) { return pMesh.has_key(key); }
	int contains_id(int meshId) { for (int idx = 0; idx < (int)pMesh.size(); idx++) if (pMesh[idx]->get_id() == meshId) return idx; return -1; }

	//from unique mesh id get key of entry in pMesh vector
	string key_from_meshId(int meshId) { for (int idx = 0; idx < (int)pMesh.size(); idx++) if (pMesh[idx]->get_id() == meshId) return pMesh.get_key_from_index(idx); return ""; }
	string key_from_meshIdx(int meshIdx) { return pMesh.get_key_from_index(meshIdx); }

	int size(void) { return (int)pMesh.size(); }

	//---------------------------------------------------------IMPORTANT CONTROL METHODS : SuperMeshControl.cpp

	//call whenever changes are made to meshes (e.g. dimensions, cellsize, modules added/deleted)
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	//couple ferromagnetic meshes to any touching dipole meshes, setting interface cell values and flags
	void CoupleToDipoles(void);

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

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

#if COMPILECUDA == 1
	void AdvanceTimeCUDA(void);
	void ComputeFieldsCUDA(void);
#endif

	//----------------------------------- ODE SOLVER CONTROL  : SuperMesh.cpp

	//Reset all ODE solvers in meshes with on ODE (e.g. ferromagnetic meshes with LLG)
	void ResetODE(void) { odeSolver.Reset(); }

	//set new stage for ODE solvers in the ferromagnetic meshes
	void NewStageODE(void) { odeSolver.NewStage(); }

	//set the ode and evaluation method. Any new ODE in a ferromagnetic mesh will use these settings.
	BError SetODE(ODE_ setOde, EVAL_ evalMethod);

	//set the time step for the magnetisation solver
	void SetTimeStep(double dT) { odeSolver.SetdT(dT); }

	//set parameters for adaptive time step control
	void SetAdaptiveTimeStepCtrl(double err_fail, double err_high, double err_low, double dT_incr, double dT_min, double dT_max) { odeSolver.SetAdaptiveTimeStepCtrl(err_fail, err_high, err_low, dT_incr, dT_min, dT_max); }

	//is the current time step fully finished? - most evaluation schemes need multiple sub-steps
	bool CurrentTimeStepSolved(void) { return odeSolver.TimeStepSolved(); }

	//check if ODE solver needs spin accumulation solved
	bool SolveSpinCurrent(void) { return odeSolver.SolveSpinCurrent(); }

	void SetMoveMeshTrigger(bool status, string meshName = "");

	void SetMoveMeshAntisymmetric(bool antisymmetric) { odeSolver.SetMoveMeshAntisymmetric(antisymmetric); }
	void SetMoveMeshThreshold(double threshold) { odeSolver.SetMoveMeshThreshold(threshold); }

	//---ODE Getters

	//get set ODE and evaluation method
	void QueryODE(ODE_ &setODE, EVAL_ &evalMethod) { odeSolver.QueryODE(setODE, evalMethod); }
	void QueryODE(ODE_ &setODE) { odeSolver.QueryODE(setODE); }

	int GetIteration(void) { return odeSolver.GetIteration(); }
	int GetStageIteration(void) { return odeSolver.GetStageIteration(); }

	double GetTime(void) { return odeSolver.GetTime(); }
	double GetStageTime(void) { return odeSolver.GetStageTime(); }
	double GetTimeStep(void) { return odeSolver.GetTimeStep(); }
	double Get_mxh(void) { return odeSolver.Get_mxh(); }
	double Get_dmdt(void) { return odeSolver.Get_dmdt(); }

	DBL3 Get_AStepRelErrCtrl(void) { return odeSolver.Get_AStepRelErrCtrl(); }
	DBL3 Get_AStepdTCtrl(void) { return odeSolver.Get_AStepdTCtrl(); }

	bool IsMovingMeshSet(void) { return odeSolver.IsMovingMeshSet(); }
	int GetId_of_MoveMeshTrigger(void) { return odeSolver.GetId_of_MoveMeshTrigger(); }
	double Get_dwshift(void) { return odeSolver.Get_dwshift(); }

	bool MoveMeshAntisymmetric(void) { return odeSolver.MoveMeshAntisymmetric(); }
	double MoveMeshThreshold(void) { return odeSolver.MoveMeshThreshold(); }

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

	//set mesh rect for named mesh (any if applicable any other dependent meshes) and update dependent save data rects by calling the provided function.
	BError SetMeshRect(string meshName, Rect meshRect, std::function<void(string, Rect)> save_data_updater);

	//--------------------------------------------------------- MESH HANDLING - SHAPES : SuperMeshMeshes.cpp

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

	//--------------------------------------------------------- MESH HANDLING - SETTINGS : SuperMeshMeshes.cpp

	BError SetMagnetisationAngle(string meshName, double polar, double azim);
	BError SetMagnetisationAngle_Rect(string meshName, double polar, double azim, Rect rectangle);

	//Invert magnetisation direction in given mesh (must be ferromagnetic)
	BError SetInvertedMagnetisation(string meshName);

	//longitudinal and transverse are the components specified as string literals : "-z", "-y", "-x", "x", "y", "z"
	BError SetMagnetisationDomainWall(string meshName, string longitudinal, string transverse, double width, double position);

	//Set Neel skyrmion in given mesh with chirality, diameter and relative x-y plane position
	BError SetSkyrmion(string meshName, int orientation, int chirality, double diameter, DBL2 position);
	
	//Set Bloch skyrmion in given mesh with chirality, diameter and relative x-y plane position
	BError SetSkyrmionBloch(string meshName, int orientation, int chirality, double diameter, DBL2 position);

	BError SetField(string meshName, DBL3 field_cartesian);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - transverse wall
	BError PrepareMovingMesh(string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - bloch wall
	BError PrepareMovingMesh_Bloch(string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - neel wall
	BError PrepareMovingMesh_Neel(string meshName);

	//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - skyrmion
	BError PrepareMovingMesh_Skyrmion(string meshName);

	//set ferromagnetic mesh roughness refinement if Roughness module enabled in given mesh
	BError SetMeshRoughnessRefinement(string meshName, INT3 refine);

	//--------------------------------------------------------- MESH PARAMETERS : SuperMeshParams.cpp

	//these set parameter values and temperature dependence in the indicated mesh - call through these since it's important to call UpdateConfiguration also
	BError set_meshparam_value(string meshName, string paramHandle, string value_text);
	
	//get named parameter value from given mesh. Set value as a string in value_text, without units
	BError get_meshparam_value(string meshName, string paramHandle, string& value_text);

	//temperature dependence

	BError set_meshparam_formula(string meshName, string paramHandle, string formulaName, vector<double> coefficients);
	
	BError set_meshparam_tscaling_array(string meshName, string paramHandle, vector<double>& temp, vector<double>& scaling);

	//clear parameters temperature dependence in given mesh (all meshes if empty string)
	BError clear_meshparam_temp(string meshName);
	
	//spatial dependence

	//clear parameters spatial dependence (variation) in given mesh (all meshes if empty string)
	BError clear_meshparam_variation(string meshName);

	//clear parameter spatial dependence (variation) in given mesh for named parameter only
	BError clear_meshparam_variation(string meshName, string paramHandle);

	//set parameter to display in given mesh when ParamVar spatial variation display is enabled
	BError set_meshparamvar_display(string meshName, string paramHandle);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
	BError set_meshparam_var(string meshName, string paramHandle, string generatorHandle, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader);

	//others

	BError SetBaseTemperature(string meshName, double Temperature);

	//ambient and alpha boundary coefficient for Robin boundary conditions - set in Heat module if active
	BError SetAmbientTemperature(string meshName, double T_ambient);
	BError SetAlphaHeatBoundary(string meshName, double alpha_boundary);
	//insulating mesh sides for heat equation (Neumann boundary conditions). literal can be "x", "-x", "y", "-y", "z", "-z"
	BError SetInsulatingSides(string meshName, string literal, bool status);

	//set Curie temperature/atomic moment as Bohr magneton multiple for named mesh or all meshes (if meshName is the supermesh handle)
	//this is for the actually set Tc value
	BError SetCurieTemperature(string meshName, double T_Curie);
	//this is for the indicative material Tc value
	BError SetCurieTemperatureMaterial(string meshName, double T_Curie_material);
	BError SetAtomicMagneticMoment(string meshName, double atomic_moment);

	//copy all parameters from another Mesh
	BError copy_mesh_parameters(string meshName_from, string meshName_to);

	//----------------------------------- MODULES CONTROL : SuperMeshModules.cpp

	//Add module to given mesh (delegate implementation to AddModule in the referenced mesh), checking for super-mesh module clashes
	BError AddModule(string meshName, MOD_ moduleId);

	//Delete module from given mesh (delegate implementation to DelModule in the referenced mesh), checking for super-mesh module clashes
	BError DelModule(string meshName, MOD_ moduleId);

	bool IsSuperMeshModuleSet(MOD_ moduleId) { return pSMod.is_ID_set(moduleId); }

	//--------------------------------------------------------- GET PROPERTIES / VALUES

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
	double GetTotalEnergy(void);

	//--------------------------------------------------------- GET MODULE SPECIFIC PROPERTIES

	//status for multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
	int Get_Multilayered_Convolution_Status(void)
	{
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return CallModuleMethod(&SDemag::Get_Multilayered_Convolution_Status);
		}
		else return -1;
	}

	//status for force 2d multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
	int Get_2D_Multilayered_Convolution_Status(void)
	{
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return CallModuleMethod(&SDemag::Get_2D_Multilayered_Convolution_Status);
		}
		else return -1;
	}

	//status for use default n_common for multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
	int Use_Default_n_Status(void)
	{
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return CallModuleMethod(&SDemag::Use_Default_n_Status);
		}
		else return -1;
	}

	//get n_common for multilayerd convolution (if return = SZ3() : N/A)
	SZ3 Get_n_common(void)
	{
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return CallModuleMethod(&SDemag::Get_n_common);
		}
		else return SZ3();
	}

	//--------------------------------------------------------- VALUE SETTERS : SuperMesh.cpp

	BError SetFMSMeshCellsize(DBL3 h_fm_);
	BError SetESMeshCellsize(DBL3 h_e_);

	void Set_Scale_Rects(bool status) { scale_rects = status; }

	void Set_Coupled_To_Dipoles(bool status) { coupled_dipoles = status; CoupleToDipoles(); }

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS : SuperMeshDisplay.cpp

	vector<PhysQ> FetchOnScreenPhysicalQuantity(double detail_level = 0.0);

	int GetDisplayedPhysicalQuantity(void) { return displayedPhysicalQuantity; }

	BError SetDisplayedPhysicalQuantity(string meshName, int displayedPhysicalQuantity_);

	//--------------------------------------------------------- MODULE METHODS TEMPLATED CALLERS

	//IMPORTANT NOTE: These templated callers work by trying to cast the Modules* (Module is the abstract base class) to a derived implementation type (i.e. to Owner*) - dynamic_cast results in nullptr if couldn't cast.
	//				  Thus it is important that Owner is actually the implementation and not the base class Modules. If runThisMethod points to a method in Modules, e.g. GetEnergy, don't use the template deduction mechanism!
	//				  e.g. if you use CallModuleMethod(&STransport::GetEnergy), then Owner will not be STransport but Modules, so dynamic_cast will succeed on first attempt which may not be the right module!
	//				  In this case explicitly specify the template parameters as : CallModuleMethod<double, STransport>(&STransport::GetEnergy)

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