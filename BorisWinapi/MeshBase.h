#pragma once

#include <omp.h>

#include "CompileFlags.h"
#include "ErrorHandler.h"

#include "MeshDefs.h"
#include "PhysQRep.h"
#include "SimSharedData.h"

#include "Modules.h"
#include "ParametersDefs.h"

#include "MeshParamsBase.h"

#if COMPILECUDA == 1
#include "MeshBaseCUDA.h"
#endif

using namespace std;

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

protected:

	//--------Auxiliary

	int OmpThreads;

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

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

	//--------Modules

	//effective field modules currently assigned to this mesh
	vector_lut<Modules*> pMod;

	//--------Configuration

	//if this mesh can participate in multilayered demag convolution, you have the option of excluding it (e.g. antiferromagnetic mesh) - by default all meshes with magnetic computation enabled are included.
	bool exclude_from_multiconvdemag = false;

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

private:

protected:

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

	//----------------------------------- OTHER CONTROL METHODS

	//used by move mesh algorithm : shift mesh quantities (e.g. magnetisation) by the given shift (metric units) value within this mesh. The shift is along the x-axis direction only (+ve or -ve).
	virtual void MoveMesh(double x_shift) = 0;

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	virtual double CheckMoveMesh(void) = 0;

	//set PBC for required VECs : should only be called from a demag module
	virtual BError Set_Magnetic_PBC(INT3 pbc_images) = 0;

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

	//----------------------------------- PARAMETERS CONTROL/INFO : MeshBaseParamsControl.cpp

	//set/get mesh base temperature; by default any text equation dependence will be cleared unless indicated specifically not to (e.g. called when setting base temperature value after evaluating the text equation)
	virtual void SetBaseTemperature(double Temperature, bool clear_equation = true) = 0;
	double GetBaseTemperature(void) { return base_temperature; }

	//set text equation for base temperature : when iterating, the base temperature will be evaluated and set using the text equation
	BError SetBaseTemperatureEquation(string equation_string, int step);
	void UpdateTEquationUserConstants(void);

	//parameters spatial variation setters

	//update all mesh parameters spatial variation (if needed)
	bool update_meshparam_var(void);

	//update text equations for mesh parameters with user constants, mesh dimensions, Curie temperature, base temperature
	bool update_all_meshparam_equations(void);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
	BError set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader);

	//others	

	//copy all parameters from another Mesh
	virtual BError copy_mesh_parameters(MeshBase& copy_this) = 0;

	//each parameter has a certain type (PARAMTYPE_) - return cellsize in this mesh associated with this parameter type (e.g. ferromagnetic, electric, thermal, elastic)
	DBL3 get_paramtype_cellsize(PARAM_ paramID);

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	//Return quantity currently set to display on screen.
	//Detail_level : the number of displayed elements in each mesh is obtained by dividing its rectangle by this cubic cell (to the nearest integer in each dimension).
	//detail_level is only needed here when CUDA is enabled so we can fetch data from gpu memory to a coarser array.
	//getBackground : if true then use displayedBackgroundPhysicalQuantity instead of displayedPhysicalQuantity
	virtual PhysQ FetchOnScreenPhysicalQuantity(double detail_level = 0.0, bool getBackground = false) = 0;

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	virtual BError SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType) = 0;

	//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
	virtual void PrepareDisplayedMeshValue(void) = 0;

	//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
	//Return an Any as the displayed quantity could be either a scalar or a vector.
	virtual Any GetDisplayedMeshValue(DBL3 abs_pos) = 0;

	//return average value for currently displayed mesh quantity in the given relative rectangle
	virtual Any GetAverageDisplayedMeshValue(Rect rel_rect) = 0;

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

	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS : MeshBase.cpp

	int get_id(void) { return meshId; }
	void swap_ids(MeshBase* pMesh_swap);

	MESH_ GetMeshType(void) { return (MESH_)meshType; }

	Rect GetMeshRect(void) { return meshRect; }
	DBL3 GetMeshDimensions(void) { return meshRect.size(); }
	DBL3 GetOrigin(void) { return meshRect.s; }

	//set new mesh rectangle and adjust cellsizes to default values if needed - the cellsizes can be manually set after changing the mesh rectangle
	virtual BError SetMeshRect(Rect meshRect_) = 0;

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
	virtual bool Is_PBC_x(void) = 0;
	virtual bool Is_PBC_y(void) = 0;
	virtual bool Is_PBC_z(void) = 0;

	//check if this mesh has an exchange module enabled (surf exchange doesn't count)
	bool ExchangeComputation_Enabled(void) { return IsModuleSet(MOD_EXCHANGE) || IsModuleSet(MOD_DMEXCHANGE) || IsModuleSet(MOD_IDMEXCHANGE); }

	//is there a demag-type module set for this mesh? (SDemag not included as this is a SuperMesh module)
	virtual bool Is_Demag_Enabled(void) = 0;

	//----------------------------------- VALUE GETTERS

	//get energy value for given module or one of its exclusive modules (if none active return 0); call it with MOD_ALL to return total energy density in this mesh; 
	double GetEnergyDensity(MOD_ moduleType);

	//As above, but get energy in specified rectangle only, where applicable
	double GetEnergyDensity(MOD_ moduleType, Rect& avRect);

	//get maximum exchange energy density modulus over specified rectangle
	virtual double Get_Max_Exchange_EnergyDensity(Rect& rectangle) = 0;

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	virtual double GetTopologicalCharge(Rect rectangle = Rect()) = 0;

	virtual DBL3 GetAverageChargeCurrentDensity(Rect rectangle = Rect()) = 0;

	virtual DBL3 GetAverageSpinCurrentX(Rect rectangle = Rect()) = 0;
	virtual DBL3 GetAverageSpinCurrentY(Rect rectangle = Rect()) = 0;
	virtual DBL3 GetAverageSpinCurrentZ(Rect rectangle = Rect()) = 0;

	virtual double GetAverageElectricalPotential(Rect rectangle = Rect()) = 0;
	virtual DBL3 GetAverageSpinAccumulation(Rect rectangle = Rect()) = 0;

	virtual double GetAverageElectricalConductivity(Rect rectangle = Rect()) = 0;

	//get base temperature or average temperature (if Temp enabled)
	virtual double GetAverageTemperature(Rect rectangle = Rect()) = 0;
	virtual double GetAverageLatticeTemperature(Rect rectangle = Rect()) = 0;

	//----------------------------------- OTHER CALCULATION METHODS

	//compute exchange energy spatial variation and have it available to display in Cust_S
	virtual void Compute_Exchange(void) = 0;

	//compute topological charge density spatial dependence and have it available to display in Cust_S
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	virtual void Compute_TopoChargeDensity(void) = 0;

	//----------------------------------- OTHER MESH SHAPE CONTROL

	virtual BError copy_mesh_data(MeshBase& copy_this) = 0;

	//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
	virtual BError applymask(double zDepth_m, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader) = 0;

	//set cells to empty in given rectangle (delete by setting entries to zero). The rectangle is relative to this mesh.
	virtual BError delrect(Rect rectangle) = 0;

	//set cells to non-empty in given box
	virtual BError setrect(Rect rectangle) = 0;

	//roughen mesh sides perpendicular to a named axis (axis = "x", "y", "z") to given depth (same units as h) with prng instantiated with given seed.
	virtual BError RoughenMeshSides(string axis, double depth, unsigned seed) = 0;

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	virtual BError RoughenMeshSurfaces_Jagged(double depth, double spacing, unsigned seed, string sides) = 0;

	//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	virtual BError GenerateGrains2D(double spacing, unsigned seed) = 0;

	//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	virtual BError GenerateGrains3D(double spacing, unsigned seed) = 0;

	//----------------------------------- METHODS REDEFINED IN SOME IMPLEMENTATIONS (virtual here)

	//use virtual to allow calling using base pointer - if called on "incorrect" mesh then nothing happens (the versions defined here used instead)

	virtual bool GetMoveMeshTrigger(void) { return false; }
	virtual void SetMoveMeshTrigger(bool status) {}

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in output data, not with ShowData or with data box)
	virtual DBL2 Get_skyshift(Rect skyRect) { return DBL2(); }

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in output data, not with ShowData or with data box), as well as diameters along x and y directions.
	virtual DBL4 Get_skypos_diameters(Rect skyRect) { return DBL4(); }

	//couple this mesh to touching dipoles by setting skip cells as required : used for domain wall moving mesh algorithms
	virtual void CoupleToDipoles(bool status) {}

	//----------------------------------- MESH QUANTITIES CONTROL

	//this method is also used by the dipole mesh where it does something else - sets the dipole direction
	virtual void SetMagAngle(double polar, double azim, Rect rectangle = Rect()) {}

	//Invert magnetisation direction in given mesh (must be magnetic)
	virtual void SetInvertedMag(bool x, bool y, bool z) {}

	//Mirror magnetisation in given axis (literal x, y, or z) in given mesh (must be magnetic)
	virtual void SetMirroredMag(string axis) {}

	//Set random magentisation distribution in given mesh (must be magnetic)
	virtual void SetRandomMag(void) {}

	//set a domain wall with given width (metric units) at position within mesh (metric units). 
	//Longitudinal and transverse are magnetisation componets as: 1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
	virtual void SetMagDomainWall(int longitudinal, int transverse, double width, double position) {}

	//set Neel skyrmion with given orientation (core is up: 1, core is down: -1), chirality (1 for towards centre, -1 away from it) in given rectangle (relative to mesh), calculated in the x-y plane
	virtual void SetSkyrmion(int orientation, int chirality, Rect skyrmion_rect) {}

	//set Bloch skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
	virtual void SetSkyrmionBloch(int orientation, int chirality, Rect skyrmion_rect) {}

	//set M from given data VEC (0 values mean empty points) -> stretch data to M dimensions if needed.
	virtual void SetMagFromData(VEC<DBL3>& data, const Rect& dstRect = Rect()) {}

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
