#pragma once

#include <omp.h>

#include "ErrorHandler.h"

#include "PhysQRep.h"
#include "SimSharedData.h"

#include "MeshDefs.h"

#include "MeshParams.h"

#include "DiffEq.h"

#include "Exchange.h"
#include "DMExchange.h"
#include "iDMExchange.h"
#include "SurfExchange.h"
#include "Demag.h"
#include "Demag_N.h"
#include "SDemag_Demag.h"
#include "Zeeman.h"
#include "Anisotropy.h"
#include "AnisotropyCubi.h"
#include "Transport.h"
#include "Heat.h"
#include "SOTField.h"
#include "Roughness.h"

#if COMPILECUDA == 1
#include "MeshCUDA.h"
#endif

using namespace std;

class SuperMesh;

/////////////////////////////////////////////////////////////////////
//
//abstract base Mesh class - implement various types of meshes using this base

class Mesh : 
	public SimulationSharedData,
	public MeshParams
{
	//List modules which can be used by more than one type of mesh - these modules will hold a pointer to this base class directly
	friend Transport;
	friend Heat;

#if COMPILECUDA == 1
	friend MeshCUDA;

	friend TransportCUDA;
	friend HeatCUDA;
#endif

protected:

	int OmpThreads;

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//the type of mesh
	int meshType;

	//unique mesh id generated when the object is created (meshIdCounter simply increments every time a new mesh is created and sets a new meshId)
	static int meshIdCounter;

	//a unique mesh id number generated when the object is created, using meshIdCounter.
	int meshId = ++meshIdCounter;
	   
	//current quantity displayed on screen
	int displayedPhysicalQuantity;

	//if displaying a parameter spatial variation, this value will say which parameter (a value from PARAM_ enum)
	int displayedParamVar = PARAM_GREL;

	//effective field modules currently assigned to this mesh
	vector_lut<Modules*> pMod;

	//pointer to the supermesh holding this mesh (some modules in this mesh will need to see other meshes).
	SuperMesh* pSMesh;

public:

#if COMPILECUDA == 1
	//the CUDA version of this Mesh
	MeshCUDA* pMeshCUDA = nullptr;
#endif

	//the mesh rectangle in meters : this defines the mesh position in the world (in the super-mesh). cellsizes must divide the mesh rectangle, giving the number of cells in each dimension
	Rect meshRect;

	//-----Ferromagnetic properties

	//number of cells (n.x, n.y, n.z)
	SZ3 n = SZ3(1);

	//cellsizes (h.x, h.y, h.z)
	DBL3 h = DBL3(5e-9);

	//Magnetization using double floating point precision
	VEC_VC<DBL3> M;

	//effective field - sum total field of all the added modules
	VEC<DBL3> Heff;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//number of cells for electrical properties
	SZ3 n_e = SZ3(1);

	//cellsize for electrical properties
	DBL3 h_e = DBL3(5e-9);

	//electrical potential - on n_e, h_e mesh
	VEC_VC<double> V;
	
	//electrical conductivity - on n_e, h_e mesh
	VEC_VC<double> elC;
	
	//electrical current density - on n_e, h_e mesh
	VEC_VC<DBL3> Jc;

	//spin accumulation - on n_e, h_e mesh
	VEC_VC<DBL3> S;

	//-----Thermal conduction properties

	//number of cells for thermal properties
	SZ3 n_t = SZ3(1);

	//cellsize for thermal properties
	DBL3 h_t = DBL3(5e-9);

	//temperature calculated by Heat module
	VEC_VC<double> Temp;

	//-----Mechanical properties (TO DO)

	//number of cells for mechanical properties
	SZ3 n_m = SZ3(1);

	//cellsize for mechanical properties
	DBL3 h_m = DBL3(5e-9);

private:

	//When changing the mesh shape then some of the primary VEC_VC quantities held in Mesh need to shaped (NF_EMPTY flags set inside VEC_VC)
	//Which quantities are shaped and in which order is decided in this method. All public methods which change shape must define code in a lambda (run_this) then invoke this method with any required parameters for the lambda
	template  <typename Lambda, typename ... PType>
	BError change_mesh_shape(Lambda& run_this, PType& ... params);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (AUXILIARY) (MeshParamsControl.h)

	//UPDATER M CORSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_mcoarse_spatial(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	void update_parameters_mcoarse_spatial(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_mcoarse_full(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_mcoarse_full(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//UPDATER E CORSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_ecoarse_spatial(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	void update_parameters_ecoarse_spatial(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_ecoarse_full(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_ecoarse_full(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//UPDATER T CORSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_tcoarse_spatial(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	void update_parameters_tcoarse_spatial(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_tcoarse_full(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_tcoarse_full(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value);

public:
	
	Mesh(MESH_ meshType, SuperMesh *pSMesh_);

	virtual ~Mesh();

	//obtain error_on_create from this mesh, as well as any set modules - return first error found
	BError Error_On_Create(void);

	//----------------------------------- MODULES CONTROL : MeshModules.cpp

	//Add module to list of set modules, also deleting any exclusive modules to this one
	//If you set force_add = true then duplicate modules are allowed : in this case you must keep track of them
	BError AddModule(MOD_ moduleId, bool force_add = false);
	//delete module
	void DelModule(MOD_ moduleId);

	bool IsModuleSet(MOD_ moduleId) { return pMod.is_ID_set(moduleId); }
	Modules* GetModule(MOD_ moduleId) { if (IsModuleSet(moduleId)) return pMod(moduleId); else return nullptr; }

#if COMPILECUDA == 1
	ModulesCUDA* GetCUDAModule(MOD_ moduleId) { if (IsModuleSet(moduleId)) return pMod(moduleId)->GetCUDAModule(); else return nullptr; }
#endif

	//---- INITIALIZE MODULES -> Call Modules::Initialize

	//initialize all modules prior to computations starting
	BError InitializeAllModules(void);

#if COMPILECUDA == 1
	BError InitializeAllModulesCUDA(void);
#endif

	//---- UPDATE MODULES -> Call Modules::UpdateField

	//update computational state of all modules in this mesh
	void UpdateModules(void);

#if COMPILECUDA == 1
	void UpdateModulesCUDA(void);
#endif

	//----------------------------------- PARAMETERS CONTROL/INFO : MeshParamsControl.cpp

	//set/get mesh base temperature
	void SetBaseTemperature(double Temperature);
	double GetBaseTemperature(void) { return base_temperature; }

	//parameters spatial variation setters
	
	//update all mesh parameters spatial variation (if needed)
	bool update_meshparam_var(void);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
	BError set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader);

	//others	

	//copy all parameters from another Mesh
	BError copy_mesh_parameters(Mesh& copy_this);

	//each parameter has a certain type (PARAMTYPE_) - return cellsize in this mesh associated with this parameter type (e.g. ferromagnetic, electric, thermal, elastic)
	DBL3 get_paramtype_cellsize(PARAM_ paramID);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (MeshParamsControl.h)
	
	//SPATIAL DEPENDENCE ONLY - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_spatial(const DBL3& position, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_spatial(const DBL3& position, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_full(const DBL3& position, const double& Temperature, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_full(const DBL3& position, const double& Temperature, MatP<PType, SType>& matp, PType& matp_value);

	//UPDATER M CORSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_mcoarse(int mcell_idx, MeshParam_List& ... params);

	//UPDATER E CORSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_ecoarse(int ecell_idx, MeshParam_List& ... params);

	//UPDATER T CORSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_tcoarse(int tcell_idx, MeshParam_List& ... params);

	//UPDATER POSITION KNOWN - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_atposition(const DBL3& position, MeshParam_List& ... params);

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS : MeshControl.cpp

	//call when the mesh dimensions have changed - sets every quantity to the right dimensions
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) = 0;

	//at the start of each iteration the mesh might have to be prepared (e.g. state flags set)
	virtual void PrepareNewIteration(void) = 0;

#if COMPILECUDA == 1
	virtual void PrepareNewIterationCUDA(void) = 0;
#endif

	//used by move mesh algorithm : shift mesh quantities (e.g. magnetisation) by the given shift (metric units) value within this mesh. The shift is along the x-axis direction only (+ve or -ve).
	void MoveMesh(double x_shift);

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	virtual double CheckMoveMesh(void) = 0;

	//---- CUDA SWITCH

	//switch CUDA state on/off
	virtual BError SwitchCUDAState(bool cudaState) = 0;

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS : MeshDisplay.cpp

	//Return quantity currently set to display on screen.
	//Detail_level : the number of displayed elements in each mesh is obtained by dividing its rectangle by this cubic cell (to the nearest integer in each dimension).
	//detail_level is only needed here when CUDA is enabled so we can fetch data from gpu memory to a coarser array.
	PhysQ FetchOnScreenPhysicalQuantity(double detail_level = 0.0);

	int GetDisplayedPhysicalQuantity(void) { return displayedPhysicalQuantity; }

	int GetDisplayedParamVar(void) { return displayedParamVar; }

	void SetDisplayedPhysicalQuantity(int displayedPhysicalQuantity_) { displayedPhysicalQuantity = displayedPhysicalQuantity_; }

	void SetDisplayedParamVar(int displayedParamVar_) { displayedParamVar = displayedParamVar_; }

	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS : Mesh.cpp

	int get_id(void) { return meshId; }
	MESH_ GetMeshType(void) { return (MESH_)meshType; }

	//set new mesh rectangle and adjust cellsizes to default values if needed - the cellsizes can be manually set after changing the mesh rectangle
	BError SetMeshRect(Rect meshRect_);

	Rect GetMeshRect(void) { return meshRect; }
	DBL3 GetMeshDimensions(void) { return meshRect.size(); }
	DBL3 GetOrigin(void) { return meshRect.s; }

	//ferromagnetic properties
	SZ3 GetMeshSize(void) { return n; }
	//electrical conduction properties
	SZ3 GetMeshESize(void) { return n_e; }
	//thermal conduction properties
	SZ3 GetMeshTSize(void) { return n_t; }

	//ferromagnetic properties
	BError SetMeshCellsize(DBL3 h_);
	DBL3 GetMeshCellsize(void) { return h; }

	//electrical conduction properties
	BError SetMeshECellsize(DBL3 h_e_);
	DBL3 GetMeshECellsize(void) { return h_e; }

	//thermal conduction properties
	BError SetMeshTCellsize(DBL3 h_t_);
	DBL3 GetMeshTCellsize(void) { return h_t; }

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//magnetization dynamics computation enabled (check Heff not M - M is not empty for dipole meshes)
	bool MComputation_Enabled(void) { return Heff.linear_size(); }

	//slighly more general than MComputation_Enabled - this means the mesh can have magnetic properties but doesn't have an ODE set, e.g. dipole mesh
	bool Magnetisation_Enabled(void) { return M.linear_size(); }

	//electrical conduction computation enabled
	bool EComputation_Enabled(void) { return V.linear_size(); }

	//thermal conduction computation enabled
	bool TComputation_Enabled(void) { return Temp.linear_size(); }

	//check if interface conductance is enabled (for spin transport solver)
	bool GInterface_Enabled(void) { return (DBL2(Gi).norm() > 0); }

	//----------------------------------- VALUE GETTERS

	//get energy value for given module or one of its exclusive modules (if none active return 0)
	double GetEnergy(MOD_ moduleType);

	//get average magnetisation in given rectangle (entire mesh if none specified)
	DBL3 GetAverageMagnetisation(Rect rectangle = Rect());

	DBL3 GetAverageChargeCurrentDensity(Rect rectangle = Rect());

	DBL3 GetAverageSpinCurrentX(Rect rectangle = Rect());
	DBL3 GetAverageSpinCurrentY(Rect rectangle = Rect());
	DBL3 GetAverageSpinCurrentZ(Rect rectangle = Rect());

	double GetAverageElectricalPotential(Rect rectangle = Rect());
	DBL3 GetAverageSpinAccumulation(Rect rectangle = Rect());

	double GetAverageElectricalConductivity(Rect rectangle = Rect());

	//get base temperature or average temperature (if Temp enabled)
	double GetAverageTemperature(Rect rectangle = Rect());

	//get Curie temperature (the set value)
	double GetCurieTemperature(void) { return T_Curie; }

	//get Curie temperature for the material (the indicative value)
	double GetCurieTemperatureMaterial(void) { return T_Curie_material; }

	//----------------------------------- QUANTITY GETTERS

	//returns M on the cpu, thus transfers M from gpu to cpu before returning if cuda enabled
	VEC_VC<DBL3>& Get_M(void);

	//----------------------------------- OTHER MESH SHAPE CONTROL

	BError copy_mesh_data(Mesh& copy_this);

	//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
	BError applymask(double zDepth_m, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader);

	//set cells to empty in given rectangle (delete by setting entries to zero). The rectangle is relative to this mesh.
	BError delrect(Rect rectangle);

	//set cells to non-empty in given box
	BError setrect(Rect rectangle);

	//roughen mesh sides perpendicular to a named axis (axis = "x", "y", "z") to given depth (same units as h) with prng instantiated with given seed.
	BError RoughenMeshSides(string axis, double depth, unsigned seed);

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	BError RoughenMeshSurfaces_Jagged(double depth, double spacing, unsigned seed, string sides);

	//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains2D(double spacing, unsigned seed);

	//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains3D(double spacing, unsigned seed);

	//----------------------------------- METHODS REDEFINED IN SOME IMPLEMENTATIONS (virtual here)

	//use virtual to allow calling using base pointer - if called on "incorrect" mesh then nothing happens (the versions defined here used instead)

	virtual bool GetMoveMeshTrigger(void) { return false; }
	virtual void SetMoveMeshTrigger(bool status) {}

	//Curie temperature for ferromagnetic meshes. Calling this forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten.
	virtual void SetCurieTemperature(double Tc) {}
	//this just sets the indicative material Tc value
	virtual void SetCurieTemperatureMaterial(double Tc_material) {}

	//atomic moment (as multiple of Bohr magneton) for ferromagnetic meshes. Calling this forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten.
	virtual void SetAtomicMoment(double atomic_moment_ub) {}
	virtual double GetAtomicMoment(void) { return 0.0; }

	//get rate of change of magnetization (overloaded by Ferromagentic meshes)
	virtual DBL3 dMdt(int idx) { return DBL3(); }

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in output data, not with ShowData or with data box)
	virtual DBL2 Get_skyshift(Rect skyRect) { return DBL2(); }

	//----------------------------------- FERROMAGNETIC MESH QUANTITIES CONTROL : Mesh_Ferromagnetic_Control.cpp

	//this method is also used by the dipole mesh where it does something else - sets the dipole direction
	virtual void SetMagnetisationAngle(double polar, double azim, Rect rectangle = Rect()) {}

	//Invert magnetisation direction in given mesh (must be ferromagnetic)
	virtual void SetInvertedMagnetisation(void) {}

	//set a domain wall with given width (metric units) at position within mesh (metric units). 
	//Longitudinal and transverse are magnetisation componets as: 1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
	virtual void SetMagnetisationDomainWall(int longitudinal, int transverse, double width, double position) {}

	//set Neel skyrmion with given orientation (core is up: 1, core is down: -1), chirality (1 for towards centre, -1 away from it) in given rectangle (relative to mesh), calculated in the x-y plane
	virtual void SetSkyrmion(int orientation, int chirality, Rect skyrmion_rect) {}

	//set Bloch skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
	virtual void SetSkyrmionBloch(int orientation, int chirality, Rect skyrmion_rect) {}

	//set M from given data VEC (0 values mean empty points) -> stretch data to M dimensions if needed.
	virtual void SetMagnetisationFromData(VEC<DBL3>& data) {}

	//----------------------------------- MODULE METHODS TEMPLATED CALLERS

	//IMPORTANT NOTE: read note for CallModuleMethod in SuperMesh, repeated here:
	//IMPORTANT NOTE: These templated callers work by trying to cast the Modules* (Module is the abstract base class) to a derived implementation type (i.e. to Owner*) - dynamic_cast results in nullptr if couldn't cast.
	//				  Thus it is important that Owner is actually the implementation and not the base class Modules. If runThisMethod points to a method in Modules, e.g. GetEnergy, don't use the template deduction mechanism!
	//				  e.g. if you use CallModuleMethod(&STransport::GetEnergy), then Owner will not be STransport but Modules, so dynamic_cast will succeed on first attempt which may not be the right module!
	//				  In this case explicitly specify the template parameters as : CallModuleMethod<double, STransport>(&STransport::GetEnergy)

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
};

//!!!NOTES!!!
//When adding new VECs to the list here remember to also modify : 1) Mesh::copy_mesh_data and 2) MeshCUDA::copy_shapes_from_cpu
//Ideally there would be only one function where the VECs are explicitly listed to make maintainence and updates easier, but the 2 cases above are not so easy.
//It's possible for them to make use of change_mesh_shape using nested lambdas but you need a procedure to check if a VEC from one mesh is the same as the VEC from another (e.g. elC and copy_this.elC)
//When you have time you could do this, it will probably need an additional enum and a vector to check their types, don't think it's possible to do it without some sort of additional info (e.g. checking template type is not enough).
//e.g. see MeshParams::copy_parameters
template  <typename Lambda, typename ... PType>
BError Mesh::change_mesh_shape(Lambda& run_this, PType& ... params)
{
	//When changing the mesh shape then some of the primary VEC_VC quantities held in Mesh need to be shaped (NF_EMPTY flags set inside VEC_VC)
	//Which quantities are shaped and in which order is decided in this method. All public methods which change shape must define code in a lambda (run_this) then invoke this method with any required parameters for the lambda

	//Primary quantities are : M, elC, Temp

	BError error(__FUNCTION__);

	//1. shape magnetization
	if (M.linear_size()) {

		//if Roughness module is enabled then apply shape via the Roughness module instead
		if (IsModuleSet(MOD_ROUGHNESS)) {

			error = reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->change_mesh_shape(run_this, params...);
		}
		else error = run_this(M, DBL3(-Ms, 0, 0), params...);
	}
	
	//2. shape electrical conductivity
	if (elC.linear_size()) error = run_this(elC, elecCond, params...);
	
	//3. shape temperature
	if (Temp.linear_size()) error = run_this(Temp, base_temperature, params...);

#if COMPILECUDA == 1
	if (!error && pMeshCUDA) {

		error = pMeshCUDA->copy_shapes_from_cpu();
	}
#endif

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}
	
	return error;
}
