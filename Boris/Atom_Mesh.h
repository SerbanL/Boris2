#pragma once

#include <atomic>

#include "Boris_Enums_Defs.h"

#include "MeshBase.h"

#include "Atom_MeshParams.h"
#include "Atom_DiffEq.h"

#if COMPILECUDA == 1
#include "Atom_MeshCUDA.h"
#endif

class SuperMesh;

//Monte-Carlo Algorithm : minimum allowed cone angle
#define MONTECARLO_CONEANGLEDEG_MIN		1.0
//Monte-Carlo Algorithm : maximum allowed cone angle
#define MONTECARLO_CONEANGLEDEG_MAX		180.0
//Monte-Carlo Algorithm : change in cone angle per step
#define MONTECARLO_CONEANGLEDEG_DELTA	1.0
//Monte-Carlo Algorithm : target aceptance probability (vary cone angle to reach this)
#define MONTECARLO_TARGETACCEPTANCE		0.5

/////////////////////////////////////////////////////////////////////
//
//abstract base Atomic Mesh class - implement various types of atomic meshes using this base

class Atom_Mesh :
	public MeshBase,
	public Atom_MeshParams
{

#if COMPILECUDA == 1
	friend Atom_MeshCUDA;
#endif

protected:

	// MONTE-CARLO DATA

	//random number generator - used by Monte Carlo methods
	BorisRand prng;

	//Monte-Carlo current cone angle (vary to reach MONTECARLO_TARGETACCEPTANCE)
	double mc_cone_angledeg = 30.0;

	//last Monte-Carlo step acceptance probability (save it so we can read it out)
	double mc_acceptance_rate = 0.0;

	//use parallel Monte-Carlo algorithms?
	bool mc_parallel = true;

	// Constrained MONTE-CARLO DATA

	//use constrained Monte-Carlo?
	bool mc_constrain = false;

public:

#if COMPILECUDA == 1
	//the CUDA version of this Mesh
	Atom_MeshCUDA* paMeshCUDA = nullptr;
#endif

	//-----Demagnetizing field macrocell : compute demagnetizing field or macroscopic dipole-dipole interaction using this cellsize

	//number of cells (n.x, n.y, n.z)
	SZ3 n_dm = SZ3(1);

	//cellsize; doesn't have to be an integer number of atomistic cells
	DBL3 h_dm = DBL3(5e-9);

	//-----Magnetic properties

	//Atomic meshes which inherit this base class will have their defined crystal structure
	//The different types of crystal structures which can be simulated (cubic, bcc, fcc, hcp) can be though of, for computational purposes, as interleaved cubic (or rectangular in general) sub-lattices.
	//Thus cubic : one lattice
	//bcc : two sub-lattices
	//fcc : four sub-lattices
	//hcp : four rectangular sub-lattices
	//This determines the number of VECs used in each atomic mesh, labelled as 1, 2, 3, 4 for the required sub-lattice.

	//Atomic moments in units of Bohr magneton using double floating point precision (first sub-lattice : used in cubic, bcc, fcc, hcp)
	VEC_VC<DBL3> M1;

	//effective field (units of A/m : easier to integrate with micromagnetic meshes in a multiscale simulation this way) - sum total field of all the added modules (first sub-lattice : used in cubic, bcc, fcc, hcp)
	VEC<DBL3> Heff1;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//In Meshbase

	//-----Thermal conduction properties

	//In Meshbase

	//-----Mechanical properties

	//In Meshbase

private:

	//When changing the mesh shape then some of the primary VEC_VC quantities held in Mesh need to be shaped (NF_EMPTY flags set inside VEC_VC)
	//Which quantities are shaped and in which order is decided in this method. All public methods which change shape must define code in a lambda (run_this) then invoke this method with any required parameters for the lambda
	template  <typename Lambda, typename ... PType>
	BError change_mesh_shape(Lambda& run_this, PType& ... params);

	//calls pSMesh->UpdateConfiguration : use it in this file to avoid use of incomplete type pSMesh
	BError SuperMesh_UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (AUXILIARY) (Atom_MeshParamsControl.h)

	//SPATIAL DEPENDENCE ONLY - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_spatial(const DBL3& position, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_spatial(const DBL3& position, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_full(const DBL3& position, const double& Temperature, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_full(const DBL3& position, const double& Temperature, MatP<PType, SType>& matp, PType& matp_value);

	//UPDATER M COARSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_mcoarse_spatial(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	void update_parameters_mcoarse_spatial(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_mcoarse_full(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_mcoarse_full(int mcell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//UPDATER E COARSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_ecoarse_spatial(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	void update_parameters_ecoarse_spatial(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_ecoarse_full(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_ecoarse_full(int ecell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//UPDATER T COARSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_tcoarse_spatial(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	void update_parameters_tcoarse_spatial(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	void update_parameters_tcoarse_full(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	void update_parameters_tcoarse_full(int tcell_idx, MatP<PType, SType>& matp, PType& matp_value);

public:

	//------------------------CTOR/DTOR

	Atom_Mesh(MESH_ meshType, SuperMesh *pSMesh_);

	virtual ~Atom_Mesh();

	//----------------------------------- OTHER CONTROL METHODS : Atom_MeshControl.cpp

	//used by move mesh algorithm : shift mesh quantities (e.g. atomic moments) by the given shift (metric units) value within this mesh. The shift is along the x-axis direction only (+ve or -ve).
	void MoveMesh(double x_shift);

	//set PBC for required VECs : should only be called from a demag module
	BError Set_Magnetic_PBC(INT3 pbc_images);
	INT3 Get_Magnetic_PBC(void) { return INT3(M1.is_pbc_x(), M1.is_pbc_y(), M1.is_pbc_z()); }

	void Set_MonteCarlo_Serial(bool status) { mc_parallel = !status; }
	bool Get_MonteCarlo_Serial(void) { return !mc_parallel; }

	bool Get_MonteCarlo_Constrained(void) { return mc_constrain; }

	virtual void Set_MonteCarlo_Constrained(DBL3 cmc_n_) = 0;
	virtual DBL3 Get_MonteCarlo_Constrained_Direction(void) = 0;

	//----------------------------------- MODULES CONTROL (implement MeshBase) : Atom_MeshModules.cpp

	//Add module to list of set modules, also deleting any exclusive modules to this one
	//If you set force_add = true then duplicate modules are allowed : in this case you must keep track of them
	BError AddModule(MOD_ moduleID, bool force_add = false);

	//---- UPDATE MODULES -> Call Modules::UpdateField

	//update MOD_TRANSPORT module only if set
	void UpdateTransportSolver(void);

#if COMPILECUDA == 1
	void UpdateTransportSolverCUDA(void);
#endif

	//----------------------------------- PARAMETERS CONTROL/INFO : Atom_MeshParamsControl.cpp

	//set/get mesh base temperature; by default any text equation dependence will be cleared unless indicated specifically not to (e.g. called when setting base temperature value after evaluating the text equation)
	void SetBaseTemperature(double Temperature, bool clear_equation = true);

	//others	

	//copy all parameters from another Mesh
	BError copy_mesh_parameters(MeshBase& copy_this);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (Atom_MeshParamsControl.h)

	//UPDATER M COARSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_mcoarse(int mcell_idx, MeshParam_List& ... params);

	//UPDATER E COARSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_ecoarse(int ecell_idx, MeshParam_List& ... params);

	//UPDATER T COARSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_tcoarse(int tcell_idx, MeshParam_List& ... params);

	//UPDATER POSITION KNOWN - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	void update_parameters_atposition(const DBL3& position, MeshParam_List& ... params);

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS : Atom_MeshDisplay.cpp

	//Return quantity currently set to display on screen.
	//Detail_level : the number of displayed elements in each mesh is obtained by dividing its rectangle by this cubic cell (to the nearest integer in each dimension).
	//detail_level is only needed here when CUDA is enabled so we can fetch data from gpu memory to a coarser array.
	//getBackground : if true then use displayedBackgroundPhysicalQuantity instead of displayedPhysicalQuantity
	PhysQ FetchOnScreenPhysicalQuantity(double detail_level = 0.0, bool getBackground = false);

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType);

	//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
	void PrepareDisplayedMeshValue(void);

	//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
	//Return an Any as the displayed quantity could be either a scalar or a vector.
	Any GetDisplayedMeshValue(DBL3 abs_pos);

	//return average value for currently displayed mesh quantity in the given relative rectangle
	Any GetAverageDisplayedMeshValue(Rect rel_rect);

	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS : Atom_MeshDimensions.cpp

	//set new mesh rectangle and adjust cellsizes to default values if needed - the cellsizes can be manually set after changing the mesh rectangle
	//there's the option of adjusting the number of cells (adjust_num_cells) so set limits are not exceeded (useful when creating a new large mesh which might exceed memory availability with too small a cellsize)
	BError SetMeshRect(Rect meshRect_, bool adjust_num_cells);

	//ferromagnetic properties
	BError SetMeshCellsize(DBL3 h_);

	//electrical conduction properties
	BError SetMeshECellsize(DBL3 h_e_);

	//thermal conduction properties
	BError SetMeshTCellsize(DBL3 h_t_);

	//mechanical properties
	BError SetMeshMCellsize(DBL3 h_m_);

	//Set demagnetizing field evaluation macrocell size for atomistic meshes; may need to be adjusted if it doesn't result in an integer number of cells
	BError Set_Demag_Cellsize(DBL3 h_dm_);
	DBL3 Get_Demag_Cellsize(void) { return h_dm; }

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//magnetization dynamics computation enabled (check Beff not M)
	bool MComputation_Enabled(void) { return Heff1.linear_size(); }

	//slighly more general than MComputation_Enabled - this means the mesh can have magnetic properties but doesn't have an ODE set
	bool Magnetism_Enabled(void) { return M1.linear_size(); }

	//electrical conduction computation enabled
	bool EComputation_Enabled(void) { return V.linear_size(); }

	//thermal conduction computation enabled
	bool TComputation_Enabled(void) { return Temp.linear_size(); }

	//mechanical computation enabled
	bool MechComputation_Enabled(void) { return strain_diag.linear_size(); }

	//check if interface conductance is enabled (for spin transport solver)
	//TO DO
	bool GInterface_Enabled(void)
	{
		return false;
		//return (DBL2(Gmix.get0()).norm() > 0); 
	}

	//are periodic boundary conditions set for atomic moments?
	bool Is_PBC_x(void) { return M1.is_pbc_x(); }
	bool Is_PBC_y(void) { return M1.is_pbc_y(); }
	bool Is_PBC_z(void) { return M1.is_pbc_z(); }

	//is there a demag-type module set for this mesh? (SDemag not included as this is a SuperMesh module)
	bool Is_Demag_Enabled(void) { return IsModuleSet(MOD_DEMAG) || IsModuleSet(MOD_ATOM_DIPOLEDIPOLE); }

	//----------------------------------- VALUE GETTERS : Atom_MeshGetData.cpp

	//------Specific to Atom_Mesh

	//get average magnetization in given rectangle (entire mesh if none specified)
	virtual DBL3 GetAverageMoment(Rect rectangle = Rect()) = 0;
	
	//Average square of components
	virtual double GetAverageXMomentSq(Rect rectangle = Rect()) = 0;
	virtual double GetAverageYMomentSq(Rect rectangle = Rect()) = 0;
	virtual double GetAverageZMomentSq(Rect rectangle = Rect()) = 0;

	//get moment magnitude min-max in given rectangle (entire mesh if none specified)
	virtual DBL2 GetMomentMinMax(Rect rectangle = Rect()) = 0;

	//get moment component min-max in given rectangle (entire mesh if none specified)
	virtual DBL2 GetMomentXMinMax(Rect rectangle = Rect()) = 0;
	virtual DBL2 GetMomentYMinMax(Rect rectangle = Rect()) = 0;
	virtual DBL2 GetMomentZMinMax(Rect rectangle = Rect()) = 0;

	DBL2 Get_MonteCarlo_Params(void) { return DBL2(mc_cone_angledeg, mc_acceptance_rate); }

	//------Implementing MeshBase

	DBL3 GetAverageChargeCurrentDensity(Rect rectangle = Rect());

	DBL3 GetAverageSpinCurrentX(Rect rectangle = Rect());
	DBL3 GetAverageSpinCurrentY(Rect rectangle = Rect());
	DBL3 GetAverageSpinCurrentZ(Rect rectangle = Rect());

	double GetAverageElectricalPotential(Rect rectangle = Rect());
	DBL3 GetAverageSpinAccumulation(Rect rectangle = Rect());

	double GetAverageElectricalConductivity(Rect rectangle = Rect());

	//get base temperature or average temperature (if Temp enabled)
	double GetAverageTemperature(Rect rectangle = Rect());
	double GetAverageLatticeTemperature(Rect rectangle = Rect());

	//----------------------------------- QUANTITY GETTERS : Atom_MeshGetQuantities.cpp

	//returns M1 on the cpu, thus transfers M1 from gpu to cpu before returning if cuda enabled
	VEC_VC<DBL3>& Get_M1(void);

	//----------------------------------- VALUE GETTERS : Atom_MeshCompute.cpp

	//get maximum exchange energy density modulus over specified rectangle
	double Get_Max_Exchange_EnergyDensity(Rect& rectangle);

	//return phase transition temperature (K) based on formula Tc = J*e*z/3kB
	virtual double Show_Transition_Temperature(void) = 0;

	//return saturation magnetization (A/m) based on formula Ms = mu_s*n/a^3
	virtual double Show_Ms(void) = 0;

	//return exchange stiffness (J/m) based on formula A = J*n/2a
	virtual double Show_A(void) = 0;

	//return uniaxial anisotropy constant (J/m^3) based on formula K = k*n/a^3
	virtual double Show_Ku(void) = 0;

	//----------------------------------- OTHER CALCULATION METHODS : Atom_MeshCompute.cpp

	//compute exchange energy spatial variation and have it available to display in Cust_S
	void Compute_Exchange(void);

	//----------------------------------- OTHER MESH SHAPE CONTROL : Atom_MeshShape.cpp

	BError copy_mesh_data(MeshBase& copy_this);

	//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
	BError applymask(double zDepth_m, std::string fileName, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader);

	//set cells to empty in given rectangle (delete by setting entries to zero). The rectangle is relative to this mesh.
	BError delrect(Rect rectangle);

	//set cells to non-empty in given box
	BError setrect(Rect rectangle);

	//roughen mesh sides (side = "x", "y", "z", "-x", "-y", or "-z") to given depth (same units as h) with prng instantiated with given seed.
	BError RoughenMeshSides(std::string side, double depth, int seed);

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	BError RoughenMeshSurfaces_Jagged(double depth, double spacing, int seed, std::string sides);

	//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains2D(double spacing, int seed);

	//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
	BError GenerateGrains3D(double spacing, int seed);

	//Advanced mesh shape control methods

	//Disk with dimensions (x, y diameters, thickness), centre position relative to mesh, rotation angles, number of repetitions along x, y, z (1, 1, 1 for no repetitions), displacement values if repetitions used
	//Method: or (add shape) / not (delete shape) / xor (add and delete overlaps)
	BError shape_disk(MeshShape shape);

	//rectangle shape
	BError shape_rect(MeshShape shape);

	//triangle shape
	BError shape_triangle(MeshShape shape);

	//prolate ellipsoid
	BError shape_ellipsoid(MeshShape shape);

	//pyramid
	BError shape_pyramid(MeshShape shape);

	//tetrahedron
	BError shape_tetrahedron(MeshShape shape);

	//cone
	BError shape_cone(MeshShape shape);

	//torus
	BError shape_torus(MeshShape shape);

	//general shape setting function, can set composite shape using combination of the above elementary shapes
	BError shape_set(std::vector<MeshShape> shapes);

	//----------------------------------- METHODS REDEFINED IN SOME IMPLEMENTATIONS (virtual here - with exceptions)
};

//!!!NOTES!!!
//When adding new VECs to the list here remember to also modify : 1) Atom_Mesh::copy_mesh_data and 2) Atom_MeshCUDA::copy_shapes_from_cpu
//Ideally there would be only one function where the VECs are explicitly listed to make maintenance and updates easier, but the 2 cases above are not so easy.
//It's possible for them to make use of change_mesh_shape using nested lambdas but you need a procedure to check if a VEC from one mesh is the same as the VEC from another (e.g. elC and copy_this.elC)
//When you have time you could do this, it will probably need an additional enum and a vector to check their types, don't think it's possible to do it without some sort of additional info (e.g. checking template type is not enough).
//e.g. see Atom_MeshParams::copy_parameters
template  <typename Lambda, typename ... PType>
BError Atom_Mesh::change_mesh_shape(Lambda& run_this, PType& ... params)
{
	//When changing the mesh shape then some of the primary VEC_VC quantities held in Atom_Mesh need to be shaped (NF_EMPTY flags set inside VEC_VC)
	//Which quantities are shaped and in which order is decided in this method. All public methods which change shape must define code in a lambda (run_this) then invoke this method with any required parameters for the lambda

	//Primary quantities are : M1 (M2, M3, M4)

	BError error(__FUNCTION__);

#if COMPILECUDA == 1
	if (paMeshCUDA) {

		error = paMeshCUDA->copy_shapes_to_cpu();
	}
#endif

	//Atomic moments
	if (!shape_change_individual || (shape_change_individual && (displayedPhysicalQuantity == MESHDISPLAY_MOMENT)))
	{
		//1. shape moments
		if (M1.linear_size()) {

			error = run_this(M1, DBL3(-mu_s.get_current(), 0, 0), params...);
		}
	}

	//Electrical Conductivity
	if (!shape_change_individual || (shape_change_individual && displayedPhysicalQuantity == MESHDISPLAY_ELCOND)) {

		//2. shape electrical conductivity
		if (elC.linear_size()) error = run_this(elC, elecCond, params...);
	}

	//Temperature
	if (!shape_change_individual || (shape_change_individual && displayedPhysicalQuantity == MESHDISPLAY_TEMPERATURE)) {

		//3. shape temperature
		if (Temp.linear_size()) error = run_this(Temp, base_temperature, params...);
	}

	//Mechanical Displacement
	if (!shape_change_individual || (shape_change_individual && displayedPhysicalQuantity == MESHDISPLAY_UDISP)) {

		//4. shape mechanical displacement
		if (u_disp.linear_size()) error = run_this(u_disp, DBL3(), params...);
	}

#if COMPILECUDA == 1
	if (!error && paMeshCUDA) {

		error = paMeshCUDA->copy_shapes_from_cpu();
	}
#endif

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = SuperMesh_UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}
	
	return error;
}

