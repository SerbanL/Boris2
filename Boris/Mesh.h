#pragma once

#include "MeshBase.h"

#include "MeshParams.h"
#include "DiffEq.h"

class SuperMesh;

#include "Roughness.h"

/////////////////////////////////////////////////////////////////////
//
// abstract base Mesh class for micromagnetic meshes
// implement various types of micromagnetic meshes using this base

class Mesh : 
	public MeshBase,
	public MeshParams
{

#if COMPILECUDA == 1
	friend MeshCUDA;
#endif

private:

protected:

public:

#if COMPILECUDA == 1
	//the CUDA version of this Mesh
	MeshCUDA* pMeshCUDA = nullptr;
#endif

	//-----Ferromagnetic properties

	//Magnetization using double floating point precision
	VEC_VC<DBL3> M;

	//Additional magnetization used for antiferromagnetic meshes with 2 sub-lattice local approximation; exactly same dimensions as M
	VEC_VC<DBL3> M2;

	//effective field - sum total field of all the added modules
	VEC<DBL3> Heff;

	//Additional effective field used for antiferromagnetic meshes with 2 sub-lattice local approximation; exactly same dimensions as Heff
	VEC<DBL3> Heff2;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//In Meshbase

	//-----Thermal conduction properties

	//In Meshbase

	//-----Stochastic cellsize (VECs held in DiffEq)

	//number of cells for stochastic VECs
	SZ3 n_s = SZ3(1);

	//cellsize for stochastic VECs
	DBL3 h_s = DBL3(5e-9);

	//link stochastic cellsize to magnetic cellsize by default (set this to false if you want to control h_s independently)
	bool link_stochastic = true;

	//-----Mechanical properties

	//In Meshbase

private:

	//When changing the mesh shape then some of the primary VEC_VC quantities held in Mesh need to shaped (NF_EMPTY flags set inside VEC_VC)
	//Which quantities are shaped and in which order is decided in this method. All public methods which change shape must define code in a lambda (run_this) then invoke this method with any required parameters for the lambda
	template  <typename Lambda, typename ... PType>
	BError change_mesh_shape(Lambda& run_this, PType& ... params);

	//calls pSMesh->UpdateConfiguration : use it in this file to avoid use of incomplete type pSMesh
	BError SuperMesh_UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (AUXILIARY) (MeshParamsControl.h)

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

	Mesh(MESH_ meshType, SuperMesh *pSMesh_);

	virtual ~Mesh();

	//----------------------------------- OTHER CONTROL METHODS : MeshControl.cpp

	//used by move mesh algorithm : shift mesh quantities (e.g. magnetization) by the given shift (metric units) value within this mesh. The shift is along the x-axis direction only (+ve or -ve).
	void MoveMesh(double x_shift);

	//set PBC for required VECs : should only be called from a demag module
	BError Set_Magnetic_PBC(INT3 pbc_images);
	INT3 Get_Magnetic_PBC(void) { return INT3(M.is_pbc_x(), M.is_pbc_y(), M.is_pbc_z()); }

	//----------------------------------- MODULES CONTROL (implement MeshBase) : MeshModules.cpp

	//Add module to list of set modules, also deleting any exclusive modules to this one
	//If you set force_add = true then duplicate modules are allowed : in this case you must keep track of them
	BError AddModule(MOD_ moduleID, bool force_add = false);

	//---- UPDATE MODULES -> Call Modules::UpdateField

	//update MOD_TRANSPORT module only if set
	void UpdateTransportSolver(void);

#if COMPILECUDA == 1
	void UpdateTransportSolverCUDA(void);
#endif

	//----------------------------------- PARAMETERS CONTROL/INFO : MeshParamsControl.cpp

	//set/get mesh base temperature; by default any text equation dependence will be cleared unless indicated specifically not to (e.g. called when setting base temperature value after evaluating the text equation)
	void SetBaseTemperature(double Temperature, bool clear_equation = true);

	//others	

	//copy all parameters from another Mesh
	BError copy_mesh_parameters(MeshBase& copy_this);

	//this just sets the indicative material Tc value
	void SetCurieTemperatureMaterial(double Tc_material) { T_Curie_material = Tc_material; }

	//set tensorial anisotropy terms
	BError set_tensorial_anisotropy(std::vector<DBL4> Kt);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (MeshParamsControl.h)

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

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS : MeshDisplay.cpp

	//Return quantity currently set to display on screen.
	//Detail_level : the number of displayed elements in each mesh is obtained by dividing its rectangle by this cubic cell (to the nearest integer in each dimension).
	//detail_level is only needed here when CUDA is enabled so we can fetch data from gpu memory to a coarser array.
	//getBackground : if true then use displayedBackgroundPhysicalQuantity instead of displayedPhysicalQuantity
	PhysQ FetchOnScreenPhysicalQuantity(double detail_level = 0.0, bool getBackground = false);

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType);

	//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
	//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
	//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
	void GetPhysicalQuantityProfile(DBL3 start, DBL3 end, double step, DBL3 stencil, std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, bool do_average, bool read_average);

	//return average value for currently displayed mesh quantity in the given relative rectangle
	Any GetAverageDisplayedMeshValue(Rect rel_rect, std::vector<MeshShape> shapes = {});

	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS : MeshDimensions.cpp

	//set new mesh rectangle and adjust cellsizes to default values if needed - the cellsizes can be manually set after changing the mesh rectangle
	//there's the option of adjusting the number of cells (adjust_num_cells) so set limits are not exceeded (useful when creating a new large mesh which might exceed memory availability with too small a cellsize)
	BError SetMeshRect(Rect meshRect_, bool adjust_num_cells);

	//ferromagnetic properties
	BError SetMeshCellsize(DBL3 h_);

	//electrical conduction properties
	BError SetMeshECellsize(DBL3 h_e_);

	//thermal conduction properties
	BError SetMeshTCellsize(DBL3 h_t_);

	//stochastic properties
	BError SetLinkStochastic(bool link_stochastic_);
	bool GetLinkStochastic(void) { return link_stochastic; }
	
	BError SetMeshSCellsize(DBL3 h_s_, bool link_stochastic_ = true);
	DBL3 GetMeshSCellsize(void) { return h_s; }

	//mechanical properties
	BError SetMeshMCellsize(DBL3 h_m_);

	double Get_NonEmpty_Magnetic_Volume(void) { return M.get_nonempty_cells() * M.h.dim(); }

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//magnetization dynamics computation enabled (check Heff not M - M is not empty for dipole meshes)
	bool MComputation_Enabled(void) { return Heff.linear_size(); }

	//slighly more general than MComputation_Enabled - this means the mesh can have magnetic properties but doesn't necessarily have an ODE set, e.g. dipole mesh
	bool Magnetism_Enabled(void) { return M.linear_size(); }

	//electrical conduction computation enabled
	bool EComputation_Enabled(void) { return V.linear_size(); }

	//thermal conduction computation enabled
	bool TComputation_Enabled(void) { return Temp.linear_size(); }

	//mechanical computation enabled
	bool MechComputation_Enabled(void) { return strain_diag.linear_size(); }

	//check if interface conductance is enabled (for spin transport solver)
	bool GInterface_Enabled(void) { return (DBL2(Gmix.get0()).norm() > 0); }

	//are periodic boundary conditions set for magnetization?
	int Is_PBC_x(void) { return M.is_pbc_x(); }
	int Is_PBC_y(void) { return M.is_pbc_y(); }
	int Is_PBC_z(void) { return M.is_pbc_z(); }

	//is there a demag-type module set for this mesh? (SDemag not included as this is a SuperMesh module)
	bool Is_Demag_Enabled(void) { return IsModuleSet(MOD_DEMAG); }

	//----------------------------------- VALUE GETTERS : MeshGetData.cpp

	//------Specific to Mesh

	//get average magnetization in given rectangle (entire mesh if none specified)
	//1: ferromagnetic or antiferromagnetic sub-lattice A
	//2: antiferromagnetic only, sub-lattice B
	DBL3 GetAverageMagnetization(Rect rectangle = Rect());
	DBL3 GetAveragemagnetization2(Rect rectangle = Rect());

	//get magnetization magnitude min-max in given rectangle (entire mesh if none specified)
	DBL2 GetmagnetizationMinMax(Rect rectangle = Rect());

	//get magnetization component min-max in given rectangle (entire mesh if none specified)
	DBL2 GetmagnetizationXMinMax(Rect rectangle = Rect());
	DBL2 GetmagnetizationYMinMax(Rect rectangle = Rect());
	DBL2 GetmagnetizationZMinMax(Rect rectangle = Rect());

	//get Curie temperature (the set value)
	double GetCurieTemperature(void) { return T_Curie; }

	//get Curie temperature for the material (the indicative value)
	double GetCurieTemperatureMaterial(void) { return T_Curie_material; }

	double GetAtomicMoment(void) { return atomic_moment; }
	DBL2 GetAtomicMoment_AFM(void) { return atomic_moment_AFM; }

	DBL4 GetTcCoupling(void) { DBL2 tau_intra = tau_ii; DBL2 tau_inter = tau_ij; return DBL4(tau_intra.i, tau_intra.j, tau_inter.i, tau_inter.j); }

	//------Implementing MeshBase

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	double GetTopologicalCharge(Rect rectangle = Rect());

	double GetAverageElectricalPotential(Rect rectangle = Rect());
	DBL3 GetAverageSpinAccumulation(Rect rectangle = Rect());

	double GetAverageElectricalConductivity(Rect rectangle = Rect());

	//get base temperature or average temperature (if Temp enabled)
	double GetAverageTemperature(Rect rectangle = Rect());
	double GetAverageLatticeTemperature(Rect rectangle = Rect());

	//----------------------------------- QUANTITY GETTERS : MeshGetQuantities.cpp

	//returns M on the cpu, thus transfers M from gpu to cpu before returning if cuda enabled
	VEC_VC<DBL3>& Get_M(void);
	VEC_VC<DBL3>& Get_M2(void);

	//returns charge current on the cpu, assuming transport module is enabled
	VEC_VC<DBL3>& Get_Jc(void);

	//returns bulk self-consistent spin torque on the cpu, assuming transport module is enabled and spin solver is enabled
	VEC<DBL3>& Get_SpinTorque(void);

	//returns interfacial self-consistent spin torque on the cpu, assuming transport module is enabled and spin solver is enabled
	VEC<DBL3>& Get_InterfacialSpinTorque(void);

	//----------------------------------- VALUE GETTERS : MeshCompute.cpp

	//----------------------------------- OTHER CALCULATION METHODS : MeshCompute.cpp

	//compute topological charge density spatial dependence and have it available to display in Cust_S
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void);

	//----------------------------------- OTHER MESH SHAPE CONTROL : MeshShape.cpp

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

	//use virtual to allow calling using base pointer - if called on "incorrect" mesh then nothing happens (the versions defined here used instead)

	//Curie temperature for ferromagnetic meshes. Calling this forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten.
	virtual void SetCurieTemperature(double Tc, bool set_default_dependences) {}

	//atomic moment (as multiple of Bohr magneton) for ferromagnetic meshes. Calling this forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten.
	//Use DBL2 to allow 2-sublattice model as well. For single lattice just get the first component of the DBL2.
	virtual void SetAtomicMoment(DBL2 atomic_moment_ub) {}

	//set tau_ii and tau_ij values (applicable for 2-sublattice model only)
	virtual void SetTcCoupling(DBL2 tau_intra, DBL2 tau_inter) {}
	virtual void SetTcCoupling_Intra(DBL2 tau) {}
	virtual void SetTcCoupling_Inter(DBL2 tau) {}

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_..._ODEControl.cpp

	//get rate of change of magnetization (overloaded by Ferromagnetic meshes)
	virtual DBL3 dMdt(int idx) { return DBL3(); }
};

//!!!NOTES!!!
//When adding new VECs to the list here remember to also modify : 1) Mesh::copy_mesh_data and 2) MeshCUDA::copy_shapes_from_cpu
//Ideally there would be only one function where the VECs are explicitly listed to make maintenance and updates easier, but the 2 cases above are not so easy.
//It's possible for them to make use of change_mesh_shape using nested lambdas but you need a procedure to check if a VEC from one mesh is the same as the VEC from another (e.g. elC and copy_this.elC)
//When you have time you could do this, it will probably need an additional enum and a vector to check their types, don't think it's possible to do it without some sort of additional info (e.g. checking template type is not enough).
//e.g. see MeshParams::copy_parameters
template  <typename Lambda, typename ... PType>
BError Mesh::change_mesh_shape(Lambda& run_this, PType& ... params)
{
	//When changing the mesh shape then some of the primary VEC_VC quantities held in Mesh need to be shaped (NF_EMPTY flags set inside VEC_VC)
	//Which quantities are shaped and in which order is decided in this method. All public methods which change shape must define code in a lambda (run_this) then invoke this method with any required parameters for the lambda

	//Primary quantities are : M, elC, Temp, u_disp

	BError error(__FUNCTION__);

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		error = pMeshCUDA->copy_shapes_to_cpu();
	}
#endif

	//magnetization
	if (!shape_change_individual || 
		(shape_change_individual && (displayedPhysicalQuantity == MESHDISPLAY_MAGNETIZATION || displayedPhysicalQuantity == MESHDISPLAY_MAGNETIZATION2 || displayedPhysicalQuantity == MESHDISPLAY_MAGNETIZATION12)))
	{
		//1. shape magnetization
		if (M.linear_size()) {

			//if Roughness module is enabled then apply shape via the Roughness module instead
			if (IsModuleSet(MOD_ROUGHNESS)) {

				error = dynamic_cast<Roughness*>(pMod(MOD_ROUGHNESS))->change_mesh_shape(run_this, params...);
			}
			else {

				if (M2.linear_size()) {

					error = run_this(M, DBL3(-Ms_AFM.get_current().i, 0, 0), params...);
					error = run_this(M2, DBL3(Ms_AFM.get_current().j, 0, 0), params...);
				}
				else error = run_this(M, DBL3(-Ms.get_current(), 0, 0), params...);
			}
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
	if (!error && pMeshCUDA) {

		error = pMeshCUDA->copy_shapes_from_cpu();
	}
#endif

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = SuperMesh_UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}
	
	return error;
}

