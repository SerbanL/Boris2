#pragma once

#include "Atom_Mesh.h"

#include "Atom_DiffEqCubic.h"

#include "SkyrmionTrack.h"
#include "DWRunTimeFit.h"

#include "Atom_Demag.h"
#include "Atom_Demag_N.h"
#include "Atom_DipoleDipole.h"
#include "Atom_Zeeman.h"
#include "Atom_Exchange.h"
#include "Atom_DMExchange.h"
#include "Atom_iDMExchange.h"
#include "Atom_viDMExchange.h"
#include "Atom_SurfExchange.h"
#include "Atom_MOptical.h"
#include "Atom_Anisotropy.h"
#include "Atom_AnisotropyCubi.h"
#include "Atom_AnisotropyBiaxial.h"
#include "Atom_AnisotropyTensorial.h"
#include "Atom_Transport.h"
#include "Atom_Heat.h"
#include "Atom_SOTField.h"
#include "Atom_STField.h"

#if COMPILECUDA == 1
#include "Atom_Mesh_CubicCUDA.h"
#endif

#ifdef MESH_COMPILATION_ATOM_CUBIC

/////////////////////////////////////////////////////////////////////
//
//Cubic Atomistic Mesh

class SuperMesh;

class Atom_Mesh_Cubic :
	public Atom_Mesh,
	public ProgramState<Atom_Mesh_Cubic,
	std::tuple<
	//Mesh members
	int, int, int, 
	int, int, int, int, int, int,
	Rect, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3,
	VEC_VC<DBL3>, 
	VEC_VC<double>, VEC_VC<DBL3>, VEC_VC<DBL3>, VEC_VC<double>,
	VEC_VC<double>, VEC_VC<double>,
	vector_lut<Modules*>, 
	bool,
	//Members in this derived class
	bool, SkyrmionTrack, bool,
	double, double, bool, bool, bool, DBL3,
	//Material Parameters
	MatP<double, double>, MatP<double, double>, MatP<DBL2, double>,
	MatP<double, double>, MatP<double, double>, MatP<DBL3, DBL3>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL3, DBL3>, MatP<DBL3, DBL3>, MatP<DBL3, DBL3>,
	std::vector<DBL4>,
	MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<DBL2, double>, MatP<DBL3, DBL3>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<DBL2, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>,
	double, TEquation<double>,
	MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>
	>,
	//Module Implementations
	std::tuple<
	Atom_Demag_N, Atom_Demag, Atom_DipoleDipole, 
	Atom_Zeeman, Atom_MOptical,
	Atom_Exchange, Atom_DMExchange, Atom_iDMExchange, Atom_viDMExchange, Atom_SurfExchange,
	Atom_Anisotropy_Uniaxial, Atom_Anisotropy_Cubic, Atom_Anisotropy_Biaxial, Atom_Anisotropy_Tensorial,
	Atom_Transport,
	Atom_Heat,
	Atom_SOTField, Atom_STField> >
{
#if COMPILECUDA == 1
	friend Atom_Mesh_CubicCUDA;
#endif

private:

	//The set ODE, associated with this cubic atomistic mesh (the ODE type and evaluation method is controlled from SuperMesh)
	Atom_DifferentialEquationCubic meshODE;

	//is this mesh used to trigger mesh movement? i.e. the CheckMoveMesh method should only be used if this flag is set
	bool move_mesh_trigger = false;

	//object used to track one or more skyrmions in this mesh
	SkyrmionTrack skyShift;

	//domain wall run-time position and width fitting
	DWPosWidth dwPos;

	//direct exchange coupling to neighboring meshes?
	//If true this is applicable for this mesh only for cells at contacts with other magnetic meshes 
	bool exchange_couple_to_meshes = false;

	//spin-wave factor for Tc estimation
	const double spinwave_factor = 0.723;
	
	//number of nearest neighbors
	const int coordination_number = 6;

	//number of atomic moments per cell
	const int atoms_per_cell = 1;

	// MONTE-CARLO DATA

	//vector of Monte-Carlo algorithm indices, used for shuffling spin picking order (used for serial MC algorithms)
	std::vector<unsigned> mc_indices;

	// Constrained MONTE-CARLO DATA

	//used for parallel constrained MC algorithm
	//Need to shuffle the mc indices every iteration (this is important, can't have the same spin pairings every step)
	//define as double, unsigned pair so we can use a sort-based shuffle algorithm, with sort from std::sort : generate random doubles, then sort them, with the mc indices moved in the same way, so result is shuffling.
	//std::sort can be executed as a parallel algorithm with C++17 stl, and random doubles also generated in parallel, so this is a fully parallel shuffle algorithm.
	//TO DO : Direct parallel shuffling is possible but a bit of a pain - probably best to use a bijective hash function to generate random permutations but need to look into this carefully. Probably not worth the effort for CPU code.
	std::vector<std::pair<double, unsigned>> mc_indices_red, mc_indices_black;

private:

	//Take a Monte Carlo step in this atomistic mesh : these functions implement the actual algorithms
	void Iterate_MonteCarlo_Serial_Classic(void);
	void Iterate_MonteCarlo_Serial_Constrained(void);
	void Iterate_MonteCarlo_Parallel_Classic(void);
	void Iterate_MonteCarlo_Parallel_Constrained(void);

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	Atom_Mesh_Cubic(SuperMesh *pSMesh_);

	Atom_Mesh_Cubic(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_);

	~Atom_Mesh_Cubic();

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError SwitchCUDAState(bool cudaState);

	//called at the start of each iteration
	void PrepareNewIteration(void) { if (!pMod.is_ID_set(MOD_ZEEMAN)) Heff1.set(DBL3(0)); }

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) { if (paMeshCUDA && !pMod.is_ID_set(MOD_ZEEMAN)) paMeshCUDA->Heff1()->set(n.dim(), cuReal3()); }
#endif

	//Take a Monte Carlo step in this atomistic mesh
	void Iterate_MonteCarlo(double acceptance_rate);

#if COMPILECUDA == 1
	//Take a Monte Carlo step in this atomistic mesh
	void Iterate_MonteCarloCUDA(double acceptance_rate);
#endif

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	double CheckMoveMesh(void);

	//couple this mesh to touching dipoles by setting skip cells as required : used for domain wall moving mesh algorithms
	void CoupleToDipoles(bool status);

	//----------------------------------- MOVING MESH TRIGGER FLAG

	bool GetMoveMeshTrigger(void) { return move_mesh_trigger; }
	void SetMoveMeshTrigger(bool status) { move_mesh_trigger = status; }

	//----------------------------------- ODE METHODS : Atom_Mesh_Cubic_ODEControl.cpp

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_dmdt(Rect avRect);

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_mxdmdt(Rect avRect);

	//----------------------------------- MESH QUANTITIES CONTROL : Atom_Mesh_Cubic_Control.cpp

	//this method is also used by the dipole mesh where it does something else - sets the dipole direction
	void SetMagAngle(double polar, double azim, Rect rectangle = Rect());

	//set magnetization angle only in given shape
	void SetMagAngle_Shape(double polar, double azim, std::vector<MeshShape> shapes);

	//Set magnetization angle in solid object only containing given relative position uniformly using polar coordinates
	void SetMagAngle_Object(double polar, double azim, DBL3 position);

	//Flower state magnetization
	void SetMagFlower(int direction, DBL3 centre, double radius, double thickness);

	//Onion state magnetization
	void SetMagOnion(int direction, DBL3 centre, double radius1, double radius2, double thickness);

	//Crosstie state magnetization
	void SetMagCrosstie(int direction, DBL3 centre, double radius, double thickness);

	//Invert magnetization direction in given mesh (must be magnetic)
	void SetInvertedMag(bool x, bool y, bool z);

	//Mirror magnetization in given axis (literal x, y, or z) in given mesh (must be magnetic)
	void SetMirroredMag(std::string axis);

	//Set random magentisation distribution in given mesh (must be magnetic)
	void SetRandomMag(int seed);
	void SetRandomXYMag(int seed);

	//set a domain wall with given width (metric units) at position within mesh (metric units). 
	//Longitudinal and transverse are magnetization componets as: 1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
	void SetMagDomainWall(int longitudinal, int transverse, double width, double position);

	//set Neel skyrmion with given orientation (core is up: 1, core is down: -1), chirality (1 for towards centre, -1 away from it) in given rectangle (relative to mesh), calculated in the x-y plane
	void SetSkyrmion(int orientation, int chirality, Rect skyrmion_rect);

	//set Bloch skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
	void SetSkyrmionBloch(int orientation, int chirality, Rect skyrmion_rect);

	//set M from given data VEC (0 values mean empty points) -> stretch data to M dimensions if needed.
	void SetMagFromData(VEC<DBL3>& data, const Rect& dstRect = Rect());

	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS

	double Get_NonEmpty_Magnetic_Volume(void) { return M1.get_nonempty_cells() * M1.h.dim(); }

	//----------------------------------- OVERLOAD MESH VIRTUAL METHODS

	//set/get exchange_couple_to_meshes status flag
	void SetMeshExchangeCoupling(bool status) { exchange_couple_to_meshes = status; }
	bool GetMeshExchangeCoupling(void) { return exchange_couple_to_meshes; }

	//Fit domain wall along the x direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
	DBL2 FitDomainWall_X(Rect rectangle);
	//Fit domain wall along the y direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
	DBL2 FitDomainWall_Y(Rect rectangle);
	//Fit domain wall along the z direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
	DBL2 FitDomainWall_Z(Rect rectangle);

	//compute magnitude histogram data
	//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first. 
	//output probabilities in histogram_p, corresponding to values set in histogram_x min, min + bin, ..., max, where bin = (max - min) / (num_bins - 1)
	//if macrocell_dims is not INT3(1) then first average in macrocells containing given number of individual mesh cells, then obtain histogram
	bool Get_Histogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims);

	//As for Get_Histogram, but use thermal averaging in each macrocell
	bool Get_ThAvHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims);

	//angular deviation histogram computed from ndir unit vector direction. If ndir not given (DBL3()), then angular deviation computed from average magnetization direction
	bool Get_AngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir);

	//As for Get_AngHistogram, but use thermal averaging in each macrocell
	bool Get_ThAvAngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir);

	//calculate thermodynamic average of magnetization
	DBL3 GetThermodynamicAverageMagnetization(Rect rectangle);

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in data box or output data, not with ShowData)
	//the rectangle must use relative coordinates
	DBL2 Get_skyshift(Rect skyRect)
	{
#if COMPILECUDA == 1
		if (paMeshCUDA) return skyShift.Get_skyshiftCUDA(n.dim(), h, meshRect, paMeshCUDA->M1, skyRect);
#endif
		return skyShift.Get_skyshift(M1, skyRect);
	}

	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in output data, not with ShowData or with data box), as well as diameters along x and y directions.
	DBL4 Get_skypos_diameters(Rect skyRect)
	{
#if COMPILECUDA == 1
		if (paMeshCUDA) return skyShift.Get_skypos_diametersCUDA(n.dim(), h, meshRect, paMeshCUDA->M1, skyRect);
#endif
		return skyShift.Get_skypos_diameters(M1, skyRect);
	}

	//----------------------------------- OTHER CONTROL METHODS : implement pure virtual Atom_Mesh methods

	//----------------------------------- OTHER CALCULATION METHODS : Atom_Mesh_Cubic_Compute.cpp

	//compute topological charge density spatial dependence and have it available to display in Cust_S
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void);

	//return phase transition temperature (K) based on formula Tc = J*e*z/3kB
	double Show_Transition_Temperature(void);

	//return saturation magnetization (A/m) based on formula Ms = mu_s*n/a^3
	double Show_Ms(void);

	//return exchange stiffness (J/m) based on formula A = J*n/2a
	double Show_A(void);

	//return uniaxial anisotropy constant (J/m^3) based on formula K = k*n/a^3
	double Show_Ku(void);

	//----------------------------------- GETTERS : Atom_Mesh_Cubic_GetData.cpp

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	double GetTopologicalCharge(Rect rectangle = Rect());

	//get average magnetic moment in given rectangle (entire mesh if none specified)
	DBL3 GetAverageMoment(Rect rectangle = Rect());

	//get moment magnitude min-max in given rectangle (entire mesh if none specified)
	DBL2 GetMomentMinMax(Rect rectangle = Rect());

	//get moment component min-max in given rectangle (entire mesh if none specified)
	DBL2 GetMomentXMinMax(Rect rectangle = Rect());
	DBL2 GetMomentYMinMax(Rect rectangle = Rect());
	DBL2 GetMomentZMinMax(Rect rectangle = Rect());

#if COMPILECUDA == 1
	//get reference to stored differential equation object (meshODE)
	Atom_DifferentialEquationCubic& Get_DifferentialEquation(void) { return meshODE; }
#endif
};

#else

class Atom_Mesh_Cubic :
	public Atom_Mesh
{

private:

	//The set ODE, associated with this cubic atomistic mesh (the ODE type and evaluation method is controlled from SuperMesh)
	Atom_DifferentialEquationCubic meshODE;

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	Atom_Mesh_Cubic(SuperMesh *pSMesh_) :
		Atom_Mesh(MESH_ATOM_CUBIC, pSMesh_),
		meshODE(this)
	{}

	Atom_Mesh_Cubic(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
		Atom_Mesh(MESH_ATOM_CUBIC, pSMesh_),
		meshODE(this)
	{}

	~Atom_Mesh_Cubic() {}

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError SwitchCUDAState(bool cudaState) { return BError(); }

	//called at the start of each iteration
	void PrepareNewIteration(void) {}

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) {}
#endif

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	double CheckMoveMesh(void) { return 0.0; }

	//couple this mesh to touching dipoles by setting skip cells as required : used for domain wall moving mesh algorithms
	void CoupleToDipoles(bool status) {}

	//----------------------------------- MOVING MESH TRIGGER FLAG

	bool GetMoveMeshTrigger(void) { return false; }
	void SetMoveMeshTrigger(bool status) {}


	//----------------------------------- MESH INFO AND SIZE GET/SET METHODS

	double Get_NonEmpty_Magnetic_Volume(void) { return 0.0; }

	//----------------------------------- OVERLOAD MESH VIRTUAL METHODS

	//set/get exchange_couple_to_meshes status flag
	void SetMeshExchangeCoupling(bool status) {}
	bool GetMeshExchangeCoupling(void) { return false; }

	//----------------------------------- OTHER CONTROL METHODS : implement pure virtual Atom_Mesh methods

	//----------------------------------- OTHER CALCULATION METHODS : Atom_Mesh_Cubic_Compute.cpp

	//compute topological charge density spatial dependence and have it available to display in Cust_S
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void) {}

	//----------------------------------- GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	double GetTopologicalCharge(Rect rectangle = Rect()) { return 0.0; }

	//get average magnetic moment in given rectangle (entire mesh if none specified)
	DBL3 GetAverageMoment(Rect rectangle = Rect()) { return DBL3(); }

	//get moment magnitude min-max in given rectangle (entire mesh if none specified)
	DBL2 GetMomentMinMax(Rect rectangle = Rect()) { return DBL2(); }

	//get moment component min-max in given rectangle (entire mesh if none specified)
	DBL2 GetMomentXMinMax(Rect rectangle = Rect()) { return DBL2(); }
	DBL2 GetMomentYMinMax(Rect rectangle = Rect()) { return DBL2(); }
	DBL2 GetMomentZMinMax(Rect rectangle = Rect()) { return DBL2(); }

	//----------------------------------- VALUE GETTERS

	//return phase transition temperature (K) based on formula Tc = J*e*z/3kB
	double Show_Transition_Temperature(void) { return 0.0; }

	//return saturation magnetization (A/m) based on formula Ms = mu_s*n/a^3
	double Show_Ms(void) { return 0.0; }

	//return exchange stiffness (J/m) based on formula A = J*n/2a
	double Show_A(void) { return 0.0; }

	//return uniaxial anisotropy constant (J/m^3) based on formula K = k*n/a^3
	double Show_Ku(void) { return 0.0; }

#if COMPILECUDA == 1
	//get reference to stored differential equation object (meshODE)
	Atom_DifferentialEquationCubic& Get_DifferentialEquation(void) { return meshODE; }
#endif
};

#endif