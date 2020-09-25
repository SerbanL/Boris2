#pragma once

#include "Atom_Mesh.h"

#include "Atom_DiffEqCubic.h"

#include "Atom_Demag.h"
#include "Atom_Demag_N.h"
#include "Atom_DipoleDipole.h"
#include "Atom_Zeeman.h"
#include "Atom_Exchange.h"
#include "Atom_DMExchange.h"
#include "Atom_iDMExchange.h"
#include "Atom_MOptical.h"
#include "Atom_Anisotropy.h"
#include "Atom_AnisotropyCubi.h"
#include "Atom_Heat.h"

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
	tuple<
	//Mesh members
	int, int, int, 
	int, int, int, int, 
	Rect, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3,
	VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<double>,
	vector_lut<Modules*>, 
	bool,
	//Members in this derived class
	bool, bool,
	double, double, bool, bool, DBL3,
	//Material Parameters
	MatP<double, double>, MatP<double, double>, MatP<DBL2, double>,
	MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<DBL3, DBL3>, MatP<DBL3, DBL3>,
	MatP<double, double>, MatP<double, double>,
	MatP<double, double>,
	double, TEquation<double>,
	MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>
	>,
	//Module Implementations
	tuple<Atom_Demag_N, Atom_Demag, Atom_DipoleDipole, Atom_Zeeman, Atom_Exchange, Atom_DMExchange, Atom_iDMExchange, Atom_MOptical, Atom_Anisotropy_Uniaxial, Atom_Anisotropy_Cubic, Atom_Heat> >
{
#if COMPILECUDA == 1
	friend Atom_Mesh_CubicCUDA;
#endif

private:

	//The set ODE, associated with this cubic atomistic mesh (the ODE type and evaluation method is controlled from SuperMesh)
	Atom_DifferentialEquationCubic meshODE;

	//is this mesh used to trigger mesh movement? i.e. the CheckMoveMesh method should only be used if this flag is set
	bool move_mesh_trigger = false;

	//direct exchange coupling to neighboring meshes?
	//If true this is applicable for this mesh only for cells at contacts with other magnetic meshes 
	bool exchange_couple_to_meshes = false;

	//spin-wave factor for Tc estimation
	const double spinwave_factor = 0.71;
	
	//number of nearest neighbors
	const int coordination_number = 6;

	//number of atomic moments per cell
	const int atoms_per_cell = 1;

	// MONTE-CARLO DATA

	//vector of Monte-Carlo algorithm indices, used for shuffling spin picking order (used for serial MC algorithms)
	std::vector<unsigned> mc_indices;

	//used for parallel constrained MC algorithm, where we use the Fisher-Yates algorithm to shuffle spin indices so spins are paired randomly (this must surely be important, can't have the same spin pairings every step).
	std::vector<unsigned> mc_indices_red, mc_indices_black;

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

	//----------------------------------- MESH QUANTITIES CONTROL : Atom_Mesh_Cubic_Control.cpp

	//this method is also used by the dipole mesh where it does something else - sets the dipole direction
	void SetMagAngle(double polar, double azim, Rect rectangle = Rect());

	//Invert magnetisation direction in given mesh (must be magnetic)
	void SetInvertedMag(bool x, bool y, bool z);

	//Mirror magnetisation in given axis (literal x, y, or z) in given mesh (must be magnetic)
	void SetMirroredMag(string axis);

	//Set random magentisation distribution in given mesh (must be magnetic)
	void SetRandomMag(void);

	//set a domain wall with given width (metric units) at position within mesh (metric units). 
	//Longitudinal and transverse are magnetisation componets as: 1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
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

	//----------------------------------- OTHER CALCULATION METHODS : Atom_Mesh_Cubic_Compute.cpp

	//compute topological charge density spatial dependence and have it available to display in Cust_S
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void);

	//return phase transition temperature (K) based on formula Tc = J*e*z/3kB
	double Show_Transition_Temperature(void);

	//return saturation magnetisation (A/m) based on formula Ms = mu_s*n/a^3
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

	//Average square of components
	double GetAverageXMomentSq(Rect rectangle = Rect());
	double GetAverageYMomentSq(Rect rectangle = Rect());
	double GetAverageZMomentSq(Rect rectangle = Rect());

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

	//return saturation magnetisation (A/m) based on formula Ms = mu_s*n/a^3
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