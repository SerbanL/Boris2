#pragma once

#include "Mesh.h"

#include "SkyrmionTrack.h"

#if COMPILECUDA == 1
#include "Mesh_FerromagneticCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////
//
//Ferromagnetic Material Mesh

class SuperMesh;

class FMesh :
	public Mesh,
	public ProgramState<FMesh,
	//Mesh members
	tuple<int, int, int, int, int, Rect, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<double>, vector_lut<Modules*>,
	//Members in this derived class
	bool, SkyrmionTrack,
	//Material Parameters
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL3, DBL3>, MatP<DBL3, DBL3>, MatP<double, double>, MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<DBL2, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>,
	double, double, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>>,
	//Module Implementations
	tuple<Demag_N, Demag, SDemag_Demag, Exch_6ngbr_Neu, DMExchange, iDMExchange, SurfExchange, Zeeman, Anisotropy_Uniaxial, Anisotropy_Cubic, Transport, Heat, SOTField, Roughness> >
{
	//List modules which can only be used in this type of mesh (these modules will hold a pointer to this mesh directly, not the base Mesh class)
	friend Demag_N;
	friend Demag;
	friend SDemag_Demag;
	friend Exch_6ngbr_Neu;
	friend DMExchange;
	friend iDMExchange;
	friend SurfExchange;
	friend Zeeman;
	friend Anisotropy_Uniaxial;
	friend Anisotropy_Cubic;
	friend SOTField;
	friend Roughness;

	friend DifferentialEquation;

#if COMPILECUDA == 1
	friend FMeshCUDA;
#endif

private:

	//The set ODE, associated with this ferromagnetic mesh (the ODE type and evaluation method is controlled from SuperMesh)
	DifferentialEquation meshODE;

	//is this mesh used to trigger mesh movement? i.e. the CheckMoveMesh method should only be used if this flag is set
	bool move_mesh_trigger = false;

	//object used to track one or more skyrmions in this mesh
	SkyrmionTrack skyShift;

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	FMesh(SuperMesh *pSMesh_);

	FMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_);

	~FMesh();

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when the mesh dimensions have changed - sets every quantity to the right dimensions
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError SwitchCUDAState(bool cudaState);

	//called at the start of each iteration
	void PrepareNewIteration(void) { Heff.set(DBL3(0)); }

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) { if (pMeshCUDA) pMeshCUDA->Heff()->set(n.dim(), cuReal3()); }
#endif

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	double CheckMoveMesh(void);

	//couple this ferromagnetic mesh to any touching dipole meshes, setting interface cell values and flags
	void CoupleToDipoles(bool status);

	//----------------------------------- MOVING MESH TRIGGER FLAG

	bool GetMoveMeshTrigger(void) { return move_mesh_trigger; }
	void SetMoveMeshTrigger(bool status) { move_mesh_trigger = status; }

	//----------------------------------- FERROMAGNETIC MESH QUANTITIES CONTROL : Mesh_Ferromagnetic_Control.cpp

	//this method is also used by the dipole mesh where it does something else - sets the dipole direction
	void SetMagnetisationAngle(double polar, double azim, Rect rectangle = Rect());

	//Invert magnetisation direction in given mesh (must be ferromagnetic)
	void SetInvertedMagnetisation(void);

	//Set random magentisation distribution in given mesh (must be ferromagnetic)
	void SetRandomMagnetisation(void);

	//set a domain wall with given width (metric units) at position within mesh (metric units). 
	//Longitudinal and transverse are magnetisation componets as: 1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
	void SetMagnetisationDomainWall(int longitudinal, int transverse, double width, double position);

	//set Neel skyrmion with given orientation (core is up: 1, core is down: -1), chirality (1 for towards centre, -1 away from it) in given rectangle (relative to mesh), calculated in the x-y plane
	void SetSkyrmion(int orientation, int chirality, Rect skyrmion_rect);

	//set Bloch skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
	void SetSkyrmionBloch(int orientation, int chirality, Rect skyrmion_rect);

	//set M from given data VEC (0 values mean empty points) -> stretch data to M dimensions if needed.
	void SetMagnetisationFromData(VEC<DBL3>& data);

	//set periodic boundary conditions for magnetization
	void Set_PBC_X(int pbc_x);
	void Set_PBC_Y(int pbc_y);
	void Set_PBC_Z(int pbc_z);

	//----------------------------------- OVERLOAD MESH VIRTUAL METHODS

	//Curie temperature for ferromagnetic meshes. Calling this forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten.
	void SetCurieTemperature(double Tc);
	//this just sets the indicative material Tc value
	void SetCurieTemperatureMaterial(double Tc_material) { T_Curie_material = Tc_material; }

	//atomic moment (as multiple of Bohr magneton) for ferromagnetic meshes. Calling this forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten.
	void SetAtomicMoment(double atomic_moment_ub);
	double GetAtomicMoment(void) { return atomic_moment; }

	DBL3 dMdt(int idx) { return meshODE.dMdt(idx); }
	
	//get skyrmion shift for a skyrmion initially in the given rectangle (works only with data in data box or output data, not with ShowData)
	//the rectangle must use relative coordinates
	DBL2 Get_skyshift(Rect skyRect) 
	{ 
#if COMPILECUDA == 1
		if(pMeshCUDA) return skyShift.Get_skyshiftCUDA(n.dim(), h, pMeshCUDA->M, skyRect);
#endif
		return skyShift.Get_skyshift(M, skyRect); 
	}

	//----------------------------------- GETTERS

#if COMPILECUDA == 1
	//get reference to stored differential equation object (meshODE)
	DifferentialEquation& Get_DifferentialEquation(void) { return meshODE; }
#endif
};