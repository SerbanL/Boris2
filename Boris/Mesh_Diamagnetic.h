#pragma once

#include "Mesh.h"

#if COMPILECUDA == 1
#include "Mesh_DiamagneticCUDA.h"
#endif

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "DiffEqDM.h"

#include "SurfExchange.h"
#include "Demag.h"
#include "SDemag_Demag.h"
#include "Zeeman.h"
#include "MOptical.h"
#include "Transport.h"
#include "Heat.h"

/////////////////////////////////////////////////////////////////////
//
//A diamagnetic mesh : a small negative temperature-independent susceptibility

class SuperMesh;

class DiaMesh :
	public Mesh,
	public ProgramState<DiaMesh,
	std::tuple<
	
	//Mesh members
	int, int, int, int, int, int, int, Rect, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<double>, VEC_VC<double>, VEC_VC<DBL3>, VEC_VC<DBL3>, VEC_VC<DBL3>, vector_lut<Modules*>, bool,
	//Members in this derived class
	//Material Parameters
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, 
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<DBL2, double>,
	double, TEquation<double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>
	>,
	//Module Implementations
	std::tuple<Demag, SDemag_Demag, SurfExchange, Zeeman, MOptical, Transport, Heat> >
{

#if COMPILECUDA == 1
	friend DiaMeshCUDA;
#endif

private:

	//The set ODE, associated with this diamagnetic mesh (the ODE type and evaluation method is controlled from SuperMesh)
	DifferentialEquationDM meshODE;

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	DiaMesh(SuperMesh *pSMesh_);

	DiaMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_);

	~DiaMesh();

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError SwitchCUDAState(bool cudaState);

	//called at the start of each iteration
	void PrepareNewIteration(void) { if (!pMod.is_ID_set(MOD_ZEEMAN)) Heff.set(DBL3(0)); }

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) { if (pMeshCUDA && !pMod.is_ID_set(MOD_ZEEMAN)) pMeshCUDA->Heff()->set(n.dim(), cuReal3()); }
#endif

	double CheckMoveMesh(void) { return 0.0; }
};

#else

class DiaMesh :
	public Mesh
{

private:

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	DiaMesh(SuperMesh *pSMesh_) :
		Mesh(MESH_DIAMAGNETIC, pSMesh_)
	{}

	DiaMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
		Mesh(MESH_DIAMAGNETIC, pSMesh_)
	{}

	~DiaMesh() {}

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

	double CheckMoveMesh(void) { return 0.0; }
};

#endif
