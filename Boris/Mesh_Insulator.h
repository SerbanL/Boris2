#pragma once

#include "Mesh.h"

#if COMPILECUDA == 1
#include "Mesh_InsulatorCUDA.h"
#endif

#include "Heat.h"
#include "TMR.h"

#ifdef MESH_COMPILATION_INSULATOR

/////////////////////////////////////////////////////////////////////
//
//Electrical Insulator Material Mesh

class SuperMesh;

class InsulatorMesh :
	public Mesh,
	public ProgramState<InsulatorMesh,
	std::tuple<
	//Mesh members
	int, int, int, 
	int, int, int, int, int, int,
	Rect, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3,
	VEC_VC<double>, VEC_VC<DBL3>, VEC_VC<DBL3>, VEC_VC<double>,
	VEC_VC<double>, 
	vector_lut<Modules*>,
	//Members in this derived class
	int,
	//Material Parameters
	double, TEquation<double>, 
	MatP<DBL2, double>, MatP<DBL2, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>
	>,
	//Module Implementations
	std::tuple<Heat, TMR> >
{

#if COMPILECUDA == 1
	friend InsulatorMeshCUDA;
#endif

private:

	//formula selection for TMR
	int TMR_type = TMR_SLONCZEWSKI;

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	InsulatorMesh(SuperMesh *pSMesh_);

	InsulatorMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_);

	~InsulatorMesh() {}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//----------------------------------- PROPERTIES

	void SetTMRType(TMR_ type) { TMR_type = type; }
	TMR_ GetTMRType(void) { return (TMR_)TMR_type; }

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError SwitchCUDAState(bool cudaState);

	//called at the start of each iteration
	void PrepareNewIteration(void) { }

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) {}
#endif
};

#else

class InsulatorMesh :
	public Mesh
{

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	InsulatorMesh(SuperMesh *pSMesh_) :
		Mesh(MESH_INSULATOR, pSMesh_)
	{}

	InsulatorMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
		Mesh(MESH_INSULATOR, pSMesh_)
	{}

	~InsulatorMesh() {}

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError SwitchCUDAState(bool cudaState) { return BError(); }

	//called at the start of each iteration
	void PrepareNewIteration(void) { }

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) {}
#endif
};

#endif
