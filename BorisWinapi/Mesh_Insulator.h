#pragma once

#include "Mesh.h"

#if COMPILECUDA == 1
#include "Mesh_InsulatorCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////
//
//Electrical Insulator Material Mesh

class SuperMesh;

class InsulatorMesh :
	public Mesh,
	public ProgramState<InsulatorMesh,
	tuple<
	//Mesh members
	int, int, int, int, int, Rect, SZ3, DBL3, VEC_VC<double>, vector_lut<Modules*>,
	//Members in this derived class

	//Material Parameters
	double, MatP<double, double>, MatP<double, double>, MatP<double, double>
	>,
	//Module Implementations
	tuple<Heat> >
{

#if COMPILECUDA == 1
	friend InsulatorMeshCUDA;
#endif

private:

public:

private:

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	InsulatorMesh(SuperMesh *pSMesh_);

	InsulatorMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_);

	~InsulatorMesh() {}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when the mesh dimensions have changed - sets every quantity to the right dimensions
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError SwitchCUDAState(bool cudaState);

	//called at the start of each iteration
	void PrepareNewIteration(void) { }

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) {}
#endif

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	double CheckMoveMesh(void) { return 0.0; }
};