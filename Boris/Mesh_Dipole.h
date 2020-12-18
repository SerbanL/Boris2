#pragma once

#include "Mesh.h"

#if COMPILECUDA == 1
#include "Mesh_DipoleCUDA.h"
#endif

#ifdef MESH_COMPILATION_DIPOLE

#include "Transport.h"
#include "Heat.h"

/////////////////////////////////////////////////////////////////////
//
//A simple magnetic dipole

class SuperMesh;

class DipoleMesh :
	public Mesh,
	public ProgramState<DipoleMesh,
	std::tuple<
	//Mesh members
	int, int, int, int, int, int, int, Rect, SZ3, DBL3, SZ3, DBL3, SZ3, DBL3, VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<DBL3>, VEC_VC<double>, VEC_VC<double>, vector_lut<Modules*>,
	//Members in this derived clas
	//Material Parameters
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<DBL2, double>, double, TEquation<double>, double, MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>>,
	std::tuple<Transport, Heat> >
{

#if COMPILECUDA == 1
	friend DipoleMeshCUDA;
#endif

private:

	//dipoles are used to generate stray fields, to be used in ferromagnetic meshes. 
	//When a dipole value changes (direction or strength) then the stray field generated must be recalculated - this flag is checked by the StrayField module where the recalculation is done.
	bool recalculateStrayField = true;

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	DipoleMesh(SuperMesh *pSMesh_);

	DipoleMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_);

	~DipoleMesh() {}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//----------------------------------- INITIALIZATION

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError SwitchCUDAState(bool cudaState);

	//called at the start of each iteration
	void PrepareNewIteration(void) {}

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) {}
#endif

	double CheckMoveMesh(void) { return 0.0; }

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole (also setting recalculateStrayField flag)
	void Reset_Mdipole(void);

	void SetMagAngle(double polar, double azim, Rect rectangle = Rect());

	void Reset_recalculateStrayField(void) { recalculateStrayField = false; }

	//also need to set magnetization temperature dependence in this mesh when a Curie temperature is set - don't need to store Curie temperature here, just set temperature dependence for Ms
	void SetCurieTemperature(double Tc, bool set_default_dependences);

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool Check_recalculateStrayField(void);
};

#else

class DipoleMesh :
	public Mesh
{

private:

public:

	//constructor taking only a SuperMesh pointer (SuperMesh is the owner) only needed for loading : all required values will be set by LoadObjectState method in ProgramState
	DipoleMesh(SuperMesh *pSMesh_) :
		Mesh(MESH_DIPOLE, pSMesh_)
	{}

	DipoleMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
		Mesh(MESH_DIPOLE, pSMesh_)
	{}

	~DipoleMesh() {}

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

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole (also setting recalculateStrayField flag)
	void Reset_Mdipole(void) {}

	void Reset_recalculateStrayField(void) {}

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool Check_recalculateStrayField(void) { return false; }
};

#endif
