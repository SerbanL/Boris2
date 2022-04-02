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
	DBL3, DBL3, DBL3, double,
	//Material Parameters
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<DBL2, double>, MatP<DBL2, double>, double, TEquation<double>, double, MatP<double, double>, MatP<double, double>,
	MatP<double, double>, MatP<double, double>, MatP<double, double>, MatP<double, double>>,
	std::tuple<Transport, Heat> >
{

#if COMPILECUDA == 1
	friend DipoleMeshCUDA;
#endif

private:

	//dipole velocity (e.g. for implementing moving dipoles/track scanning etc.)
	DBL3 dipole_velocity = DBL3();
	//when a dipole is shifted due to non-zero velocity, if displacement is below the cliping distance (dipole_shift_clip), then add it to the shift_debt
	//when dipole_shift_debt exceeds clipping distance, then perform shift and reduce debt. This is to limit excessive shifting by small factors, since for non-zero velocity shifting is performed every iteration.
	DBL3 dipole_shift_debt = DBL3();
	DBL3 dipole_shift_clip = DBL3();
	//time at last attempt to shift dipole (i.e. previous iteration - keep a copy of it here, as we can't be sure what the ODE solver is doing, e.g. backtracking etc.)
	double dipole_last_time = 0.0;

	//dipoles are used to generate stray fields, to be used in ferromagnetic meshes. 
	//When a dipole value changes (direction or strength) then the stray field generated must be recalculated - this flag is checked by the StrayField module where the recalculation is done.
	bool recalculateStrayField = true;
	bool strayField_recalculated = false;

private:

	//at the start of each iteration see if we need to implement a moving dipole
	void Dipole_Shifting_Algorithm(void);

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
	void PrepareNewIteration(void);

#if COMPILECUDA == 1
	void PrepareNewIterationCUDA(void) { PrepareNewIteration(); }
#endif

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole (also setting recalculateStrayField flag)
	void Reset_Mdipole(void);

	void SetMagAngle(double polar, double azim, Rect rectangle = Rect());

	//also need to set magnetization temperature dependence in this mesh when a Curie temperature is set - don't need to store Curie temperature here, just set temperature dependence for Ms
	void SetCurieTemperature(double Tc, bool set_default_dependences);

	//shift dipole mesh rectangle by given amount
	void Shift_Dipole(DBL3 shift);

	//methods for dipole shifting algorithm
	void Set_Dipole_Velocity(DBL3 velocity, DBL3 clipping) { dipole_velocity = velocity; dipole_shift_clip = clipping; }
	DBL3 Get_Dipole_Velocity(void) { return dipole_velocity; }
	DBL3 Get_Dipole_Clipping(void) { return dipole_shift_clip; }

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool CheckRecalculateStrayField(void);
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

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole (also setting recalculateStrayField flag)
	void Reset_Mdipole(void) {}

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool CheckRecalculateStrayField(void) { return false; }
};

#endif
