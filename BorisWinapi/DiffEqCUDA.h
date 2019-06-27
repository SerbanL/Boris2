#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"
#include "ManagedDiffEqCUDA.h"

//DiffEq_IterateCUDA.cpp : coordinates call sequences and time steps for all evaluation methods (DONE)

//DiffEq_EvalsCUDA.cu : launchers and evaluation kernels for all evaluation types and equations (DONE) - also contains some auxiliary functions

//DiffEq_EquationsCUDA.h : All __device__ non-stochastic equation (DONE)
//DiffEq_SEquationsCUDA.h : All __device__ stochastic equation (DONE)
//DiffEq_SEquationsCUDA.cu : random field generation (DONE)

//NOTES : renormalize at last step for all equations except all LLB versions

class ODECommon;
class DifferentialEquation;
class FMesh;
class FMeshCUDA;

class ODECommonCUDA 
{
	friend ODECommon;

private:

	static ODECommon *pODE;

protected:

	//these need to be pointers, not cu_obj directly : we only want to make the cuda objects when ODECommonCUDA is made (cuda switched on), not at the start of the program - what if cuda not available on the system?
	static cu_obj<cuReal>* pdT;
	static cu_obj<cuReal>* pdT_last;

	static cu_obj<cuReal>* pmxh;
	static cu_obj<cuReal3>* pmxh_av;
	static cu_obj<size_t>* pavpoints;

	static cu_obj<cuReal>* plte;

	static cu_obj<bool>* prenormalize;

	static cu_obj<bool>* psolve_spin_current;

	static cu_obj<int>* psetODE;

	static cu_obj<bool>* palternator;

private:

	//---------------------------------------- SET-UP METHODS : DiffEqCUDA.cpp

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void Zero_reduction_values(void);
	void mxhav_to_mxh(void);

	//set all cuda values here from their cpu values held in ODECommon
	void SyncODEValues(void);

	//set specific cuda values (used often)
	void Sync_dT(void);
	void Sync_dT_last(void);
	void Sync_alternator(void);

public:

	ODECommonCUDA(void) {}
	ODECommonCUDA(ODECommon *pODE_);

	virtual ~ODECommonCUDA();

	//---------------------------------------- GET METHODS

	cuReal Get_mxh(void) { return pmxh->to_cpu(); }
	cuReal Get_lte(void) { return plte->to_cpu(); }

};

class DifferentialEquationCUDA :
	public ODECommonCUDA
{
	friend ODECommon;
	friend DifferentialEquation;

	friend ManagedDiffEqCUDA;

private:

	//ManagedDiffEqCUDA holds pointers to data in DifferentialEquationCUDA and ODECommonCUDA in an object in gpu memory.
	//pass cuDiffEq to a cuda kernel then all gpu data held here in cu_obj objects can be accessed in device code.
	//Initialize ManagedDiffEqCUDA with all the pointers you need then forget about it - no book-keeping required.
	cu_obj<ManagedDiffEqCUDA> cuDiffEq;

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//Used to save starting magnetization - all evaluation methods do this, even when not needed by the method itself, so we can calculate dM/dt when needed.
	cu_obj<cuVEC<cuReal3>> sM1;

	//Used for RK4 (0, 1, 2); ABM (0, 1)
	cu_obj<cuVEC<cuReal3>> sEval0, sEval1, sEval2;

	//Additional for use with RKF45
	cu_obj<cuVEC<cuReal3>> sEval3, sEval4;

	//Thermal field and torques, enabled only for the stochastic equations
	cu_obj<cuVEC<cuReal3>> H_Thermal, Torque_Thermal;

	//pseudo-random number generator for use in cuda kernels
	cu_obj<cuBorisRand> prng;

	//pointer to mesh with this ODE set
	FMesh *pMesh;

	//pointer to CUDA version of mesh with this ODE set
	FMeshCUDA *pMeshCUDA;

	//pointer to DifferentialEquation holding this CUDA object
	DifferentialEquation *pmeshODE;

private:

	//---------------------------------------- SET-UP METHODS  : DiffEqCUDA.cpp and DiffEqCUDA.cu

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	BError AllocateMemory(void);

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	void CleanupMemory(void);

	void SetODEMethodPointers(void);

	//---------------------------------------- EQUATIONS : these are defined as __device__ methods in ManagedDiffEqCUDA

	//---------------------------------------- SOLVER KERNEL LAUNCHERS generic : DiffEq_EvalsCUDA.cu

	//Euler evaluation of ODE
	void RunEuler(void);

	//Trapezoidal Euler evaluation of ODE
	void RunTEuler(int step);

	//RK4 evaluation of ODE
	void RunRK4(int step);

	//ABM
	void RunABM(int stepr);
	void RunABMTEuler(int step);

	//RKF45
	void RunRKF45(int step);

	//---------------------------------------- SOLVER KERNEL LAUNCHERS with in-lined LLG (faster) : DiffEq_EvalsLLGCUDA.cu

	//Euler evaluation of ODE
	void RunEuler_LLG(void);

	//Trapezoidal Euler evaluation of ODE
	void RunTEuler_LLG(int step);

	//RK4 evaluation of ODE
	void RunRK4_LLG(int step);

	//ABM
	void RunABM_LLG(int stepr);
	void RunABMTEuler_LLG(int step);

	//RKF45
	void RunRKF45_LLG(int step);

	//---------------------------------------- OTHERS : DiffEq_EvalsCUDA.cu

	//Restore magnetisation after a failed step for adaptive time-step methods
	void RestoreMagnetisation(void);

	//---------------------------------------- OTHER CALCULATION METHODS : DiffEq_SEquationsCUDA.cu

	//called when using stochastic equations
	void GenerateThermalField(void);
	void GenerateThermalField_and_Torque(void);

public:

	DifferentialEquationCUDA(DifferentialEquation *pmeshODE);

	~DifferentialEquationCUDA() { CleanupMemory(); }

	BError Error_On_Create(void) { return error_on_create; }

	//---------------------------------------- SET-UP METHODS : DiffEqCUDA.cpp

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	//----------------------------------- GETTERS

	//get reference to stored managed cuda differential equation object (cuDiffEq)
	cu_obj<ManagedDiffEqCUDA>& Get_ManagedDiffEqCUDA(void) { return cuDiffEq; }
};

#endif
