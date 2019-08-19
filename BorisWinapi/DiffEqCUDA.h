#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"
#include "ManagedDiffEqCUDA.h"

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

	//-----------------------------------Time step

	//these need to be pointers, not cu_obj directly : we only want to make the cuda objects when ODECommonCUDA is made (cuda switched on), not at the start of the program - what if cuda not available on the system?
	static cu_obj<cuReal>* pdT;
	static cu_obj<cuReal>* pdT_last;

	//-----------------------------------Primary data

	static cu_obj<cuReal>* pmxh;
	static cu_obj<cuReal3>* pmxh_av;
	static cu_obj<size_t>* pavpoints;

	static cu_obj<cuReal>* pdmdt;
	static cu_obj<cuReal3>* pdmdt_av;
	static cu_obj<size_t>* pavpoints2;

	static cu_obj<cuReal>* plte;

	//-----------------------------------Evaluation method modifiers

	static cu_obj<bool>* prenormalize;

	//-----------------------------------Properties flags

	static cu_obj<bool>* psolve_spin_current;

	//-----------------------------------Equation and Evaluation method values

	static cu_obj<int>* psetODE;

	//-----------------------------------Special values

	static cu_obj<bool>* palternator;

	//-----------------------------------Steepest Descent Solver

	//quantities used to calculate Barzilai-Borwein stepsizes across multiple meshes
	//Accumulate values in these quantities, then obtain stepsizes as:
	//step1 = delta_M_sq / delta_M_dot_delta_G
	//step2 = delta_M_dot_delta_G / delta_G_sq
	static cu_obj<cuReal>* pdelta_M_sq;
	static cu_obj<cuReal>* pdelta_G_sq;
	static cu_obj<cuReal>* pdelta_M_dot_delta_G;

private:

	//---------------------------------------- SET-UP METHODS : DiffEqCUDA.cpp, DiffEq_EvalsCUDA.cu

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	//zero all main reduction values : mxh, dmdt, lte
	void Zero_reduction_values(void);
	void Zero_mxh_lte_values(void);
	void Zero_dmdt_lte_values(void);
	void Zero_lte_value(void);
	
	void mxhav_to_mxh(void);
	void dmdtav_to_dmdt(void);

	//set all cuda values here from their cpu values held in ODECommon
	void SyncODEValues(void);

	//set specific cuda values (used often)
	void Sync_dT(void);
	void Sync_dT_last(void);
	void Sync_alternator(void);

	//specific to SD solver
	void Zero_SD_Solver_BB_Values(void);
	void Get_SD_Solver_BB_Values(double* pdelta_M_sq_cpu, double* pdelta_G_sq_cpu, double* pdelta_M_dot_delta_G_cpu)
	{
		*pdelta_M_sq_cpu = pdelta_M_sq->to_cpu();
		*pdelta_G_sq_cpu = pdelta_G_sq->to_cpu();
		*pdelta_M_dot_delta_G_cpu = pdelta_M_dot_delta_G->to_cpu();
	}

public:

	ODECommonCUDA(void) {}
	ODECommonCUDA(ODECommon *pODE_);

	virtual ~ODECommonCUDA();

	//---------------------------------------- GET METHODS

	cuReal Get_mxh(void) { return pmxh->to_cpu(); }
	cuReal Get_dmdt(void) { return pdmdt->to_cpu(); }
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

	//evalution scratch spaces
	cu_obj<cuVEC<cuReal3>> sEval0, sEval1, sEval2, sEval3, sEval4, sEval5;

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

#ifdef ODE_EVAL_EULER
	//Euler evaluation of ODE
	void RunEuler(bool calculate_mxh, bool calculate_dmdt);
#endif

#ifdef ODE_EVAL_TEULER
	//Trapezoidal Euler evaluation of ODE
	void RunTEuler(int step, bool calculate_mxh, bool calculate_dmdt);
#endif

#ifdef ODE_EVAL_AHEUN
	//Adaptive Huen evaluation of ODE
	void RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt);
#endif

#ifdef ODE_EVAL_ABM
	//ABM
	void RunABM(int step, bool calculate_mxh, bool calculate_dmdt);
	void RunABMTEuler(int step);
#endif

#ifdef ODE_EVAL_RK23
	//RK23 (Bogacki-Shampine)
	void RunRK23_Step0_NoAdvance(bool calculate_mxh = false);
	void RunRK23(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_RK4
	//RK4 evaluation of ODE
	void RunRK4(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_RKF
	//RKF45
	void RunRKF45(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_RKCK
	//RKCK45
	void RunRKCK45(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_RKDP
	//RKDP54
	void RunRKDP54_Step0_NoAdvance(bool calculate_mxh = false);
	void RunRKDP54(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_SD
	//0. prime the SD solver
	void RunSD_Start(void);
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	void RunSD_BB(void);
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	void RunSD_Advance(bool calculate_mxh, bool calculate_dmdt);
#endif

	//---------------------------------------- SOLVER KERNEL LAUNCHERS with in-lined LLG (faster) : DiffEq_EvalsLLGCUDA.cu

	//Euler evaluation of ODE
	void RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt);

	//Trapezoidal Euler evaluation of ODE
	void RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt);

	//Adaptive Huen evaluation of ODE
	void RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt);

	//RK4 evaluation of ODE
	void RunRK4_LLG(int step, bool calculate_mxh = false, bool calculate_dmdt = false);

	//ABM
	void RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt);
	void RunABMTEuler_LLG(int step);

	//RKF45
	void RunRKF45_LLG(int step, bool calculate_mxh = false, bool calculate_dmdt = false);

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
