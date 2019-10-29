#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"

#include "DiffEq_CommonCUDA.h"

//NOTES : renormalize at last step for all equations except all LLB versions

class ODECommon;
class DifferentialEquation;
class DifferentialEquationFM;
class DifferentialEquationAFM;
class ManagedDiffEqFMCUDA;
class ManagedDiffEqAFMCUDA;
class Mesh;
class MeshCUDA;

class DifferentialEquationCUDA :
	public ODECommonCUDA
{
	friend ODECommon;
	friend DifferentialEquationFM;
	friend DifferentialEquationAFM;
	friend ManagedDiffEqFMCUDA;
	friend ManagedDiffEqAFMCUDA;

protected:

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
	Mesh *pMesh;

	//pointer to CUDA version of mesh with this ODE set
	MeshCUDA *pMeshCUDA;

	//pointer to DifferentialEquation holding this CUDA object
	DifferentialEquation *pmeshODE;

protected:

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	virtual BError AllocateMemory(void) = 0;

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	virtual void CleanupMemory(void) = 0;

	virtual void SetODEMethodPointers(void) = 0;

	//---------------------------------------- EQUATIONS : these are defined as __device__ methods in ManagedDiffEqFMCUDA

	//---------------------------------------- SOLVER KERNEL LAUNCHERS generic : DiffEq_EvalsCUDA.cu

#ifdef ODE_EVAL_EULER
	//Euler evaluation of ODE
	virtual void RunEuler(bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

#ifdef ODE_EVAL_TEULER
	//Trapezoidal Euler evaluation of ODE
	virtual void RunTEuler(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

#ifdef ODE_EVAL_AHEUN
	//Adaptive Huen evaluation of ODE
	virtual void RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

#ifdef ODE_EVAL_ABM
	//ABM
	virtual void RunABM(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
	virtual void RunABMTEuler(int step) = 0;
#endif

#ifdef ODE_EVAL_RK23
	//RK23 (Bogacki-Shampine)
	virtual void RunRK23_Step0_NoAdvance(bool calculate_mxh = false) = 0;
	virtual void RunRK23(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_RK4
	//RK4 evaluation of ODE
	virtual void RunRK4(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_RKF
	//RKF45
	virtual void RunRKF45(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_RKCK
	//RKCK45
	virtual void RunRKCK45(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_RKDP
	//RKDP54
	virtual void RunRKDP54_Step0_NoAdvance(bool calculate_mxh = false) = 0;
	virtual void RunRKDP54(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_SD
	//0. prime the SD solver
	virtual void RunSD_Start(void) = 0;
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	virtual void RunSD_BB(void) = 0;
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	virtual void RunSD_Advance(bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

	//---------------------------------------- SOLVER KERNEL LAUNCHERS with in-lined LLG (faster) : DiffEq_EvalsLLGCUDA.cu

	//Euler evaluation of ODE
	virtual void RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt) = 0;

	//Trapezoidal Euler evaluation of ODE
	virtual void RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt) = 0;

	//Adaptive Huen evaluation of ODE
	virtual void RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt) = 0;

	//RK4 evaluation of ODE
	virtual void RunRK4_LLG(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;

	//ABM
	virtual void RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
	virtual void RunABMTEuler_LLG(int step) = 0;

	//RKF45
	virtual void RunRKF45_LLG(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;

	//---------------------------------------- OTHERS : DiffEq_EvalsCUDA.cu

	//Restore magnetisation after a failed step for adaptive time-step methods
	virtual void RestoreMagnetisation(void) = 0;

	//---------------------------------------- OTHER CALCULATION METHODS : DiffEq_SEquationsCUDA.cu

	//called when using stochastic equations
	virtual void GenerateThermalField(void) = 0;
	virtual void GenerateThermalField_and_Torque(void) = 0;

public:

	DifferentialEquationCUDA(DifferentialEquation *pmeshODE);

	virtual ~DifferentialEquationCUDA() {}

	BError Error_On_Create(void) { return error_on_create; }

	//---------------------------------------- SET-UP METHODS

	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) = 0;
};

#endif
