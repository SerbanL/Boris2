#pragma once

#include "DiffEqCUDA.h"
#include "ManagedDiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

class DifferentialEquationDMCUDA :
	public DifferentialEquationCUDA
{
private:

	//ManagedDiffEqDMCUDA holds pointers to data in DifferentialEquationCUDA and ODECommonCUDA in an object in gpu memory.
	//pass cuDiffEq to a cuda kernel then all gpu data held here in cu_obj objects can be accessed in device code.
	//Initialize ManagedDiffEqDMCUDA with all the pointers you need then forget about it - no book-keeping required.
	cu_obj<ManagedDiffEqDMCUDA> cuDiffEq;

protected:

	//---------------------------------------- OTHER CALCULATION METHODS

	//called when using stochastic equations
	void GenerateThermalField_CUDA(cu_obj<cuBReal>& deltaT) {}
	void GenerateThermalField_and_Torque_CUDA(cu_obj<cuBReal>& deltaT) {}

public:

	DifferentialEquationDMCUDA(DifferentialEquation *pmeshODE);
	~DifferentialEquationDMCUDA();

	//---------------------------------------- EQUATIONS : these are defined as __device__ methods in ManagedDiffEqDMCUDA

	//---------------------------------------- SOLVER KERNEL LAUNCHERS generic

#ifdef ODE_EVAL_COMPILATION_EULER
	//Euler evaluation of ODE
	void RunEuler(bool calculate_mxh, bool calculate_dmdt);
#endif

#ifdef ODE_EVAL_COMPILATION_TEULER
	//Trapezoidal Euler evaluation of ODE
	void RunTEuler(int step, bool calculate_mxh, bool calculate_dmdt);
#endif

#ifdef ODE_EVAL_COMPILATION_AHEUN
	//Adaptive Huen evaluation of ODE
	void RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt);
#endif

#ifdef ODE_EVAL_COMPILATION_ABM
	//ABM
	void RunABM(int step, bool calculate_mxh, bool calculate_dmdt);
	void RunABMTEuler(int step);
#endif

#ifdef ODE_EVAL_COMPILATION_RK23
	//RK23 (Bogacki-Shampine)
	void RunRK23_Step0_NoAdvance(bool calculate_mxh = false);
	void RunRK23(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_COMPILATION_RK4
	//RK4 evaluation of ODE
	void RunRK4(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_COMPILATION_RKF
	//RKF45
	void RunRKF45(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_COMPILATION_RKCK
	//RKCK45
	void RunRKCK45(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_COMPILATION_RKDP
	//RKDP54
	void RunRKDP54_Step0_NoAdvance(bool calculate_mxh = false);
	void RunRKDP54(int step, bool calculate_mxh = false, bool calculate_dmdt = false);
#endif

#ifdef ODE_EVAL_COMPILATION_SD
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

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	BError AllocateMemory(bool copy_from_cpu = false);

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	void CleanupMemory(bool copy_to_cpu = false);

	void SetODEMethodPointers(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//---------------------------------------- OTHERS

	//Restore magnetisation after a failed step for adaptive time-step methods
	void RestoreMagnetisation(void);

	//----------------------------------- GETTERS

	//get reference to stored managed cuda differential equation object (cuDiffEq)
	cu_obj<ManagedDiffEqDMCUDA>& Get_ManagedDiffEqCUDA(void) { return cuDiffEq; }
};

#else

class DifferentialEquationDMCUDA :
	public DifferentialEquationCUDA
{

private:

	cu_obj<ManagedDiffEqDMCUDA> cuDiffEq;

protected:

	//---------------------------------------- OTHER CALCULATION METHODS

	//called when using stochastic equations
	void GenerateThermalField_CUDA(cu_obj<cuBReal>& deltaT) {}
	void GenerateThermalField_and_Torque_CUDA(cu_obj<cuBReal>& deltaT) {}

public:

	DifferentialEquationDMCUDA(DifferentialEquation *pmeshODE) : DifferentialEquationCUDA() {}
	~DifferentialEquationDMCUDA() {}

	//---------------------------------------- EQUATIONS : these are defined as __device__ methods in ManagedDiffEqAFMCUDA

	//---------------------------------------- SOLVER KERNEL LAUNCHERS generic

#ifdef ODE_EVAL_COMPILATION_EULER
	//Euler evaluation of ODE
	void RunEuler(bool calculate_mxh, bool calculate_dmdt) {}
#endif

#ifdef ODE_EVAL_COMPILATION_TEULER
	//Trapezoidal Euler evaluation of ODE
	void RunTEuler(int step, bool calculate_mxh, bool calculate_dmdt) {}
#endif

#ifdef ODE_EVAL_COMPILATION_AHEUN
	//Adaptive Huen evaluation of ODE
	void RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt) {}
#endif

#ifdef ODE_EVAL_COMPILATION_ABM
	//ABM
	void RunABM(int step, bool calculate_mxh, bool calculate_dmdt) {}
	void RunABMTEuler(int step);
#endif

#ifdef ODE_EVAL_COMPILATION_RK23
	//RK23 (Bogacki-Shampine)
	void RunRK23_Step0_NoAdvance(bool calculate_mxh = false) {}
	void RunRK23(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RK4
	//RK4 evaluation of ODE
	void RunRK4(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RKF
	//RKF45
	void RunRKF45(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RKCK
	//RKCK45
	void RunRKCK45(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RKDP
	//RKDP54
	void RunRKDP54_Step0_NoAdvance(bool calculate_mxh = false) {}
	void RunRKDP54(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}
#endif

#ifdef ODE_EVAL_COMPILATION_SD
	//0. prime the SD solver
	void RunSD_Start(void) {}
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	void RunSD_BB(void) {}
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	void RunSD_Advance(bool calculate_mxh, bool calculate_dmdt) {}
#endif

	//---------------------------------------- SOLVER KERNEL LAUNCHERS with in-lined LLG (faster) : DiffEq_EvalsLLGCUDA.cu

	//Euler evaluation of ODE
	void RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt) {}

	//Trapezoidal Euler evaluation of ODE
	void RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt) {}

	//Adaptive Huen evaluation of ODE
	void RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt) {}

	//RK4 evaluation of ODE
	void RunRK4_LLG(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}

	//ABM
	void RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt) {}
	void RunABMTEuler_LLG(int step) {}

	//RKF45
	void RunRKF45_LLG(int step, bool calculate_mxh = false, bool calculate_dmdt = false) {}

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	BError AllocateMemory(bool copy_from_cpu = false) { return BError(); }

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	void CleanupMemory(bool copy_to_cpu = false) {}

	void SetODEMethodPointers(void) {}

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//---------------------------------------- OTHERS

	//Restore magnetisation after a failed step for adaptive time-step methods
	void RestoreMagnetisation(void) {}

	//----------------------------------- GETTERS

	//get reference to stored managed cuda differential equation object (cuDiffEq)
	cu_obj<ManagedDiffEqDMCUDA>& Get_ManagedDiffEqCUDA(void) { return cuDiffEq; }
};

#endif
#endif
