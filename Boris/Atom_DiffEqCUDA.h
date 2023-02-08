#pragma once

#include "CompileFlags.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "Atom_MeshParamsCUDA.h"

#include "Atom_DiffEq_CommonCUDA.h"

class Atom_ODECommon;
class ODECommon_Base;
class Atom_DifferentialEquation;
class Atom_DifferentialEquationCubic;
class ManagedAtom_DiffEqCubicCUDA;
class Atom_Mesh;
class Atom_MeshCUDA;

class Atom_DifferentialEquationCUDA :
	public Atom_ODECommonCUDA
{
	friend Atom_ODECommon;
	friend ODECommon_Base;
	friend Atom_DifferentialEquationCubic;
	friend ManagedAtom_DiffEqCubicCUDA;

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//Used to save starting moment - all evaluation methods do this, even when not needed by the method itself, so we can calculate dM/dt when needed.
	cu_obj<cuVEC<cuReal3>> sM1;

	//evalution scratch spaces
	cu_obj<cuVEC<cuReal3>> sEval0, sEval1, sEval2, sEval3, sEval4, sEval5, sEval6;

	//Thermal field, enabled only for the stochastic equations
	cu_obj<cuVEC<cuReal3>> H_Thermal;

	//pseudo-random number generator for use in cuda kernels
	cu_obj<cuBorisRand<>> prng;

	//deltaT to use for stochastic field generation : set before generating stochastic fields if not linked to dT
	cu_obj<cuBReal> deltaTstoch;

	//pointer to mesh with this ODE set
	Atom_Mesh *paMesh = nullptr;

	//pointer to CUDA version of mesh with this ODE set
	Atom_MeshCUDA *paMeshCUDA;

	//pointer to DifferentialEquation holding this CUDA object
	Atom_DifferentialEquation *pameshODE;

protected:

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	virtual BError AllocateMemory(bool copy_from_cpu = false) = 0;

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	virtual void CleanupMemory(bool copy_to_cpu = false) = 0;

	virtual void SetODEMethodPointers(void) = 0;

	//---------------------------------------- EQUATIONS : these are defined as __device__ methods in ManagedAtom_DiffEq...CUDA

	//---------------------------------------- SOLVER KERNEL LAUNCHERS generic : Atom_DiffEq_EvalsCUDA.cu

#ifdef ODE_EVAL_COMPILATION_EULER
	//Euler evaluation of ODE
	virtual void RunEuler(bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_TEULER
	//Trapezoidal Euler evaluation of ODE
	virtual void RunTEuler(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_AHEUN
	//Adaptive Huen evaluation of ODE
	virtual void RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_ABM
	//ABM
	virtual void RunABM(int step, bool calculate_mxh, bool calculate_dmdt) = 0;
	virtual void RunABMTEuler(int step) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RK23
	//RK23 (Bogacki-Shampine)
	virtual void RunRK23_Step0_NoAdvance(bool calculate_mxh = false) = 0;
	virtual void RunRK23(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RK4
	//RK4 evaluation of ODE
	virtual void RunRK4(int step, bool calculate_mxh = false, bool calculate_dmdt = false, bool stochastic = false) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKF45
	//RKF45
	virtual void RunRKF45(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKF56
	//RKF56
	virtual void RunRKF56(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKCK
	//RKCK45
	virtual void RunRKCK45(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKDP
	//RKDP54
	virtual void RunRKDP54_Step0_NoAdvance(bool calculate_mxh = false) = 0;
	virtual void RunRKDP54(int step, bool calculate_mxh = false, bool calculate_dmdt = false) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_SD
	//0. prime the SD solver
	virtual void RunSD_Start(void) = 0;
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	virtual void RunSD_BB(void) = 0;
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	virtual void RunSD_Advance(bool calculate_mxh, bool calculate_dmdt) = 0;
#endif

	//---------------------------------------- SOLVER KERNEL LAUNCHERS with in-lined LLG (faster) : Atom_DiffEq_EvalsLLGCUDA.cu

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

	//Restore moments after a failed step for adaptive time-step methods
	virtual void RestoreMoments(void) = 0;

	//Save current moments in sM VECs (e.g. useful to reset dM / dt calculation)
	virtual void SaveMoments(void) = 0;

	//renormalize vectors to set moment length value (which could have a spatial variation)
	virtual void RenormalizeMoments(void) = 0;

	//---------------------------------------- OTHER CALCULATION METHODS : Atom_DiffEq_SEquationsCUDA.cu

	//called when using stochastic equations
	void GenerateThermalField(void);

	//called when using stochastic equations : these hold kernel launchers
	virtual void GenerateThermalField_CUDA(cu_obj<cuBReal>& deltaT) = 0;

public:

	Atom_DifferentialEquationCUDA(void) :
		Atom_ODECommonCUDA()
	{}

	Atom_DifferentialEquationCUDA(Atom_DifferentialEquation *pameshODE);

	virtual ~Atom_DifferentialEquationCUDA() {}

	BError Error_On_Create(void) { return error_on_create; }

	//---------------------------------------- SET-UP METHODS

	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) = 0;
};

#endif
