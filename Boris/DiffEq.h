#pragma once

#include "BorisLib.h"

#include "ErrorHandler.h"
#include "CompileFlags.h"

#include "Boris_Enums_Defs.h"

#include "DiffEq_Defs.h"

#include "DiffEq_Common.h"

#if COMPILECUDA == 1
#include "DiffEqFMCUDA.h"
#include "DiffEqAFMCUDA.h"
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Abstract base class for DifferentialEquation type objects, e.g. for ferromagnetic or antiferromagnetic meshes.

class DifferentialEquation : 
	public ODECommon 
{
	friend ODECommon;
	friend ODECommon_Base;

#if COMPILECUDA == 1
	friend DifferentialEquationCUDA;
	friend DifferentialEquationFMCUDA;
	friend DifferentialEquationAFMCUDA;
#endif

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//reductions
	OmpReduction<double> mxh_reduction;
	OmpReduction<DBL3> mxh_av_reduction;
	OmpReduction<double> dmdt_reduction;
	OmpReduction<DBL3> dmdt_av_reduction;
	OmpReduction<double> lte_reduction;

	//Used to save starting magnetization - all evaluation methods do this, even when not needed by the method itself, so we can calculate dM/dt when needed.
	VEC<DBL3> sM1;

	//evalution scratch spaces
	VEC<DBL3> sEval0, sEval1, sEval2, sEval3, sEval4, sEval5, sEval6;

	//Thermal field and torques, enabled only for the stochastic equations
	VEC<DBL3> H_Thermal, Torque_Thermal;

	//random number generator
	BorisRand prng;

	Mesh *pMesh = nullptr;

	//unique odeId generated when a new entry is made in the pODE vector : used to delete it in the destructor.
	INT2 odeId;

#if COMPILECUDA == 1
	//When deleting DifferentialEquation object (e.g. deleting mesh) with CUDA switched on, we still need to cleanup by deleting pmeshODECUDA.
	//This will invoke the DifferentialEquationCUDA destructor which attempts to transfer data to cpu
	//The problem is, some of this data can be held in a derived class of DifferentialEquation, whose destructor was already executed - so we shouldn't attempt to copy data there if called_from_destructor flag is true
	//called_from_destructor is false by default; it will only be set to true in the destructor, and thereafter the object is gone, thus called_from_destructor is false during the lifetime of this object
	bool called_from_destructor = false;
	DifferentialEquationCUDA *pmeshODECUDA = nullptr;
#endif

protected:

	//---------------------------------------- SOLVER METHODS : DiffEq_Evals.cpp

#ifdef ODE_EVAL_COMPILATION_EULER
	//Euler evaluation of ODE
	virtual void RunEuler_withReductions(void) = 0;
	virtual void RunEuler(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_TEULER
	//Trapezoidal Euler evaluation of ODE
	virtual void RunTEuler_Step0_withReductions(void) = 0;
	virtual void RunTEuler_Step0(void) = 0;
	virtual void RunTEuler_Step1_withReductions(void) = 0;
	virtual void RunTEuler_Step1(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_AHEUN
	//Adaptive Heun evaluation of ODE
	virtual void RunAHeun_Step0_withReductions(void) = 0;
	virtual void RunAHeun_Step0(void) = 0;
	virtual void RunAHeun_Step1_withReductions(void) = 0;
	virtual void RunAHeun_Step1(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_ABM
	//ABM
	virtual void RunABM_Predictor_withReductions(void) = 0;
	virtual void RunABM_Predictor(void) = 0;
	virtual void RunABM_Corrector_withReductions(void) = 0;
	virtual void RunABM_Corrector(void) = 0;
	virtual void RunABM_TEuler0(void) = 0;
	virtual void RunABM_TEuler1(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RK23
	//RK23
	virtual void RunRK23_Step0_withReductions(void) = 0;
	virtual void RunRK23_Step0(void) = 0;
	virtual void RunRK23_Step0_Advance(void) = 0;
	virtual void RunRK23_Step1(void) = 0;
	virtual void RunRK23_Step2_withReductions(void) = 0;
	virtual void RunRK23_Step2(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RK4
	//RK4 evaluation of ODE
	virtual void RunRK4_Step0_withReductions(void) = 0;
	virtual void RunRK4_Step0(void) = 0;
	virtual void RunRK4_Step1(void) = 0;
	virtual void RunRK4_Step2(void) = 0;
	virtual void RunRK4_Step3_withReductions(void) = 0;
	virtual void RunRK4_Step3(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKF45
	//RKF45
	virtual void RunRKF45_Step0_withReductions(void) = 0;
	virtual void RunRKF45_Step0(void) = 0;
	virtual void RunRKF45_Step1(void) = 0;
	virtual void RunRKF45_Step2(void) = 0;
	virtual void RunRKF45_Step3(void) = 0;
	virtual void RunRKF45_Step4(void) = 0;
	virtual void RunRKF45_Step5_withReductions(void) = 0;
	virtual void RunRKF45_Step5(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKF56
	//RKF45
	virtual void RunRKF56_Step0_withReductions(void) = 0;
	virtual void RunRKF56_Step0(void) = 0;
	virtual void RunRKF56_Step1(void) = 0;
	virtual void RunRKF56_Step2(void) = 0;
	virtual void RunRKF56_Step3(void) = 0;
	virtual void RunRKF56_Step4(void) = 0;
	virtual void RunRKF56_Step5(void) = 0;
	virtual void RunRKF56_Step6(void) = 0;
	virtual void RunRKF56_Step7_withReductions(void) = 0;
	virtual void RunRKF56_Step7(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKCK
	//RKCK45
	virtual void RunRKCK45_Step0_withReductions(void) = 0;
	virtual void RunRKCK45_Step0(void) = 0;
	virtual void RunRKCK45_Step1(void) = 0;
	virtual void RunRKCK45_Step2(void) = 0;
	virtual void RunRKCK45_Step3(void) = 0;
	virtual void RunRKCK45_Step4(void) = 0;
	virtual void RunRKCK45_Step5_withReductions(void) = 0;
	virtual void RunRKCK45_Step5(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_RKDP
	//RKDP54
	virtual void RunRKDP54_Step0_withReductions(void) = 0;
	virtual void RunRKDP54_Step0(void) = 0;
	virtual void RunRKDP54_Step0_Advance(void) = 0;
	virtual void RunRKDP54_Step1(void) = 0;
	virtual void RunRKDP54_Step2(void) = 0;
	virtual void RunRKDP54_Step3(void) = 0;
	virtual void RunRKDP54_Step4(void) = 0;
	virtual void RunRKDP54_Step5_withReductions(void) = 0;
	virtual void RunRKDP54_Step5(void) = 0;
#endif

#ifdef ODE_EVAL_COMPILATION_SD
	//0. prime the SD solver
	virtual void RunSD_Start(void) = 0;
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	virtual void RunSD_BB(void) = 0;
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	virtual void RunSD_Advance_withReductions(void) = 0;
	virtual void RunSD_Advance(void) = 0;
#endif

	//---------------------------------------- OTHERS

	//Restore magnetization after a failed step for adaptive time-step methods
	virtual void RestoreMagnetization(void) = 0;

	//Save current magnetization in sM VECs (e.g. useful to reset dM / dt calculation)
	virtual void SaveMagnetization(void) = 0;

	//renormalize vectors to set magnetization length value (which could have a spatial variation)
	virtual void RenormalizeMagnetization(void) = 0;

	//---------------------------------------- OTHER CALCULATION METHODS

	//called when using stochastic equations
	virtual void GenerateThermalField(void) = 0;
	virtual void GenerateThermalField_and_Torque(void) = 0;

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	virtual BError AllocateMemory(void) = 0;

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	virtual void CleanupMemory(void) = 0;

	//---------------------------------------- EQUATIONS : DiffEq_Equations.cpp and DiffEq_SEquations.cpp

	//Landau-Lifshitz-Gilbert equation
	virtual DBL3 LLG(int idx) = 0;

	//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
	virtual DBL3 LLGStatic(int idx) = 0;

	//Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	virtual DBL3 LLGSTT(int idx) = 0;

	//Landau-Lifshitz-Bloch equation
	virtual DBL3 LLB(int idx) = 0;

	//Landau-Lifshitz-Bloch equation with Zhang-Li STT
	virtual DBL3 LLBSTT(int idx) = 0;

	//Stochastic Landau-Lifshitz-Gilbert equation
	virtual DBL3 SLLG(int idx) = 0;

	//Stochastic Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	virtual DBL3 SLLGSTT(int idx) = 0;

	//Stochastic Landau-Lifshitz-Bloch equation
	virtual DBL3 SLLB(int idx) = 0;

	//Stochastic Landau-Lifshitz-Bloch equation with Zhang-Li STT
	virtual DBL3 SLLBSTT(int idx) = 0;

public:  //public methods

	DifferentialEquation(void) :
		ODECommon(true),
		prng(GetSystemTickCount())
	{}

	DifferentialEquation(Mesh *pMesh);
	virtual ~DifferentialEquation();

	BError Error_On_Create(void) { return error_on_create; }

	//---------------------------------------- SET-UP METHODS

	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) = 0;

	//switch CUDA state on/off
	virtual BError SwitchCUDAState(bool cudaState) = 0;

	//---------------------------------------- GETTERS

	//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
	//Note, this is available even when dM_dt VEC is disabled.
	virtual DBL3 dMdt(int idx) = 0;

#if COMPILECUDA == 1
	//get stored cuda differential equation pointer (pmeshODECUDA)
	DifferentialEquationCUDA* Get_DifferentialEquationCUDA_ptr(void) { return pmeshODECUDA; }
#endif
};
