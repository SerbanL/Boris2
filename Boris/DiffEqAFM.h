#pragma once

#include "DiffEq.h"

class AFMesh;

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#if COMPILECUDA == 1
#include "DiffEqAFMCUDA.h"
#endif

class DifferentialEquationAFM :
	public DifferentialEquation
{

#if COMPILECUDA == 1
	friend DifferentialEquationAFMCUDA;
#endif

private:

	//Extra scratch spaces for AFM equations

	int OmpThreads;

	//Used to save starting magnetization - all evaluation methods do this, even when not needed by the method itself, so we can calculate dM/dt when needed.
	VEC<DBL3> sM1_2;

	//evalution scratch spaces
	VEC<DBL3> sEval0_2, sEval1_2, sEval2_2, sEval3_2, sEval4_2, sEval5_2;

	//Thermal field and torques, enabled only for the stochastic equations
	VEC<DBL3> H_Thermal_2, Torque_Thermal_2;

	//When calling an equation to evaluate it only the sub-lattice A value is returned.
	//The additional sub-lattice B value is set here so we can read it
	vector<DBL3> Equation_Eval_2;

public:

	DifferentialEquationAFM(AFMesh *pMesh);
	~DifferentialEquationAFM();

	//---------------------------------------- SOLVER METHODS : DiffEq_Evals.cpp

#ifdef ODE_EVAL_COMPILATION_EULER
	//Euler evaluation of ODE
	void RunEuler_withReductions(void);
	void RunEuler(void);
#endif

#ifdef ODE_EVAL_COMPILATION_TEULER
	//Trapezoidal Euler evaluation of ODE
	void RunTEuler_Step0_withReductions(void);
	void RunTEuler_Step0(void);
	void RunTEuler_Step1_withReductions(void);
	void RunTEuler_Step1(void);
#endif

#ifdef ODE_EVAL_COMPILATION_AHEUN
	//Adaptive Heun evaluation of ODE
	void RunAHeun_Step0_withReductions(void);
	void RunAHeun_Step0(void);
	void RunAHeun_Step1_withReductions(void);
	void RunAHeun_Step1(void);
#endif

#ifdef ODE_EVAL_COMPILATION_ABM
	//ABM
	void RunABM_Predictor_withReductions(void);
	void RunABM_Predictor(void);
	void RunABM_Corrector_withReductions(void);
	void RunABM_Corrector(void);
	void RunABM_TEuler0(void);
	void RunABM_TEuler1(void);
#endif

#ifdef ODE_EVAL_COMPILATION_RK23
	//RK23
	void RunRK23_Step0_withReductions(void);
	void RunRK23_Step0(void);
	void RunRK23_Step0_Advance(void);
	void RunRK23_Step1(void);
	void RunRK23_Step2_withReductions(void);
	void RunRK23_Step2(void);
#endif

#ifdef ODE_EVAL_COMPILATION_RK4
	//RK4 evaluation of ODE
	void RunRK4_Step0_withReductions(void);
	void RunRK4_Step0(void);
	void RunRK4_Step1(void);
	void RunRK4_Step2(void);
	void RunRK4_Step3_withReductions(void);
	void RunRK4_Step3(void);
#endif

#ifdef ODE_EVAL_COMPILATION_RKF
	//RKF45
	void RunRKF45_Step0_withReductions(void);
	void RunRKF45_Step0(void);
	void RunRKF45_Step1(void);
	void RunRKF45_Step2(void);
	void RunRKF45_Step3(void);
	void RunRKF45_Step4(void);
	void RunRKF45_Step5_withReductions(void);
	void RunRKF45_Step5(void);
#endif

#ifdef ODE_EVAL_COMPILATION_RKCK
	//RKCK45
	void RunRKCK45_Step0_withReductions(void);
	void RunRKCK45_Step0(void);
	void RunRKCK45_Step1(void);
	void RunRKCK45_Step2(void);
	void RunRKCK45_Step3(void);
	void RunRKCK45_Step4(void);
	void RunRKCK45_Step5_withReductions(void);
	void RunRKCK45_Step5(void);
#endif

#ifdef ODE_EVAL_COMPILATION_RKDP
	//RKDP54
	void RunRKDP54_Step0_withReductions(void);
	void RunRKDP54_Step0(void);
	void RunRKDP54_Step0_Advance(void);
	void RunRKDP54_Step1(void);
	void RunRKDP54_Step2(void);
	void RunRKDP54_Step3(void);
	void RunRKDP54_Step4(void);
	void RunRKDP54_Step5_withReductions(void);
	void RunRKDP54_Step5(void);
#endif

#ifdef ODE_EVAL_COMPILATION_SD
	//0. prime the SD solver
	void RunSD_Start(void);
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	void RunSD_BB(void);
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	void RunSD_Advance_withReductions(void);
	void RunSD_Advance(void);
#endif

	//---------------------------------------- EQUATIONS : DiffEq_Equations.cpp and DiffEq_SEquations.cpp

	//Landau-Lifshitz-Gilbert equation
	DBL3 LLG(int idx);

	//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
	DBL3 LLGStatic(int idx);

	//Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	DBL3 LLGSTT(int idx);

	//Landau-Lifshitz-Bloch equation
	DBL3 LLB(int idx);

	//Landau-Lifshitz-Bloch equation with Zhang-Li STT
	DBL3 LLBSTT(int idx);

	//Stochastic Landau-Lifshitz-Gilbert equation
	DBL3 SLLG(int idx);

	//Stochastic Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	DBL3 SLLGSTT(int idx);

	//Stochastic Landau-Lifshitz-Bloch equation
	DBL3 SLLB(int idx);

	//Stochastic Landau-Lifshitz-Bloch equation with Zhang-Li STT
	DBL3 SLLBSTT(int idx);

	//---------------------------------------- OTHERS : DiffEqFM.cpp

	void Restoremagnetization(void);

	//---------------------------------------- OTHER CALCULATION METHODS : DiffEqFM_SEquations.cpp

	//called when using stochastic equations
	void GenerateThermalField(void);
	void GenerateThermalField_and_Torque(void);

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	BError AllocateMemory(void);

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	void CleanupMemory(void);

	//---------------------------------------- SET-UP METHODS

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

	//---------------------------------------- GETTERS

	//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
	DBL3 dMdt(int idx);
};

#else

class DifferentialEquationAFM :
	public DifferentialEquation
{

private:

public:

	DifferentialEquationAFM(AFMesh *pMesh) : DifferentialEquation() {}
	~DifferentialEquationAFM() {}

	//---------------------------------------- SOLVER METHODS : DiffEq_Evals.cpp

#ifdef ODE_EVAL_COMPILATION_EULER
	//Euler evaluation of ODE
	void RunEuler_withReductions(void) {}
	void RunEuler(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_TEULER
	//Trapezoidal Euler evaluation of ODE
	void RunTEuler_Step0_withReductions(void) {}
	void RunTEuler_Step0(void) {}
	void RunTEuler_Step1_withReductions(void) {}
	void RunTEuler_Step1(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_AHEUN
	//Adaptive Heun evaluation of ODE
	void RunAHeun_Step0_withReductions(void) {}
	void RunAHeun_Step0(void) {}
	void RunAHeun_Step1_withReductions(void) {}
	void RunAHeun_Step1(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_ABM
	//ABM
	void RunABM_Predictor_withReductions(void) {}
	void RunABM_Predictor(void) {}
	void RunABM_Corrector_withReductions(void) {}
	void RunABM_Corrector(void) {}
	void RunABM_TEuler0(void) {}
	void RunABM_TEuler1(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RK23
	//RK23
	void RunRK23_Step0_withReductions(void) {}
	void RunRK23_Step0(void) {}
	void RunRK23_Step0_Advance(void) {}
	void RunRK23_Step1(void) {}
	void RunRK23_Step2_withReductions(void) {}
	void RunRK23_Step2(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RK4
	//RK4 evaluation of ODE
	void RunRK4_Step0_withReductions(void) {}
	void RunRK4_Step0(void) {}
	void RunRK4_Step1(void) {}
	void RunRK4_Step2(void) {}
	void RunRK4_Step3_withReductions(void) {}
	void RunRK4_Step3(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RKF
	//RKF45
	void RunRKF45_Step0_withReductions(void) {}
	void RunRKF45_Step0(void) {}
	void RunRKF45_Step1(void) {}
	void RunRKF45_Step2(void) {}
	void RunRKF45_Step3(void) {}
	void RunRKF45_Step4(void) {}
	void RunRKF45_Step5_withReductions(void) {}
	void RunRKF45_Step5(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RKCK
	//RKCK45
	void RunRKCK45_Step0_withReductions(void) {}
	void RunRKCK45_Step0(void) {}
	void RunRKCK45_Step1(void) {}
	void RunRKCK45_Step2(void) {}
	void RunRKCK45_Step3(void) {}
	void RunRKCK45_Step4(void) {}
	void RunRKCK45_Step5_withReductions(void) {}
	void RunRKCK45_Step5(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_RKDP
	//RKDP54
	void RunRKDP54_Step0_withReductions(void) {}
	void RunRKDP54_Step0(void) {}
	void RunRKDP54_Step0_Advance(void) {}
	void RunRKDP54_Step1(void) {}
	void RunRKDP54_Step2(void) {}
	void RunRKDP54_Step3(void) {}
	void RunRKDP54_Step4(void) {}
	void RunRKDP54_Step5_withReductions(void) {}
	void RunRKDP54_Step5(void) {}
#endif

#ifdef ODE_EVAL_COMPILATION_SD
	//0. prime the SD solver
	void RunSD_Start(void) {}
	//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
	//must reset the static delta_... quantities before running these across all meshes
	void RunSD_BB(void) {}
	//2. Set stepsize -> done in the Iterate method
	//3. set new magnetization vectors
	void RunSD_Advance_withReductions(void) {}
	void RunSD_Advance(void) {}
#endif

	//---------------------------------------- EQUATIONS : DiffEq_Equations.cpp and DiffEq_SEquations.cpp

	//Landau-Lifshitz-Gilbert equation
	DBL3 LLG(int idx) { return DBL3(); }

	//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
	DBL3 LLGStatic(int idx) { return DBL3(); }

	//Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	DBL3 LLGSTT(int idx) { return DBL3(); }

	//Landau-Lifshitz-Bloch equation
	DBL3 LLB(int idx) { return DBL3(); }

	//Landau-Lifshitz-Bloch equation with Zhang-Li STT
	DBL3 LLBSTT(int idx) { return DBL3(); }

	//Stochastic Landau-Lifshitz-Gilbert equation
	DBL3 SLLG(int idx) { return DBL3(); }

	//Stochastic Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	DBL3 SLLGSTT(int idx) { return DBL3(); }

	//Stochastic Landau-Lifshitz-Bloch equation
	DBL3 SLLB(int idx) { return DBL3(); }

	//Stochastic Landau-Lifshitz-Bloch equation with Zhang-Li STT
	DBL3 SLLBSTT(int idx) { return DBL3(); }

	//---------------------------------------- OTHERS : DiffEqFM.cpp

	void Restoremagnetization(void) {}

	//---------------------------------------- OTHER CALCULATION METHODS : DiffEqFM_SEquations.cpp

	//called when using stochastic equations
	void GenerateThermalField(void) {}
	void GenerateThermalField_and_Torque(void) {}

	//---------------------------------------- SET-UP METHODS

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	BError AllocateMemory(void) { return BError(); }

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	void CleanupMemory(void) {}

	//---------------------------------------- SET-UP METHODS

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState) { return BError(); }

	//---------------------------------------- GETTERS

	//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
	DBL3 dMdt(int idx) { return DBL3(); }
};

#endif
