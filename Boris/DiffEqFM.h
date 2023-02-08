#pragma once

#include "DiffEq.h"

class FMesh;

#ifdef MESH_COMPILATION_FERROMAGNETIC

#if COMPILECUDA == 1
#include "DiffEqFMCUDA.h"
#endif

class DifferentialEquationFM :
	public DifferentialEquation
{
private:

public:

	DifferentialEquationFM(FMesh *pMesh);
	~DifferentialEquationFM();

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

#ifdef ODE_EVAL_COMPILATION_RKF45
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

#ifdef ODE_EVAL_COMPILATION_RKF56
	//RKF56
	void RunRKF56_Step0_withReductions(void);
	void RunRKF56_Step0(void);
	void RunRKF56_Step1(void);
	void RunRKF56_Step2(void);
	void RunRKF56_Step3(void);
	void RunRKF56_Step4(void);
	void RunRKF56_Step5(void);
	void RunRKF56_Step6(void);
	void RunRKF56_Step7_withReductions(void);
	void RunRKF56_Step7(void);
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

	void RestoreMagnetization(void);

	//Save current magnetization in sM VECs (e.g. useful to reset dM / dt calculation)
	void SaveMagnetization(void);

	//renormalize vectors to set magnetization length value (which could have a spatial variation)
	void RenormalizeMagnetization(void);

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

class DifferentialEquationFM :
	public DifferentialEquation
{
private:

public:

	DifferentialEquationFM(FMesh *pMesh) : DifferentialEquation() {}
	~DifferentialEquationFM() {}

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

#ifdef ODE_EVAL_COMPILATION_RKF45
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

#ifdef ODE_EVAL_COMPILATION_RKF56
	//RKF56
	void RunRKF56_Step0_withReductions(void) {}
	void RunRKF56_Step0(void) {}
	void RunRKF56_Step1(void) {}
	void RunRKF56_Step2(void) {}
	void RunRKF56_Step3(void) {}
	void RunRKF56_Step4(void) {}
	void RunRKF56_Step5(void) {}
	void RunRKF56_Step6(void) {}
	void RunRKF56_Step7_withReductions(void) {}
	void RunRKF56_Step7(void) {}
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

	void RestoreMagnetization(void) {}

	//Save current magnetization in sM VECs (e.g. useful to reset dM / dt calculation)
	void SaveMagnetization(void) {}

	//renormalize vectors to set magnetization length value (which could have a spatial variation)
	void RenormalizeMagnetization(void) {}

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