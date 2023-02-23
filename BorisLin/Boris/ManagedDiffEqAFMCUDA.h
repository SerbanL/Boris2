#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"
#include "ErrorHandler.h"

class DifferentialEquationAFMCUDA;

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "Mesh_AntiFerromagneticCUDA.h"

//This holds pointers to managed objects in DiffEqCUDA : set and forget. They are available for use in cuda kernels by passing a cu_obj-managed object ManagedDiffEqAFMCUDA

class ManagedDiffEqAFMCUDA {

	typedef cuReal3(ManagedDiffEqAFMCUDA::*pODE_t)(int, cuReal3&);

public:

	//Pointers to data in ODECommonCUDA

	//-----------------------------------Time step

	cuBReal* pdT;
	cuBReal* pdT_last;
	
	//-----------------------------------Primary data

	cuBReal* pmxh;
	cuReal3* pmxh_av;
	size_t* pavpoints;

	cuBReal* pdmdt;
	cuReal3* pdmdt_av;
	size_t* pavpoints2;
	
	cuBReal* plte;
	
	//-----------------------------------Evaluation method modifiers

	bool* prenormalize;
	
	//-----------------------------------Properties flags

	bool* psolve_spin_current;
	
	//-----------------------------------Equation and Evaluation method values

	int* psetODE;
	
	//-----------------------------------Special values

	bool* palternator;

	//-----------------------------------Steepest Descent Solver

	//quantities used to calculate Barzilai-Borwein stepsizes across multiple meshes
	//Accumulate values in these quantities, then obtain stepsizes as:
	//step1 = delta_m_sq / delta_m_dot_delta_G
	//step2 = delta_m_dot_delta_G / delta_G_sq
	cuBReal* pdelta_M_sq;
	cuBReal* pdelta_G_sq;
	cuBReal* pdelta_M_dot_delta_G;

	cuBReal* pdelta_M2_sq;
	cuBReal* pdelta_G2_sq;
	cuBReal* pdelta_M2_dot_delta_G2;

	//-----------------------------------Pointers to data in DifferentialEquationCUDA

	//Used for Trapezoidal Euler, RK4, ABM
	cuVEC<cuReal3>* psM1;
	cuVEC<cuReal3>* psM1_2;

	//scratch spaces for evaluations
	cuVEC<cuReal3>* psEval0;
	cuVEC<cuReal3>* psEval1;
	cuVEC<cuReal3>* psEval2;
	cuVEC<cuReal3>* psEval3;
	cuVEC<cuReal3>* psEval4;
	cuVEC<cuReal3>* psEval5;
	cuVEC<cuReal3>* psEval6;

	cuVEC<cuReal3>* psEval0_2;
	cuVEC<cuReal3>* psEval1_2;
	cuVEC<cuReal3>* psEval2_2;
	cuVEC<cuReal3>* psEval3_2;
	cuVEC<cuReal3>* psEval4_2;
	cuVEC<cuReal3>* psEval5_2;
	cuVEC<cuReal3>* psEval6_2;

	//Thermal field and torques, enabled only for the stochastic equations
	cuVEC<cuReal3> *pH_Thermal, *pH_Thermal_2;
	cuVEC<cuReal3> *pTorque_Thermal, *pTorque_Thermal_2;

	//Managed cuda mesh pointer so all mesh data can be accessed in device code
	ManagedMeshCUDA* pcuMesh;

	//pointer to device methods ODEs
	pODE_t pODEFunc;

public:

	//---------------------------------------- CONSTRUCTION

	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	__host__ BError set_pointers(DifferentialEquationAFMCUDA* pDiffEqCUDA);

	//---------------------------------------- EQUATIONS : DiffEq_EquationsCUDA.h and DiffEq_SEquationsCUDA.h
	
	//Landau-Lifshitz-Gilbert equation
	__device__ cuReal3 LLG(int idx, cuReal3& value_B);

	//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
	__device__ cuReal3 LLGStatic(int idx, cuReal3& value_B);

	//Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	__device__ cuReal3 LLGSTT(int idx, cuReal3& value_B);

	//Landau-Lifshitz-Bloch equation
	__device__ cuReal3 LLB(int idx, cuReal3& value_B);
	
	//Landau-Lifshitz-Bloch equation with Zhang-Li STT
	__device__ cuReal3 LLBSTT(int idx, cuReal3& value_B);
	
	//Stochastic Landau-Lifshitz-Gilbert equation
	__device__ cuReal3 SLLG(int idx, cuReal3& value_B);

	//Stochastic Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	__device__ cuReal3 SLLGSTT(int idx, cuReal3& value_B);

	//Stochastic Landau-Lifshitz-Bloch equation
	__device__ cuReal3 SLLB(int idx, cuReal3& value_B);

	//Stochastic Landau-Lifshitz-Bloch equation with Zhang-Li STT
	__device__ cuReal3 SLLBSTT(int idx, cuReal3& value_B);

	//---------------------------------------- GETTERS

	//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
	__device__ cuReal3 dMdt(int idx) { return ((*(pcuMesh->pM))[idx] - (*psM1)[idx]) / *pdT_last; }
	__device__ cuReal3 dMdt2(int idx) { return ((*(pcuMesh->pM2))[idx] - (*psM1_2)[idx]) / *pdT_last; }
};

#else

class ManagedDiffEqAFMCUDA {

public:

public:

	//---------------------------------------- CONSTRUCTION

	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	__host__ BError set_pointers(DifferentialEquationAFMCUDA* pDiffEqCUDA) { return BError(); }
};

#endif
#endif