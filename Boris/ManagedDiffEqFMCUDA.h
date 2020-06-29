#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"
#include "ErrorHandler.h"

class DifferentialEquationFMCUDA;

#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "Mesh_FerromagneticCUDA.h"

//This holds pointers to managed objects in DiffEqCUDA : set and forget. They are available for use in cuda kernels by passing a cu_obj-managed object ManagedDiffEqFMCUDA

class ManagedDiffEqFMCUDA {

	typedef cuReal3(ManagedDiffEqFMCUDA::*pODE_t)(int);

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
	//step1 = delta_M_sq / delta_M_dot_delta_G
	//step2 = delta_M_dot_delta_G / delta_G_sq
	cuBReal* pdelta_M_sq;
	cuBReal* pdelta_G_sq;
	cuBReal* pdelta_M_dot_delta_G;

	//-----------------------------------Pointers to data in DifferentialEquationCUDA

	//Used for Trapezoidal Euler, RK4, ABM
	cuVEC<cuReal3>* psM1;

	//scratch spaces for evaluations
	cuVEC<cuReal3>* psEval0;
	cuVEC<cuReal3>* psEval1;
	cuVEC<cuReal3>* psEval2;
	cuVEC<cuReal3>* psEval3;
	cuVEC<cuReal3>* psEval4;
	cuVEC<cuReal3>* psEval5;

	//Thermal field and torques, enabled only for the stochastic equations
	cuVEC<cuReal3>* pH_Thermal;
	cuVEC<cuReal3>* pTorque_Thermal;

	//Managed cuda mesh pointer so all mesh data can be accessed in device code
	ManagedMeshCUDA* pcuMesh;

	//pointer to device methods ODEs
	pODE_t pODEFunc;

public:

	//---------------------------------------- CONSTRUCTION

	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	__host__ BError set_pointers(DifferentialEquationFMCUDA* pDiffEqCUDA);

	//---------------------------------------- EQUATIONS : DiffEq_EquationsCUDA.h and DiffEq_SEquationsCUDA.h
	
	//Landau-Lifshitz-Gilbert equation
	__device__ cuReal3 LLG(int idx);

	//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
	__device__ cuReal3 LLGStatic(int idx);

	//Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	__device__ cuReal3 LLGSTT(int idx);

	//Landau-Lifshitz-Bloch equation
	__device__ cuReal3 LLB(int idx);
	
	//Landau-Lifshitz-Bloch equation with Zhang-Li STT
	__device__ cuReal3 LLBSTT(int idx);
	
	//Stochastic Landau-Lifshitz-Gilbert equation
	__device__ cuReal3 SLLG(int idx);

	//Stochastic Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	__device__ cuReal3 SLLGSTT(int idx);

	//Stochastic Landau-Lifshitz-Bloch equation
	__device__ cuReal3 SLLB(int idx);

	//Stochastic Landau-Lifshitz-Bloch equation with Zhang-Li STT
	__device__ cuReal3 SLLBSTT(int idx);

	//---------------------------------------- GETTERS

	//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
	__device__ cuReal3 dMdt(int idx) { return ((*(pcuMesh->pM))[idx] - (*psM1)[idx]) / *pdT_last; }
};

#else

class ManagedDiffEqFMCUDA {

public:

public:

	//---------------------------------------- CONSTRUCTION

	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	__host__ BError set_pointers(DifferentialEquationFMCUDA* pDiffEqCUDA) { return BError(); }
};

#endif
#endif