#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

class ODECommonCUDA;

//This holds pointers to managed objects in ODECommonCUDA : set and forget. They are available for use in cuda kernels by passing a cu_obj-managed object ManagedDiffEq_CommonCUDA

class ManagedDiffEq_CommonCUDA {

public:

	//Pointers to data in ODECommonCUDA

	//-----------------------------------Time and Stage Time

	cuBReal* ptime;
	cuBReal* pstagetime;

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

	cuBReal* pdelta_M2_sq;
	cuBReal* pdelta_G2_sq;
	cuBReal* pdelta_M2_dot_delta_G2;

public:

	//---------------------------------------- CONSTRUCTION

	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	__host__ BError set_pointers(ODECommonCUDA* pDiffEqCUDA);
};

#endif
