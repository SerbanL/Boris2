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

	//-----------------------------------Primary Data

	cuBReal* ptime;
	cuBReal* pstagetime;

	//-----------------------------------Time step

	cuBReal* pdT;
	cuBReal* pdT_last;

	//-----------------------------------Equation

	int* psetODE;

	//-----------------------------------mxh and dmdt

	cuBReal* pmxh;
	cuReal3* pmxh_av;
	size_t* pavpoints;

	cuBReal* pdmdt;
	cuReal3* pdmdt_av;
	size_t* pavpoints2;

	//----------------------------------Adaptive time step control

	cuBReal* plte;

	//-----------------------------------Special evaluation values

	bool* palternator;

	//-----------------------------------Evaluation method modifiers

	bool* prenormalize;

	//-----------------------------------Special Properties

	bool* psolve_spin_current;

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

public:

	//---------------------------------------- CONSTRUCTION

	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	__host__ BError set_pointers(ODECommonCUDA* pDiffEqCUDA);
};

#endif
