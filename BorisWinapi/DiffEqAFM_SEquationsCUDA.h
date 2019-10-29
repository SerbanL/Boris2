#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

//Defines stochastic equations

#include "BorisCUDALib.h"

#include "DiffEqAFMCUDA.h"
#include "ManagedDiffEqAFMCUDA.h"

#include "DiffEq_Defs.h"

#include "Funcs_Math_base.h" //includes constant values

#include "MeshParamsControlCUDA.h"

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::SLLG(int idx, cuReal3& value_B)
{
	return cuReal3();
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::SLLGSTT(int idx, cuReal3& value_B)
{
	return cuReal3();
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::SLLB(int idx, cuReal3& value_B)
{
	return cuReal3();
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::SLLBSTT(int idx, cuReal3& value_B)
{
	return cuReal3();
}

#endif