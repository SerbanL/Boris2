#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

//Defines stochastic equations

#include "BorisCUDALib.h"

#include "DiffEqDMCUDA.h"
#include "ManagedDiffEqDMCUDA.h"

#include "DiffEq_Defs.h"

#include "Funcs_Math_base.h" //includes constant values

#include "MeshParamsControlCUDA.h"

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::SLLG(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::SLLGSTT(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::SLLB(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::SLLBSTT(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

#endif
#endif