#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

//Defines non-stochastic equations

#include "BorisCUDALib.h"

#include "ManagedDiffEqDMCUDA.h"

#include "DiffEq_Defs.h"

#include "Funcs_Math_base.h" //includes constant values

#include "MeshParamsControlCUDA.h"

//----------------------------------------- EQUATIONS

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::LLG(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

__device__ cuReal3 ManagedDiffEqDMCUDA::LLGStatic(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::LLGSTT(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::LLB(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqDMCUDA::LLBSTT(int idx)
{
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal susrel = *pcuMesh->psusrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->psusrel, susrel);

	return susrel * Heff[idx];
}

#endif
#endif
