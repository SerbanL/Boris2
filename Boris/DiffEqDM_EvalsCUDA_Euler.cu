#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_EULER
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void DifferentialEquationDMCUDA::RunEuler(bool calculate_mxh, bool calculate_dmdt)
{
	RunEuler_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
}

#endif
#endif
#endif