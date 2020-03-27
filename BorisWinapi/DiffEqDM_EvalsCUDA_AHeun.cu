#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_AHEUN
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunAHeun_Step0_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

__global__ void RunAHeun_Step1_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//TRAPEZOIDAL EULER

void DifferentialEquationDMCUDA::RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		RunAHeun_Step0_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunAHeun_Step1_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

#endif
#endif
#endif