#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKDP
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKDP54

__global__ void RunRKDP54_Step0_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRKDP54_Step0_Advance_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

__global__ void RunRKDP54_Step1_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRKDP54_Step2_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRKDP54_Step3_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRKDP54_Step4_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRKDP54_Step5_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

//RUNGE KUTTA DORMAND-PRINCE

void DifferentialEquationDMCUDA::RunRKDP54_Step0_NoAdvance(bool calculate_mxh)
{
	RunRKDP54_Step0_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
}

void DifferentialEquationDMCUDA::RunRKDP54(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		RunRKDP54_Step0_Advance_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 1:

		RunRKDP54_Step1_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRKDP54_Step2_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRKDP54_Step3_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 4:

		RunRKDP54_Step4_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 5:

		RunRKDP54_Step5_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;
	}
}

#endif
#endif
#endif