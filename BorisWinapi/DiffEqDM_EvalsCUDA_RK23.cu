#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_RK23
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RK23

__global__ void RunRK23_Step0_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRK23_Step0_Advance_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

__global__ void RunRK23_Step1_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);
		}
	}
}

__global__ void RunRK23_Step2_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

//RUNGE KUTTA 23 (Bogacki - Shampine) (2nd order adaptive step with FSAL, 3rd order evaluation)

void DifferentialEquationDMCUDA::RunRK23_Step0_NoAdvance(bool calculate_mxh)
{
	RunRK23_Step0_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
}

void DifferentialEquationDMCUDA::RunRK23(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		RunRK23_Step0_Advance_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 1:

		RunRK23_Step1_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRK23_Step2_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;
	}
}

#endif
#endif
#endif