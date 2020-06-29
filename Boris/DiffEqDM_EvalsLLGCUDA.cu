#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers. LLG equation in-lined for faster evaluation

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunTEuler_Step0_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunTEuler_Step1_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

//----------------------------------------- EVALUATIONS : Adaptive Heun

//Step0 same as for TEuler

__global__ void RunAHeun_Step1_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRK4_Step1_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRK4_Step2_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRK4_Step3_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

//----------------------------------------- EVALUATIONS : ABM

__global__ void RunABM_Predictor_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunABM_Corrector_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunABMTEuler_Step0_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunABMTEuler_Step1_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF45_Step0_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRKF45_Step1_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRKF45_Step2_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRKF45_Step3_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRKF45_Step4_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

__global__ void RunRKF45_Step5_LLG_Kernel(ManagedDiffEqDMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal susrel = *cuMesh.psusrel;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.psusrel, susrel);

			//Set M from diamagnetic susceptibility
			(*cuMesh.pM)[idx] = susrel * Heff[idx];
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void DifferentialEquationDMCUDA::RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt)
{
	RunEuler_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
}

//TRAPEZOIDAL EULER

void DifferentialEquationDMCUDA::RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		RunTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//ADAPTIVE HEUN

void DifferentialEquationDMCUDA::RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		RunTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunAHeun_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//RUNGE KUTTA 4th order

void DifferentialEquationDMCUDA::RunRK4_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		RunRK4_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 1:

		RunRK4_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRK4_Step2_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRK4_Step3_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;
	}
}

//Adams-Bashforth-Moulton 2nd order

void DifferentialEquationDMCUDA::RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		RunABM_Predictor_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABM_Corrector_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//Adams-Bashforth-Moulton 2nd order priming using Trapezoidal Euler

void DifferentialEquationDMCUDA::RunABMTEuler_LLG(int step)
{
	if (step == 0) {

		RunABMTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABMTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//RUNGE KUTTA FEHLBERG

void DifferentialEquationDMCUDA::RunRKF45_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		RunRKF45_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 1:

		RunRKF45_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRKF45_Step2_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRKF45_Step3_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 4:

		RunRKF45_Step4_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 5:

		RunRKF45_Step5_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;
	}
}

#endif
#endif