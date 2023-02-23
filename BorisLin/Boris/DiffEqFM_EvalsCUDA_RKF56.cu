#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKF56
#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKF56

__global__ void RunRKF56_Step0_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal mxh = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//obtain maximum normalized torque term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				mxh = cu_GetMagnitude((*cuMesh.pM)[idx] ^ (*cuMesh.pHeff)[idx]) / (Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now estimate magnetization using RKF first step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 6);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRKF56_Step0_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now estimate magnetization using RKF first step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 6);
			}
		}
	}
}

__global__ void RunRKF56_Step1_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval1)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RKF midle step 1
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (4 * (*cuDiffEq.psEval0)[idx] + 16 * (*cuDiffEq.psEval1)[idx]) * dT / 75;
		}
	}
}

__global__ void RunRKF56_Step2_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval2)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RKF midle step 2
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (5 * (*cuDiffEq.psEval0)[idx] / 6 - 8 * (*cuDiffEq.psEval1)[idx] / 3 + 5 * (*cuDiffEq.psEval2)[idx] / 2) * dT;
		}
	}
}

__global__ void RunRKF56_Step3_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval3)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RKF midle step 3
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (-8 * (*cuDiffEq.psEval0)[idx] / 5 + 144 * (*cuDiffEq.psEval1)[idx] / 25 - 4 * (*cuDiffEq.psEval2)[idx] + 16 * (*cuDiffEq.psEval3)[idx] / 25) * dT;
		}
	}
}

__global__ void RunRKF56_Step4_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval4)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RKF midle step 4
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (361 * (*cuDiffEq.psEval0)[idx] / 320 - 18 * (*cuDiffEq.psEval1)[idx] / 5 + 407 * (*cuDiffEq.psEval2)[idx] / 128 - 11 * (*cuDiffEq.psEval3)[idx] / 80 + 55 * (*cuDiffEq.psEval4)[idx] / 128) * dT;
		}
	}
}

__global__ void RunRKF56_Step5_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval5)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RKF midle step 4
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (-11 * (*cuDiffEq.psEval0)[idx] / 640 + 11 * (*cuDiffEq.psEval2)[idx] / 256 - 11 * (*cuDiffEq.psEval3)[idx] / 160 + 11 * (*cuDiffEq.psEval4)[idx] / 256) * dT;
		}
	}
}

__global__ void RunRKF56_Step6_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval6)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RKF midle step 4
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (93 * (*cuDiffEq.psEval0)[idx] / 640 - 18 * (*cuDiffEq.psEval1)[idx] / 5 + 803 * (*cuDiffEq.psEval2)[idx] / 256 - 11 * (*cuDiffEq.psEval3)[idx] / 160 + 99 * (*cuDiffEq.psEval4)[idx] / 256 + (*cuDiffEq.psEval6)[idx]) * dT;
		}
	}
}

__global__ void RunRKF56_Step7_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//5th order evaluation
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (31 * (*cuDiffEq.psEval0)[idx] / 384 + 1125 * (*cuDiffEq.psEval2)[idx] / 2816 + 9 * (*cuDiffEq.psEval3)[idx] / 32 + 125 * (*cuDiffEq.psEval4)[idx] / 768 + 5 * (*cuDiffEq.psEval5)[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				cuReal3 lte_diff = 5 * ((*cuDiffEq.psEval0)[idx] + (*cuDiffEq.psEval5)[idx] - (*cuDiffEq.psEval6)[idx] - rhs) * dT / 66;

				if (*cuDiffEq.prenormalize) {

					cuBReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(lte_diff) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuBReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				(*cuMesh.pM)[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunRKF56_Step7_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//5th order evaluation
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (31 * (*cuDiffEq.psEval0)[idx] / 384 + 1125 * (*cuDiffEq.psEval2)[idx] / 2816 + 9 * (*cuDiffEq.psEval3)[idx] / 32 + 125 * (*cuDiffEq.psEval4)[idx] / 768 + 5 * (*cuDiffEq.psEval5)[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				cuReal3 lte_diff = 5 * ((*cuDiffEq.psEval0)[idx] + (*cuDiffEq.psEval5)[idx] - (*cuDiffEq.psEval6)[idx] - rhs) * dT / 66;

				if (*cuDiffEq.prenormalize) {

					cuBReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(lte_diff) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuBReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				(*cuMesh.pM)[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA FEHLBERG 5(6)

void DifferentialEquationFMCUDA::RunRKF56(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKF56_Step0_withReductions_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF56_Step0_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;

	case 1:

		RunRKF56_Step1_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRKF56_Step2_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRKF56_Step3_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 4:

		RunRKF56_Step4_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 5:

		RunRKF56_Step5_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 6:

		RunRKF56_Step6_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 7:

		if (calculate_dmdt) {

			RunRKF56_Step7_withReductions_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF56_Step7_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

#endif
#endif
#endif