#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKF
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF45_Step0_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal mxh = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
			(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//obtain maximum normalized torque term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				mxh = cu_GetMagnitude((*cuMesh.pM)[idx] ^ (*cuMesh.pHeff)[idx]) / (Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval0_2)[idx]);

				//Now estimate magnetization using RKF first step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
				(*cuMesh.pM2)[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 4);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRKF45_Step0_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
			(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval0_2)[idx]);

				//Now estimate magnetization using RKF first step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
				(*cuMesh.pM2)[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 4);
			}
		}
	}
}

__global__ void RunRKF45_Step1_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval1)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval1_2)[idx]);

			//Now estimate magnetization using RKF midle step 1
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (3 * (*cuDiffEq.psEval0)[idx] + 9 * (*cuDiffEq.psEval1)[idx]) * dT / 32;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (3 * (*cuDiffEq.psEval0_2)[idx] + 9 * (*cuDiffEq.psEval1_2)[idx]) * dT / 32;
		}
	}
}

__global__ void RunRKF45_Step2_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval2)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval2_2)[idx]);

			//Now estimate magnetization using RKF midle step 2
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (1932 * (*cuDiffEq.psEval0)[idx] - 7200 * (*cuDiffEq.psEval1)[idx] + 7296 * (*cuDiffEq.psEval2)[idx]) * dT / 2197;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (1932 * (*cuDiffEq.psEval0_2)[idx] - 7200 * (*cuDiffEq.psEval1_2)[idx] + 7296 * (*cuDiffEq.psEval2_2)[idx]) * dT / 2197;
		}
	}
}

__global__ void RunRKF45_Step3_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval3)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval3_2)[idx]);

			//Now estimate magnetization using RKF midle step 3
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (439 * (*cuDiffEq.psEval0)[idx] / 216 - 8 * (*cuDiffEq.psEval1)[idx] + 3680 * (*cuDiffEq.psEval2)[idx] / 513 - 845 * (*cuDiffEq.psEval3)[idx] / 4104) * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (439 * (*cuDiffEq.psEval0_2)[idx] / 216 - 8 * (*cuDiffEq.psEval1_2)[idx] + 3680 * (*cuDiffEq.psEval2_2)[idx] / 513 - 845 * (*cuDiffEq.psEval3_2)[idx] / 4104) * dT;
		}
	}
}

__global__ void RunRKF45_Step4_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval4)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval4_2)[idx]);

			//Now estimate magnetization using RKF midle step 4
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (-8 * (*cuDiffEq.psEval0)[idx] / 27 + 2 * (*cuDiffEq.psEval1)[idx] - 3544 * (*cuDiffEq.psEval2)[idx] / 2565 + 1859 * (*cuDiffEq.psEval3)[idx] / 4104 - 11 * (*cuDiffEq.psEval4)[idx] / 40) * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (-8 * (*cuDiffEq.psEval0_2)[idx] / 27 + 2 * (*cuDiffEq.psEval1_2)[idx] - 3544 * (*cuDiffEq.psEval2_2)[idx] / 2565 + 1859 * (*cuDiffEq.psEval3_2)[idx] / 4104 - 11 * (*cuDiffEq.psEval4_2)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKF45_Step5_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs2;
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, rhs2);

				//4th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;
				(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (16 * (*cuDiffEq.psEval0_2)[idx] / 135 + 6656 * (*cuDiffEq.psEval2_2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3_2)[idx] / 56430 - 9 * (*cuDiffEq.psEval4_2)[idx] / 50 + 2 * rhs2 / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
					(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
					(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - prediction) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunRKF45_Step5_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs2;
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, rhs2);

				//4th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;
				(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (16 * (*cuDiffEq.psEval0_2)[idx] / 135 + 6656 * (*cuDiffEq.psEval2_2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3_2)[idx] / 56430 - 9 * (*cuDiffEq.psEval4_2)[idx] / 50 + 2 * rhs2 / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
					(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
					(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - prediction) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA FEHLBERG

void DifferentialEquationAFMCUDA::RunRKF45(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKF45_Step0_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF45_Step0_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;

	case 1:

		RunRKF45_Step1_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRKF45_Step2_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRKF45_Step3_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 4:

		RunRKF45_Step4_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 5:

		if (calculate_dmdt) {

			RunRKF45_Step5_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF45_Step5_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

#endif
#endif
#endif