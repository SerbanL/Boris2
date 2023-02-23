#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKF45
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
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (2 * dT / 9);
				(*cuMesh.pM2)[idx] += (*cuDiffEq.psEval0_2)[idx] * (2 * dT / 9);
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
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (2 * dT / 9);
				(*cuMesh.pM2)[idx] += (*cuDiffEq.psEval0_2)[idx] * (2 * dT / 9);
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
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] / 12 + (*cuDiffEq.psEval1)[idx] / 4) * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] / 12 + (*cuDiffEq.psEval1_2)[idx] / 4) * dT;
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
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (69 * (*cuDiffEq.psEval0)[idx] / 128 - 243 * (*cuDiffEq.psEval1)[idx] / 128 + 135 * (*cuDiffEq.psEval2)[idx] / 64) * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (69 * (*cuDiffEq.psEval0_2)[idx] / 128 - 243 * (*cuDiffEq.psEval1_2)[idx] / 128 + 135 * (*cuDiffEq.psEval2_2)[idx] / 64) * dT;
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
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (-17 * (*cuDiffEq.psEval0)[idx] / 12 + 27 * (*cuDiffEq.psEval1)[idx] / 4 - 27 * (*cuDiffEq.psEval2)[idx] / 5 + 16 * (*cuDiffEq.psEval3)[idx] / 15) * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (-17 * (*cuDiffEq.psEval0_2)[idx] / 12 + 27 * (*cuDiffEq.psEval1_2)[idx] / 4 - 27 * (*cuDiffEq.psEval2_2)[idx] / 5 + 16 * (*cuDiffEq.psEval3_2)[idx] / 15) * dT;
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
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (65 * (*cuDiffEq.psEval0)[idx] / 432 - 5 * (*cuDiffEq.psEval1)[idx] / 16 + 13 * (*cuDiffEq.psEval2)[idx] / 16 + 4 * (*cuDiffEq.psEval3)[idx] / 27 + 5 * (*cuDiffEq.psEval4)[idx] / 144) * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (65 * (*cuDiffEq.psEval0_2)[idx] / 432 - 5 * (*cuDiffEq.psEval1_2)[idx] / 16 + 13 * (*cuDiffEq.psEval2_2)[idx] / 16 + 4 * (*cuDiffEq.psEval3_2)[idx] / 27 + 5 * (*cuDiffEq.psEval4_2)[idx] / 144) * dT;
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
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] / 9 + 9 * (*cuDiffEq.psEval2)[idx] / 20 + 16 * (*cuDiffEq.psEval3)[idx] / 45 + (*cuDiffEq.psEval4)[idx] / 12) * dT;
				(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] / 9 + 9 * (*cuDiffEq.psEval2_2)[idx] / 20 + 16 * (*cuDiffEq.psEval3_2)[idx] / 45 + (*cuDiffEq.psEval4_2)[idx] / 12) * dT;

				//5th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (47 * (*cuDiffEq.psEval0)[idx] / 450 + 12 * (*cuDiffEq.psEval2)[idx] / 25 + 32 * (*cuDiffEq.psEval3)[idx] / 225 + 1 * (*cuDiffEq.psEval4)[idx] / 30 + 6 * rhs / 25) * dT;

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
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] / 9 + 9 * (*cuDiffEq.psEval2)[idx] / 20 + 16 * (*cuDiffEq.psEval3)[idx] / 45 + (*cuDiffEq.psEval4)[idx] / 12) * dT;
				(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] / 9 + 9 * (*cuDiffEq.psEval2_2)[idx] / 20 + 16 * (*cuDiffEq.psEval3_2)[idx] / 45 + (*cuDiffEq.psEval4_2)[idx] / 12) * dT;

				//5th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (47 * (*cuDiffEq.psEval0)[idx] / 450 + 12 * (*cuDiffEq.psEval2)[idx] / 25 + 32 * (*cuDiffEq.psEval3)[idx] / 225 + 1 * (*cuDiffEq.psEval4)[idx] / 30 + 6 * rhs / 25) * dT;

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

			RunRKF45_Step0_withReductions_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF45_Step0_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;

	case 1:

		RunRKF45_Step1_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRKF45_Step2_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRKF45_Step3_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 4:

		RunRKF45_Step4_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 5:

		if (calculate_dmdt) {

			RunRKF45_Step5_withReductions_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF45_Step5_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

#endif
#endif
#endif