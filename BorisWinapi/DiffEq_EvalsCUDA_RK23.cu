#include "DiffEqCUDA.h"
#include "MeshParamsControlCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_RK23

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RK23

__global__ void RunRK23_Step0_withReductions_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//obtain maximum normalized torque term
				cuReal Mnorm = (*cuMesh.pM)[idx].norm();
				cuReal mxh = cu_GetMagnitude((*cuMesh.pM)[idx] ^ (*cuMesh.pHeff)[idx]) / (Mnorm * Mnorm);

				//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
				if (cuMesh.pgrel->get0()) {

					reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
				}

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//2nd order evaluation for adaptive step
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (7 * (*cuDiffEq.psEval0)[idx] / 24 + 1 * (*cuDiffEq.psEval1)[idx] / 4 + 1 * (*cuDiffEq.psEval2)[idx] / 3 + 1 * rhs / 8) * dT;

				//Save current magnetization for later use
				(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

				//local truncation error (between predicted and corrected)
				cuReal lte = cu_GetMagnitude((*cuMesh.pM)[idx] - prediction) / Mnorm;
				reduction_max(0, 1, &lte, *cuDiffEq.plte);

				//save evaluation for later use
				(*cuDiffEq.psEval0)[idx] = rhs;

				//Now estimate magnetization using RK23 first step
				(*cuMesh.pM)[idx] += rhs * (dT / 2);
			}
		}
	}
}

__global__ void RunRK23_Step0_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//2nd order evaluation for adaptive step
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (7 * (*cuDiffEq.psEval0)[idx] / 24 + 1 * (*cuDiffEq.psEval1)[idx] / 4 + 1 * (*cuDiffEq.psEval2)[idx] / 3 + 1 * rhs / 8) * dT;

				//Save current magnetization for later use
				(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

				//local truncation error (between predicted and corrected)
				cuReal lte = cu_GetMagnitude((*cuMesh.pM)[idx] - prediction) / (*cuMesh.pM)[idx].norm();
				reduction_max(0, 1, &lte, *cuDiffEq.plte);

				//save evaluation for later use
				(*cuDiffEq.psEval0)[idx] = rhs;

				//Now estimate magnetization using RK23 first step
				(*cuMesh.pM)[idx] += rhs * (dT / 2);
			}
		}
	}
}

__global__ void RunRK23_Step1_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuDiffEq.psEval1)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

			//Now estimate magnetization using RK23 midle step 1
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + 3 * (*cuDiffEq.psEval1)[idx] * dT / 4;
		}
	}
}

__global__ void RunRK23_Step2_withReductions_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval2)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now calculate 3rd order evaluation
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (2 * (*cuDiffEq.psEval0)[idx] / 9 + 1 * (*cuDiffEq.psEval1)[idx] / 3 + 4 * (*cuDiffEq.psEval2)[idx] / 9) * dT;

				if (*cuDiffEq.prenormalize) {

					cuReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuReal Mnorm = (*cuMesh.pM)[idx].norm();
				cuReal dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuReal)GAMMA * Mnorm * Mnorm);

				//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
				if (cuMesh.pgrel->get0()) {

					reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
				}
			}
			else {

				cuReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				(*cuMesh.pM)[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}
}

__global__ void RunRK23_Step2_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval2)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now calculate 3rd order evaluation
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (2 * (*cuDiffEq.psEval0)[idx] / 9 + 1 * (*cuDiffEq.psEval1)[idx] / 3 + 4 * (*cuDiffEq.psEval2)[idx] / 9) * dT;

				if (*cuDiffEq.prenormalize) {

					cuReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}
			}
			else {

				cuReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				(*cuMesh.pM)[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA 23 (Bogacki - Shampine) (2nd order adaptive step with FSAL, 3rd order evaluation)

void DifferentialEquationCUDA::RunRK23(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRK23_Step0_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRK23_Step0_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;

	case 1:

		RunRK23_Step1_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		if (calculate_dmdt) {

			RunRK23_Step2_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRK23_Step2_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

#endif
#endif