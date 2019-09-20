#include "DiffEqCUDA.h"
#include "MeshParamsControlCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_TEULER

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunAHeun_Step0_withReductions_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuReal3 mxh = cuReal3();

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//obtain average normalized torque term
			cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
			
			if (cuDiffEq.pH_Thermal->linear_size()) {

				mxh = ((*cuMesh.pM)[idx] ^ ((*cuMesh.pHeff)[idx] + (*cuDiffEq.pH_Thermal)[idx])) / (Mnorm * Mnorm);
			}
			else mxh = ((*cuMesh.pM)[idx] ^ (*cuMesh.pHeff)[idx]) / (Mnorm * Mnorm);

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now estimate magnetization for the next time step
				(*cuMesh.pM)[idx] += rhs * dT;
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints);
	}
}

__global__ void RunAHeun_Step0_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now estimate magnetization for the next time step
				(*cuMesh.pM)[idx] += rhs * dT;
			}
		}
	}
}

__global__ void RunAHeun_Step1_withReductions_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuReal3 dmdt = cuReal3();
	cuBReal lte = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//Now estimate magnetization using the second trapezoidal Euler step formula
				(*cuMesh.pM)[idx] = ((*cuDiffEq.psM1)[idx] + (*cuMesh.pM)[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					cuBReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - saveM) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuBReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				(*cuMesh.pM)[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt (and mxh) if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunAHeun_Step1_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//Now estimate magnetization using the second trapezoidal Euler step formula
				(*cuMesh.pM)[idx] = ((*cuDiffEq.psM1)[idx] + (*cuMesh.pM)[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					cuBReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - saveM) / (*cuMesh.pM)[idx].norm();
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

//TRAPEZOIDAL EULER

void DifferentialEquationCUDA::RunAHeun(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunAHeun_Step0_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunAHeun_Step0_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunAHeun_Step1_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunAHeun_Step1_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
}

#endif
#endif