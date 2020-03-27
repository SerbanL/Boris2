#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_ABM
#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : ABM

__global__ void RunABM_Predictor_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal mxh = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//obtain maximum normalized torque term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				mxh = cu_GetMagnitude((*cuMesh.pM)[idx] ^ (*cuMesh.pHeff)[idx]) / (Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuDiffEq.palternator) {

					(*cuMesh.pM)[idx] += dT * (3 * rhs - (*cuDiffEq.psEval0)[idx]) / 2;
					(*cuDiffEq.psEval1)[idx] = rhs;
				}
				else {

					(*cuMesh.pM)[idx] += dT * (3 * rhs - (*cuDiffEq.psEval1)[idx]) / 2;
					(*cuDiffEq.psEval0)[idx] = rhs;
				}
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunABM_Predictor_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuDiffEq.palternator) {

					(*cuMesh.pM)[idx] += dT * (3 * rhs - (*cuDiffEq.psEval0)[idx]) / 2;
					(*cuDiffEq.psEval1)[idx] = rhs;
				}
				else {

					(*cuMesh.pM)[idx] += dT * (3 * rhs - (*cuDiffEq.psEval1)[idx]) / 2;
					(*cuDiffEq.psEval0)[idx] = rhs;
				}
			}
		}
	}
}

__global__ void RunABM_Corrector_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuDiffEq.palternator) {

					(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval1)[idx]) / 2;
				}
				else {

					(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval0)[idx]) / 2;
				}

				if (*cuDiffEq.prenormalize) {

					cuBReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

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

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunABM_Corrector_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuDiffEq.palternator) {

					(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval1)[idx]) / 2;
				}
				else {

					(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval0)[idx]) / 2;
				}

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

__global__ void RunABMTEuler_Step0_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now estimate magnetization for the next time step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * dT;
			}
		}
	}
}

__global__ void RunABMTEuler_Step1_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx);

				//Now estimate magnetization using the second trapezoidal Euler step equation
				(*cuMesh.pM)[idx] = ((*cuDiffEq.psM1)[idx] + (*cuMesh.pM)[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					cuBReal Ms = *cuMesh.pMs;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
					(*cuMesh.pM)[idx].renormalize(Ms);
				}
			}
			else {

				cuBReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				(*cuMesh.pM)[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//Adams-Bashforth-Moulton 2nd order

void DifferentialEquationFMCUDA::RunABM(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunABM_Predictor_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunABM_Predictor_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunABM_Corrector_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunABM_Corrector_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
}

//Adams-Bashforth-Moulton 2nd order priming using Trapezoidal Euler

void DifferentialEquationFMCUDA::RunABMTEuler(int step)
{
	if (step == 0) {

		RunABMTEuler_Step0_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABMTEuler_Step1_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

#endif
#endif
#endif