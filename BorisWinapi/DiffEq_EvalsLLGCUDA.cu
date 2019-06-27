#include "DiffEqCUDA.h"
#include "MeshParamsControlCUDA.h"

#if COMPILECUDA == 1

//defines evaluation methods kernel launchers. LLG equation in-lined for faster evaluation

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtained average normalized torque term
				cuReal3 mxh = (M[idx] ^ Heff[idx]) / (Ms * Ms);

				//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
				if (cuMesh.pgrel->get0()) {

					reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunTEuler_Step0_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
			}
		}
	}
}

__global__ void RunTEuler_Step1_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step formula
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtained average normalized torque term
				cuReal3 mxh = (M[idx] ^ Heff[idx]) / (Ms * Ms);

				//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
				if (cuMesh.pgrel->get0()) {

					reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				
				(*cuDiffEq.psEval0)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using RK4 midle step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK4_Step1_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			
			(*cuDiffEq.psEval1)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RK4 midle step
			M[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval1)[idx] * (dT / 2);
		}
	}
}

__global__ void RunRK4_Step2_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			
			(*cuDiffEq.psEval2)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RK4 last step
			M[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval2)[idx] * dT;
		}
	}
}

__global__ void RunRK4_Step3_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal _mxh = 0.0;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using previous RK4 evaluations
				M[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtained maximum normalized torque term
				_mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Ms * Ms);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &_mxh, *cuDiffEq.pmxh);
	}
}

//----------------------------------------- EVALUATIONS : ABM

__global__ void RunABM_Predictor_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuDiffEq.palternator) {

					M[idx] += dT * (3 * rhs - (*cuDiffEq.psEval0)[idx]) / 2;
					(*cuDiffEq.psEval1)[idx] = rhs;
				}
				else {

					M[idx] += dT * (3 * rhs - (*cuDiffEq.psEval1)[idx]) / 2;
					(*cuDiffEq.psEval0)[idx] = rhs;
				}
			}
		}
	}
}

__global__ void RunABM_Corrector_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal _mxh = 0.0;
	cuReal _lte = 0.0;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = M[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuDiffEq.palternator) {

					M[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval1)[idx]) / 2;
				}
				else {

					M[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval0)[idx]) / 2;
				}

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				_lte = cu_GetMagnitude(M[idx] - saveM) / Ms;

				//obtained maximum normalized torque term
				_mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Ms * Ms);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &_mxh, *cuDiffEq.pmxh);
		reduction_max(0, 1, &_lte, *cuDiffEq.plte);
	}
}

__global__ void RunABMTEuler_Step0_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				(*cuDiffEq.psEval0)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += (*cuDiffEq.psEval0)[idx] * dT;
			}
		}
	}
}

__global__ void RunABMTEuler_Step1_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step formula
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF45_Step0_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				(*cuDiffEq.psEval0)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using RKF first step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
			}
		}
	}
}

__global__ void RunRKF45_Step1_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval1)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 1
			M[idx] = (*cuDiffEq.psM1)[idx] + (3 * (*cuDiffEq.psEval0)[idx] + 9 * (*cuDiffEq.psEval1)[idx]) * dT / 32;
		}
	}
}

__global__ void RunRKF45_Step2_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval2)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 2
			M[idx] = (*cuDiffEq.psM1)[idx] + (1932 * (*cuDiffEq.psEval0)[idx] - 7200 * (*cuDiffEq.psEval1)[idx] + 7296 * (*cuDiffEq.psEval2)[idx]) * dT / 2197;
		}
	}
}

__global__ void RunRKF45_Step3_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval3)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 3
			M[idx] = (*cuDiffEq.psM1)[idx] + (439 * (*cuDiffEq.psEval0)[idx] / 216 - 8 * (*cuDiffEq.psEval1)[idx] + 3680 * (*cuDiffEq.psEval2)[idx] / 513 - 845 * (*cuDiffEq.psEval3)[idx] / 4104) * dT;
		}
	}
}

__global__ void RunRKF45_Step4_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval4)[idx] = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 4
			M[idx] = (*cuDiffEq.psM1)[idx] + (-8 * (*cuDiffEq.psEval0)[idx] / 27 + 2 * (*cuDiffEq.psEval1)[idx] - 3544 * (*cuDiffEq.psEval2)[idx] / 2565 + 1859 * (*cuDiffEq.psEval3)[idx] / 4104 - 11 * (*cuDiffEq.psEval4)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKF45_Step5_LLG_Kernel(ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal _mxh = 0.0;
	cuReal _lte = 0.0;

	cuReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuReal Ms = *cuMesh.pMs;
	cuReal alpha = *cuMesh.palpha;
	cuReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//RKF45 : 4th order predictor for adaptive time step
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//Now calculate 5th order evaluation
				M[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				_lte = cu_GetMagnitude(M[idx] - prediction) / Ms;

				//obtained maximum normalized torque term
				_mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Ms * Ms);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &_mxh, *cuDiffEq.pmxh);
		reduction_max(0, 1, &_lte, *cuDiffEq.plte);
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void DifferentialEquationCUDA::RunEuler_LLG(void)
{
	RunEuler_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
}

//TRAPEZOIDAL EULER

void DifferentialEquationCUDA::RunTEuler_LLG(int step)
{
	if (step == 0) {

		RunTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//RUNGE KUTTA 4th order

void DifferentialEquationCUDA::RunRK4_LLG(int step)
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

void DifferentialEquationCUDA::RunABM_LLG(int step)
{
	if (step == 0) {

		RunABM_Predictor_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABM_Corrector_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//Adams-Bashforth-Moulton 2nd order priming using Trapezoidal Euler

void DifferentialEquationCUDA::RunABMTEuler_LLG(int step)
{
	if (step == 0) {

		RunABMTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABMTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//RUNGE KUTTA FEHLBERG 4th order predictor with 5th order evaluator

void DifferentialEquationCUDA::RunRKF45_LLG(int step)
{
	switch (step) {

	case 0:

		RunRKF45_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 1:

		RunRKF45_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRKF45_Step2_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		RunRKF45_Step3_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 4:

		RunRKF45_Step4_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 5:

		RunRKF45_Step5_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);

		break;
	}
}

#endif