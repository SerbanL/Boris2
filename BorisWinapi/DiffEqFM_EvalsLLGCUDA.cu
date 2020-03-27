#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers. LLG equation in-lined for faster evaluation

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuReal3 mxh = cuReal3();
	cuReal3 dmdt = cuReal3();

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//obtain average normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = (M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}
				
				//obtain maximum normalized dmdt term
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt (and mxh) if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints);
		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2);
	}
}

__global__ void RunEuler_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
				
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

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunTEuler_Step0_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuReal3 mxh = cuReal3();

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//obtain average normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = (M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints);
	}
}

__global__ void RunTEuler_Step0_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
			}
		}
	}
}

__global__ void RunTEuler_Step1_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuReal3 dmdt = cuReal3();

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2);
	}
}

__global__ void RunTEuler_Step1_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step equation
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

//----------------------------------------- EVALUATIONS : Adaptive Heun

//Step0 same as for TEuler

__global__ void RunAHeun_Step1_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal lte = 0.0;
	cuReal3 dmdt = cuReal3();

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - saveM) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunAHeun_Step1_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - saveM) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal mxh = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);

				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//Now estimate magnetization using RK4 midle step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRK4_Step0_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				
				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using RK4 midle step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK4_Step1_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			
			(*cuDiffEq.psEval1)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RK4 midle step
			M[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval1)[idx] * (dT / 2);
		}
	}
}

__global__ void RunRK4_Step2_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			
			(*cuDiffEq.psEval2)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RK4 last step
			M[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval2)[idx] * dT;
		}
	}
}

__global__ void RunRK4_Step3_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal dmdt = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using previous RK4 evaluations
				M[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}
}

__global__ void RunRK4_Step3_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using previous RK4 evaluations
				M[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);

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

//----------------------------------------- EVALUATIONS : ABM

__global__ void RunABM_Predictor_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal mxh = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

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

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunABM_Predictor_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

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

__global__ void RunABM_Corrector_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

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

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);		

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - saveM) / Ms;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunABM_Corrector_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

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
				lte = cu_GetMagnitude(M[idx] - saveM) / Ms;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunABMTEuler_Step0_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += (*cuDiffEq.psEval0)[idx] * dT;
			}
		}
	}
}

__global__ void RunABMTEuler_Step1_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step equation
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

__global__ void RunRKF45_Step0_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal mxh = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//Now estimate magnetization using RKF first step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRKF45_Step0_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//Now estimate magnetization using RKF first step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
			}
		}
	}
}

__global__ void RunRKF45_Step1_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval1)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 1
			M[idx] = (*cuDiffEq.psM1)[idx] + (3 * (*cuDiffEq.psEval0)[idx] + 9 * (*cuDiffEq.psEval1)[idx]) * dT / 32;
		}
	}
}

__global__ void RunRKF45_Step2_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval2)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 2
			M[idx] = (*cuDiffEq.psM1)[idx] + (1932 * (*cuDiffEq.psEval0)[idx] - 7200 * (*cuDiffEq.psEval1)[idx] + 7296 * (*cuDiffEq.psEval2)[idx]) * dT / 2197;
		}
	}
}

__global__ void RunRKF45_Step3_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval3)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 3
			M[idx] = (*cuDiffEq.psM1)[idx] + (439 * (*cuDiffEq.psEval0)[idx] / 216 - 8 * (*cuDiffEq.psEval1)[idx] + 3680 * (*cuDiffEq.psEval2)[idx] / 513 - 845 * (*cuDiffEq.psEval3)[idx] / 4104) * dT;
		}
	}
}

__global__ void RunRKF45_Step4_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
			(*cuDiffEq.psEval4)[idx] = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

			//Now estimate magnetization using RKF midle step 4
			M[idx] = (*cuDiffEq.psM1)[idx] + (-8 * (*cuDiffEq.psEval0)[idx] / 27 + 2 * (*cuDiffEq.psEval1)[idx] - 3544 * (*cuDiffEq.psEval2)[idx] / 2565 + 1859 * (*cuDiffEq.psEval3)[idx] / 4104 - 11 * (*cuDiffEq.psEval4)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKF45_Step5_LLG_withReductions_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//4th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				M[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - prediction) / Ms;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunRKF45_Step5_LLG_Kernel(ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuBReal Ms = *cuMesh.pMs;
	cuBReal alpha = *cuMesh.palpha;
	cuBReal grel = *cuMesh.pgrel;

	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.palpha, alpha, *cuMesh.pgrel, grel);
				cuReal3 rhs = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ Heff[idx]) + alpha * ((M[idx] / Ms) ^ (M[idx] ^ Heff[idx])));

				//4th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				M[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - prediction) / Ms;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);
				M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void DifferentialEquationFMCUDA::RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh || calculate_dmdt) {

		RunEuler_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunEuler_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//TRAPEZOIDAL EULER

void DifferentialEquationFMCUDA::RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunTEuler_Step0_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunTEuler_Step1_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
}

//ADAPTIVE HEUN

void DifferentialEquationFMCUDA::RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunTEuler_Step0_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunAHeun_Step1_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunAHeun_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
}

//RUNGE KUTTA 4th order

void DifferentialEquationFMCUDA::RunRK4_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRK4_Step0_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRK4_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;

	case 1:

		RunRK4_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRK4_Step2_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		if (calculate_dmdt) {

			RunRK4_Step3_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRK4_Step3_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

//Adams-Bashforth-Moulton 2nd order

void DifferentialEquationFMCUDA::RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunABM_Predictor_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunABM_Predictor_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunABM_Corrector_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunABM_Corrector_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
	}
}

//Adams-Bashforth-Moulton 2nd order priming using Trapezoidal Euler

void DifferentialEquationFMCUDA::RunABMTEuler_LLG(int step)
{
	if (step == 0) {

		RunABMTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABMTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//RUNGE KUTTA FEHLBERG

void DifferentialEquationFMCUDA::RunRKF45_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKF45_Step0_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF45_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

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

		if (calculate_dmdt) {

			RunRKF45_Step5_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRKF45_Step5_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

#endif
#endif