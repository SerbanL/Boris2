#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers. LLG equation in-lined for faster evaluation

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuReal3 mxh = cuReal3();
	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//obtain average normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = (M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);
				include_in_average = true;

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
				M2[idx] += rhs2 * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}
				
				//obtain maximum normalized dmdt term
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt (and mxh) if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints, include_in_average);
		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2, include_in_average);
	}
}

__global__ void RunEuler_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
				M2[idx] += rhs2 * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunTEuler_Step0_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuReal3 mxh = cuReal3();
	bool include_in_average = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//obtain average normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = (M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);
				include_in_average = true;

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
				M2[idx] += rhs2 * dT;
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints, include_in_average);
	}
}

__global__ void RunTEuler_Step0_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += rhs * dT;
				M2[idx] += rhs2 * dT;
			}
		}
	}
}

__global__ void RunTEuler_Step1_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;
				M2[idx] = ((*cuDiffEq.psM1_2)[idx] + M2[idx] + rhs2 * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
				include_in_average = true;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2, include_in_average);
	}
}

__global__ void RunTEuler_Step1_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;
				M2[idx] = ((*cuDiffEq.psM1_2)[idx] + M2[idx] + rhs2 * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : Adaptive Heun

//Step0 same as for TEuler

__global__ void RunAHeun_Step1_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal lte = 0.0;
	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;
				M2[idx] = ((*cuDiffEq.psM1_2)[idx] + M2[idx] + rhs2 * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
				include_in_average = true;

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - saveM) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2, include_in_average);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunAHeun_Step1_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = (*cuMesh.pM)[idx];

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;
				M2[idx] = ((*cuDiffEq.psM1_2)[idx] + M2[idx] + rhs2 * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuMesh.pM)[idx] - saveM) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal mxh = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				(*cuDiffEq.psEval0_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//Now estimate magnetization using RK4 midle step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
				M2[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 2);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRK4_Step0_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);
				
				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				(*cuDiffEq.psEval0_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using RK4 midle step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
				M2[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK4_Step1_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);
			
			(*cuDiffEq.psEval1)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
			(*cuDiffEq.psEval1_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

			//Now estimate magnetization using RK4 midle step
			M[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval1)[idx] * (dT / 2);
			M2[idx] = (*cuDiffEq.psM1_2)[idx] + (*cuDiffEq.psEval1_2)[idx] * (dT / 2);
		}
	}
}

__global__ void RunRK4_Step2_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);
			
			(*cuDiffEq.psEval2)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
			(*cuDiffEq.psEval2_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

			//Now estimate magnetization using RK4 last step
			M[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval2)[idx] * dT;
			M2[idx] = (*cuDiffEq.psM1_2)[idx] + (*cuDiffEq.psEval2_2)[idx] * dT;

		}
	}
}

__global__ void RunRK4_Step3_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal dmdt = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using previous RK4 evaluations
				M[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);
				M2[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] + 2 * (*cuDiffEq.psEval1_2)[idx] + 2 * (*cuDiffEq.psEval2_2)[idx] + rhs2) * (dT / 6);

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}
}

__global__ void RunRK4_Step3_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);
				
				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using previous RK4 evaluations
				M[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);
				M2[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] + 2 * (*cuDiffEq.psEval1_2)[idx] + 2 * (*cuDiffEq.psEval2_2)[idx] + rhs2) * (dT / 6);

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : ABM

__global__ void RunABM_Predictor_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal mxh = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuDiffEq.palternator) {

					M[idx] += dT * (3 * rhs - (*cuDiffEq.psEval0)[idx]) / 2;
					M2[idx] += dT * (3 * rhs2 - (*cuDiffEq.psEval0_2)[idx]) / 2;

					(*cuDiffEq.psEval1)[idx] = rhs;
					(*cuDiffEq.psEval1_2)[idx] = rhs2;
				}
				else {

					M[idx] += dT * (3 * rhs - (*cuDiffEq.psEval1)[idx]) / 2;
					M2[idx] += dT * (3 * rhs2 - (*cuDiffEq.psEval1_2)[idx]) / 2;

					(*cuDiffEq.psEval0)[idx] = rhs;
					(*cuDiffEq.psEval0_2)[idx] = rhs2;
				}
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunABM_Predictor_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuDiffEq.palternator) {

					M[idx] += dT * (3 * rhs - (*cuDiffEq.psEval0)[idx]) / 2;
					M2[idx] += dT * (3 * rhs2 - (*cuDiffEq.psEval0_2)[idx]) / 2;

					(*cuDiffEq.psEval1)[idx] = rhs;
					(*cuDiffEq.psEval1_2)[idx] = rhs2;
				}
				else {

					M[idx] += dT * (3 * rhs - (*cuDiffEq.psEval1)[idx]) / 2;
					M2[idx] += dT * (3 * rhs2 - (*cuDiffEq.psEval1_2)[idx]) / 2;

					(*cuDiffEq.psEval0)[idx] = rhs;
					(*cuDiffEq.psEval0_2)[idx] = rhs2;
				}
			}
		}
	}
}

__global__ void RunABM_Corrector_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = M[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuDiffEq.palternator) {

					M[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval1)[idx]) / 2;
					M2[idx] = (*cuDiffEq.psM1_2)[idx] + dT * (rhs2 + (*cuDiffEq.psEval1_2)[idx]) / 2;
				}
				else {

					M[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval0)[idx]) / 2;
					M2[idx] = (*cuDiffEq.psM1_2)[idx] + dT * (rhs2 + (*cuDiffEq.psEval0_2)[idx]) / 2;
				}

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);		

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - saveM) / Mnorm;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunABM_Corrector_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//First save predicted magnetization for lte calculation
				cuReal3 saveM = M[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuDiffEq.palternator) {

					M[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval1)[idx]) / 2;
					M2[idx] = (*cuDiffEq.psM1_2)[idx] + dT * (rhs2 + (*cuDiffEq.psEval1_2)[idx]) / 2;
				}
				else {

					M[idx] = (*cuDiffEq.psM1)[idx] + dT * (rhs + (*cuDiffEq.psEval0)[idx]) / 2;
					M2[idx] = (*cuDiffEq.psM1_2)[idx] + dT * (rhs2 + (*cuDiffEq.psEval0_2)[idx]) / 2;
				}

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - saveM) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunABMTEuler_Step0_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				(*cuDiffEq.psEval0_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization for the next time step
				M[idx] += (*cuDiffEq.psEval0)[idx] * dT;
				M2[idx] += (*cuDiffEq.psEval0_2)[idx] * dT;
			}
		}
	}
}

__global__ void RunABMTEuler_Step1_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using the second trapezoidal Euler step equation
				M[idx] = ((*cuDiffEq.psM1)[idx] + M[idx] + rhs * dT) / 2;
				M2[idx] = ((*cuDiffEq.psM1_2)[idx] + M2[idx] + rhs2 * dT) / 2;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF45_Step0_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal mxh = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				(*cuDiffEq.psEval0_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M[idx].norm();
				mxh = cu_GetMagnitude(M[idx] ^ Heff[idx]) / (Mnorm * Mnorm);

				//Now estimate magnetization using RKF first step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
				M2[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 4);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRKF45_Step0_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = M[idx];
			(*cuDiffEq.psM1_2)[idx] = M2[idx];

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				(*cuDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				(*cuDiffEq.psEval0_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//Now estimate magnetization using RKF first step
				M[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 4);
				M2[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 4);
			}
		}
	}
}

__global__ void RunRKF45_Step1_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

			(*cuDiffEq.psEval1)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
			(*cuDiffEq.psEval1_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

			//Now estimate magnetization using RKF midle step 1
			M[idx] = (*cuDiffEq.psM1)[idx] + (3 * (*cuDiffEq.psEval0)[idx] + 9 * (*cuDiffEq.psEval1)[idx]) * dT / 32;
			M2[idx] = (*cuDiffEq.psM1_2)[idx] + (3 * (*cuDiffEq.psEval0_2)[idx] + 9 * (*cuDiffEq.psEval1_2)[idx]) * dT / 32;
		}
	}
}

__global__ void RunRKF45_Step2_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

			(*cuDiffEq.psEval2)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
			(*cuDiffEq.psEval2_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

			//Now estimate magnetization using RKF midle step 2
			M[idx] = (*cuDiffEq.psM1)[idx] + (1932 * (*cuDiffEq.psEval0)[idx] - 7200 * (*cuDiffEq.psEval1)[idx] + 7296 * (*cuDiffEq.psEval2)[idx]) * dT / 2197;
			M2[idx] = (*cuDiffEq.psM1_2)[idx] + (1932 * (*cuDiffEq.psEval0_2)[idx] - 7200 * (*cuDiffEq.psEval1_2)[idx] + 7296 * (*cuDiffEq.psEval2_2)[idx]) * dT / 2197;
		}
	}
}

__global__ void RunRKF45_Step3_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

			(*cuDiffEq.psEval3)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
			(*cuDiffEq.psEval3_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

			//Now estimate magnetization using RKF midle step 3
			M[idx] = (*cuDiffEq.psM1)[idx] + (439 * (*cuDiffEq.psEval0)[idx] / 216 - 8 * (*cuDiffEq.psEval1)[idx] + 3680 * (*cuDiffEq.psEval2)[idx] / 513 - 845 * (*cuDiffEq.psEval3)[idx] / 4104) * dT;
			M2[idx] = (*cuDiffEq.psM1_2)[idx] + (439 * (*cuDiffEq.psEval0_2)[idx] / 216 - 8 * (*cuDiffEq.psEval1_2)[idx] + 3680 * (*cuDiffEq.psEval2_2)[idx] / 513 - 845 * (*cuDiffEq.psEval3_2)[idx] / 4104) * dT;
		}
	}
}

__global__ void RunRKF45_Step4_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && !M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

			(*cuDiffEq.psEval4)[idx] = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
			(*cuDiffEq.psEval4_2)[idx] = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

			//Now estimate magnetization using RKF midle step 4
			M[idx] = (*cuDiffEq.psM1)[idx] + (-8 * (*cuDiffEq.psEval0)[idx] / 27 + 2 * (*cuDiffEq.psEval1)[idx] - 3544 * (*cuDiffEq.psEval2)[idx] / 2565 + 1859 * (*cuDiffEq.psEval3)[idx] / 4104 - 11 * (*cuDiffEq.psEval4)[idx] / 40) * dT;
			M2[idx] = (*cuDiffEq.psM1_2)[idx] + (-8 * (*cuDiffEq.psEval0_2)[idx] / 27 + 2 * (*cuDiffEq.psEval1_2)[idx] - 3544 * (*cuDiffEq.psEval2_2)[idx] / 2565 + 1859 * (*cuDiffEq.psEval3_2)[idx] / 4104 - 11 * (*cuDiffEq.psEval4_2)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKF45_Step5_LLG_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//4th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				M[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;
				M2[idx] = (*cuDiffEq.psM1_2)[idx] + (16 * (*cuDiffEq.psEval0_2)[idx] / 135 + 6656 * (*cuDiffEq.psEval2_2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3_2)[idx] / 56430 - 9 * (*cuDiffEq.psEval4_2)[idx] / 50 + 2 * rhs2 / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - prediction) / Mnorm;
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

__global__ void RunRKF45_Step5_LLG_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
	cuReal2 alpha_AFM = *cuMesh.palpha_AFM;
	cuReal2 grel_AFM = *cuMesh.pgrel_AFM;

	cuBReal lte = 0.0;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			if (!M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.palpha_AFM, alpha_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 rhs = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
				cuReal3 rhs2 = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

				//4th order evaluation
				cuReal3 prediction = (*cuDiffEq.psM1)[idx] + (25 * (*cuDiffEq.psEval0)[idx] / 216 + 1408 * (*cuDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuDiffEq.psEval3)[idx] / 4101 - (*cuDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				M[idx] = (*cuDiffEq.psM1)[idx] + (16 * (*cuDiffEq.psEval0)[idx] / 135 + 6656 * (*cuDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3)[idx] / 56430 - 9 * (*cuDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;
				M2[idx] = (*cuDiffEq.psM1_2)[idx] + (16 * (*cuDiffEq.psEval0_2)[idx] / 135 + 6656 * (*cuDiffEq.psEval2_2)[idx] / 12825 + 28561 * (*cuDiffEq.psEval3_2)[idx] / 56430 - 9 * (*cuDiffEq.psEval4_2)[idx] / 50 + 2 * rhs2 / 55) * dT;

				if (*cuDiffEq.prenormalize) {

					M[idx].renormalize(Ms_AFM.i);
					M2[idx].renormalize(Ms_AFM.j);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M[idx] - prediction) / (*cuMesh.pM)[idx].norm();
			}
			else {

				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				M[idx].renormalize(Ms_AFM.i);
				M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	reduction_max(0, 1, &lte, *cuDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void DifferentialEquationAFMCUDA::RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh || calculate_dmdt) {

		RunEuler_LLG_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunEuler_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//TRAPEZOIDAL EULER

void DifferentialEquationAFMCUDA::RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
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

void DifferentialEquationAFMCUDA::RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
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

void DifferentialEquationAFMCUDA::RunRK4_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
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

void DifferentialEquationAFMCUDA::RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
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

void DifferentialEquationAFMCUDA::RunABMTEuler_LLG(int step)
{
	if (step == 0) {

		RunABMTEuler_Step0_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunABMTEuler_Step1_LLG_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >(cuDiffEq, pMeshCUDA->cuMesh);
	}
}

//RUNGE KUTTA FEHLBERG

void DifferentialEquationAFMCUDA::RunRKF45_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
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
