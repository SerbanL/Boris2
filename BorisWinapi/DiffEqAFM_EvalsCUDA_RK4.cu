#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_RK4
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
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

				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval0_2)[idx]);

				//Now estimate magnetization using RK4 midle step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
				(*cuMesh.pM2)[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 2);
			}
		}
	}

	//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	}
}

__global__ void RunRK4_Step0_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization for later use
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
			(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				(*cuDiffEq.psEval0)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval0_2)[idx]);

				//Now estimate magnetization using RK4 midle step
				(*cuMesh.pM)[idx] += (*cuDiffEq.psEval0)[idx] * (dT / 2);
				(*cuMesh.pM2)[idx] += (*cuDiffEq.psEval0_2)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK4_Step1_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			(*cuDiffEq.psEval1)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval1_2)[idx]);

			//Now estimate magnetization using RK4 midle step
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval1)[idx] * (dT / 2);
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (*cuDiffEq.psEval1_2)[idx] * (dT / 2);
		}
	}
}

__global__ void RunRK4_Step2_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			(*cuDiffEq.psEval2)[idx] = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, (*cuDiffEq.psEval2_2)[idx]);

			//Now estimate magnetization using RK4 last step
			(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + (*cuDiffEq.psEval2)[idx] * dT;
			(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + (*cuDiffEq.psEval2_2)[idx] * dT;
		}
	}
}

__global__ void RunRK4_Step3_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal dmdt = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs2;
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, rhs2);

				//Now estimate magnetization using previous RK4 evaluations
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);
				(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] + 2 * (*cuDiffEq.psEval1_2)[idx] + 2 * (*cuDiffEq.psEval2_2)[idx] + rhs2) * (dT / 6);

				if (*cuDiffEq.prenormalize) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
					(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
					(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
	}
}

__global__ void RunRK4_Step3_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs2;
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, rhs2);

				//Now estimate magnetization using previous RK4 evaluations
				(*cuMesh.pM)[idx] = (*cuDiffEq.psM1)[idx] + ((*cuDiffEq.psEval0)[idx] + 2 * (*cuDiffEq.psEval1)[idx] + 2 * (*cuDiffEq.psEval2)[idx] + rhs) * (dT / 6);
				(*cuMesh.pM2)[idx] = (*cuDiffEq.psM1_2)[idx] + ((*cuDiffEq.psEval0_2)[idx] + 2 * (*cuDiffEq.psEval1_2)[idx] + 2 * (*cuDiffEq.psEval2_2)[idx] + rhs2) * (dT / 6);

				if (*cuDiffEq.prenormalize) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
					(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
					(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA 4th order

void DifferentialEquationAFMCUDA::RunRK4(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRK4_Step0_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRK4_Step0_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;

	case 1:

		RunRK4_Step1_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 2:

		RunRK4_Step2_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);

		break;

	case 3:

		if (calculate_dmdt) {

			RunRK4_Step3_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}
		else {

			RunRK4_Step3_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
		}

		break;
	}
}

#endif
#endif
#endif