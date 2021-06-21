#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKDP
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKDP54

__global__ void RunRKDP54_Step0_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal mxh = 0.0;
	cuBReal lte = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//obtain maximum normalized torque term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				mxh = cu_GetMagnitude((*cuaMesh.pM1)[idx] ^ (*cuaMesh.pHeff1)[idx]) / (conversion * Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now calculate 5th order evaluation for adaptive time step -> FSAL property (a full pass required for this to be valid)
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (5179 * (*cuaDiffEq.psEval0)[idx] / 57600 + 7571 * (*cuaDiffEq.psEval2)[idx] / 16695 + 393 * (*cuaDiffEq.psEval3)[idx] / 640 - 92097 * (*cuaDiffEq.psEval4)[idx] / 339200 + 187 * (*cuaDiffEq.psEval5)[idx] / 2100 + rhs / 40) * dT;

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - prediction) / Mnorm;

				//save evaluation for later use
				(*cuaDiffEq.psEval0)[idx] = rhs;
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRKDP54_Step0_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now calculate 5th order evaluation for adaptive time step -> FSAL property (a full pass required for this to be valid)
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (5179 * (*cuaDiffEq.psEval0)[idx] / 57600 + 7571 * (*cuaDiffEq.psEval2)[idx] / 16695 + 393 * (*cuaDiffEq.psEval3)[idx] / 640 - 92097 * (*cuaDiffEq.psEval4)[idx] / 339200 + 187 * (*cuaDiffEq.psEval5)[idx] / 2100 + rhs / 40) * dT;

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - prediction) / (*cuaMesh.pM1)[idx].norm();

				//save evaluation for later use
				(*cuaDiffEq.psEval0)[idx] = rhs;
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRKDP54_Step0_Advance_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//Save current moment for later use
				(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

				//Now estimate moment using RKDP first step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 5);
			}
		}
	}
}

__global__ void RunRKDP54_Step1_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval1)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKDP midle step 1
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (3 * (*cuaDiffEq.psEval0)[idx] / 40 + 9 * (*cuaDiffEq.psEval1)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKDP54_Step2_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKDP midle step 2
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (44 * (*cuaDiffEq.psEval0)[idx] / 45 - 56 * (*cuaDiffEq.psEval1)[idx] / 15 + 32 * (*cuaDiffEq.psEval2)[idx] / 9) * dT;
		}
	}
}

__global__ void RunRKDP54_Step3_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval3)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKDP midle step 3
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (19372 * (*cuaDiffEq.psEval0)[idx] / 6561 - 25360 * (*cuaDiffEq.psEval1)[idx] / 2187 + 64448 * (*cuaDiffEq.psEval2)[idx] / 6561 - 212 * (*cuaDiffEq.psEval3)[idx] / 729) * dT;
		}
	}
}

__global__ void RunRKDP54_Step4_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval4)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKDP midle step 4
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (9017 * (*cuaDiffEq.psEval0)[idx] / 3168 - 355 * (*cuaDiffEq.psEval1)[idx] / 33 + 46732 * (*cuaDiffEq.psEval2)[idx] / 5247 + 49 * (*cuaDiffEq.psEval3)[idx] / 176 - 5103 * (*cuaDiffEq.psEval4)[idx] / 18656) * dT;
		}
	}
}

__global__ void RunRKDP54_Step5_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal dmdt = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuaDiffEq.psEval5)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//RKDP54 : 5th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (35 * (*cuaDiffEq.psEval0)[idx] / 384 + 500 * (*cuaDiffEq.psEval2)[idx] / 1113 + 125 * (*cuaDiffEq.psEval3)[idx] / 192 - 2187 * (*cuaDiffEq.psEval4)[idx] / 6784 + 11 * (*cuaDiffEq.psEval5)[idx] / 84) * dT;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
}

__global__ void RunRKDP54_Step5_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuaDiffEq.psEval5)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//RKDP54 : 5th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (35 * (*cuaDiffEq.psEval0)[idx] / 384 + 500 * (*cuaDiffEq.psEval2)[idx] / 1113 + 125 * (*cuaDiffEq.psEval3)[idx] / 192 - 2187 * (*cuaDiffEq.psEval4)[idx] / 6784 + 11 * (*cuaDiffEq.psEval5)[idx] / 84) * dT;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}
			}
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA DORMAND-PRINCE

void Atom_DifferentialEquationCubicCUDA::RunRKDP54_Step0_NoAdvance(bool calculate_mxh)
{
	if (calculate_mxh) {

		RunRKDP54_Step0_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else {

		RunRKDP54_Step0_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
}

void Atom_DifferentialEquationCubicCUDA::RunRKDP54(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		RunRKDP54_Step0_Advance_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 1:

		RunRKDP54_Step1_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRKDP54_Step2_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		RunRKDP54_Step3_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 4:

		RunRKDP54_Step4_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 5:

		if (calculate_dmdt) {

			RunRKDP54_Step5_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKDP54_Step5_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif
#endif