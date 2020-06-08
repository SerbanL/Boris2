#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_RKCK
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKCK45

__global__ void RunRKCK45_Step0_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal mxh = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//obtain maximum normalized torque term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				mxh = cu_GetMagnitude((*cuaMesh.pM1)[idx] ^ (*cuaMesh.pHeff1)[idx]) / (conversion * Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				(*cuaDiffEq.psEval0)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using RKCK first step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 5);
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunRKCK45_Step0_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuaDiffEq.psEval0)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using RKCK first step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 5);
			}
		}
	}
}

__global__ void RunRKCK45_Step1_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval1)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKCK midle step 1
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (3 * (*cuaDiffEq.psEval0)[idx] / 40 + 9 * (*cuaDiffEq.psEval1)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKCK45_Step2_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKCK midle step 2
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (3 * (*cuaDiffEq.psEval0)[idx] / 10 - 9 * (*cuaDiffEq.psEval1)[idx] / 10 + 6 * (*cuaDiffEq.psEval2)[idx] / 5) * dT;
		}
	}
}

__global__ void RunRKCK45_Step3_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval3)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKCK midle step 3
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (-11 * (*cuaDiffEq.psEval0)[idx] / 54 + 5 * (*cuaDiffEq.psEval1)[idx] / 2 - 70 * (*cuaDiffEq.psEval2)[idx] / 27 + 35 * (*cuaDiffEq.psEval3)[idx] / 27) * dT;
		}
	}
}

__global__ void RunRKCK45_Step4_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval4)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 4
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (1631 * (*cuaDiffEq.psEval0)[idx] / 55296 + 175 * (*cuaDiffEq.psEval1)[idx] / 512 + 575 * (*cuaDiffEq.psEval2)[idx] / 13824 + 44275 * (*cuaDiffEq.psEval3)[idx] / 110592 + 253 * (*cuaDiffEq.psEval4)[idx] / 4096) * dT;
		}
	}
}

__global__ void RunRKCK45_Step5_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal dmdt = 0.0;
	cuBReal lte = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//RKCK45 : 4th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (2825 * (*cuaDiffEq.psEval0)[idx] / 27648 + 18575 * (*cuaDiffEq.psEval2)[idx] / 48384 + 13525 * (*cuaDiffEq.psEval3)[idx] / 55296 + 277 * (*cuaDiffEq.psEval4)[idx] / 14336 + rhs / 4) * dT;

				//Now calculate 5th order evaluation for adaptive time step
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (37 * (*cuaDiffEq.psEval0)[idx] / 378 + 250 * (*cuaDiffEq.psEval2)[idx] / 621 + 125 * (*cuaDiffEq.psEval3)[idx] / 594 + 512 * rhs / 1771) * dT;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - prediction) / (*cuaMesh.pM1)[idx].norm();
			}
		}
	}

	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRKCK45_Step5_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//RKCK45 : 4th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (2825 * (*cuaDiffEq.psEval0)[idx] / 27648 + 18575 * (*cuaDiffEq.psEval2)[idx] / 48384 + 13525 * (*cuaDiffEq.psEval3)[idx] / 55296 + 277 * (*cuaDiffEq.psEval4)[idx] / 14336 + rhs / 4) * dT;

				//Now calculate 5th order evaluation for adaptive time step
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (37 * (*cuaDiffEq.psEval0)[idx] / 378 + 250 * (*cuaDiffEq.psEval2)[idx] / 621 + 125 * (*cuaDiffEq.psEval3)[idx] / 594 + 512 * rhs / 1771) * dT;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - prediction) / (*cuaMesh.pM1)[idx].norm();
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA CASH-KARP

void Atom_DifferentialEquationCubicCUDA::RunRKCK45(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKCK45_Step0_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKCK45_Step0_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;

	case 1:

		RunRKCK45_Step1_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRKCK45_Step2_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		RunRKCK45_Step3_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 4:

		RunRKCK45_Step4_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 5:

		if (calculate_dmdt) {

			RunRKCK45_Step5_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKCK45_Step5_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif
#endif