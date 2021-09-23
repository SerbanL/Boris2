#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKF56
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF56_Step0_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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

				//Now estimate moment using RKF first step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 6);
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunRKF56_Step0_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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

				//Now estimate moment using RKF first step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 6);
			}
		}
	}
}

__global__ void RunRKF56_Step1_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval1)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 1
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (4 * (*cuaDiffEq.psEval0)[idx] + 16 * (*cuaDiffEq.psEval1)[idx]) * dT / 75;
		}
	}
}

__global__ void RunRKF56_Step2_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 2
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (5 * (*cuaDiffEq.psEval0)[idx] / 6 - 8 * (*cuaDiffEq.psEval1)[idx] / 3 + 5 * (*cuaDiffEq.psEval2)[idx] / 2) * dT;
		}
	}
}

__global__ void RunRKF56_Step3_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval3)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 3
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (-8 * (*cuaDiffEq.psEval0)[idx] / 5 + 144 * (*cuaDiffEq.psEval1)[idx] / 25 - 4 * (*cuaDiffEq.psEval2)[idx] + 16 * (*cuaDiffEq.psEval3)[idx] / 25) * dT;
		}
	}
}

__global__ void RunRKF56_Step4_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval4)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 4
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (361 * (*cuaDiffEq.psEval0)[idx] / 320 - 18 * (*cuaDiffEq.psEval1)[idx] / 5 + 407 * (*cuaDiffEq.psEval2)[idx] / 128 - 11 * (*cuaDiffEq.psEval3)[idx] / 80 + 55 * (*cuaDiffEq.psEval4)[idx] / 128) * dT;
		}
	}
}

__global__ void RunRKF56_Step5_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval5)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (-11 * (*cuaDiffEq.psEval0)[idx] / 640 + 11 * (*cuaDiffEq.psEval2)[idx] / 256 - 11 * (*cuaDiffEq.psEval3)[idx] / 160 + 11 * (*cuaDiffEq.psEval4)[idx] / 256) * dT;
		}
	}
}

__global__ void RunRKF56_Step6_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval6)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (93 * (*cuaDiffEq.psEval0)[idx] / 640 - 18 * (*cuaDiffEq.psEval1)[idx] / 5 + 803 * (*cuaDiffEq.psEval2)[idx] / 256 - 11 * (*cuaDiffEq.psEval3)[idx] / 160 + 99 * (*cuaDiffEq.psEval4)[idx] / 256 + (*cuaDiffEq.psEval6)[idx]) * dT;
		}
	}
}

__global__ void RunRKF56_Step7_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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

				//5th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (31 * (*cuaDiffEq.psEval0)[idx] / 384 + 1125 * (*cuaDiffEq.psEval2)[idx] / 2816 + 9 * (*cuaDiffEq.psEval3)[idx] / 32 + 125 * (*cuaDiffEq.psEval4)[idx] / 768 + 5 * (*cuaDiffEq.psEval5)[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				cuReal3 lte_diff = 5 * ((*cuaDiffEq.psEval0)[idx] + (*cuaDiffEq.psEval5)[idx] - (*cuaDiffEq.psEval6)[idx] - rhs) * dT / 66;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(lte_diff) / (*cuaMesh.pM1)[idx].norm();
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRKF56_Step7_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//5th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (31 * (*cuaDiffEq.psEval0)[idx] / 384 + 1125 * (*cuaDiffEq.psEval2)[idx] / 2816 + 9 * (*cuaDiffEq.psEval3)[idx] / 32 + 125 * (*cuaDiffEq.psEval4)[idx] / 768 + 5 * (*cuaDiffEq.psEval5)[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				cuReal3 lte_diff = 5 * ((*cuaDiffEq.psEval0)[idx] + (*cuaDiffEq.psEval5)[idx] - (*cuaDiffEq.psEval6)[idx] - rhs) * dT / 66;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(lte_diff) / (*cuaMesh.pM1)[idx].norm();
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//RUNGE KUTTA FEHLBERG 5(6)

void Atom_DifferentialEquationCubicCUDA::RunRKF56(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKF56_Step0_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKF56_Step0_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;

	case 1:

		RunRKF56_Step1_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRKF56_Step2_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		RunRKF56_Step3_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 4:

		RunRKF56_Step4_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 5:

		RunRKF56_Step5_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 6:

		RunRKF56_Step6_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 7:

		if (calculate_dmdt) {

			RunRKF56_Step7_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKF56_Step7_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif
#endif