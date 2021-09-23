#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RKF45
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF45_Step0_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (2 * dT / 9);
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunRKF45_Step0_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (2 * dT / 9);
			}
		}
	}
}

__global__ void RunRKF45_Step1_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval1)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 1
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] / 12 + (*cuaDiffEq.psEval1)[idx] / 4) * dT;
		}
	}
}

__global__ void RunRKF45_Step2_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 2
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (69 * (*cuaDiffEq.psEval0)[idx] / 128 - 243 * (*cuaDiffEq.psEval1)[idx] / 128 + 135 * (*cuaDiffEq.psEval2)[idx] / 64) * dT;
		}
	}
}

__global__ void RunRKF45_Step3_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval3)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 3
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (-17 * (*cuaDiffEq.psEval0)[idx] / 12 + 27 * (*cuaDiffEq.psEval1)[idx] / 4 - 27 * (*cuaDiffEq.psEval2)[idx] / 5 + 16 * (*cuaDiffEq.psEval3)[idx] / 15) * dT;
		}
	}
}

__global__ void RunRKF45_Step4_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval4)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RKF midle step 4
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (65 * (*cuaDiffEq.psEval0)[idx] / 432 - 5 * (*cuaDiffEq.psEval1)[idx] / 16 + 13 * (*cuaDiffEq.psEval2)[idx] / 16 + 4 * (*cuaDiffEq.psEval3)[idx] / 27 + 5 * (*cuaDiffEq.psEval4)[idx] / 144) * dT;
		}
	}
}

__global__ void RunRKF45_Step5_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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

				//4th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] / 9 + 9 * (*cuaDiffEq.psEval2)[idx] / 20 + 16 * (*cuaDiffEq.psEval3)[idx] / 45 + (*cuaDiffEq.psEval4)[idx] / 12) * dT;

				//5th order evaluation
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (47 * (*cuaDiffEq.psEval0)[idx] / 450 + 12 * (*cuaDiffEq.psEval2)[idx] / 25 + 32 * (*cuaDiffEq.psEval3)[idx] / 225 + 1 * (*cuaDiffEq.psEval4)[idx] / 30 + 6 * rhs / 25) * dT;

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

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRKF45_Step5_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//4th order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] / 9 + 9 * (*cuaDiffEq.psEval2)[idx] / 20 + 16 * (*cuaDiffEq.psEval3)[idx] / 45 + (*cuaDiffEq.psEval4)[idx] / 12) * dT;

				//5th order evaluation
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (47 * (*cuaDiffEq.psEval0)[idx] / 450 + 12 * (*cuaDiffEq.psEval2)[idx] / 25 + 32 * (*cuaDiffEq.psEval3)[idx] / 225 + 1 * (*cuaDiffEq.psEval4)[idx] / 30 + 6 * rhs / 25) * dT;
				
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

//RUNGE KUTTA FEHLBERG

void Atom_DifferentialEquationCubicCUDA::RunRKF45(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKF45_Step0_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKF45_Step0_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;

	case 1:

		RunRKF45_Step1_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRKF45_Step2_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		RunRKF45_Step3_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 4:

		RunRKF45_Step4_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 5:

		if (calculate_dmdt) {

			RunRKF45_Step5_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKF45_Step5_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif
#endif