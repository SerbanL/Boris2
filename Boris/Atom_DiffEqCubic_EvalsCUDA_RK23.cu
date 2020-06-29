#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RK23
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RK23

__global__ void RunRK23_Step0_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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

				//2nd order evaluation for adaptive step
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (7 * (*cuaDiffEq.psEval0)[idx] / 24 + 1 * (*cuaDiffEq.psEval1)[idx] / 4 + 1 * (*cuaDiffEq.psEval2)[idx] / 3 + 1 * rhs / 8) * dT;

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - prediction) / Mnorm;

				//save evaluation for later use
				(*cuaDiffEq.psEval0)[idx] = rhs;
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRK23_Step0_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal lte = 0.0;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//2nd order evaluation for adaptive step
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (7 * (*cuaDiffEq.psEval0)[idx] / 24 + 1 * (*cuaDiffEq.psEval1)[idx] / 4 + 1 * (*cuaDiffEq.psEval2)[idx] / 3 + 1 * rhs / 8) * dT;

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - prediction) / (*cuaMesh.pM1)[idx].norm();

				//save evaluation for later use
				(*cuaDiffEq.psEval0)[idx] = rhs;
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRK23_Step0_Advance_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//Save current moment for later use
				(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

				//Now estimate moment using RK23 first step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK23_Step1_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			(*cuaDiffEq.psEval1)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RK23 midle step 1
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + 3 * (*cuaDiffEq.psEval1)[idx] * dT / 4;
		}
	}
}

__global__ void RunRK23_Step2_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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
				(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now calculate 3rd order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (2 * (*cuaDiffEq.psEval0)[idx] / 9 + 1 * (*cuaDiffEq.psEval1)[idx] / 3 + 4 * (*cuaDiffEq.psEval2)[idx] / 9) * dT;

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

	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
}

__global__ void RunRK23_Step2_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now calculate 3rd order evaluation
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (2 * (*cuaDiffEq.psEval0)[idx] / 9 + 1 * (*cuaDiffEq.psEval1)[idx] / 3 + 4 * (*cuaDiffEq.psEval2)[idx] / 9) * dT;

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

//RUNGE KUTTA 23 (Bogacki - Shampine) (2nd order adaptive step with FSAL, 3rd order evaluation)

void Atom_DifferentialEquationCubicCUDA::RunRK23_Step0_NoAdvance(bool calculate_mxh)
{
	if (calculate_mxh) {

		RunRK23_Step0_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else {

		RunRK23_Step0_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
}

void Atom_DifferentialEquationCubicCUDA::RunRK23(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		RunRK23_Step0_Advance_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 1:

		RunRK23_Step1_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		if (calculate_dmdt) {

			RunRK23_Step2_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRK23_Step2_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif
#endif