#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_RK4
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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

				(*cuaDiffEq.psEval0)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using RK4 midle step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunRK4_Step0_withAverageReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuReal3 mxh = cuReal3();
	bool include_in_average = false;

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

				(*cuaDiffEq.psEval0)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using RK4 midle step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_avg(0, 1, &mxh, *cuaDiffEq.pmxh_av, *cuaDiffEq.pavpoints, include_in_average);
}

__global__ void RunRK4_Step0_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				(*cuaDiffEq.psEval0)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using RK4 midle step
				(*cuaMesh.pM1)[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK4_Step1_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			(*cuaDiffEq.psEval1)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RK4 midle step
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (*cuaDiffEq.psEval1)[idx] * (dT / 2);
		}
	}
}

__global__ void RunRK4_Step2_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx) && !cuaMesh.pM1->is_skipcell(idx)) {

			(*cuaDiffEq.psEval2)[idx] = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

			//Now estimate moment using RK4 last step
			(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + (*cuaDiffEq.psEval2)[idx] * dT;
		}
	}
}

__global__ void RunRK4_Step3_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
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
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using previous RK4 evaluations
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] + 2 * (*cuaDiffEq.psEval1)[idx] + 2 * (*cuaDiffEq.psEval2)[idx] + rhs) * (dT / 6);

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

__global__ void RunRK4_Step3_withAverageReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using previous RK4 evaluations
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] + 2 * (*cuaDiffEq.psEval1)[idx] + 2 * (*cuaDiffEq.psEval2)[idx] + rhs) * (dT / 6);

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = ((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
				include_in_average = true;
			}
		}
	}

	if (cuaMesh.pgrel->get0()) reduction_avg(0, 1, &dmdt, *cuaDiffEq.pdmdt_av, *cuaDiffEq.pavpoints2, include_in_average);
}

__global__ void RunRK4_Step3_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment using previous RK4 evaluations
				(*cuaMesh.pM1)[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] + 2 * (*cuaDiffEq.psEval1)[idx] + 2 * (*cuaDiffEq.psEval2)[idx] + rhs) * (dT / 6);

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

//RUNGE KUTTA 4th order

void Atom_DifferentialEquationCubicCUDA::RunRK4(int step, bool calculate_mxh, bool calculate_dmdt, bool stochastic)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			if (stochastic) RunRK4_Step0_withAverageReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
			else RunRK4_Step0_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRK4_Step0_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;

	case 1:

		RunRK4_Step1_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRK4_Step2_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		if (calculate_dmdt) {

			if (stochastic) RunRK4_Step3_withAverageReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
			else RunRK4_Step3_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRK4_Step3_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif
#endif