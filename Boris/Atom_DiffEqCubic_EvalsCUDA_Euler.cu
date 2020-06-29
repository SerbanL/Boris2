#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_EULER
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_Kernel_withReductions(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuReal3 mxh = cuReal3();
	cuReal3 dmdt = cuReal3();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			//Save current moment
			(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//obtain average normalized torque term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				mxh = ((*cuaMesh.pM1)[idx] ^ (*cuaMesh.pHeff1)[idx]) / (conversion * Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment for the next time step
				(*cuaMesh.pM1)[idx] += rhs * dT;

				if (*cuaDiffEq.prenormalize) {

					cuBReal mu_s = *cuaMesh.pmu_s;
					cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);
					(*cuaMesh.pM1)[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				dmdt = ((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
			}
		}
	}

	reduction_avg(0, 1, &mxh, *cuaDiffEq.pmxh_av, *cuaDiffEq.pavpoints);
	reduction_avg(0, 1, &dmdt, *cuaDiffEq.pdmdt_av, *cuaDiffEq.pavpoints2);
}

__global__ void RunEuler_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			//Save current moment
			(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs = (cuaDiffEq.*(cuaDiffEq.pODEFunc))(idx);

				//Now estimate moment for the next time step
				(*cuaMesh.pM1)[idx] += rhs * dT;

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

//EULER

void Atom_DifferentialEquationCubicCUDA::RunEuler(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh || calculate_dmdt) {

		RunEuler_Kernel_withReductions <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else {

		RunEuler_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
}

#endif
#endif
#endif