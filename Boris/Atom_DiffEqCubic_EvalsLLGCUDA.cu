#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers. LLG equation in-lined for faster evaluation

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuReal3 mxh = cuReal3();
	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//obtain average normalized torque term
				cuBReal Mnorm = M1[idx].norm();
				mxh = (M1[idx] ^ Heff1[idx]) / (conversion * Mnorm * Mnorm);

				//Now estimate moment for the next time step
				M1[idx] += rhs * dT;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}
				
				//obtain maximum normalized dmdt term
				dmdt = ((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
				include_in_average = true;
			}
		}
	}

	reduction_avg(0, 1, &mxh, *cuaDiffEq.pmxh_av, *cuaDiffEq.pavpoints, include_in_average);
	reduction_avg(0, 1, &dmdt, *cuaDiffEq.pdmdt_av, *cuaDiffEq.pavpoints2, include_in_average);
}

__global__ void RunEuler_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment for the next time step
				M1[idx] += rhs * dT;
				
				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : Trapezoidal Euler

__global__ void RunTEuler_Step0_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuReal3 mxh = cuReal3();
	bool include_in_average = false;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for the next step
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//obtain average normalized torque term
				cuBReal Mnorm = M1[idx].norm();
				mxh = (M1[idx] ^ Heff1[idx]) / (conversion * Mnorm * Mnorm);
				include_in_average = true;

				//Now estimate moment for the next time step
				M1[idx] += rhs * dT;
			}
		}
	}

	reduction_avg(0, 1, &mxh, *cuaDiffEq.pmxh_av, *cuaDiffEq.pavpoints, include_in_average);
}

__global__ void RunTEuler_Step0_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for the next step
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment for the next time step
				M1[idx] += rhs * dT;
			}
		}
	}
}

__global__ void RunTEuler_Step1_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using the second trapezoidal Euler step equation
				M1[idx] = ((*cuaDiffEq.psM1)[idx] + M1[idx] + rhs * dT) / 2;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = ((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
				include_in_average = true;
			}
		}
	}

	reduction_avg(0, 1, &dmdt, *cuaDiffEq.pdmdt_av, *cuaDiffEq.pavpoints2, include_in_average);
}

__global__ void RunTEuler_Step1_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using the second trapezoidal Euler step equation
				M1[idx] = ((*cuaDiffEq.psM1)[idx] + M1[idx] + rhs * dT) / 2;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : Adaptive Heun

//Step0 same as for TEuler

__global__ void RunAHeun_Step1_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuReal3 dmdt = cuReal3();
	bool include_in_average = false;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	cuBReal lte = 0.0;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//First save predicted moment for lte calculation
				cuReal3 saveM = (*cuaMesh.pM1)[idx];

				//Now estimate moment using the second trapezoidal Euler step equation
				M1[idx] = ((*cuaDiffEq.psM1)[idx] + M1[idx] + rhs * dT) / 2;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = ((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
				include_in_average = true;

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - saveM) / (*cuaMesh.pM1)[idx].norm();
			}
		}
	}

	reduction_avg(0, 1, &dmdt, *cuaDiffEq.pdmdt_av, *cuaDiffEq.pavpoints2, include_in_average);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunAHeun_Step1_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal lte = 0.0;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//First save predicted moment for lte calculation
				cuReal3 saveM = (*cuaMesh.pM1)[idx];

				//Now estimate moment using the second trapezoidal Euler step equation
				M1[idx] = ((*cuaDiffEq.psM1)[idx] + M1[idx] + rhs * dT) / 2;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude((*cuaMesh.pM1)[idx] - saveM) / (*cuaMesh.pM1)[idx].norm();
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

//----------------------------------------- EVALUATIONS : RK4

__global__ void RunRK4_Step0_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal mxh = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);

				(*cuaDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M1[idx].norm();
				mxh = cu_GetMagnitude(M1[idx] ^ Heff1[idx]) / (conversion * Mnorm * Mnorm);

				//Now estimate moment using RK4 midle step
				M1[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunRK4_Step0_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				
				(*cuaDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using RK4 midle step
				M1[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 2);
			}
		}
	}
}

__global__ void RunRK4_Step1_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && !M1.is_skipcell(idx)) {

			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
			
			(*cuaDiffEq.psEval1)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

			//Now estimate moment using RK4 midle step
			M1[idx] = (*cuaDiffEq.psM1)[idx] + (*cuaDiffEq.psEval1)[idx] * (dT / 2);
		}
	}
}

__global__ void RunRK4_Step2_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && !M1.is_skipcell(idx)) {

			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
			
			(*cuaDiffEq.psEval2)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

			//Now estimate moment using RK4 last step
			M1[idx] = (*cuaDiffEq.psM1)[idx] + (*cuaDiffEq.psEval2)[idx] * dT;
		}
	}
}

__global__ void RunRK4_Step3_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal dmdt = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);

				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using previous RK4 evaluations
				M1[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] + 2 * (*cuaDiffEq.psEval1)[idx] + 2 * (*cuaDiffEq.psEval2)[idx] + rhs) * (dT / 6);

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
			}
		}
	}

	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
}

__global__ void RunRK4_Step3_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using previous RK4 evaluations
				M1[idx] = (*cuaDiffEq.psM1)[idx] + ((*cuaDiffEq.psEval0)[idx] + 2 * (*cuaDiffEq.psEval1)[idx] + 2 * (*cuaDiffEq.psEval2)[idx] + rhs) * (dT / 6);

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : ABM

__global__ void RunABM_Predictor_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal mxh = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for the next step
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M1[idx].norm();
				mxh = cu_GetMagnitude(M1[idx] ^ Heff1[idx]) / (conversion * Mnorm * Mnorm);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuaDiffEq.palternator) {

					M1[idx] += dT * (3 * rhs - (*cuaDiffEq.psEval0)[idx]) / 2;
					(*cuaDiffEq.psEval1)[idx] = rhs;
				}
				else {

					M1[idx] += dT * (3 * rhs - (*cuaDiffEq.psEval1)[idx]) / 2;
					(*cuaDiffEq.psEval0)[idx] = rhs;
				}
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunABM_Predictor_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for the next step
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (*cuaDiffEq.palternator) {

					M1[idx] += dT * (3 * rhs - (*cuaDiffEq.psEval0)[idx]) / 2;
					(*cuaDiffEq.psEval1)[idx] = rhs;
				}
				else {

					M1[idx] += dT * (3 * rhs - (*cuaDiffEq.psEval1)[idx]) / 2;
					(*cuaDiffEq.psEval0)[idx] = rhs;
				}
			}
		}
	}
}

__global__ void RunABM_Corrector_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal dmdt = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	cuBReal lte = 0.0;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//First save predicted moment for lte calculation
				cuReal3 saveM = M1[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuaDiffEq.palternator) {

					M1[idx] = (*cuaDiffEq.psM1)[idx] + dT * (rhs + (*cuaDiffEq.psEval1)[idx]) / 2;
				}
				else {

					M1[idx] = (*cuaDiffEq.psM1)[idx] + dT * (rhs + (*cuaDiffEq.psEval0)[idx]) / 2;
				}

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M1[idx] - saveM) / mu_s;
			}
		}
	}

	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunABM_Corrector_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal lte = 0.0;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//First save predicted moment for lte calculation
				cuReal3 saveM = M1[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (*cuaDiffEq.palternator) {

					M1[idx] = (*cuaDiffEq.psM1)[idx] + dT * (rhs + (*cuaDiffEq.psEval1)[idx]) / 2;
				}
				else {

					M1[idx] = (*cuaDiffEq.psM1)[idx] + dT * (rhs + (*cuaDiffEq.psEval0)[idx]) / 2;
				}

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M1[idx] - saveM) / mu_s;
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunABMTEuler_Step0_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for the next step
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				(*cuaDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment for the next time step
				M1[idx] += (*cuaDiffEq.psEval0)[idx] * dT;
			}
		}
	}
}

__global__ void RunABMTEuler_Step1_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using the second trapezoidal Euler step equation
				M1[idx] = ((*cuaDiffEq.psM1)[idx] + M1[idx] + rhs * dT) / 2;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}
			}
		}
	}
}

//----------------------------------------- EVALUATIONS : RKF45

__global__ void RunRKF45_Step0_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal mxh = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				(*cuaDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//obtain maximum normalized torque term
				cuBReal Mnorm = M1[idx].norm();
				mxh = cu_GetMagnitude(M1[idx] ^ Heff1[idx]) / (conversion * Mnorm * Mnorm);

				//Now estimate moment using RKF first step
				M1[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 4);
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunRKF45_Step0_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//Save current moment for later use
			(*cuaDiffEq.psM1)[idx] = M1[idx];

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				(*cuaDiffEq.psEval0)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//Now estimate moment using RKF first step
				M1[idx] += (*cuaDiffEq.psEval0)[idx] * (dT / 4);
			}
		}
	}
}

__global__ void RunRKF45_Step1_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && !M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
			(*cuaDiffEq.psEval1)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

			//Now estimate moment using RKF midle step 1
			M1[idx] = (*cuaDiffEq.psM1)[idx] + (3 * (*cuaDiffEq.psEval0)[idx] + 9 * (*cuaDiffEq.psEval1)[idx]) * dT / 32;
		}
	}
}

__global__ void RunRKF45_Step2_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && !M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
			(*cuaDiffEq.psEval2)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

			//Now estimate moment using RKF midle step 2
			M1[idx] = (*cuaDiffEq.psM1)[idx] + (1932 * (*cuaDiffEq.psEval0)[idx] - 7200 * (*cuaDiffEq.psEval1)[idx] + 7296 * (*cuaDiffEq.psEval2)[idx]) * dT / 2197;
		}
	}
}

__global__ void RunRKF45_Step3_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && !M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
			(*cuaDiffEq.psEval3)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

			//Now estimate moment using RKF midle step 3
			M1[idx] = (*cuaDiffEq.psM1)[idx] + (439 * (*cuaDiffEq.psEval0)[idx] / 216 - 8 * (*cuaDiffEq.psEval1)[idx] + 3680 * (*cuaDiffEq.psEval2)[idx] / 513 - 845 * (*cuaDiffEq.psEval3)[idx] / 4104) * dT;
		}
	}
}

__global__ void RunRKF45_Step4_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && !M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
			(*cuaDiffEq.psEval4)[idx] = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

			//Now estimate moment using RKF midle step 4
			M1[idx] = (*cuaDiffEq.psM1)[idx] + (-8 * (*cuaDiffEq.psEval0)[idx] / 27 + 2 * (*cuaDiffEq.psEval1)[idx] - 3544 * (*cuaDiffEq.psEval2)[idx] / 2565 + 1859 * (*cuaDiffEq.psEval3)[idx] / 4104 - 11 * (*cuaDiffEq.psEval4)[idx] / 40) * dT;
		}
	}
}

__global__ void RunRKF45_Step5_LLG_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal dmdt = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / M1.h.dim();

	cuBReal lte = 0.0;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//4th order evaluation
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (25 * (*cuaDiffEq.psEval0)[idx] / 216 + 1408 * (*cuaDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuaDiffEq.psEval3)[idx] / 4101 - (*cuaDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				M1[idx] = (*cuaDiffEq.psM1)[idx] + (16 * (*cuaDiffEq.psEval0)[idx] / 135 + 6656 * (*cuaDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuaDiffEq.psEval3)[idx] / 56430 - 9 * (*cuaDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//obtain maximum normalized dmdt term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M1[idx] - prediction) / mu_s;
			}
		}
	}

	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

__global__ void RunRKF45_Step5_LLG_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	cuBReal mu_s = *cuaMesh.pmu_s;
	cuBReal alpha = *cuaMesh.palpha;

	cuBReal lte = 0.0;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			if (!M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.palpha, alpha);
				cuReal3 rhs = (-(cuBReal)GAMMA / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

				//4th order evaluation
				cuReal3 prediction = (*cuaDiffEq.psM1)[idx] + (25 * (*cuaDiffEq.psEval0)[idx] / 216 + 1408 * (*cuaDiffEq.psEval2)[idx] / 2565 + 2197 * (*cuaDiffEq.psEval3)[idx] / 4101 - (*cuaDiffEq.psEval4)[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				M1[idx] = (*cuaDiffEq.psM1)[idx] + (16 * (*cuaDiffEq.psEval0)[idx] / 135 + 6656 * (*cuaDiffEq.psEval2)[idx] / 12825 + 28561 * (*cuaDiffEq.psEval3)[idx] / 56430 - 9 * (*cuaDiffEq.psEval4)[idx] / 50 + 2 * rhs / 55) * dT;

				if (*cuaDiffEq.prenormalize) {

					M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				lte = cu_GetMagnitude(M1[idx] - prediction) / mu_s;
			}
		}
	}

	reduction_max(0, 1, &lte, *cuaDiffEq.plte);
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void Atom_DifferentialEquationCubicCUDA::RunEuler_LLG(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh || calculate_dmdt) {

		RunEuler_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else {

		RunEuler_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
}

//TRAPEZOIDAL EULER

void Atom_DifferentialEquationCubicCUDA::RunTEuler_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunTEuler_Step0_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunTEuler_Step0_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunTEuler_Step1_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunTEuler_Step1_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
	}
}

//ADAPTIVE HEUN

void Atom_DifferentialEquationCubicCUDA::RunAHeun_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunTEuler_Step0_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunTEuler_Step0_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunAHeun_Step1_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunAHeun_Step1_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
	}
}

//RUNGE KUTTA 4th order

void Atom_DifferentialEquationCubicCUDA::RunRK4_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRK4_Step0_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRK4_Step0_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;

	case 1:

		RunRK4_Step1_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRK4_Step2_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		if (calculate_dmdt) {

			RunRK4_Step3_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRK4_Step3_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

//Adams-Bashforth-Moulton 2nd order

void Atom_DifferentialEquationCubicCUDA::RunABM_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	if (step == 0) {

		if (calculate_mxh) {

			RunABM_Predictor_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunABM_Predictor_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
	}
	else {

		if (calculate_dmdt) {

			RunABM_Corrector_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunABM_Corrector_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
	}
}

//Adams-Bashforth-Moulton 2nd order priming using Trapezoidal Euler

void Atom_DifferentialEquationCubicCUDA::RunABMTEuler_LLG(int step)
{
	if (step == 0) {

		RunABMTEuler_Step0_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>(cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else {

		RunABMTEuler_Step1_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>(cuaDiffEq, paMeshCUDA->cuaMesh);
	}
}

//RUNGE KUTTA FEHLBERG

void Atom_DifferentialEquationCubicCUDA::RunRKF45_LLG(int step, bool calculate_mxh, bool calculate_dmdt)
{
	switch (step) {

	case 0:

		if (calculate_mxh) {

			RunRKF45_Step0_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKF45_Step0_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;

	case 1:

		RunRKF45_Step1_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 2:

		RunRKF45_Step2_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 3:

		RunRKF45_Step3_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 4:

		RunRKF45_Step4_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);

		break;

	case 5:

		if (calculate_dmdt) {

			RunRKF45_Step5_LLG_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}
		else {

			RunRKF45_Step5_LLG_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
		}

		break;
	}
}

#endif
#endif