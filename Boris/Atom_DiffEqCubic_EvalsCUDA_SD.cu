#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_SD
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: SD Solver

__global__ void RunSD_Start_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				/////////////////////////

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuBReal grel = *cuaMesh.pgrel;
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pgrel, grel);

				/////////////////////////

				cuReal3 m = (*cuaMesh.pM1)[idx] / mu_s;
				cuReal3 H = (*cuaMesh.pHeff1)[idx];

				//calculate M cross Heff (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				cuReal3 MxHeff = ((cuBReal)GAMMA / 2) * ((*cuaMesh.pM1)[idx] ^ H);

				/////////////////////////

				//current torque value G = m x (M x H)
				cuReal3 G = m ^ MxHeff;

				//save calculated torque for next time
				(*cuaDiffEq.psEval0)[idx] = G;

				//save current M for next time
				(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

				/////////////////////////

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel / 4.0;

				cuReal3 mxH = m ^ H;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));

				//set new M
				(*cuaMesh.pM1)[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuaMesh.pM1)[idx].renormalize(mu_s);

				/////////////////////////
			}
		}
	}
}

__global__ void RunSD_BB_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal _delta_M_sq = 0.0;
	cuBReal _delta_G_sq = 0.0;
	cuBReal _delta_M_dot_delta_G = 0.0;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				/////////////////////////

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);

				cuReal3 m = (*cuaMesh.pM1)[idx] / mu_s;
				cuReal3 H = (*cuaMesh.pHeff1)[idx];

				//calculate M cross Heff (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				cuReal3 MxHeff = ((cuBReal)GAMMA / 2) * ((*cuaMesh.pM1)[idx] ^ H);

				/////////////////////////

				//current torque value G = m x (M x H)
				cuReal3 G = m ^ MxHeff;

				//change in torque
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				cuReal3 delta_G = (G - (*cuaDiffEq.psEval0)[idx]) / 1e6;

				//save calculated torque for next time
				(*cuaDiffEq.psEval0)[idx] = G;

				/////////////////////////

				//change in m
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				cuReal3 delta_M = ((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / 1e6;

				//save current M for next time
				(*cuaDiffEq.psM1)[idx] = (*cuaMesh.pM1)[idx];

				/////////////////////////

				//calculate num and denom for the two Barzilai-Borwein stepsize solutions (see Journal of Numerical Analysis (1988) 8, 141-148) so we can find new stepsize
				_delta_M_sq = delta_M * delta_M;
				_delta_G_sq = delta_G * delta_G;
				_delta_M_dot_delta_G = delta_M * delta_G;
			}
		}
	}

	//the delta_... quantities in which we accumulate are cuBReal (thus likely single precision depending on compilation flag SINGLEPRECISION)
	//Bear in mind this could potentially cause catastrophic loss of precision for large simulations even though we normalized to 1e6*1e6 - unlikely however so leave them like this for now.
	//e.g. when I used these quantities as not normalized the stepsize calculation was all wrong even for medium sized meshes.
	//If you ever have problems then this is the likely culprit
	reduction_sum(0, 1, &_delta_M_sq, *cuaDiffEq.pdelta_M_sq);
	reduction_sum(0, 1, &_delta_G_sq, *cuaDiffEq.pdelta_G_sq);
	reduction_sum(0, 1, &_delta_M_dot_delta_G, *cuaDiffEq.pdelta_M_dot_delta_G);
}

__global__ void RunSD_Advance_withReductions_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal mxh = 0.0;
	cuBReal dmdt = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				/////////////////////////

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuBReal grel = *cuaMesh.pgrel;
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pgrel, grel);

				cuReal3 m = (*cuaMesh.pM1)[idx] / mu_s;
				cuReal3 H = (*cuaMesh.pHeff1)[idx];

				//obtained maximum normalized torque term
				cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
				if (cuIsNZ(grel)) mxh = cu_GetMagnitude(m ^ H) / (conversion * Mnorm);

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel / 4.0;

				cuReal3 mxH = m ^ H;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));

				//set new M
				(*cuaMesh.pM1)[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuaMesh.pM1)[idx].renormalize(mu_s);

				//obtain maximum normalized dmdt term
				if (cuIsNZ(grel)) dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
}

__global__ void RunSD_Advance_withReduction_mxh_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal mxh = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				/////////////////////////

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuBReal grel = *cuaMesh.pgrel;
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pgrel, grel);

				cuReal3 m = (*cuaMesh.pM1)[idx] / mu_s;
				cuReal3 H = (*cuaMesh.pHeff1)[idx];

				//obtained maximum normalized torque term
				if (cuIsNZ(grel)) mxh = cu_GetMagnitude(m ^ H) / (conversion * (*cuaMesh.pM1)[idx].norm());

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel / 4.0;

				cuReal3 mxH = m ^ H;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));

				//set new M
				(*cuaMesh.pM1)[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuaMesh.pM1)[idx].renormalize(mu_s);
			}
		}
	}

	
	reduction_max(0, 1, &mxh, *cuaDiffEq.pmxh);
}

__global__ void RunSD_Advance_withReduction_dmdt_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	cuBReal dmdt = 0.0;

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	cuBReal conversion = (cuBReal)MUB / cuaMesh.pM1->h.dim();

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				/////////////////////////

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuBReal grel = *cuaMesh.pgrel;
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pgrel, grel);

				cuReal3 m = (*cuaMesh.pM1)[idx] / mu_s;
				cuReal3 H = (*cuaMesh.pHeff1)[idx];

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel / 4.0;

				cuReal3 mxH = m ^ H;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));

				//set new M
				(*cuaMesh.pM1)[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuaMesh.pM1)[idx].renormalize(mu_s);

				//obtain maximum normalized dmdt term
				if (cuIsNZ(grel)) {

					cuBReal Mnorm = (*cuaMesh.pM1)[idx].norm();
					dmdt = cu_GetMagnitude((*cuaMesh.pM1)[idx] - (*cuaDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * conversion * Mnorm * Mnorm);
				}
			}
		}
	}

	reduction_max(0, 1, &dmdt, *cuaDiffEq.pdmdt);
}

__global__ void RunSD_Advance_Kernel(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuaDiffEq.pdT;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			if (!cuaMesh.pM1->is_skipcell(idx)) {

				/////////////////////////

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuBReal grel = *cuaMesh.pgrel;
				cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pgrel, grel);

				cuReal3 m = (*cuaMesh.pM1)[idx] / mu_s;
				cuReal3 H = (*cuaMesh.pHeff1)[idx];

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel / 4.0;

				cuReal3 mxH = m ^ H;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));

				//set new M
				(*cuaMesh.pM1)[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuaMesh.pM1)[idx].renormalize(mu_s);
			}
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//SD Solver

void Atom_DifferentialEquationCubicCUDA::RunSD_Start(void)
{
	RunSD_Start_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
}

void Atom_DifferentialEquationCubicCUDA::RunSD_BB(void)
{
	RunSD_BB_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
}

void Atom_DifferentialEquationCubicCUDA::RunSD_Advance(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh && calculate_dmdt) {

		RunSD_Advance_withReductions_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else if (calculate_mxh) {

		RunSD_Advance_withReduction_mxh_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else if (calculate_dmdt) {

		RunSD_Advance_withReduction_dmdt_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
	else {

		RunSD_Advance_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaDiffEq, paMeshCUDA->cuaMesh);
	}
}

#endif
#endif
#endif