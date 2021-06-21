#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef ODE_EVAL_COMPILATION_SD
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: SD Solver

__global__ void RunSD_Start_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				/////////////////////////

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuReal2 grel_AFM = *cuMesh.pgrel_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				/////////////////////////

				cuReal3 m = (*cuMesh.pM)[idx] / Ms_AFM.i;
				cuReal3 H = (*cuMesh.pHeff)[idx];

				cuReal3 m2 = (*cuMesh.pM2)[idx] / Ms_AFM.j;
				cuReal3 H2 = (*cuMesh.pHeff2)[idx];

				//calculate M cross Heff (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				cuReal3 MxHeff = ((cuBReal)GAMMA / 2) * ((*cuMesh.pM)[idx] ^ H);
				cuReal3 MxHeff2 = ((cuBReal)GAMMA / 2) * ((*cuMesh.pM2)[idx] ^ H2);

				/////////////////////////

				//current torque value G = m x (M x H)
				cuReal3 G = m ^ MxHeff;
				cuReal3 G2 = m2 ^ MxHeff2;

				//save calculated torque for next time
				(*cuDiffEq.psEval0)[idx] = G;
				(*cuDiffEq.psEval0_2)[idx] = G2;

				//save current M for next time
				(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
				(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

				/////////////////////////

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel_AFM.i / 4.0;

				cuReal3 mxH = m ^ H;
				cuReal3 mxH2 = m2 ^ H2;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));
				m2 = ((1 - s*s*(mxH2*mxH2)) * m2 - 2*s*(m2 ^ mxH2)) / (1 + s*s*(mxH2*mxH2));

				//set new M
				(*cuMesh.pM)[idx] = m * Ms_AFM.i;
				(*cuMesh.pM2)[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);

				/////////////////////////
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

__global__ void RunSD_BB_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal _delta_M_sq = 0.0;
	cuBReal _delta_G_sq = 0.0;
	cuBReal _delta_M_dot_delta_G = 0.0;

	cuBReal _delta_M2_sq = 0.0;
	cuBReal _delta_G2_sq = 0.0;
	cuBReal _delta_M2_dot_delta_G2 = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				/////////////////////////

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);

				cuReal3 m = (*cuMesh.pM)[idx] / Ms_AFM.i;
				cuReal3 H = (*cuMesh.pHeff)[idx];

				cuReal3 m2 = (*cuMesh.pM2)[idx] / Ms_AFM.j;
				cuReal3 H2 = (*cuMesh.pHeff2)[idx];

				//calculate M cross Heff (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				cuReal3 MxHeff = ((cuBReal)GAMMA / 2) * ((*cuMesh.pM)[idx] ^ H);
				cuReal3 MxHeff2 = ((cuBReal)GAMMA / 2) * ((*cuMesh.pM2)[idx] ^ H2);

				/////////////////////////

				//current torque value G = m x (M x H)
				cuReal3 G = m ^ MxHeff;
				cuReal3 G2 = m2 ^ MxHeff2;

				//change in torque
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				cuReal3 delta_G = (G - (*cuDiffEq.psEval0)[idx]) / 1e6;
				cuReal3 delta_G2 = (G2 - (*cuDiffEq.psEval0_2)[idx]) / 1e6;

				//save calculated torque for next time
				(*cuDiffEq.psEval0)[idx] = G;
				(*cuDiffEq.psEval0_2)[idx] = G2;

				/////////////////////////

				//change in m
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				cuReal3 delta_M = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / 1e6;
				cuReal3 delta_M2 = ((*cuMesh.pM2)[idx] - (*cuDiffEq.psM1_2)[idx]) / 1e6;

				//save current M for next time
				(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
				(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

				/////////////////////////

				//calculate num and denom for the two Barzilai-Borwein stepsize solutions (see Journal of Numerical Analysis (1988) 8, 141-148) so we can find new stepsize
				_delta_M_sq = delta_M * delta_M;
				_delta_G_sq = delta_G * delta_G;
				_delta_M_dot_delta_G = delta_M * delta_G;

				_delta_M2_sq = delta_M2 * delta_M2;
				_delta_G2_sq = delta_G2 * delta_G2;
				_delta_M2_dot_delta_G2 = delta_M2 * delta_G2;
			}
		}
	}

	//the delta_... quantities in which we accumulate are cuBReal (thus likely single precision depending on compilation flag SINGLEPRECISION)
	//Bear in mind this could potentially cause catastrophic loss of precision for large simulations even though we normalized to 1e6*1e6 - unlikely however so leave them like this for now.
	//e.g. when I used these quantities as not normalized the stepsize calculation was all wrong even for medium sized meshes.
	//If you ever have problems then this is the likely culprit
	reduction_sum(0, 1, &_delta_M_sq, *cuDiffEq.pdelta_M_sq);
	reduction_sum(0, 1, &_delta_G_sq, *cuDiffEq.pdelta_G_sq);
	reduction_sum(0, 1, &_delta_M_dot_delta_G, *cuDiffEq.pdelta_M_dot_delta_G);

	reduction_sum(0, 1, &_delta_M2_sq, *cuDiffEq.pdelta_M2_sq);
	reduction_sum(0, 1, &_delta_G2_sq, *cuDiffEq.pdelta_G2_sq);
	reduction_sum(0, 1, &_delta_M2_dot_delta_G2, *cuDiffEq.pdelta_M2_dot_delta_G2);
}

__global__ void RunSD_Advance_withReductions_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal mxh = 0.0;
	cuBReal dmdt = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				/////////////////////////

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuReal2 grel_AFM = *cuMesh.pgrel_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 m = (*cuMesh.pM)[idx] / Ms_AFM.i;
				cuReal3 H = (*cuMesh.pHeff)[idx];

				cuReal3 m2 = (*cuMesh.pM2)[idx] / Ms_AFM.j;
				cuReal3 H2 = (*cuMesh.pHeff2)[idx];

				//obtained maximum normalized torque term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				if (cuIsNZ(grel_AFM.i)) mxh = cu_GetMagnitude(m ^ H) / Mnorm;

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel_AFM.i / 4.0;

				cuReal3 mxH = m ^ H;
				cuReal3 mxH2 = m2 ^ H2;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));
				m2 = ((1 - s*s*(mxH2*mxH2)) * m2 - 2*s*(m2 ^ mxH2)) / (1 + s*s*(mxH2*mxH2));

				//set new M
				(*cuMesh.pM)[idx] = m * Ms_AFM.i;
				(*cuMesh.pM2)[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);

				//obtain maximum normalized dmdt term
				if (cuIsNZ(grel_AFM.i)) dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * grel_AFM.i * Mnorm * Mnorm);
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
	reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
}

__global__ void RunSD_Advance_withReduction_mxh_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal mxh = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				/////////////////////////

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuReal2 grel_AFM = *cuMesh.pgrel_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 m = (*cuMesh.pM)[idx] / Ms_AFM.i;
				cuReal3 H = (*cuMesh.pHeff)[idx];

				cuReal3 m2 = (*cuMesh.pM2)[idx] / Ms_AFM.j;
				cuReal3 H2 = (*cuMesh.pHeff2)[idx];

				//obtained maximum normalized torque term
				if (cuIsNZ(grel_AFM.i)) mxh = cu_GetMagnitude(m ^ H) / (*cuMesh.pM)[idx].norm();

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel_AFM.i / 4.0;

				cuReal3 mxH = m ^ H;
				cuReal3 mxH2 = m2 ^ H2;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));
				m2 = ((1 - s*s*(mxH2*mxH2)) * m2 - 2*s*(m2 ^ mxH2)) / (1 + s*s*(mxH2*mxH2));

				//set new M
				(*cuMesh.pM)[idx] = m * Ms_AFM.i;
				(*cuMesh.pM2)[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	reduction_max(0, 1, &mxh, *cuDiffEq.pmxh);
}

__global__ void RunSD_Advance_withReduction_dmdt_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuBReal dmdt = 0.0;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				/////////////////////////

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuReal2 grel_AFM = *cuMesh.pgrel_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 m = (*cuMesh.pM)[idx] / Ms_AFM.i;
				cuReal3 H = (*cuMesh.pHeff)[idx];

				cuReal3 m2 = (*cuMesh.pM2)[idx] / Ms_AFM.j;
				cuReal3 H2 = (*cuMesh.pHeff2)[idx];

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel_AFM.i / 4.0;

				cuReal3 mxH = m ^ H;
				cuReal3 mxH2 = m2 ^ H2;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));
				m2 = ((1 - s*s*(mxH2*mxH2)) * m2 - 2*s*(m2 ^ mxH2)) / (1 + s*s*(mxH2*mxH2));

				//set new M
				(*cuMesh.pM)[idx] = m * Ms_AFM.i;
				(*cuMesh.pM2)[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);

				//obtain maximum normalized dmdt term
				if (cuIsNZ(grel_AFM.i)) {

					cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
					dmdt = cu_GetMagnitude((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * grel_AFM.i * Mnorm * Mnorm);
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

	reduction_max(0, 1, &dmdt, *cuDiffEq.pdmdt);
}

__global__ void RunSD_Advance_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			if (!cuMesh.pM->is_skipcell(idx)) {

				/////////////////////////

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuReal2 grel_AFM = *cuMesh.pgrel_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pgrel_AFM, grel_AFM);

				cuReal3 m = (*cuMesh.pM)[idx] / Ms_AFM.i;
				cuReal3 H = (*cuMesh.pHeff)[idx];

				cuReal3 m2 = (*cuMesh.pM2)[idx] / Ms_AFM.j;
				cuReal3 H2 = (*cuMesh.pHeff2)[idx];

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				cuBReal s = dT * (cuBReal)GAMMA * grel_AFM.i / 4.0;

				cuReal3 mxH = m ^ H;
				cuReal3 mxH2 = m2 ^ H2;
				m = ((1 - s*s*(mxH*mxH)) * m - 2*s*(m ^ mxH)) / (1 + s*s*(mxH*mxH));
				m2 = ((1 - s*s*(mxH2*mxH2)) * m2 - 2*s*(m2 ^ mxH2)) / (1 + s*s*(mxH2*mxH2));

				//set new M
				(*cuMesh.pM)[idx] = m * Ms_AFM.i;
				(*cuMesh.pM2)[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
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

//SD Solver

void DifferentialEquationAFMCUDA::RunSD_Start(void)
{
	RunSD_Start_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
}

void DifferentialEquationAFMCUDA::RunSD_BB(void)
{
	RunSD_BB_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuDiffEq, pMeshCUDA->cuMesh);
}

void DifferentialEquationAFMCUDA::RunSD_Advance(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh && calculate_dmdt) {

		RunSD_Advance_withReductions_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else if (calculate_mxh) {

		RunSD_Advance_withReduction_mxh_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else if (calculate_dmdt) {

		RunSD_Advance_withReduction_dmdt_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunSD_Advance_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

#endif
#endif
#endif