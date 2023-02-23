#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

//Defines non-stochastic equations

#include "BorisCUDALib.h"

#include "ManagedAtom_DiffEqCubicCUDA.h"

#include "DiffEq_Defs.h"

#include "Funcs_Math_base.h" //includes constant values

#include "Atom_MeshParamsControlCUDA.h"

//----------------------------------------- EQUATIONS

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedAtom_DiffEqCubicCUDA::LLG(int idx)
{
	//gamma_e is the electron gyromagnetic ratio, which is negative

	//LLG in explicit form : dm/dt = [gamma_e/(1+alpha^2)] * [m*B + alpha * m*(m*B) / mu_s]
	//B is effective field in Tesla; prefer to include the effective field in A/m as it integrates easier in a multiscale simulation; thus gamma = |gamma_e| * mu0
	//m is atomic moment in units of muB : has magnitude mu_s (muB)

	cuVEC_VC<cuReal3>& M1 = *pcuaMesh->pM1;
	cuVEC<cuReal3>& Heff1 = *pcuaMesh->pHeff1;

	cuBReal mu_s = *pcuaMesh->pmu_s;
	cuBReal alpha = *pcuaMesh->palpha;
	cuBReal grel = *pcuaMesh->pgrel;
	pcuaMesh->update_parameters_mcoarse(idx, *pcuaMesh->pmu_s, mu_s, *pcuaMesh->palpha, alpha, *pcuaMesh->pgrel, grel);

	return (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));
}

//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
__device__ cuReal3 ManagedAtom_DiffEqCubicCUDA::LLGStatic(int idx)
{
	cuVEC_VC<cuReal3>& M1 = *pcuaMesh->pM1;
	cuVEC<cuReal3>& Heff1 = *pcuaMesh->pHeff1;

	cuBReal mu_s = *pcuaMesh->pmu_s;
	cuBReal grel = *pcuaMesh->pgrel;
	pcuaMesh->update_parameters_mcoarse(idx, *pcuaMesh->pmu_s, mu_s, *pcuaMesh->pgrel, grel);

	return (-(cuBReal)GAMMA * grel / 2) * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx]));
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedAtom_DiffEqCubicCUDA::LLGSTT(int idx)
{
	//Currently same as LLG : TO DO

	//gamma_e is the electron gyromagnetic ratio, which is negative

	//LLG in explicit form : dm/dt = [gamma_e/(1+alpha^2)] * [m*B + alpha * m*(m*B) / mu_s]
	//B is effective field in Tesla; prefer to include the effective field in A/m as it integrates easier in a multiscale simulation; thus gamma = |gamma_e| * mu0
	//m is atomic moment in units of muB : has magnitude mu_s (muB)

	cuVEC_VC<cuReal3>& M1 = *pcuaMesh->pM1;
	cuVEC<cuReal3>& Heff1 = *pcuaMesh->pHeff1;
	cuVEC_VC<cuReal3>& E = *pcuaMesh->pE;
	cuVEC_VC<cuBReal>& elC = *pcuaMesh->pelC;

	cuBReal mu_s = *pcuaMesh->pmu_s;
	cuBReal alpha = *pcuaMesh->palpha;
	cuBReal grel = *pcuaMesh->pgrel;
	cuBReal P = *pcuaMesh->pP;
	cuBReal beta = *pcuaMesh->pbeta;
	pcuaMesh->update_parameters_mcoarse(idx, *pcuaMesh->pmu_s, mu_s, *pcuaMesh->palpha, alpha, *pcuaMesh->pgrel, grel, *pcuaMesh->pP, P, *pcuaMesh->pbeta, beta);

	cuReal3 LLGSTT_Eval = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));

	if (E.linear_size()) {

		cuSZ3 n = M1.n;
		cuReal3 h = M1.h;

		cuReal33 grad_M1 = M1.grad_neu(idx);

		cuReal3 position = M1.cellidx_to_position(idx);

		cuBReal conv = h.dim() / MUB;
		cuReal3 u = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E * conv) / (mu_s * (1 + beta * beta));

		cuReal3 u_dot_del_M1 = (u.x * grad_M1.x) + (u.y * grad_M1.y) + (u.z * grad_M1.z);

		LLGSTT_Eval += (((1 + alpha * beta) * u_dot_del_M1) - ((beta - alpha) * ((M1[idx] / mu_s) ^ u_dot_del_M1))) / (1 + alpha * alpha);
	}

	return LLGSTT_Eval;
}

#endif
#endif
