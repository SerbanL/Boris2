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

	cuBReal mu_s = *pcuaMesh->pmu_s;
	cuBReal alpha = *pcuaMesh->palpha;
	cuBReal grel = *pcuaMesh->pgrel;
	pcuaMesh->update_parameters_mcoarse(idx, *pcuaMesh->pmu_s, mu_s, *pcuaMesh->palpha, alpha, *pcuaMesh->pgrel, grel);

	return (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M1[idx] ^ Heff1[idx]) + alpha * ((M1[idx] / mu_s) ^ (M1[idx] ^ Heff1[idx])));
}

#endif
#endif
