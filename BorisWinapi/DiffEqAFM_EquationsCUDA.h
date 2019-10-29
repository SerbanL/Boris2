#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

//Defines non-stochastic equations

#include "BorisCUDALib.h"

#include "ManagedDiffEqAFMCUDA.h"

#include "DiffEq_Defs.h"

#include "Funcs_Math_base.h" //includes constant values

#include "MeshParamsControlCUDA.h"

//----------------------------------------- EQUATIONS

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLG(int idx, cuReal3& value_B)
{
	//gamma = -mu0 * gamma_e = mu0 * g e / 2m_e = 2.212761569e5 m/As

	//LLG in explicit form : dm/dt = [mu0*gamma_e/(1+alpha^2)] * [m*H + alpha * m*(m*H)]

	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuVEC_VC<cuReal3>& M2 = *pcuMesh->pM2;
	cuVEC<cuReal3>& Heff2 = *pcuMesh->pHeff2;

	cuReal2 Ms_AFM = *pcuMesh->pMs_AFM;
	cuReal2 alpha_AFM = *pcuMesh->palpha_AFM;
	cuReal2 grel_AFM = *pcuMesh->pgrel_AFM;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs_AFM, Ms_AFM, *pcuMesh->palpha_AFM, alpha_AFM, *pcuMesh->pgrel_AFM, grel_AFM);

	//sub-lattice B value so we can read it after
	value_B = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

	//return the sub-lattice A value as normal
	return (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
}

//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
__device__ cuReal3 ManagedDiffEqAFMCUDA::LLGStatic(int idx, cuReal3& value_B)
{
	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuVEC_VC<cuReal3>& M2 = *pcuMesh->pM2;
	cuVEC<cuReal3>& Heff2 = *pcuMesh->pHeff2;

	cuReal2 Ms_AFM = *pcuMesh->pMs_AFM;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs_AFM, Ms_AFM);

	//sub-lattice B value so we can read it after
	value_B = (-(cuBReal)GAMMA / 2) * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx]));

	//return the sub-lattice A value as normal
	return (-(cuBReal)GAMMA / 2) * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx]));
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLGSTT(int idx, cuReal3& value_B)
{
	return cuReal3();
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLB(int idx, cuReal3& value_B)
{
	return cuReal3();
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLBSTT(int idx, cuReal3& value_B)
{
	return cuReal3();
}

#endif
