#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

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
	cuReal2 grel_AFM = *pcuMesh->pgrel_AFM;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs_AFM, Ms_AFM, *pcuMesh->pgrel_AFM, grel_AFM);

	//sub-lattice B value so we can read it after
	value_B = (-(cuBReal)GAMMA * grel_AFM.j / 2) * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx]));

	//return the sub-lattice A value as normal
	return (-(cuBReal)GAMMA * grel_AFM.i / 2) * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx]));
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLGSTT(int idx, cuReal3& value_B)
{
	//gmub_2e is -hbar * gamma_e / 2e = g mu_b / 2e)

	// LLG with STT in explicit form : dm/dt = [mu0*gamma_e/(1+alpha^2)] * [m*H + alpha * m*(m*H)] + (1+alpha*beta)/((1+alpha^2)*(1+beta^2)) * (u.del)m - (beta - alpha)/(1+alpha^2) * m * (u.del) m
	// where u = j * P g mu_b / 2e Ms = -(hbar * gamma_e * P / 2 *e * Ms) * j, j is the current density = conductivity * E (A/m^2)

	// STT is Zhang-Li equationtion (not Thiaville, the velocity used by Thiaville needs to be divided by (1+beta^2) to obtain Zhang-Li, also Thiaville's EPL paper has wrong STT signs!!)

	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;
	cuVEC_VC<cuReal3>& M2 = *pcuMesh->pM2;
	cuVEC<cuReal3>& Heff2 = *pcuMesh->pHeff2;
	cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
	cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

	cuReal2 Ms_AFM = *pcuMesh->pMs_AFM;
	cuReal2 alpha_AFM = *pcuMesh->palpha_AFM;
	cuReal2 grel_AFM = *pcuMesh->pgrel_AFM;
	cuBReal P = *pcuMesh->pP;
	cuBReal beta = *pcuMesh->pbeta;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs_AFM, Ms_AFM, *pcuMesh->palpha_AFM, alpha_AFM, *pcuMesh->pgrel_AFM, grel_AFM, *pcuMesh->pP, P, *pcuMesh->pbeta, beta);

	cuReal3 LLGSTT_Eval_A = (-(cuBReal)GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((M[idx] ^ Heff[idx]) + alpha_AFM.i * ((M[idx] / Ms_AFM.i) ^ (M[idx] ^ Heff[idx])));
	cuReal3 LLGSTT_Eval_B = (-(cuBReal)GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((M2[idx] ^ Heff2[idx]) + alpha_AFM.j * ((M2[idx] / Ms_AFM.j) ^ (M2[idx] ^ Heff2[idx])));

	if (E.linear_size()) {

		cuSZ3 n = M.n;
		cuReal3 h = M.h;

		cuReal33 grad_M_A = M.grad_neu(idx);
		cuReal33 grad_M_B = M2.grad_neu(idx);

		cuReal3 position = M.cellidx_to_position(idx);

		cuReal3 u_A = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E) / (Ms_AFM.i * (1 + beta * beta));
		cuReal3 u_B = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E) / (Ms_AFM.j * (1 + beta * beta));

		cuReal3 u_dot_del_M_A = (u_A.x * grad_M_A.x) + (u_A.y * grad_M_A.y) + (u_A.z * grad_M_A.z);
		cuReal3 u_dot_del_M_B = (u_B.x * grad_M_B.x) + (u_B.y * grad_M_B.y) + (u_B.z * grad_M_B.z);

		LLGSTT_Eval_A += (((1 + alpha_AFM.i * beta) * u_dot_del_M_A) - ((beta - alpha_AFM.i) * ((M[idx] / Ms_AFM.i) ^ u_dot_del_M_A))) / (1 + alpha_AFM.i * alpha_AFM.i);
		LLGSTT_Eval_B += (((1 + alpha_AFM.j * beta) * u_dot_del_M_B) - ((beta - alpha_AFM.j) * ((M2[idx] / Ms_AFM.j) ^ u_dot_del_M_B))) / (1 + alpha_AFM.j * alpha_AFM.j);
	}

	//sub-lattice B value so we can read it after
	value_B = LLGSTT_Eval_B;

	//return the sub-lattice A value as normal
	return LLGSTT_Eval_A;
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLB(int idx, cuReal3& value_B)
{
	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC_VC<cuReal3>& M2 = *pcuMesh->pM2;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;
	cuVEC<cuReal3>& Heff2 = *pcuMesh->pHeff2;

	cuBReal T_Curie = *pcuMesh->pT_Curie;

	//cell temperature : the base temperature if uniform temperature, else get the temperature from Temp
	cuBReal Temperature;
	if (pcuMesh->pTemp->linear_size()) Temperature = (*pcuMesh->pTemp)[pcuMesh->pM->cellidx_to_position(idx)];
	else Temperature = *pcuMesh->pbase_temperature;

	//m is M / Ms0 : magnitude of M in this cell divided by the saturation magnetization at 0K.
	cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());
	cuReal2 Ms0 = pcuMesh->pMs_AFM->get0();
	cuReal2 m = Mmag / Ms0;
	cuReal2 msq = m & m;

	cuReal2 Ms = *pcuMesh->pMs_AFM;
	cuReal2 alpha = *pcuMesh->palpha_AFM;
	cuReal2 grel = *pcuMesh->pgrel_AFM;
	cuReal2 susrel = *pcuMesh->psusrel_AFM;
	cuReal2 tau_ii = *pcuMesh->ptau_ii;
	cuReal2 tau_ij = *pcuMesh->ptau_ij;
	cuReal2 mu = *pcuMesh->patomic_moment_AFM;

	cuReal2 alpha_par;

	//the longitudinal relaxation field - an effective field contribution, but only need to add it to the longitudinal relaxation term as the others involve cross products with pMesh->M[idx]
	cuReal3 Hl_1, Hl_2;

	if (Temperature < T_Curie) {

		if (Temperature > T_Curie - (cuBReal)TCURIE_EPSILON) {

			Ms = pcuMesh->pMs_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			alpha = pcuMesh->palpha_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs_AFM, Ms, *pcuMesh->palpha_AFM, alpha, *pcuMesh->pgrel_AFM, grel, *pcuMesh->psusrel_AFM, susrel);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel.i > 1.0) susrel.i = 1.0;
		if (susrel.j > 1.0) susrel.j = 1.0;

		alpha_par = 2 * (pcuMesh->palpha_AFM->get0() - alpha);

		cuReal2 me = Ms / Ms0;
		cuReal2 r = m / me;

		Hl_1 = (M[idx] / (2 * (cuBReal)MU0 * Ms0.i)) * ((1.0 - r.i * r.i) / susrel.i + (3 * tau_ij.i * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.i) * ((susrel.j / susrel.i) * (1 - r.i * r.i) + (me.j / me.i) * (1 - r.j * r.j)));
		Hl_2 = (M2[idx] / (2 * (cuBReal)MU0 * Ms0.j)) * ((1.0 - r.j * r.j) / susrel.j + (3 * tau_ij.j * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.j) * ((susrel.i / susrel.j) * (1 - r.j * r.j) + (me.i / me.j) * (1 - r.i * r.i)));
	}
	else {

		if (Temperature < T_Curie + (cuBReal)TCURIE_EPSILON) {

			alpha = pcuMesh->palpha_AFM->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel_AFM->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel_AFM->get(T_Curie + (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->palpha_AFM, alpha, *pcuMesh->pgrel_AFM, grel, *pcuMesh->psusrel_AFM, susrel);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel.i > 1.0) susrel.i = 1.0;
		if (susrel.j > 1.0) susrel.j = 1.0;

		alpha_par = alpha;

		cuReal2 me = Ms / Ms0;

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		Hl_1 = -1 * (M[idx] / ((cuBReal)MU0 * Ms0.i)) * ((1.0 / susrel.i) + (3 * tau_ij.i * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.i) * ((susrel.j / susrel.i) + 1));
		Hl_2 = -1 * (M2[idx] / ((cuBReal)MU0 * Ms0.j)) * ((1.0 / susrel.j) + (3 * tau_ij.j * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.j) * ((susrel.i / susrel.j) + 1));
	}

	//sub-lattice B value so we can read it after
	value_B = (-(cuBReal)GAMMA * grel.j * msq.j / (msq.j + alpha.j * alpha.j)) * (M2[idx] ^ Heff2[idx]) + (-(cuBReal)GAMMA * grel.j * m.j * alpha.j / (msq.j + alpha.j * alpha.j)) * ((M2[idx] / Mmag.j) ^ (M2[idx] ^ Heff2[idx])) +
		(cuBReal)GAMMA * grel.j * alpha_par.j * Ms0.j * ((M2[idx] / Mmag.j) * (Heff2[idx] + Hl_2)) * (M2[idx] / Mmag.j);

	//return the sub-lattice A value as normal
	return (-(cuBReal)GAMMA * grel.i * msq.i / (msq.i + alpha.i * alpha.i)) * (M[idx] ^ Heff[idx]) + (-(cuBReal)GAMMA * grel.i * m.i * alpha.i / (msq.i + alpha.i * alpha.i)) * ((M[idx] / Mmag.i) ^ (M[idx] ^ Heff[idx])) +
		(cuBReal)GAMMA * grel.i * alpha_par.i * Ms0.i * ((M[idx] / Mmag.i) * (Heff[idx] + Hl_1)) * (M[idx] / Mmag.i);
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqAFMCUDA::LLBSTT(int idx, cuReal3& value_B)
{
	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC_VC<cuReal3>& M2 = *pcuMesh->pM2;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;
	cuVEC<cuReal3>& Heff2 = *pcuMesh->pHeff2;
	cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
	cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

	cuBReal T_Curie = *pcuMesh->pT_Curie;
	
	cuReal3 position = pcuMesh->pM->cellidx_to_position(idx);

	//cell temperature : the base temperature if uniform temperature, else get the temperature from Temp
	cuBReal Temperature;
	if (pcuMesh->pTemp->linear_size()) Temperature = (*pcuMesh->pTemp)[position];
	else Temperature = *pcuMesh->pbase_temperature;

	//m is M / Ms0 : magnitude of M in this cell divided by the saturation magnetization at 0K.
	cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());
	cuReal2 Ms0 = pcuMesh->pMs_AFM->get0();
	cuReal2 m = Mmag / Ms0;
	cuReal2 msq = m & m;

	cuReal2 Ms = *pcuMesh->pMs_AFM;
	cuReal2 alpha = *pcuMesh->palpha_AFM;
	cuReal2 grel = *pcuMesh->pgrel_AFM;
	cuReal2 susrel = *pcuMesh->psusrel_AFM;
	cuBReal P = *pcuMesh->pP;
	cuBReal beta = *pcuMesh->pbeta;
	
	cuReal2 tau_ii = *pcuMesh->ptau_ii;
	cuReal2 tau_ij = *pcuMesh->ptau_ij;
	cuReal2 mu = *pcuMesh->patomic_moment_AFM;

	cuReal2 alpha_par;

	//the longitudinal relaxation field - an effective field contribution, but only need to add it to the longitudinal relaxation term as the others involve cross products with pMesh->M[idx]
	cuReal3 Hl_1, Hl_2;

	if (Temperature < T_Curie) {

		if (Temperature > T_Curie - (cuBReal)TCURIE_EPSILON) {

			Ms = pcuMesh->pMs_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			alpha = pcuMesh->palpha_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel_AFM->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			P = pcuMesh->pP->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			beta = pcuMesh->pbeta->get(T_Curie - (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs_AFM, Ms, *pcuMesh->palpha_AFM, alpha, *pcuMesh->pgrel_AFM, grel, *pcuMesh->psusrel_AFM, susrel, *pcuMesh->pP, P, *pcuMesh->pbeta, beta);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel.i > 1.0) susrel.i = 1.0;
		if (susrel.j > 1.0) susrel.j = 1.0;

		alpha_par = 2 * (pcuMesh->palpha_AFM->get0() - alpha);

		cuReal2 me = Ms / Ms0;
		cuReal2 r = m / me;

		Hl_1 = (M[idx] / (2 * (cuBReal)MU0 * Ms0.i)) * ((1.0 - r.i * r.i) / susrel.i + (3 * tau_ij.i * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.i) * ((susrel.j / susrel.i) * (1 - r.i * r.i) + (me.j / me.i) * (1 - r.j * r.j)));
		Hl_2 = (M2[idx] / (2 * (cuBReal)MU0 * Ms0.j)) * ((1.0 - r.j * r.j) / susrel.j + (3 * tau_ij.j * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.j) * ((susrel.i / susrel.j) * (1 - r.j * r.j) + (me.i / me.j) * (1 - r.i * r.i)));
	}
	else {

		if (Temperature < T_Curie + (cuBReal)TCURIE_EPSILON) {

			alpha = pcuMesh->palpha_AFM->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel_AFM->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel_AFM->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			P = pcuMesh->pP->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			beta = pcuMesh->pbeta->get(T_Curie + (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->palpha_AFM, alpha, *pcuMesh->pgrel_AFM, grel, *pcuMesh->psusrel_AFM, susrel, *pcuMesh->pP, P, *pcuMesh->pbeta, beta);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel.i > 1.0) susrel.i = 1.0;
		if (susrel.j > 1.0) susrel.j = 1.0;

		alpha_par = alpha;

		cuReal2 me = Ms / Ms0;

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		Hl_1 = -1 * (M[idx] / ((cuBReal)MU0 * Ms0.i)) * ((1.0 / susrel.i) + (3 * tau_ij.i * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.i) * ((susrel.j / susrel.i) + 1));
		Hl_2 = -1 * (M2[idx] / ((cuBReal)MU0 * Ms0.j)) * ((1.0 / susrel.j) + (3 * tau_ij.j * T_Curie * ((cuBReal)BOLTZMANN / (cuBReal)MUB) / mu.j) * ((susrel.i / susrel.j) + 1));
	}

	cuReal3 LLGSTT_Eval_A = (-(cuBReal)GAMMA * grel.i * msq.i / (msq.i + alpha.i * alpha.i)) * (M[idx] ^ Heff[idx]) + (-(cuBReal)GAMMA * grel.i * m.i * alpha.i / (msq.i + alpha.i * alpha.i)) * ((M[idx] / Mmag.i) ^ (M[idx] ^ Heff[idx])) +
		(cuBReal)GAMMA * grel.i * alpha_par.i * Ms0.i * ((M[idx] / Mmag.i) * (Heff[idx] + Hl_1)) * (M[idx] / Mmag.i);
	
	cuReal3 LLGSTT_Eval_B = (-(cuBReal)GAMMA * grel.j * msq.j / (msq.j + alpha.j * alpha.j)) * (M2[idx] ^ Heff2[idx]) + (-(cuBReal)GAMMA * grel.j * m.j * alpha.j / (msq.j + alpha.j * alpha.j)) * ((M2[idx] / Mmag.j) ^ (M2[idx] ^ Heff2[idx])) +
		(cuBReal)GAMMA * grel.j * alpha_par.j * Ms0.j * ((M2[idx] / Mmag.j) * (Heff2[idx] + Hl_2)) * (M2[idx] / Mmag.j);

	if (E.linear_size()) {

		cuSZ3 n = M.n;
		cuReal3 h = M.h;

		cuReal33 grad_M_A = M.grad_neu(idx);
		cuReal33 grad_M_B = M2.grad_neu(idx);

		cuReal3 u_A = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E) / (Ms.i * (1 + beta * beta));
		cuReal3 u_B = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E) / (Ms.j * (1 + beta * beta));

		cuReal3 u_dot_del_M_A = (u_A.x * grad_M_A.x) + (u_A.y * grad_M_A.y) + (u_A.z * grad_M_A.z);
		cuReal3 u_dot_del_M_B = (u_B.x * grad_M_B.x) + (u_B.y * grad_M_B.y) + (u_B.z * grad_M_B.z);

		cuReal2 alpha_perp_red = alpha / m;

		LLGSTT_Eval_A +=
			(((1 + alpha_perp_red.i * beta) * u_dot_del_M_A) -
			((beta - alpha_perp_red.i) * ((M[idx] / Mmag.i) ^ u_dot_del_M_A)) -
				(alpha_perp_red.i * (beta - alpha_perp_red.i) * (M[idx] / Mmag.i) * ((M[idx] / Mmag.i) * u_dot_del_M_A))) * msq.i / (msq.i + alpha.i * alpha.i);

		LLGSTT_Eval_B +=
			(((1 + alpha_perp_red.j * beta) * u_dot_del_M_B) -
			((beta - alpha_perp_red.j) * ((M2[idx] / Mmag.j) ^ u_dot_del_M_B)) -
				(alpha_perp_red.j * (beta - alpha_perp_red.j) * (M2[idx] / Mmag.j) * ((M2[idx] / Mmag.j) * u_dot_del_M_B))) * msq.j / (msq.j + alpha.j * alpha.j);
	}

	//sub-lattice B value so we can read it after
	value_B = LLGSTT_Eval_B;

	//return the sub-lattice A value as normal
	return LLGSTT_Eval_A;
}

#endif
#endif