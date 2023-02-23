#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

//Defines stochastic equations

#include "BorisCUDALib.h"

#include "DiffEqFMCUDA.h"
#include "ManagedDiffEqFMCUDA.h"

#include "DiffEq_Defs.h"

#include "Funcs_Math_base.h" //includes constant values

#include "MeshParamsControlCUDA.h"

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqFMCUDA::SLLG(int idx)
{
	//gamma = -mu0 * gamma_e = mu0 * g e / 2m_e = 2.212761569e5 m/As

	//LLG in explicit form : dm/dt = [mu0*gamma_e/(1+alpha^2)] * [m*H + alpha * m*(m*H)]

	//Add thermal field to damping term, remembering to include damping contribution which was not included when H_Thermal was generated

	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuBReal Ms = *pcuMesh->pMs;
	cuBReal alpha = *pcuMesh->palpha;
	cuBReal grel = *pcuMesh->pgrel;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->palpha, alpha, *pcuMesh->pgrel, grel);

	cuReal3 position = M.cellidx_to_position(idx);
	cuReal3 H_Thermal_Value = (*pH_Thermal)[position] * sqrt(alpha);

	return (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ (Heff[idx] + H_Thermal_Value)) +
		alpha * ((M[idx] / Ms) ^ (M[idx] ^ (Heff[idx] + H_Thermal_Value))));
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqFMCUDA::SLLGSTT(int idx)
{
	//gmub_2e is -hbar * gamma_e / 2e = g mu_b / 2e)

	// LLG with STT in explicit form : dm/dt = [mu0*gamma_e/(1+alpha^2)] * [m*H + alpha * m*(m*H)] + (1+alpha*beta)/((1+alpha^2)*(1+beta^2)) * (u.del)m - (beta - alpha)/(1+alpha^2) * m * (u.del) m
	// where u = j * P g mu_b / 2e Ms = -(hbar * gamma_e * P / 2 *e * Ms) * j, j is the current density = conductivity * E (A/m^2)

	// STT is Zhang-Li equationtion (not Thiaville, the velocity used by Thiaville needs to be divided by (1+beta^2) to obtain Zhang-Li, also Thiaville's EPL paper has wrong STT signs!!)

	//Add thermal field to damping term, remembering to include damping contribution which was not included when H_Thermal was generated

	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;
	cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
	cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

	cuBReal Ms = *pcuMesh->pMs;
	cuBReal alpha = *pcuMesh->palpha;
	cuBReal grel = *pcuMesh->pgrel;
	cuBReal P = *pcuMesh->pP;
	cuBReal beta = *pcuMesh->pbeta;
	pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->palpha, alpha, *pcuMesh->pgrel, grel, *pcuMesh->pP, P, *pcuMesh->pbeta, beta);

	cuReal3 position = M.cellidx_to_position(idx);
	cuReal3 H_Thermal_Value = (*pH_Thermal)[position] * sqrt(alpha);

	cuReal3 LLGSTT_Eval = (-(cuBReal)GAMMA * grel / (1 + alpha * alpha)) * ((M[idx] ^ (Heff[idx] + H_Thermal_Value)) +
		alpha * ((M[idx] / Ms) ^ (M[idx] ^ (Heff[idx] + H_Thermal_Value))));

	if (E.linear_size()) {

		cuSZ3 n = M.n;
		cuReal3 h = M.h;

		cuReal33 grad_M = M.grad_neu(idx);

		cuReal3 u = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E) / (Ms * (1 + beta * beta));

		cuReal3 u_dot_del_M = (u.x * grad_M.x) + (u.y * grad_M.y) + (u.z * grad_M.z);

		LLGSTT_Eval +=
			(((1 + alpha * beta) * u_dot_del_M) -
			((beta - alpha) * ((M[idx] / Ms) ^ u_dot_del_M))) / (1 + alpha * alpha);
	}

	return LLGSTT_Eval;
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqFMCUDA::SLLB(int idx)
{
	//gamma = -mu0 * gamma_e = mu0 * g e / 2m_e = 2.212761569e5 m/As

	//LLB in explicit form : dM/dt = [mu0*gamma_e/(1+alpha_perp_red^2)] * [M*H + alpha_perp_red * (M/|M|)*(M*H)] - mu0*gamma_e* alpha_par_red * (M.(H + Hl)) * (M/|M|)

	//alpha_perp_red = alpha / m
	//alpha_par_red = 2*(alpha0 - alpha)/m up to Tc, then alpha_par_red = alpha_perp_red above Tc, where alpha0 is the zero temperature damping and alpha is the damping at a given temperature
	//m = |M| / Ms0, where Ms0 is the zero temperature saturation magnetization
	//
	//There is a longitudinal relaxation field Hl = M * (1 - (|M|/Ms)^2) / (2*suspar), where Ms is the equilibrium magnetization (i.e. the "saturation" magnetization at the given temperature - obtained from Ms)
	//
	//Ms, suspar and alpha must have temperature dependence set. In particular:
	//alpha = alpha0 * (1 - T/3Tc) up to Tc, alpha = (2*alpha0*T/3Tc) above Tc
	//For Ms and suspar see literature (e.g. S.Lepadatu, JAP 120, 163908 (2016))

	//Add thermal field to damping term, and thermal torque to evaluation, remembering to include damping contribution which was not included when H_Thermal was generated
	
	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;

	cuReal3 position = M.cellidx_to_position(idx);

	cuBReal T_Curie = *pcuMesh->pT_Curie;

	//cell temperature : the base temperature if uniform temperature, else get the temperature from Temp
	cuBReal Temperature;
	if (pcuMesh->pTemp->linear_size()) Temperature = (*pcuMesh->pTemp)[position];
	else Temperature = *pcuMesh->pbase_temperature;

	//m is M / Ms0 : magnitude of M in this cell divided by the saturation magnetization at 0K.
	cuBReal Mnorm = M[idx].norm();
	cuBReal Ms0 = pcuMesh->pMs->get0();
	cuBReal m = Mnorm / Ms0;
	cuBReal msq = m * m;

	cuBReal Ms = *pcuMesh->pMs;
	cuBReal alpha = *pcuMesh->palpha;
	cuBReal grel = *pcuMesh->pgrel;
	cuBReal susrel = *pcuMesh->psusrel;

	cuBReal alpha_par;

	//the longitudinal relaxation field - an effective field contribution, but only need to add it to the longitudinal relaxation term as the others involve cross products with M[idx]
	cuReal3 Hl = cuReal3(0.0);

	if (Temperature <= T_Curie) {
		
		if (Temperature > T_Curie - (cuBReal)TCURIE_EPSILON) {

			Ms = pcuMesh->pMs->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			alpha = pcuMesh->palpha->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel->get(T_Curie - (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->palpha, alpha, *pcuMesh->pgrel, grel, *pcuMesh->psusrel, susrel);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel > 1.0) susrel = 1.0;

		alpha_par = 2 * (pcuMesh->palpha->get0() - alpha);

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		Hl = M[idx] * ((1 - (Mnorm / Ms) * (Mnorm / Ms)) / (2 * susrel * (cuBReal)MU0 * Ms0));
	}
	else {

		if (Temperature < T_Curie + (cuBReal)TCURIE_EPSILON) {

			alpha = pcuMesh->palpha->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel->get(T_Curie + (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->palpha, alpha, *pcuMesh->pgrel, grel, *pcuMesh->psusrel, susrel);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel > 1.0) susrel = 1.0;

		alpha_par = alpha;

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		if (Temperature < T_Curie + (cuBReal)TCURIE_EPSILON) Hl = -1.0 * (M[idx] / (susrel * (cuBReal)MU0 * Ms0));
		else Hl = -1.0 * (M[idx] / (susrel * (cuBReal)MU0 * Ms0)) * (1 + 3 * msq * T_Curie / (5 * (Temperature - T_Curie)));
	}

	cuReal3 H_Thermal_Value = (*pH_Thermal)[position] * sqrt(alpha - alpha_par) / alpha;
	cuReal3 Torque_Thermal_Value = (*pTorque_Thermal)[position] * sqrt(alpha_par);

	return (-(cuBReal)GAMMA * grel * msq / (msq + alpha * alpha)) * (M[idx] ^ Heff[idx]) + (-(cuBReal)GAMMA * grel * m * alpha / (msq + alpha * alpha)) * ((M[idx] / Mnorm) ^ (M[idx] ^ (Heff[idx] + H_Thermal_Value))) +
		(cuBReal)GAMMA * grel * alpha_par * Ms0 * ((M[idx] / Mnorm) * (Heff[idx] + Hl)) * (M[idx] / Mnorm) + Torque_Thermal_Value;
}

//------------------------------------------------------------------------------------------------------

__device__ cuReal3 ManagedDiffEqFMCUDA::SLLBSTT(int idx)
{
	//gamma = -mu0 * gamma_e = mu0 * g e / 2m_e = 2.212761569e5 m/As

	//LLB in explicit form : dM/dt = [mu0*gamma_e/(1+alpha_perp_red^2)] * [M*H + alpha_perp_red * (M/|M|)*(M*H)] - mu0*gamma_e* alpha_par_red * (M.(H + Hl)) * (M/|M|)

	//alpha_perp_red = alpha / m
	//alpha_par_red = 2*(alpha0 - alpha)/m up to Tc, then alpha_par_red = alpha_perp_red above Tc, where alpha0 is the zero temperature damping and alpha is the damping at a given temperature
	//m = |M| / Ms0, where Ms0 is the zero temperature saturation magnetization
	//
	//There is a longitudinal relaxation field Hl = M * (1 - (|M|/Ms)^2) / (2*suspar), where Ms is the equilibrium magnetization (i.e. the "saturation" magnetization at the given temperature - obtained from Ms)
	//
	//Ms, suspar and alpha must have temperature dependence set. In particular:
	//alpha = alpha0 * (1 - T/3Tc) up to Tc, alpha = (2*alpha0*T/3Tc) above Tc
	//For Ms and suspar see literature (e.g. S.Lepadatu, JAP 120, 163908 (2016))

	//on top of this we have STT contributions

	//Add thermal field to damping term, and thermal torque to evaluation, remembering to include damping contribution which was not included when H_Thermal was generated

	cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
	cuVEC<cuReal3>& Heff = *pcuMesh->pHeff;
	cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
	cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

	cuReal3 position = M.cellidx_to_position(idx);

	cuBReal T_Curie = *pcuMesh->pT_Curie;

	//cell temperature : the base temperature if uniform temperature, else get the temperature from Temp
	cuBReal Temperature;
	if (pcuMesh->pTemp->linear_size()) Temperature = (*pcuMesh->pTemp)[position];
	else Temperature = *pcuMesh->pbase_temperature;

	//m is M / Ms0 : magnitude of M in this cell divided by the saturation magnetization at 0K.
	cuBReal Mnorm = M[idx].norm();
	cuBReal Ms0 = pcuMesh->pMs->get0();
	cuBReal mM = Mnorm * Mnorm / Ms0;
	cuBReal m = Mnorm / Ms0;
	cuBReal msq = m * m;

	cuBReal Ms = *pcuMesh->pMs;
	cuBReal alpha = *pcuMesh->palpha;
	cuBReal grel = *pcuMesh->pgrel;
	cuBReal susrel = *pcuMesh->psusrel;
	cuBReal P = *pcuMesh->pP;
	cuBReal beta = *pcuMesh->pbeta;

	cuBReal alpha_par;

	//the longitudinal relaxation field - an effective field contribution, but only need to add it to the longitudinal relaxation term as the others involve cross products with M[idx]
	cuReal3 Hl = cuReal3(0.0);

	if (Temperature <= T_Curie) {

		if (Temperature > T_Curie - (cuBReal)TCURIE_EPSILON) {

			Ms = pcuMesh->pMs->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			alpha = pcuMesh->palpha->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			P = pcuMesh->pP->get(T_Curie - (cuBReal)TCURIE_EPSILON);
			beta = pcuMesh->pbeta->get(T_Curie - (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->palpha, alpha, *pcuMesh->pgrel, grel, *pcuMesh->psusrel, susrel, *pcuMesh->pP, P, *pcuMesh->pbeta, beta);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel > 1.0) susrel = 1.0;

		alpha_par = 2 * (pcuMesh->palpha->get0() - alpha);

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		Hl = M[idx] * ((1 - (Mnorm / Ms) * (Mnorm / Ms)) / (2 * susrel * (cuBReal)MU0 * Ms0));
	}
	else {

		if (Temperature < T_Curie + (cuBReal)TCURIE_EPSILON) {

			alpha = pcuMesh->palpha->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			grel = pcuMesh->pgrel->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			susrel = pcuMesh->psusrel->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			P = pcuMesh->pP->get(T_Curie + (cuBReal)TCURIE_EPSILON);
			beta = pcuMesh->pbeta->get(T_Curie + (cuBReal)TCURIE_EPSILON);
		}
		else pcuMesh->update_parameters_mcoarse(idx, *pcuMesh->palpha, alpha, *pcuMesh->pgrel, grel, *pcuMesh->psusrel, susrel, *pcuMesh->pP, P, *pcuMesh->pbeta, beta);

		//limit susrel as it is difficult to solve numerically for large values in single floating point precision
		if (susrel > 1.0) susrel = 1.0;

		alpha_par = alpha;

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		if (Temperature < T_Curie + (cuBReal)TCURIE_EPSILON) Hl = -1.0 * (M[idx] / (susrel * (cuBReal)MU0 * Ms0));
		else Hl = -1.0 * (M[idx] / (susrel * (cuBReal)MU0 * Ms0)) * (1 + 3 * msq * T_Curie / (5 * (Temperature - T_Curie)));
	}

	cuReal3 H_Thermal_Value = (*pH_Thermal)[position] * sqrt(alpha - alpha_par) / alpha;
	cuReal3 Torque_Thermal_Value = (*pTorque_Thermal)[position] * sqrt(alpha_par);

	cuReal3 LLBSTT_Eval =
		(-(cuBReal)GAMMA * grel * msq / (msq + alpha * alpha)) * (M[idx] ^ Heff[idx]) + (-(cuBReal)GAMMA * grel * m * alpha / (msq + alpha * alpha)) * ((M[idx] / Mnorm) ^ (M[idx] ^ (Heff[idx] + H_Thermal_Value))) +
		(cuBReal)GAMMA * grel * alpha_par * Ms0 * ((M[idx] / Mnorm) * (Heff[idx] + Hl)) * (M[idx] / Mnorm) + Torque_Thermal_Value;

	if (E.linear_size()) {

		cuSZ3 n = M.n;
		cuReal3 h = M.h;

		cuReal33 grad_M = M.grad_neu(idx);

		cuReal3 u = (elC[position] * E.weighted_average(position, h) * P * (cuBReal)GMUB_2E) / (Ms * (1 + beta * beta));

		cuReal3 u_dot_del_M = (u.x * grad_M.x) + (u.y * grad_M.y) + (u.z * grad_M.z);

		cuBReal alpha_perp_red = alpha / m;

		LLBSTT_Eval +=
			(((1 + alpha_perp_red * beta) * u_dot_del_M) -
			((beta - alpha_perp_red) * ((M[idx] / Mnorm) ^ u_dot_del_M)) -
				(alpha_perp_red * (beta - alpha_perp_red) * (M[idx] / Mnorm) * ((M[idx] / Mnorm) * u_dot_del_M))) * msq / (msq + alpha * alpha);
	}

	return LLBSTT_Eval;
}

#endif
#endif