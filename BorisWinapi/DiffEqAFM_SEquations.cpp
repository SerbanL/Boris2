#include "stdafx.h"
#include "DiffEqAFM.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//------------------------------------------------------------------------------------------------------ THERMAL VECs GENERATIONS

//
// Thermal field given as :
//
// Hth_magnitude = rand * SQRT( 2*alpha* kb*T / (MU0*|gamma_e|*V*MU0*Ms0*dT)) / alpha   (A/m) -> the thermal field is added only to the damping torque in LLG or in the LLB transverse damping torque.
// 
// For LLB thermal damping field, it is rand * SQRT( 2*(alpha_per - alpha_par)* kb*T / (MU0*|gamma_e|*V*MU0*Ms0*dT)) / alpha_per -> reduces to  the LLG version at T = 0K
// For LLB longitudinal thermal torque, its magnitude is rand * SQRT(2*|gamma_e|*Ms0* kb*T * alpha_par / V * dT) -> this is the sLLB-II version from PRB 85, 014433 (2012)
//
// Note : alpha_per = alpha0*(1 - T/3Tc) for T < Tc, 2*alpha0*T/3Tc for T>=Tc
//		  alpha_par = 2*alpha0*T/3Tc for T < Tc, alpha_par = alpha_per for T >= Tc
//
// So : alpha_per - alpha_par = alpha0 * (1 - T/Tc) for T < Tc, then 0 for T >= Tc
//
// kB has units of m^2kg / s^2K
// gamma_e has units of As/kg
// mu0 has units of N/A^2
//
// V is the volume of the mesh cell
//
// rand is a random factor between 0 and 1
//

void DifferentialEquationAFM::GenerateThermalField(void)
{
	//NOTE !!! Do not use separate distributions for theta and phi. I tried first for the polar angle a distribution from 0 to pi, and for azimuthal 0 to 2pi - It doesn't work, the resulting field tends to be polarized towards the left. I don't understand!!!!

	DBL2 grel = pMesh->grel_AFM.get0();

	if (IsNZ(grel.i + grel.j)) {

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n_s.dim(); idx++) {

			double Temperature;

			if (pMesh->Temp.linear_size()) {

				Temperature = pMesh->Temp[H_Thermal.cellidx_to_position(idx)];
			}
			else Temperature = pMesh->GetBaseTemperature();

			//do not include any damping here - this will be included in the stochastic equations
			double Hth_const = sqrt(2 * BOLTZMANN * Temperature / (GAMMA * grel.i * pMesh->h_s.dim() * MU0 * pMesh->Ms_AFM.get0().i * dT));
			double Hth_const_2 = sqrt(2 * BOLTZMANN * Temperature / (GAMMA * grel.j * pMesh->h_s.dim() * MU0 * pMesh->Ms_AFM.get0().j * dT));

			H_Thermal[idx] = Hth_const * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
			H_Thermal_2[idx] = Hth_const_2 * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
		}
	}
}

void DifferentialEquationAFM::GenerateThermalField_and_Torque(void)
{
	DBL2 grel = pMesh->grel_AFM.get0();

	if (IsNZ(grel.i + grel.j)) {

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n_s.dim(); idx++) {

			double Temperature;

			if (pMesh->Temp.linear_size()) {

				Temperature = pMesh->Temp[H_Thermal.cellidx_to_position(idx)];
			}
			else Temperature = pMesh->GetBaseTemperature();

			//1. Thermal Field

			//do not include any damping here - this will be included in the stochastic equations
			double Hth_const = sqrt(2 * BOLTZMANN * Temperature / (GAMMA * grel.i * pMesh->h_s.dim() * MU0 * pMesh->Ms_AFM.get0().i * dT));
			double Hth_const_2 = sqrt(2 * BOLTZMANN * Temperature / (GAMMA * grel.j * pMesh->h_s.dim() * MU0 * pMesh->Ms_AFM.get0().j * dT));

			H_Thermal[idx] = Hth_const * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
			H_Thermal_2[idx] = Hth_const_2 * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));

			//2. Thermal Torque

			//do not include any damping here - this will be included in the stochastic equations
			double Tth_const = sqrt(2 * BOLTZMANN * Temperature * GAMMA * grel.i * pMesh->Ms_AFM.get0().i / (MU0 * pMesh->h_s.dim() * dT));
			double Tth_const_2 = sqrt(2 * BOLTZMANN * Temperature * GAMMA * grel.j * pMesh->Ms_AFM.get0().j / (MU0 * pMesh->h_s.dim() * dT));

			Torque_Thermal[idx] = Tth_const * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
			Torque_Thermal_2[idx] = Tth_const_2 * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
		}
	}
}


//------------------------------------------------------------------------------------------------------ STOCHASTIC EQUATIONS

DBL3 DifferentialEquationAFM::SLLG(int idx)
{
	DBL2 Ms_AFM = pMesh->Ms_AFM;
	DBL2 alpha_AFM = pMesh->alpha_AFM;
	DBL2 grel_AFM = pMesh->grel_AFM;
	pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->alpha_AFM, alpha_AFM, pMesh->grel_AFM, grel_AFM);

	int tn = omp_get_thread_num();

	DBL3 position = pMesh->M.cellidx_to_position(idx);
	DBL3 H_Thermal_Value = H_Thermal[position] / sqrt(alpha_AFM.i);
	DBL3 H_Thermal_Value_2 = H_Thermal_2[position] / sqrt(alpha_AFM.j);

	//sub-lattice B value so we can read it after
	Equation_Eval_2[tn] = (-GAMMA * grel_AFM.j / (1 + alpha_AFM.j*alpha_AFM.j)) * ((pMesh->M2[idx] ^ pMesh->Heff2[idx]) + alpha_AFM.j * ((pMesh->M2[idx] / Ms_AFM.j) ^ (pMesh->M2[idx] ^ (pMesh->Heff2[idx] + H_Thermal_Value_2))));

	//return the sub-lattice A value as normal
	return (-GAMMA * grel_AFM.i / (1 + alpha_AFM.i*alpha_AFM.i)) * ((pMesh->M[idx] ^ pMesh->Heff[idx]) + alpha_AFM.i * ((pMesh->M[idx] / Ms_AFM.i) ^ (pMesh->M[idx] ^ (pMesh->Heff[idx] + H_Thermal_Value))));
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::SLLGSTT(int idx)
{
	DBL2 Ms_AFM = pMesh->Ms_AFM;
	DBL2 alpha_AFM = pMesh->alpha_AFM;
	DBL2 grel_AFM = pMesh->grel_AFM;
	double P = pMesh->P;
	double beta = pMesh->beta;
	pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->alpha_AFM, alpha_AFM, pMesh->grel_AFM, grel_AFM, pMesh->P, P, pMesh->beta, beta);

	int tn = omp_get_thread_num();

	DBL3 position = pMesh->M.cellidx_to_position(idx);
	DBL3 H_Thermal_Value = H_Thermal[position] / sqrt(alpha_AFM.i);
	DBL3 H_Thermal_Value_2 = H_Thermal_2[position] / sqrt(alpha_AFM.j);

	DBL3 LLGSTT_Eval_A = (-GAMMA * grel_AFM.i / (1 + alpha_AFM.i * alpha_AFM.i)) * ((pMesh->M[idx] ^ pMesh->Heff[idx]) + alpha_AFM.i * ((pMesh->M[idx] / Ms_AFM.i) ^ (pMesh->M[idx] ^ (pMesh->Heff[idx] + H_Thermal_Value))));
	DBL3 LLGSTT_Eval_B = (-GAMMA * grel_AFM.j / (1 + alpha_AFM.j * alpha_AFM.j)) * ((pMesh->M2[idx] ^ pMesh->Heff2[idx]) + alpha_AFM.j * ((pMesh->M2[idx] / Ms_AFM.j) ^ (pMesh->M2[idx] ^ (pMesh->Heff2[idx] + H_Thermal_Value_2))));

	if (pMesh->E.linear_size()) {

		DBL33 grad_M_A = pMesh->M.grad_neu(idx);
		DBL33 grad_M_B = pMesh->M2.grad_neu(idx);

		DBL3 u_A = (pMesh->elC[position] * pMesh->E.weighted_average(position, pMesh->h) * P * GMUB_2E) / (Ms_AFM.i * (1 + beta * beta));
		DBL3 u_B = (pMesh->elC[position] * pMesh->E.weighted_average(position, pMesh->h) * P * GMUB_2E) / (Ms_AFM.j * (1 + beta * beta));

		DBL3 u_dot_del_M_A = (u_A.x * grad_M_A.x) + (u_A.y * grad_M_A.y) + (u_A.z * grad_M_A.z);
		DBL3 u_dot_del_M_B = (u_B.x * grad_M_B.x) + (u_B.y * grad_M_B.y) + (u_B.z * grad_M_B.z);

		LLGSTT_Eval_A +=
			(((1 + alpha_AFM.i * alpha_AFM.i) * u_dot_del_M_A) -
			((beta - alpha_AFM.i) * ((pMesh->M[idx] / Ms_AFM.i) ^ u_dot_del_M_A))) / (1 + alpha_AFM.i * alpha_AFM.i);

		LLGSTT_Eval_B +=
			(((1 + alpha_AFM.j * alpha_AFM.j) * u_dot_del_M_B) -
			((beta - alpha_AFM.j) * ((pMesh->M2[idx] / Ms_AFM.j) ^ u_dot_del_M_B))) / (1 + alpha_AFM.j * alpha_AFM.j);
	}

	//sub-lattice B value so we can read it after
	Equation_Eval_2[tn] = LLGSTT_Eval_B;

	//return the sub-lattice A value as normal
	return LLGSTT_Eval_A;
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::SLLB(int idx)
{
	int tn = omp_get_thread_num();

	DBL3 position = pMesh->M.cellidx_to_position(idx);

	double T_Curie = pMesh->GetCurieTemperature();

	//cell temperature : the base temperature if uniform temperature, else get the temperature from Temp
	double Temperature;
	if (pMesh->Temp.linear_size()) Temperature = pMesh->Temp[position];
	else Temperature = pMesh->base_temperature;

	//m is M / Ms0 : magnitude of M in this cell divided by the saturation magnetization at 0K.
	DBL2 M = DBL2(pMesh->M[idx].norm(), pMesh->M2[idx].norm());
	DBL2 Ms0 = pMesh->Ms_AFM.get0();
	DBL2 m = M / Ms0;
	DBL2 mM = m & M;

	DBL2 Ms = pMesh->Ms_AFM;
	DBL2 alpha = pMesh->alpha_AFM;
	DBL2 grel = pMesh->grel_AFM;
	DBL2 susrel = pMesh->susrel_AFM;

	DBL2 tau_ii = pMesh->tau_ii;
	DBL2 tau_ij = pMesh->tau_ij;
	DBL2 mu = pMesh->atomic_moment_AFM;
	
	DBL2 alpha_par;

	//the longitudinal relaxation field - an effective field contribution, but only need to add it to the longitudinal relaxation term as the others involve cross products with pMesh->M[idx]
	DBL3 Hl_1, Hl_2;

	if (Temperature < T_Curie) {

		if (Temperature > T_Curie - TCURIE_EPSILON) {

			Ms = pMesh->Ms_AFM.get(T_Curie - TCURIE_EPSILON);
			alpha = pMesh->alpha_AFM.get(T_Curie - TCURIE_EPSILON);
			grel = pMesh->grel_AFM.get(T_Curie - TCURIE_EPSILON);
			susrel = pMesh->susrel_AFM.get(T_Curie - TCURIE_EPSILON);
		}
		else pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms, pMesh->alpha_AFM, alpha, pMesh->grel_AFM, grel, pMesh->susrel_AFM, susrel);

		alpha_par = 2 * (pMesh->alpha_AFM.get0() - alpha);

		DBL2 me = Ms / Ms0;
		DBL2 r = m / me;
		double proj = pMesh->M[idx].normalized() * pMesh->M2[idx].normalized();

		Hl_1 = (pMesh->M[idx] / (2 * MU0 * Ms0.i)) * ((1.0 - r.i * r.i) / susrel.i + (3 * tau_ij.i * T_Curie * (BOLTZMANN / MUB) / mu.i) * ((susrel.j / susrel.i) * (1 - r.i * r.i) - proj * (me.j / me.i) * (1 - r.j * r.j)));
		Hl_2 = (pMesh->M2[idx] / (2 * MU0 * Ms0.j)) * ((1.0 - r.j * r.j) / susrel.j + (3 * tau_ij.j * T_Curie * (BOLTZMANN / MUB) / mu.j) * ((susrel.i / susrel.j) * (1 - r.j * r.j) - proj * (me.i / me.j) * (1 - r.i * r.i)));
	}
	else {

		if (Temperature < T_Curie + TCURIE_EPSILON) {

			alpha = pMesh->alpha_AFM.get(T_Curie + TCURIE_EPSILON);
			grel = pMesh->grel_AFM.get(T_Curie + TCURIE_EPSILON);
			susrel = pMesh->susrel_AFM.get(T_Curie + TCURIE_EPSILON);
		}
		else pMesh->update_parameters_mcoarse(idx, pMesh->alpha_AFM, alpha, pMesh->grel_AFM, grel, pMesh->susrel_AFM, susrel);

		alpha_par = alpha;

		DBL2 me = Ms / Ms0;
		double proj = pMesh->M[idx].normalized() * pMesh->M2[idx].normalized();

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		Hl_1 = -1 * (pMesh->M[idx] / (MU0 * Ms0.i)) * ((1.0 / susrel.i) + (3 * tau_ij.i * T_Curie * (BOLTZMANN / MUB) / mu.i) * ((susrel.j / susrel.i) - proj * (me.j / me.i)));
		Hl_2 = -1 * (pMesh->M2[idx] / (MU0 * Ms0.j)) * ((1.0 / susrel.j) + (3 * tau_ij.j * T_Curie * (BOLTZMANN / MUB) / mu.j) * ((susrel.i / susrel.j) - proj * (me.i / me.j)));
	}

	DBL3 H_Thermal_Value = H_Thermal[position] * sqrt(alpha.i - alpha_par.i) / alpha.i;
	DBL3 Torque_Thermal_Value = Torque_Thermal[position] * sqrt(alpha_par.i);

	DBL3 H_Thermal_Value_2 = H_Thermal_2[position] * sqrt(alpha.j - alpha_par.j) / alpha.j;
	DBL3 Torque_Thermal_Value_2 = Torque_Thermal_2[position] * sqrt(alpha_par.j);

	//sub-lattice B value so we can read it after
	Equation_Eval_2[tn] = (-GAMMA * grel.j / (1 + alpha.j * alpha.j)) * ((pMesh->M2[idx] ^ pMesh->Heff2[idx]) + alpha.j * ((pMesh->M2[idx] / mM.j) ^ (pMesh->M2[idx] ^ (pMesh->Heff2[idx] + H_Thermal_Value_2)))) +
		GAMMA * grel.j * alpha_par.j * (pMesh->M2[idx] * (pMesh->Heff2[idx] + Hl_2)) * (pMesh->M2[idx] / mM.j) + Torque_Thermal_Value_2;

	//return the sub-lattice A value as normal
	return (-GAMMA * grel.i / (1 + alpha.i * alpha.i)) * ((pMesh->M[idx] ^ pMesh->Heff[idx]) + alpha.i * ((pMesh->M[idx] / mM.i) ^ (pMesh->M[idx] ^ (pMesh->Heff[idx] + H_Thermal_Value)))) +
		GAMMA * grel.i * alpha_par.i * (pMesh->M[idx] * (pMesh->Heff[idx] + Hl_1)) * (pMesh->M[idx] / mM.i) + Torque_Thermal_Value;
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::SLLBSTT(int idx)
{
	int tn = omp_get_thread_num();

	DBL3 position = pMesh->M.cellidx_to_position(idx);

	double T_Curie = pMesh->GetCurieTemperature();

	//cell temperature : the base temperature if uniform temperature, else get the temperature from Temp
	double Temperature;
	if (pMesh->Temp.linear_size()) Temperature = pMesh->Temp[position];
	else Temperature = pMesh->base_temperature;

	//m is M / Ms0 : magnitude of M in this cell divided by the saturation magnetization at 0K.
	DBL2 M = DBL2(pMesh->M[idx].norm(), pMesh->M2[idx].norm());
	DBL2 Ms0 = pMesh->Ms_AFM.get0();
	DBL2 m = M / Ms0;
	DBL2 mM = m & M;

	DBL2 Ms = pMesh->Ms_AFM;
	DBL2 alpha = pMesh->alpha_AFM;
	DBL2 grel = pMesh->grel_AFM;
	DBL2 susrel = pMesh->susrel_AFM;
	double P = pMesh->P;
	double beta = pMesh->beta;

	DBL2 tau_ii = pMesh->tau_ii;
	DBL2 tau_ij = pMesh->tau_ij;
	DBL2 mu = pMesh->atomic_moment_AFM;

	DBL2 alpha_par;

	//the longitudinal relaxation field - an effective field contribution, but only need to add it to the longitudinal relaxation term as the others involve cross products with pMesh->M[idx]
	DBL3 Hl_1, Hl_2;

	if (Temperature < T_Curie) {

		if (Temperature > T_Curie - TCURIE_EPSILON) {

			Ms = pMesh->Ms_AFM.get(T_Curie - TCURIE_EPSILON);
			alpha = pMesh->alpha_AFM.get(T_Curie - TCURIE_EPSILON);
			grel = pMesh->grel_AFM.get(T_Curie - TCURIE_EPSILON);
			susrel = pMesh->susrel_AFM.get(T_Curie - TCURIE_EPSILON);
			P = pMesh->P.get(T_Curie - TCURIE_EPSILON);
			beta = pMesh->beta.get(T_Curie - TCURIE_EPSILON);
		}
		else pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms, pMesh->alpha_AFM, alpha, pMesh->grel_AFM, grel, pMesh->susrel_AFM, susrel, pMesh->P, P, pMesh->beta, beta);

		alpha_par = 2 * (pMesh->alpha_AFM.get0() - alpha);

		DBL2 me = Ms / Ms0;
		DBL2 r = m / me;
		double proj = pMesh->M[idx].normalized() * pMesh->M2[idx].normalized();

		Hl_1 = (pMesh->M[idx] / (2 * MU0 * Ms0.i)) * ((1.0 - r.i * r.i) / susrel.i + (3 * tau_ij.i * T_Curie * (BOLTZMANN / MUB) / mu.i) * ((susrel.j / susrel.i) * (1 - r.i * r.i) - proj * (me.j / me.i) * (1 - r.j * r.j)));
		Hl_2 = (pMesh->M2[idx] / (2 * MU0 * Ms0.j)) * ((1.0 - r.j * r.j) / susrel.j + (3 * tau_ij.j * T_Curie * (BOLTZMANN / MUB) / mu.j) * ((susrel.i / susrel.j) * (1 - r.j * r.j) - proj * (me.i / me.j) * (1 - r.i * r.i)));
	}
	else {

		if (Temperature < T_Curie + TCURIE_EPSILON) {

			alpha = pMesh->alpha_AFM.get(T_Curie + TCURIE_EPSILON);
			grel = pMesh->grel_AFM.get(T_Curie + TCURIE_EPSILON);
			susrel = pMesh->susrel_AFM.get(T_Curie + TCURIE_EPSILON);
			P = pMesh->P.get(T_Curie + TCURIE_EPSILON);
			beta = pMesh->beta.get(T_Curie + TCURIE_EPSILON);
		}
		else pMesh->update_parameters_mcoarse(idx, pMesh->alpha_AFM, alpha, pMesh->grel_AFM, grel, pMesh->susrel_AFM, susrel, pMesh->P, P, pMesh->beta, beta);

		alpha_par = alpha;

		DBL2 me = Ms / Ms0;
		double proj = pMesh->M[idx].normalized() * pMesh->M2[idx].normalized();

		//Note, the parallel susceptibility is related to susrel by : susrel = suspar / mu0Ms
		Hl_1 = -1 * (pMesh->M[idx] / (MU0 * Ms0.i)) * ((1.0 / susrel.i) + (3 * tau_ij.i * T_Curie * (BOLTZMANN / MUB) / mu.i) * ((susrel.j / susrel.i) - proj * (me.j / me.i)));
		Hl_2 = -1 * (pMesh->M2[idx] / (MU0 * Ms0.j)) * ((1.0 / susrel.j) + (3 * tau_ij.j * T_Curie * (BOLTZMANN / MUB) / mu.j) * ((susrel.i / susrel.j) - proj * (me.i / me.j)));
	}

	DBL3 H_Thermal_Value = H_Thermal[position] * sqrt(alpha.i - alpha_par.i) / alpha.i;
	DBL3 Torque_Thermal_Value = Torque_Thermal[position] * sqrt(alpha_par.i);

	DBL3 H_Thermal_Value_2 = H_Thermal_2[position] * sqrt(alpha.j - alpha_par.j) / alpha.j;
	DBL3 Torque_Thermal_Value_2 = Torque_Thermal_2[position] * sqrt(alpha_par.j);

	DBL3 LLGSTT_Eval_A = (-GAMMA * grel.i / (1 + alpha.i * alpha.i)) * ((pMesh->M[idx] ^ pMesh->Heff[idx]) + alpha.i * ((pMesh->M[idx] / mM.i) ^ (pMesh->M[idx] ^ (pMesh->Heff[idx] + H_Thermal_Value)))) +
		GAMMA * grel.i * alpha_par.i * (pMesh->M[idx] * (pMesh->Heff[idx] + Hl_1)) * (pMesh->M[idx] / mM.i) + Torque_Thermal_Value;

	DBL3 LLGSTT_Eval_B = (-GAMMA * grel.j / (1 + alpha.j * alpha.j)) * ((pMesh->M2[idx] ^ pMesh->Heff2[idx]) + alpha.j * ((pMesh->M2[idx] / mM.j) ^ (pMesh->M2[idx] ^ (pMesh->Heff2[idx] + H_Thermal_Value_2)))) +
		GAMMA * grel.j * alpha_par.j * (pMesh->M2[idx] * (pMesh->Heff2[idx] + Hl_2)) * (pMesh->M2[idx] / mM.j) + Torque_Thermal_Value_2;

	if (pMesh->E.linear_size()) {

		DBL33 grad_M_A = pMesh->M.grad_neu(idx);
		DBL33 grad_M_B = pMesh->M2.grad_neu(idx);

		DBL3 u_A = (pMesh->elC[position] * pMesh->E.weighted_average(position, pMesh->h) * P * GMUB_2E) / (Ms.i * (1 + beta * beta));
		DBL3 u_B = (pMesh->elC[position] * pMesh->E.weighted_average(position, pMesh->h) * P * GMUB_2E) / (Ms.j * (1 + beta * beta));

		DBL3 u_dot_del_M_A = (u_A.x * grad_M_A.x) + (u_A.y * grad_M_A.y) + (u_A.z * grad_M_A.z);
		DBL3 u_dot_del_M_B = (u_B.x * grad_M_B.x) + (u_B.y * grad_M_B.y) + (u_B.z * grad_M_B.z);

		DBL2 alpha_perp_red = alpha / m;

		LLGSTT_Eval_A +=
			(((1 + alpha_perp_red.i * beta) * u_dot_del_M_A) -
			((beta - alpha_perp_red.i) * ((pMesh->M[idx] / M.i) ^ u_dot_del_M_A)) -
				(alpha_perp_red.i * (beta - alpha_perp_red.i) * (pMesh->M[idx] / M.i) * ((pMesh->M[idx] / M.i) * u_dot_del_M_A))) / (1 + alpha.i * alpha.i);

		LLGSTT_Eval_B +=
			(((1 + alpha_perp_red.j * beta) * u_dot_del_M_B) -
			((beta - alpha_perp_red.j) * ((pMesh->M2[idx] / M.j) ^ u_dot_del_M_B)) -
				(alpha_perp_red.j * (beta - alpha_perp_red.j) * (pMesh->M2[idx] / M.j) * ((pMesh->M2[idx] / M.j) * u_dot_del_M_B))) / (1 + alpha.j * alpha.j);
	}

	//sub-lattice B value so we can read it after
	Equation_Eval_2[tn] = LLGSTT_Eval_B;

	//return the sub-lattice A value as normal
	return LLGSTT_Eval_A;
}

#endif
