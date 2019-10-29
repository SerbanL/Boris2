#include "stdafx.h"

#include "DiffEqAFM.h"
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
	//get a random number generator

	//prepare random number distributions : magnitude scaling from 0 to 1 for thermal field strength, direction angles for thermal field
	//NOTE !!! Do not use separate distributions for theta and phi. I tried first for the polar angle a distribution from 0 to pi, and for azimuthal 0 to 2pi - It doesn't work, the resulting field tends to be polarized towards the left. I don't understand!!!!

	double grel = pMesh->grel.get0();

	if (IsNZ(grel)) {

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

				double Temperature;

				if (pMesh->Temp.linear_size()) {

					Temperature = pMesh->Temp[pMesh->M.cellidx_to_position(idx)];
				}
				else Temperature = pMesh->GetBaseTemperature();

				//do not include any damping here - this will be included in the stochastic equations
				double Mag = prng.rand() * sqrt(2 * BOLTZMANN * Temperature / (GAMMA * grel * pMesh->h.dim() * MU0 * pMesh->Ms.get0() * dT));
				double theta = prng.rand() * TWO_PI;
				double phi = prng.rand() * TWO_PI;

				H_Thermal[idx] = Mag * DBL3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
			}
		}
	}
}

void DifferentialEquationAFM::GenerateThermalField_and_Torque(void)
{
	//get a random number generator

	//prepare random number distributions : magnitude scaling from 0 to 1 for thermal field strength, direction angles for thermal field

	double grel = pMesh->grel.get0();

	if (IsNZ(grel)) {

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

				double Temperature;

				if (pMesh->Temp.linear_size()) {

					Temperature = pMesh->Temp[pMesh->M.cellidx_to_position(idx)];
				}
				else Temperature = pMesh->GetBaseTemperature();

				//1. Thermal Field
				//do not include any damping here - this will be included in the stochastic equations
				double Mag = prng.rand() * sqrt(2 * BOLTZMANN * Temperature / (GAMMA * grel * pMesh->h.dim() * MU0 * pMesh->Ms.get0() * dT));
				double theta = prng.rand() * TWO_PI;
				double phi = prng.rand() * TWO_PI;

				H_Thermal[idx] = Mag * DBL3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

				//2. Thermal Torque
				//do not include any damping here - this will be included in the stochastic equations
				Mag = prng.rand() * sqrt(2 * BOLTZMANN * Temperature * GAMMA * grel * pMesh->Ms.get0() / (MU0 * pMesh->h.dim() * dT));
				theta = prng.rand() * TWO_PI;
				phi = prng.rand() * TWO_PI;

				Torque_Thermal[idx] = Mag * DBL3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
			}
		}
	}
}


//------------------------------------------------------------------------------------------------------ STOCHASTIC EQUATIONS

DBL3 DifferentialEquationAFM::SLLG(int idx)
{
	return DBL3();
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::SLLGSTT(int idx)
{
	return DBL3();
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::SLLB(int idx)
{
	return DBL3();
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::SLLBSTT(int idx)
{
	return DBL3();
}