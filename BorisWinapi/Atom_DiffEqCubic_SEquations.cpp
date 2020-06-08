#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#if ATOMISTIC == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//------------------------------------------------------------------------------------------------------ THERMAL VECs GENERATIONS

//
// Thermal field given as :
//
// Bth = rand * SQRT( 2*alpha* kB*T / (|gamma_e|*mu_s*dT)) (T)

void Atom_DifferentialEquationCubic::GenerateThermalField(void)
{
	//if not in linked dTstoch mode, then only generate stochastic field at a minimum of dTstoch spacing
	if (!link_dTstoch && GetTime() < time_stoch + dTstoch) return;

	double deltaT = (link_dTstoch ? dT : GetTime() - time_stoch);
	time_stoch = GetTime();

	double Temperature = paMesh->GetBaseTemperature();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->Temp.linear_size()) Temperature = paMesh->Temp[H_Thermal.cellidx_to_position(idx)];

		double mu_s = paMesh->mu_s;
		paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);

		//do not include any damping here - this will be included in the stochastic equations
		double Hth_const = sqrt(2 * BOLTZMANN * Temperature / (MUB_MU0 * GAMMA * mu_s * deltaT));
				
		H_Thermal[idx] = Hth_const * DBL3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
	}
}

//------------------------------------------------------------------------------------------------------ STOCHASTIC EQUATIONS

DBL3 Atom_DifferentialEquationCubic::SLLG(int idx)
{
	//gamma_e is the electron gyromagnetic ratio, which is negative

	//LLG in explicit form : dm/dt = [gamma_e/(1+alpha^2)] * [m*B + alpha * m*(m*B) / mu_s]
	//B is effective field in Tesla; prefer to include the effective field in A/m as it integrates easier in a multiscale simulation; thus gamma = |gamma_e| * mu0
	//m is atomic moment in units of muB : has magnitude mu_s (muB)

	double mu_s = paMesh->mu_s;
	double alpha = paMesh->alpha;
	paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->alpha, alpha);

	//H_Thermal has same dimensions as M1 in atomistic meshes
	DBL3 H_Thermal_Value = H_Thermal[idx] * sqrt(alpha);

	return (-GAMMA / (1 + alpha * alpha)) * ((paMesh->M1[idx] ^ (paMesh->Heff1[idx] + H_Thermal_Value)) + alpha * ((paMesh->M1[idx] / mu_s) ^ (paMesh->M1[idx] ^ (paMesh->Heff1[idx] + H_Thermal_Value))));
}

DBL3 Atom_DifferentialEquationCubic::SLLGSTT(int idx)
{
	//Currently same as SLLG : TO DO

	//gamma_e is the electron gyromagnetic ratio, which is negative

	//LLG in explicit form : dm/dt = [gamma_e/(1+alpha^2)] * [m*B + alpha * m*(m*B) / mu_s]
	//B is effective field in Tesla; prefer to include the effective field in A/m as it integrates easier in a multiscale simulation; thus gamma = |gamma_e| * mu0
	//m is atomic moment in units of muB : has magnitude mu_s (muB)

	double mu_s = paMesh->mu_s;
	double alpha = paMesh->alpha;
	paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->alpha, alpha);

	//H_Thermal has same dimensions as M1 in atomistic meshes
	DBL3 H_Thermal_Value = H_Thermal[idx] * sqrt(alpha);

	return (-GAMMA / (1 + alpha * alpha)) * ((paMesh->M1[idx] ^ (paMesh->Heff1[idx] + H_Thermal_Value)) + alpha * ((paMesh->M1[idx] / mu_s) ^ (paMesh->M1[idx] ^ (paMesh->Heff1[idx] + H_Thermal_Value))));
}
#endif
#endif