#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//------------------------------------------------------------------------------------------------------

DBL3 Atom_DifferentialEquationCubic::LLG(int idx)
{
	//gamma_e is the electron gyromagnetic ratio, which is negative

	//LLG in explicit form : dm/dt = [gamma_e/(1+alpha^2)] * [m*B + alpha * m*(m*B) / mu_s]
	//B is effective field in Tesla; prefer to include the effective field in A/m as it integrates easier in a multiscale simulation; thus gamma = |gamma_e| * mu0
	//m is atomic moment in units of muB : has magnitude mu_s (muB)
	
	double mu_s = paMesh->mu_s;
	double alpha = paMesh->alpha;
	double grel = paMesh->grel;
	paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->alpha, alpha, paMesh->grel, grel);

	return (-GAMMA * grel / (1 + alpha*alpha)) * ((paMesh->M1[idx] ^ paMesh->Heff1[idx]) + alpha * ((paMesh->M1[idx] / mu_s) ^ (paMesh->M1[idx] ^ paMesh->Heff1[idx])));
}

//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
DBL3 Atom_DifferentialEquationCubic::LLGStatic(int idx)
{
	double mu_s = paMesh->mu_s;
	double grel = paMesh->grel;
	paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->grel, grel);

	return (-GAMMA * grel / 2) * ((paMesh->M1[idx] / mu_s) ^ (paMesh->M1[idx] ^ paMesh->Heff1[idx]));
}

DBL3 Atom_DifferentialEquationCubic::LLGSTT(int idx)
{
	//Currently same as LLG : TO DO

	//gamma_e is the electron gyromagnetic ratio, which is negative

	//LLG in explicit form : dm/dt = [gamma_e/(1+alpha^2)] * [m*B + alpha * m*(m*B) / mu_s]
	//B is effective field in Tesla; prefer to include the effective field in A/m as it integrates easier in a multiscale simulation; thus gamma = |gamma_e| * mu0
	//m is atomic moment in units of muB : has magnitude mu_s (muB)

	double mu_s = paMesh->mu_s;
	double alpha = paMesh->alpha;
	double grel = paMesh->grel;
	paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->alpha, alpha, paMesh->grel, grel);

	return (-GAMMA * grel / (1 + alpha * alpha)) * ((paMesh->M1[idx] ^ paMesh->Heff1[idx]) + alpha * ((paMesh->M1[idx] / mu_s) ^ (paMesh->M1[idx] ^ paMesh->Heff1[idx])));
}

#endif