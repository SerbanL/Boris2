#include "stdafx.h"
#include "DiffEqAFM.h"
#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::LLG(int idx)
{
	//gamma = -mu0 * gamma_e = mu0 * g e / 2m_e = 2.212761569e5 m/As

	//LLG in explicit form : dm/dt = [mu0*gamma_e/(1+alpha^2)] * [m*H + alpha * m*(m*H)]
	
	DBL2 Ms_AFM = pMesh->Ms_AFM;
	DBL2 alpha_AFM = pMesh->alpha_AFM;
	DBL2 grel_AFM = pMesh->grel_AFM;
	pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->alpha_AFM, alpha_AFM, pMesh->grel_AFM, grel_AFM);

	int tn = omp_get_thread_num();

	//sub-lattice B value so we can read it after
	Equation_Eval_2[tn] = (-GAMMA * grel_AFM.j / (1 + alpha_AFM.j*alpha_AFM.j)) * ((pMesh->M2[idx] ^ pMesh->Heff2[idx]) + alpha_AFM.j * ((pMesh->M2[idx] / Ms_AFM.j) ^ (pMesh->M2[idx] ^ pMesh->Heff2[idx])));

	//return the sub-lattice A value as normal
	return (-GAMMA * grel_AFM.i / (1 + alpha_AFM.i*alpha_AFM.i)) * ((pMesh->M[idx] ^ pMesh->Heff[idx]) + alpha_AFM.i * ((pMesh->M[idx] / Ms_AFM.i) ^ (pMesh->M[idx] ^ pMesh->Heff[idx])));
}

//Landau-Lifshitz-Gilbert equation but with no precession term and damping set to 1 : faster relaxation for static problems
DBL3 DifferentialEquationAFM::LLGStatic(int idx)
{
	DBL2 Ms_AFM = pMesh->Ms_AFM;
	pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);

	int tn = omp_get_thread_num();

	//sub-lattice B value so we can read it after
	Equation_Eval_2[tn] = (-GAMMA / 2) * ((pMesh->M2[idx] / Ms_AFM.j) ^ (pMesh->M2[idx] ^ pMesh->Heff2[idx]));

	//return the sub-lattice A value as normal
	return (-GAMMA / 2) * ((pMesh->M[idx] / Ms_AFM.i) ^ (pMesh->M[idx] ^ pMesh->Heff[idx]));
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::LLGSTT(int idx)
{
	return DBL3();
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::LLB(int idx)
{
	return DBL3();
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationAFM::LLBSTT(int idx)
{
	return DBL3();
}