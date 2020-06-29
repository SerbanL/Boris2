#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_COMPILATION_ABM

//--------------------------------------------- ADAMS-BASHFORTH-MOULTON

void DifferentialEquationDM::RunABM_Predictor_withReductions(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			sM1[idx] = pMesh->M[idx];

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunABM_Predictor(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			sM1[idx] = pMesh->M[idx];

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunABM_Corrector_withReductions(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunABM_Corrector(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunABM_TEuler0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization in case we need to restore it due to evaluations in other meshes
			sM1[idx] = pMesh->M[idx];

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunABM_TEuler1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

#endif
#endif