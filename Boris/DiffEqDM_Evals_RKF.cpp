#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_COMPILATION_RKF

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order solution, 5th order error)

void DifferentialEquationDM::RunRKF45_Step0_withReductions(void)
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

void DifferentialEquationDM::RunRKF45_Step0(void)
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

void DifferentialEquationDM::RunRKF45_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKF45_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKF45_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKF45_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKF45_Step5_withReductions(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKF45_Step5(void)
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