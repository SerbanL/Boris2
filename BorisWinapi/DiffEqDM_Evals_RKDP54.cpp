#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_RKDP

//--------------------------------------------- RUNGE KUTTA DORMAND-PRINCE (4th order solution, 5th order error)

void DifferentialEquationDM::RunRKDP54_Step0_withReductions(void)
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

void DifferentialEquationDM::RunRKDP54_Step0(void)
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

void DifferentialEquationDM::RunRKDP54_Step0_Advance(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKDP54_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKDP54_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKDP54_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKDP54_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKDP54_Step5_withReductions(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Set M from diamagnetic susceptibility
			pMesh->M[idx] = CALLFP(this, equation)(idx);
		}
	}
}

void DifferentialEquationDM::RunRKDP54_Step5(void)
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