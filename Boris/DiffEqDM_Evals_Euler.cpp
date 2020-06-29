#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_COMPILATION_EULER

//--------------------------------------------- EULER

void DifferentialEquationDM::RunEuler_withReductions(void)
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

void DifferentialEquationDM::RunEuler(void)
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

#endif
#endif