#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//---------------------------------------- OTHERS

//Restore magnetisation after a failed step for adaptive time-step methods
void DifferentialEquation::RestoreMagnetisation(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++)
		pMesh->M[idx] = sM1[idx];
}