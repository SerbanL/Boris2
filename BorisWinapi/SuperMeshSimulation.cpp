#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- SIMULATION CONTROL

BError SuperMesh::InitializeAllModules(void)
{
	BError error(CLASS_STR(SuperMesh));

	//1. initialize individual mesh modules
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (!error) error = pMesh[idx]->InitializeAllModules();
	}

	//2. initialize super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error) error = pSMod[idx]->Initialize();
	}

	return error;
}

#if COMPILECUDA == 1
BError SuperMesh::InitializeAllModulesCUDA(void)
{
	BError error(CLASS_STR(SuperMesh));

	//1. initialize individual mesh modules
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (!error) error = pMesh[idx]->InitializeAllModulesCUDA();
	}

	//2. initialize super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error) error = pSMod[idx]->InitializeCUDA();
	}

	return error;
}
#endif

void SuperMesh::AdvanceTime(void)
{
	//moving mesh algorithm, if enabled
	if (odeSolver.IsMovingMeshSet()) odeSolver.MovingMeshAlgorithm(this);

	do {

		//prepare meshes for new iteration (typically involves setting some state flag)
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			pMesh[idx]->PrepareNewIteration();
		}

		//first update the effective fields in all the meshes (skipping any that have been calculated on the super-mesh
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			pMesh[idx]->UpdateModules();
		}

		//update effective field for super-mesh modules
		for (int idx = 0; idx < (int)pSMod.size(); idx++) {

			pSMod[idx]->UpdateField();
		}

		//advance time using the ODE evaluation method - ODE solvers are called separately in the ferromagnetic meshes. This is why the same evaluation method must be used in all the ferromagnetic meshes, with the same time step.
		odeSolver.Iterate();

	} while (!odeSolver.TimeStepSolved());
}

#if COMPILECUDA == 1
void SuperMesh::AdvanceTimeCUDA(void)
{
	//moving mesh algorithm, if enabled
	if (odeSolver.IsMovingMeshSet()) odeSolver.MovingMeshAlgorithm(this);

	do {

		//prepare meshes for new iteration (typically involves setting some state flag)
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			pMesh[idx]->PrepareNewIterationCUDA();
		}

		//first update the effective fields in all the meshes (skipping any that have been calculated on the super-mesh
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			pMesh[idx]->UpdateModulesCUDA();
		}

		//update effective field for super-mesh modules
		for (int idx = 0; idx < (int)pSMod.size(); idx++) {

			pSMod[idx]->UpdateFieldCUDA();
		}

		//advance time using the ODE evaluation method - ODE solvers are called separately in the ferromagnetic meshes. This is why the same evaluation method must be used in all the ferromagnetic meshes, with the same time step.
		odeSolver.IterateCUDA();

	} while (!odeSolver.TimeStepSolved());
}
#endif

void SuperMesh::ComputeFields(void)
{
	//prepare meshes for new iteration (typically involves setting some state flag)
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		pMesh[idx]->PrepareNewIteration();
	}

	//first update the effective fields in all the meshes (skipping any that have been calculated on the super-mesh
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		pMesh[idx]->UpdateModules();
	}

	//update effective field for super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		pSMod[idx]->UpdateField();
	}
}

#if COMPILECUDA == 1
void SuperMesh::ComputeFieldsCUDA(void)
{
	//prepare meshes for new iteration (typically involves setting some state flag)
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		pMesh[idx]->PrepareNewIterationCUDA();
	}

	//first update the effective fields in all the meshes (skipping any that have been calculated on the super-mesh
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		pMesh[idx]->UpdateModulesCUDA();
	}

	//update effective field for super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		pSMod[idx]->UpdateFieldCUDA();
	}
}
#endif