#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- SIMULATION CONTROL

BError SuperMesh::InitializeAllModules(void)
{
	BError error(CLASS_STR(SuperMesh));
	
	energy_density_weights.assign(pMesh.size(), 0.0);

	double total_nonempty_volume = 0.0;

	//1. initialize individual mesh modules
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (!error) error = pMesh[idx]->InitializeAllModules();

		total_nonempty_volume += (double)pMesh[idx]->M.get_nonempty_cells() * pMesh[idx]->M.h.dim();
	}

	//2. initialize super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error) error = pSMod[idx]->Initialize();
	}

	//3. set weights for total energy density calculation
	if (total_nonempty_volume) {

		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			energy_density_weights[idx] = (double)pMesh[idx]->M.get_nonempty_cells() * pMesh[idx]->M.h.dim() / total_nonempty_volume;
		}
	}

	return error;
}

#if COMPILECUDA == 1
BError SuperMesh::InitializeAllModulesCUDA(void)
{
	BError error(CLASS_STR(SuperMesh));

	energy_density_weights.assign(pMesh.size(), 0.0);

	double total_nonempty_volume = 0.0;

	//1. initialize individual mesh modules
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (!error) error = pMesh[idx]->InitializeAllModulesCUDA();

		total_nonempty_volume += (double)pMesh[idx]->M.get_nonempty_cells() * pMesh[idx]->M.h.dim();
	}

	//2. initialize super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error) error = pSMod[idx]->InitializeCUDA();
	}

	//3. set weights for total energy density calculation
	if (total_nonempty_volume) {

		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			energy_density_weights[idx] = (double)pMesh[idx]->M.get_nonempty_cells() * pMesh[idx]->M.h.dim() / total_nonempty_volume;
		}
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

		total_energy_density = 0.0;

		//first update the effective fields in all the meshes (skipping any that have been calculated on the super-mesh
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			total_energy_density += (pMesh[idx]->UpdateModules() * energy_density_weights[idx]);
		}

		//update effective field for super-mesh modules
		for (int idx = 0; idx < (int)pSMod.size(); idx++) {

			//super-mesh modules contribute with equal weights as sum total of individual mesh energy densities -> i.e. we don't need to apply a weight here
			total_energy_density += pSMod[idx]->UpdateField();
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

	total_energy_density = 0.0;

	//first update the effective fields in all the meshes (skipping any that have been calculated on the super-mesh
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		total_energy_density += pMesh[idx]->UpdateModules();
	}

	//update effective field for super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		total_energy_density += pSMod[idx]->UpdateField();
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

//iterate transport solver only, if available
void SuperMesh::UpdateTransportSolver(void)
{
	//first update MOD_TRANSPORT module in all the meshes
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		pMesh[idx]->UpdateTransportSolver();
	}

	//update MODS_STRANSPORT super-mesh module
	if (IsSuperMeshModuleSet(MODS_STRANSPORT)) pSMod(MODS_STRANSPORT)->UpdateField();
}

#if COMPILECUDA == 1
void SuperMesh::UpdateTransportSolverCUDA(void)
{
	//first update MOD_TRANSPORT module in all the meshes
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		pMesh[idx]->UpdateTransportSolverCUDA();
	}

	//update MODS_STRANSPORT super-mesh module
	if (IsSuperMeshModuleSet(MODS_STRANSPORT)) pSMod(MODS_STRANSPORT)->UpdateFieldCUDA();
}
#endif