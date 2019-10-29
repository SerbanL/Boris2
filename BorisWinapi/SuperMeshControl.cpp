#include "stdafx.h"
#include "SuperMesh.h"

//---------------------------------------------------------IMPORTANT CONTROL METHODS

BError SuperMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SuperMesh));
	
	//1. Ferromagnetic super-mesh - also construct the entire super-mesh rectangle

	//calculate ferromagnetic super-mesh from currently set meshes and super-mesh cellsize
	sMeshRect_fm = Rect();

	sMeshRect = Rect();

	//identify all existing ferrommagnetic meshes
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (pMesh[idx]->MComputation_Enabled()) {

			//build super-mesh rectangle to include it
			sMeshRect_fm += pMesh[idx]->GetMeshRect();
		}

		sMeshRect += pMesh[idx]->GetMeshRect();
	}

	//extract numbers of cells, and adjust cell-size so it divides the super-mesh rectangle exactly
	n_fm = round(sMeshRect_fm / h_fm);
	if (n_fm.x == 0) n_fm.x = 1;
	if (n_fm.y == 0) n_fm.y = 1;
	if (n_fm.z == 0) n_fm.z = 1;
	if (!sMeshRect_fm.IsNull()) h_fm = sMeshRect_fm / n_fm;

	//2. Electric super-mesh

	//calculate electric super-mesh from currently set meshes and super-mesh cellsize
	sMeshRect_e = Rect();

	//identify all existing ferrommagnetic meshes
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (pMesh[idx]->EComputation_Enabled()) {

			//build super-mesh rectangle to include it
			sMeshRect_e += pMesh[idx]->GetMeshRect();
		}
	}

	//extract numbers of cells, and adjust cell-size so it divides the super-mesh rectangle exactly
	n_e = round(sMeshRect_e / h_e);
	if (n_e.x <= 1) n_e.x = 2;
	if (n_e.y <= 1) n_e.y = 2;
	if (n_e.z <= 1) n_e.z = 2;
	if (!sMeshRect_e.IsNull()) h_e = sMeshRect_e / n_e;
	
	///////////////////////////////////////////////////////
	//Update configuration for meshes and their modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (!error) error = pMesh[idx]->UpdateConfiguration(cfgMessage);
	}
	
	///////////////////////////////////////////////////////
	//Update configuration for super-mesh modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error) error = pSMod[idx]->UpdateConfiguration(cfgMessage);
	}

	///////////////////////////////////////////////////////
	//Update configuration for ODECommon
	///////////////////////////////////////////////////////

	if (!error) error = odeSolver.UpdateConfiguration(cfgMessage);
	
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_SWITCHCUDASTATE)) {

		//Check if we need to couple ferromagnetic meshes to dipoles
		if (!error) CoupleToDipoles();
	}
	
#if COMPILECUDA == 1
	gpuMemFree_MB = cudaMemGetFree() / (1024 * 1024);
	gpuMemTotal_MB = cudaMemGetTotal() / (1024 * 1024);
#endif

	cpuMemFree_MB = MemGetFree() / (1024 * 1024);
	cpuMemTotal_MB = MemGetTotal() / (1024 * 1024);
	
	return error;
}

//couple ferromagnetic meshes to any touching dipole meshes, setting interface cell values and flags
void SuperMesh::CoupleToDipoles(void)
{
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

			reinterpret_cast<FMesh*>(pMesh[idx])->CoupleToDipoles(coupled_dipoles);
		}
	}
}

//switch CUDA state on/off
BError SuperMesh::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(SuperMesh));

#if COMPILECUDA == 1

	//Switch for SuperMesh

	if (cudaState) {

		if (pSMeshCUDA) delete pSMeshCUDA;

		pSMeshCUDA = new SuperMeshCUDA(this);
	}
	else {

		if (pSMeshCUDA) delete pSMeshCUDA;
		pSMeshCUDA = nullptr;
	}

	//Switch for ODECommon
	if (!error) error = odeSolver.SwitchCUDAState(cudaState);

	//Switch for Meshes and their modules
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (!error) error = pMesh[idx]->SwitchCUDAState(cudaState);
	}

	//Switch for super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error) error = pSMod[idx]->SwitchCUDAState(cudaState);
	}

	gpuMemFree_MB = cudaMemGetFree() / (1024 * 1024);
	gpuMemTotal_MB = cudaMemGetTotal() / (1024 * 1024);

	cpuMemFree_MB = MemGetFree() / (1024 * 1024);
	cpuMemTotal_MB = MemGetTotal() / (1024 * 1024);

	//make sure configuration is updated for the new mode
	error = UpdateConfiguration(UPDATECONFIG_SWITCHCUDASTATE);

#endif

	return error;
}