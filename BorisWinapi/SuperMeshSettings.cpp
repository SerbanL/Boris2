#include "stdafx.h"
#include "SuperMesh.h"

//--------Getters

//get total volume energy density
double SuperMesh::GetTotalEnergyDensity(void)
{
#if COMPILECUDA == 1

	//If CUDA enabled then obtain total energy density from all modules -> this involves potentially many individual cuBReal transfers from gpu to cpu so not ideal
	//This is fine since we don't need to get the total energy often (but will need optimizing if we need to use the total energy density every iteration in some solvers in the future).
	//If we've just switched CUDA on but didn't initialize yet we won't have the correct energy_density_weights (or even memory size!), so skip this for now.
	if (pSMeshCUDA && energy_density_weights.size() == pMesh.size()) {

		total_energy_density = 0.0;

		//collect energy density values from modules in individual meshes
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			total_energy_density += (pMesh[idx]->GetEnergyDensity(MOD_ALL) * energy_density_weights[idx]);
		}

		//collect energy density values from supermesh modules
		for (int idx = 0; idx < (int)pSMod.size(); idx++) {

			total_energy_density += pSMod[idx]->GetEnergyDensity();
		}
	}
#endif

	return total_energy_density;
}

//--------Getters for supermesh modules specific properties

//status for multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
int SuperMesh::Get_Multilayered_Convolution_Status(void)
{
	if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

		return CallModuleMethod(&SDemag::Get_Multilayered_Convolution_Status);
	}
	else return -1;
}

//status for force 2d multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
int SuperMesh::Get_2D_Multilayered_Convolution_Status(void)
{
	if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

		return CallModuleMethod(&SDemag::Get_2D_Multilayered_Convolution_Status);
	}
	else return -1;
}

//status for use default n_common for multi-layered convolution (-1 : N/A, 0 : off, 1 : on)
int SuperMesh::Use_Default_n_Status(void)
{
	if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

		return CallModuleMethod(&SDemag::Use_Default_n_Status);
	}
	else return -1;
}

//get n_common for multilayerd convolution (if return = SZ3() : N/A)
SZ3 SuperMesh::Get_n_common(void)
{
	if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

		return CallModuleMethod(&SDemag::Get_n_common);
	}
	else return SZ3();
}

//--------Setters

BError SuperMesh::SetFMSMeshCellsize(DBL3 h_fm_)
{
	BError error(__FUNCTION__);

	h_fm = h_fm_;

	error = UpdateConfiguration(UPDATECONFIG_SMESH_CELLSIZE);

	return error;
}

BError SuperMesh::SetESMeshCellsize(DBL3 h_e_)
{
	BError error(__FUNCTION__);

	h_e = h_e_;

	error = UpdateConfiguration(UPDATECONFIG_SMESH_CELLSIZE);

	return error;
}