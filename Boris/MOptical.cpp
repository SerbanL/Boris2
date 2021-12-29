#include "stdafx.h"
#include "MOptical.h"

#ifdef MODULE_COMPILATION_MOPTICAL

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "MOpticalCUDA.h"
#endif

MOptical::MOptical(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

MOptical::~MOptical()
{
}

BError MOptical::Initialize(void)
{
	BError error(CLASS_STR(MOptical));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_MOPTICAL || pMesh->IsOutputDataSet_withRect(DATA_E_MOPTICAL),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_MOPTICAL || pMesh->IsOutputDataSet_withRect(DATA_E_MOPTICAL),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError MOptical::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MOptical));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void MOptical::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError MOptical::MakeCUDAModule(void)
{
	BError error(CLASS_STR(MOptical));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new MOpticalCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double MOptical::UpdateField(void)
{
	double energy = 0;

	if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double cHmo = pMesh->cHmo;
			pMesh->update_parameters_mcoarse(idx, pMesh->cHmo, cHmo);

			//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
			pMesh->Heff[idx] += DBL3(0, 0, cHmo);
			pMesh->Heff2[idx] += DBL3(0, 0, cHmo);

			energy += (pMesh->M[idx] + pMesh->M2[idx]) * DBL3(0, 0, cHmo) / 2;

			if (Module_Heff.linear_size()) Module_Heff[idx] = DBL3(0, 0, cHmo);
			if (Module_Heff2.linear_size()) Module_Heff2[idx] = DBL3(0, 0, cHmo);
			if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * DBL3(0, 0, cHmo);
			if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * pMesh->M2[idx] * DBL3(0, 0, cHmo);
		}
	}

	else {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double cHmo = pMesh->cHmo;
			pMesh->update_parameters_mcoarse(idx, pMesh->cHmo, cHmo);

			//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
			pMesh->Heff[idx] += DBL3(0, 0, cHmo);

			energy += pMesh->M[idx] * DBL3(0, 0, cHmo);

			if (Module_Heff.linear_size()) Module_Heff[idx] = DBL3(0, 0, cHmo);
			if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * DBL3(0, 0, cHmo);
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / pMesh->M.get_nonempty_cells();
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//-------------------Energy methods

//FM Mesh
double MOptical::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double cHmo = pMesh->cHmo;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->cHmo, cHmo);

		if (Mnew != DBL3()) return -pMesh->h.dim() * (Mnew - pMesh->M[spin_index]) * MU0 * DBL3(0, 0, cHmo);
		else return -pMesh->h.dim() * pMesh->M[spin_index] * MU0 * DBL3(0, 0, cHmo);
	}
	else return 0.0;
}

//AFM mesh
DBL2 MOptical::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		double cHmo = pMesh->cHmo;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->cHmo, cHmo);

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) return -MU0 * pMesh->h.dim() * DBL2((Mnew_A - pMesh->M[spin_index]) * DBL3(0, 0, cHmo), (Mnew_B - pMesh->M2[spin_index]) * DBL3(0, 0, cHmo));
		else return -MU0 * pMesh->h.dim() * DBL2(pMesh->M[spin_index] * DBL3(0, 0, cHmo), pMesh->M2[spin_index] * DBL3(0, 0, cHmo));
	}
	else return DBL2();
}

#endif