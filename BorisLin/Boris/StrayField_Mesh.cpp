#include "stdafx.h"
#include "StrayField_Mesh.h"

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "Mesh.h"
#include "SuperMesh.h"

StrayField_Mesh::StrayField_Mesh(Mesh *pMesh_) :
	StrayField_Base(pMesh_->pSMesh),
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

StrayField_Mesh::~StrayField_Mesh()
{
}

BError StrayField_Mesh::Initialize(void)
{
	BError error(CLASS_STR(StrayField_Mesh));

	if (!initialized) {

		InitializeStrayField();
		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_STRAYFIELD_MESH || pMesh->IsOutputDataSet_withRect(DATA_E_STRAY),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_STRAYFIELD_MESH || pMesh->IsOutputDataSet_withRect(DATA_E_STRAY),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError StrayField_Mesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayField_Mesh));

	//only need to uninitialize if meshes have changed, added or deleted
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED)) {

		Uninitialize();
		dipoles.clear();

		strayField.assign(pMesh->h, pMesh->meshRect, DBL3());
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
#endif

	return error;
}

BError StrayField_Mesh::MakeCUDAModule(void)
{
	BError error(CLASS_STR(StrayField_Mesh));

#if COMPILECUDA == 1

	pModuleCUDA = new StrayField_MeshCUDA(pMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

double StrayField_Mesh::UpdateField(void)
{
	//recalculate stray field if needed (required when a dipole mesh changes, as indicated by its status flag)
	if (CheckRecalculateStrayField()) CalculateStrayField();

	double energy = 0;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL3 Hstray = strayField[idx];

				pMesh->Heff[idx] += Hstray;
				energy += pMesh->M[idx] * Hstray;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Hstray;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * (pMesh->M[idx] * Hstray);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL3 Hstray = strayField[idx];

				pMesh->Heff[idx] += Hstray;
				pMesh->Heff2[idx] += Hstray;

				energy += (pMesh->M[idx] * Hstray + pMesh->M2[idx] * Hstray) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Hstray;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Hstray;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * (pMesh->M[idx] * Hstray);
				if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * (pMesh->M2[idx] * Hstray);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// FINAL ENERGY DENSITY VALUE //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (pMesh->M.get_nonempty_cells());
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//-------------------Energy methods

//FM mesh
double StrayField_Mesh::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		if (Mnew != DBL3()) return -pMesh->h.dim() * (Mnew - pMesh->M[spin_index]) * MU0 * strayField[spin_index];
		else return -pMesh->h.dim() * pMesh->M[spin_index] * MU0 * strayField[spin_index];
	}
	else return 0.0;
}

//AFM mesh
DBL2 StrayField_Mesh::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			return -pMesh->h.dim() * DBL2((Mnew_A - pMesh->M[spin_index]) * MU0 * strayField[spin_index], (Mnew_B - pMesh->M2[spin_index]) * MU0 * strayField[spin_index]);
		}
		else return -pMesh->h.dim() * DBL2(pMesh->M[spin_index] * MU0 * strayField[spin_index], pMesh->M2[spin_index] * MU0 * strayField[spin_index]);
	}
	else return DBL2();
}

#endif