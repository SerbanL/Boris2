#include "stdafx.h"
#include "StrayField_AtomMesh.h"

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "Atom_Mesh.h"
#include "SuperMesh.h"

StrayField_AtomMesh::StrayField_AtomMesh(Atom_Mesh *paMesh_) :
	StrayField_Base(paMesh_->pSMesh),
	Modules(),
	ProgramStateNames(this, {}, {})
{
	paMesh = paMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

StrayField_AtomMesh::~StrayField_AtomMesh()
{
}

BError StrayField_AtomMesh::Initialize(void)
{
	BError error(CLASS_STR(StrayField_AtomMesh));

	if (!initialized) {

		InitializeStrayField();
		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect,
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_STRAYFIELD_MESH || paMesh->IsOutputDataSet_withRect(DATA_E_STRAY),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_STRAYFIELD_MESH || paMesh->IsOutputDataSet_withRect(DATA_E_STRAY));
	if (!error) initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError StrayField_AtomMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayField_AtomMesh));

	//only need to uninitialize if meshes have changed, added or deleted
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED)) {

		Uninitialize();
		dipoles.clear();

		strayField.assign(paMesh->h, paMesh->meshRect, DBL3());
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
#endif

	return error;
}

BError StrayField_AtomMesh::MakeCUDAModule(void)
{
	BError error(CLASS_STR(StrayField_AtomMesh));

#if COMPILECUDA == 1

	pModuleCUDA = new StrayField_AtomMeshCUDA(paMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

double StrayField_AtomMesh::UpdateField(void)
{
	//recalculate stray field if needed (required when a dipole mesh changes, as indicated by its status flag)
	if (CheckRecalculateStrayField()) CalculateStrayField();

	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			DBL3 Hstray = strayField[idx];

			paMesh->Heff1[idx] += Hstray;
			energy += -MUB_MU0 * paMesh->M1[idx] * Hstray;

			if (Module_Heff.linear_size()) Module_Heff[idx] = Hstray;
			if (Module_energy.linear_size()) Module_energy[idx] = -MUB_MU0 * (paMesh->M1[idx] * Hstray) / paMesh->M1.h.dim();
		}
	}

	//convert to energy density and return
	if (non_empty_volume) this->energy = energy / non_empty_volume;
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double StrayField_AtomMesh::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		if (Mnew != DBL3()) return -MUB_MU0 * (Mnew - paMesh->M1[spin_index]) * strayField[spin_index];
		else return -MUB_MU0 * paMesh->M1[spin_index] * strayField[spin_index];
	}
	else return 0.0;
}

#endif