#include "stdafx.h"
#include "StrayField_AtomMeshCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "SuperMesh.h"
#include "StrayField_AtomMesh.h"

#include "Atom_Mesh.h"
#include "Atom_MeshCUDA.h"

#include "Mesh_Dipole.h"
#include "ManagedMeshCUDA.h"

StrayField_AtomMeshCUDA::StrayField_AtomMeshCUDA(Atom_Mesh* paMesh_, StrayField_AtomMesh* pStrayField_) :
	StrayField_BaseCUDA(paMesh_->pSMesh),
	ModulesCUDA()
{
	paMesh = paMesh_;
	paMeshCUDA = paMesh->paMeshCUDA;

	Uninitialize();

	pStrayField = pStrayField_;

	//set from cpu version of sm_vals
	if (!strayField()->set_from_cpuvec(pStrayField->strayField)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

StrayField_AtomMeshCUDA::~StrayField_AtomMeshCUDA()
{
	if (Holder_Module_Available()) {

		//copy values back to cpu version
		strayField()->copy_to_cpuvec(pStrayField->strayField);
	}
}

BError StrayField_AtomMeshCUDA::Initialize(void)
{
	BError error(CLASS_STR(StrayField_AtomMeshCUDA));

	//no energy density contribution here
	ZeroEnergy();

	if (!initialized) {

		//initialize and calculate the stray field using the base class method
		InitializeStrayField();

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect,
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_STRAYFIELD_MESH || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_STRAY),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_STRAYFIELD_MESH || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_STRAY));
	if (!error) initialized = true;

	if (initialized) set_StrayField_AtomMeshCUDA_pointers();

	return error;
}

BError StrayField_AtomMeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayField_AtomMeshCUDA));

	//only need to uninitialize if meshes have changed, added or deleted
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED)) {

		Uninitialize();
		dipoles.clear();
		managed_dipoles.clear();

		if (!strayField()->resize(paMesh->h, paMesh->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		strayField_size = strayField()->size_cpu().dim();
	}

	return error;
}

void StrayField_AtomMeshCUDA::UpdateField(void)
{
	if (CheckRecalculateStrayField()) CalculateStrayField();

	UpdateFieldCUDA();
}

#endif

#endif