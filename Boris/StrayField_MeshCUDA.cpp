#include "stdafx.h"
#include "StrayField_MeshCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "SuperMesh.h"
#include "StrayField_Mesh.h"

#include "Mesh.h"
#include "MeshCUDA.h"

StrayField_MeshCUDA::StrayField_MeshCUDA(Mesh* pMesh_, StrayField_Mesh* pStrayField_) :
	StrayField_BaseCUDA(pMesh_->pSMesh),
	ModulesCUDA()
{
	pMesh = pMesh_;
	pMeshCUDA = pMesh->pMeshCUDA;

	Uninitialize();

	pStrayField = pStrayField_;

	//set from cpu version of sm_vals
	if (!strayField()->set_from_cpuvec(pStrayField->strayField)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

StrayField_MeshCUDA::~StrayField_MeshCUDA()
{
	if (Holder_Module_Available()) {

		//copy values back to cpu version
		strayField()->copy_to_cpuvec(pStrayField->strayField);
	}
}

BError StrayField_MeshCUDA::Initialize(void)
{
	BError error(CLASS_STR(StrayField_MeshCUDA));

	//no energy density contribution here
	ZeroEnergy();

	if (!initialized) {

		//initialize and calculate the stray field using the base class method
		InitializeStrayField();

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect,
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_STRAYFIELD_MESH || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_STRAY),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_STRAYFIELD_MESH || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_STRAY),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	if (initialized) set_StrayField_MeshCUDA_pointers();

	return error;
}

BError StrayField_MeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayField_MeshCUDA));

	//only need to uninitialize if meshes have changed, added or deleted
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED)) {

		Uninitialize();
		dipoles.clear();
		managed_dipoles.clear();

		if (!strayField()->resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		strayField_size = strayField()->size_cpu().dim();
	}

	return error;
}

void StrayField_MeshCUDA::UpdateField(void)
{
	if (CheckRecalculateStrayField()) CalculateStrayField();

	UpdateFieldCUDA();
}

#endif

#endif