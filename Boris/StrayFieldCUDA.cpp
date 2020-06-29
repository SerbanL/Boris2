#include "stdafx.h"
#include "StrayFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "SuperMesh.h"
#include "Mesh_Dipole.h"
#include "StrayField.h"

StrayFieldCUDA::StrayFieldCUDA(SuperMesh* pSMesh_, StrayField* pStrayField_) :
	ModulesCUDA()
{
	Uninitialize();

	pSMesh = pSMesh_;
	pStrayField = pStrayField_;

	//set from cpu version of sm_vals
	if (!strayField()->set_from_cpuvec(pStrayField->strayField)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

StrayFieldCUDA::~StrayFieldCUDA()
{
	if (Holder_Module_Available()) {

		//copy values back to cpu version
		strayField()->copy_to_cpuvec(pStrayField->strayField);
	}
}

BError StrayFieldCUDA::Initialize(void)
{
	BError error(CLASS_STR(StrayFieldCUDA));

	//no energy density contribution here
	ZeroEnergy();

	if (!initialized) {

		Mdipoles.clear();

		if (!strayField()->resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		strayField_size = strayField()->size_cpu().dim();

		//need to collect cpu versions to initialize the mesh transfer object, before transferring it to the gpu
		vector< VEC<DBL3>* > meshes_out_cpu;

		//array of pointers to oputput meshes (Heff) to transfer to
		cu_arr<cuVEC<cuReal3>> meshes_out;

		//identify all existing dipole meshes and output Heff meshes
		for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

			//collect cuda Heff pointers in pVal_to cu_arr
			if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->is_atomistic()) {

				meshes_out.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff.get_managed_object());
				
				meshes_out_cpu.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff));

				//for antiferromagnetic meshes we also need to add the stray field to the second sub-lattice
				if ((*pSMesh)[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					meshes_out.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff2.get_managed_object());

					meshes_out_cpu.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff2));
				}
			}

			//collect cuda M pointers in Mdipoles cu_arr
			if ((*pSMesh)[idx]->GetMeshType() == MESH_DIPOLE) {

				DipoleMeshCUDA *pDipoleCUDA = dynamic_cast<DipoleMeshCUDA*>(dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA);

				//save pointer to cuda M in cu_arr
				Mdipoles.push_back((cuVEC_VC<cuReal3>*&)pDipoleCUDA->M.get_managed_object());
			}
		}

		//Initialize the mesh transfer object.
		if (!pStrayField->strayField.Initialize_MeshTransfer(vector< VEC<DBL3>* >(), meshes_out_cpu, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//Now copy mesh transfer object to cuda version
		if (!strayField()->copy_transfer_info(meshes_out, meshes_out, pStrayField->strayField)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		strayField()->set(cuReal3(0.0));

		//here need to calculate stray field from currently set dipole "meshes" (if any - there should be at least one otherwise there's no point using this module)
		CalculateStrayField();

		initialized = true;
	}

	return error;
}

BError StrayFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayFieldCUDA));

	Uninitialize();

	Mdipoles.clear();

	if (!strayField()->resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	strayField_size = strayField()->size_cpu().dim();

	return error;
}

void StrayFieldCUDA::UpdateField(void)
{
	if (pStrayField->CheckRecalculateStrayField()) CalculateStrayField();

	//transfer to individual Heff meshes
	strayField()->transfer_out(pStrayField->strayField.size_transfer_out());
}

#endif

#endif