#include "stdafx.h"
#include "StrayFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "SuperMesh.h"
#include "StrayField.h"

#include "Mesh.h"
#include "Atom_Mesh.h"

StrayFieldCUDA::StrayFieldCUDA(SuperMesh* pSMesh_, StrayField* pStrayField_) :
	StrayField_BaseCUDA(pSMesh_),
	ModulesCUDA()
{
	Uninitialize();

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

		//need to collect cpu versions to initialize the mesh transfer object, before transferring it to the gpu
		std::vector< VEC<DBL3>* > meshes_out_cpu;

		//array of pointers to oputput meshes (Heff) to transfer to
		cu_arr<cuVEC<cuReal3>> meshes_out;

		//identify all output Heff meshes
		for (int idx = 0; idx < (int)pSMesh->size(); idx++) {

			//collect cuda Heff pointers in pVal_to cu_arr
			if ((*pSMesh)[idx]->MComputation_Enabled()) {

				if (!(*pSMesh)[idx]->is_atomistic()) {

					meshes_out.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff.get_managed_object());
					meshes_out_cpu.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff));
				}
				else {

					meshes_out.push_back((cuVEC<cuReal3>*&)dynamic_cast<Atom_Mesh*>((*pSMesh)[idx])->paMeshCUDA->Heff1.get_managed_object());
					meshes_out_cpu.push_back(&(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx])->Heff1));
				}

				//for antiferromagnetic meshes we also need to add the stray field to the second sub-lattice
				if ((*pSMesh)[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					meshes_out.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff2.get_managed_object());
					meshes_out_cpu.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff2));
				}
			}
		}

		//Initialize the mesh transfer object.
		if (!pStrayField->strayField.Initialize_MeshTransfer(std::vector< VEC<DBL3>* >(), meshes_out_cpu, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//Now copy mesh transfer object to cuda version
		if (!strayField()->copy_transfer_info(meshes_out, meshes_out, pStrayField->strayField)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		//initialize and calculate the stray field using the base class method
		InitializeStrayField();

		initialized = true;
	}

	return error;
}

BError StrayFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayFieldCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED)) {

		Uninitialize();
		dipoles.clear();
		managed_dipoles.clear();

		if (!strayField()->resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		strayField_size = strayField()->size_cpu().dim();
	}

	return error;
}

void StrayFieldCUDA::UpdateField(void)
{
	if (CheckRecalculateStrayField()) CalculateStrayField();

	//transfer to individual Heff meshes
	strayField()->transfer_out(pStrayField->strayField.size_transfer_out());
}

#endif

#endif