#include "stdafx.h"
#include "StrayField.h"

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "Mesh.h"
#include "Atom_Mesh.h"
#include "Mesh_Dipole.h"
#include "SuperMesh.h"

StrayField::StrayField(SuperMesh *pSMesh_) :
	StrayField_Base(pSMesh_),
	Modules(),
	ProgramStateNames(this, {}, {})
{
	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh_->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

StrayField::~StrayField()
{
}

BError StrayField::Initialize(void)
{
	BError error(CLASS_STR(StrayField));

	if (!initialized) {

		std::vector< VEC<DBL3>* > meshes_out;

		//identify all output Heff meshes
		for (int idx = 0; idx < (int)pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->MComputation_Enabled()) {

				if (!(*pSMesh)[idx]->is_atomistic()) meshes_out.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff));
				else meshes_out.push_back(&(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx])->Heff1));

				//for antiferromagnetic meshes we also need to add the stray field to the second sub-lattice
				if ((*pSMesh)[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					meshes_out.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff2));
				}
			}
		}

		//Initialize the mesh transfer object.
		if (!strayField.Initialize_MeshTransfer(std::vector< VEC<DBL3>* >(), meshes_out, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//initialize and calculate the stray field using the base class method
		InitializeStrayField();

		initialized = true;
	}

	return error;
}

BError StrayField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayField));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED)) {

		Uninitialize();
		dipoles.clear();

		if (!strayField.resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if(pModuleCUDA) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
#endif

	return error;
}

BError StrayField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(StrayField));

#if COMPILECUDA == 1

	pModuleCUDA = new StrayFieldCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

double StrayField::UpdateField(void)
{
	//recalculate stray field if needed (required when a dipole mesh changes, as indicated by its status flag)
	if (CheckRecalculateStrayField()) CalculateStrayField();

	//transfer stray field values to effective fields in the ferromagnetic meshes
	strayField.transfer_out();

	//not counting this to the total energy density for now
	return 0.0;
}

#endif