#include "stdafx.h"
#include "Oersted.h"

#ifdef MODULE_OERSTED

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

Oersted::Oersted(SuperMesh *pSMesh_) :
	Modules(),
	Convolution<OerstedKernel>(pSMesh_->n_e, pSMesh_->h_e),
	ProgramStateNames(this, {}, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration();

	//make sure to resize sm_Vals here as UpdateConfiguration will not do it now : the convolution dimensions have already been set when calling the Convolution constructor, so CheckDimensions check will return true
	if (!sm_Vals.resize(pSMesh->h_e, pSMesh->sMeshRect_e)) error_on_create(BERROR_OUTOFMEMORY_CRIT);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

BError Oersted::Initialize(void)
{
	BError error(CLASS_STR(Oersted));

	if (pSMesh->sMeshRect_e.IsNull() || pSMesh->sMeshRect_fm.IsNull()) return error(BERROR_INCORRECTMODCONFIG);

	//FFT Kernels are not so quick to calculate - if already initialized then we are guaranteed they are correct
	if (!initialized) {

		vector< VEC<DBL3>* > pVal_from;
		vector< VEC<DBL3>* > pVal_to;

		//identify all existing electric computation meshes to get Jc from and magnetic computation meshes to transfer Heff to
		for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

			if ((*pSMesh)[idx]->EComputation_Enabled()) {

				pVal_from.push_back(&((*pSMesh)[idx]->Jc));
			}

			if ((*pSMesh)[idx]->MComputation_Enabled()) {

				pVal_to.push_back(&((*pSMesh)[idx]->Heff));
			}
		}

		//Initialize the mesh transfer object.
		if (!sm_Vals.Initialize_MeshTransfer(pVal_from, pVal_to, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//transfer values from invidual M meshes to sm_Vals - we need this to get number of non-empty cells
		sm_Vals.transfer_in();

		non_empty_cells = sm_Vals.get_nonempty_cells();

		//avoid division by zero
		if (!non_empty_cells) non_empty_cells = 1;

		error = Calculate_Oersted_Kernels();

		if (!error) initialized = true;
	}

	oefield_computed = false;

	return error;
}

BError Oersted::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Oersted));

	//when the Oersted module is set then the super-mesh electric rectangle must be enlarged to also include the ferromagnetic super-mesh rectangle
	//the h_e cellsize is then also adjusted as needed.

	if (!pSMesh->sMeshRect_e.IsNull()) {

		//new super-mesh electric rect
		pSMesh->sMeshRect_e = pSMesh->sMeshRect_fm.get_union(pSMesh->sMeshRect_e);

		//adjust n_e and h_e
		pSMesh->n_e = round(pSMesh->sMeshRect_e / pSMesh->h_e);
		if (pSMesh->n_e.x <= 1) pSMesh->n_e.x = 2;
		if (pSMesh->n_e.y <= 1) pSMesh->n_e.y = 2;
		if (pSMesh->n_e.z <= 1) pSMesh->n_e.z = 2;
		pSMesh->h_e = pSMesh->sMeshRect_e / pSMesh->n_e;
	}

	//only need to uninitialize if n_e or h_e have changed
	if(!CheckDimensions(pSMesh->n_e, pSMesh->h_e, INT3())) {

		Uninitialize();
		error = SetDimensions(pSMesh->n_e, pSMesh->h_e);

		//First, initialize mesh transfer object
		if (!sm_Vals.resize(pSMesh->h_e, pSMesh->sMeshRect_e)) return error(BERROR_OUTOFMEMORY_NCRIT);
	}

	oefield_computed = false;

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
	}
#endif

	return error;
}

BError Oersted::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Oersted));

#if COMPILECUDA == 1

	pModuleCUDA = new OerstedCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

double Oersted::UpdateField(void)
{
	//only recalculate Oersted field if there was a significant change in current density (judged based on transport solver iterations)
	if (pSMesh->CallModuleMethod(&STransport::Transport_Recalculated) || !oefield_computed) {

		//transfer values from invidual Jc meshes to sm_Vals
		sm_Vals.transfer_in();

		Convolute(sm_Vals, sm_Vals, true);

		//transfer to individual Heff meshes
		sm_Vals.transfer_out(false);

		oefield_computed = true;
	}
	else {

		//transfer to individual Heff meshes
		sm_Vals.transfer_out();
	}

	//not counting this to the total energy density for now
	return 0.0;
}

#endif




