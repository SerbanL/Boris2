#include "stdafx.h"
#include "SDemag_Demag.h"

#ifdef MODULE_SDEMAG

#include "Mesh.h"
#include "SuperMesh.h"
#include "SDemag.h"
#include "MeshDefs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SDemag_Demag::SDemag_Demag(Mesh *pMesh_) :
	Modules(),
	Convolution<DemagKernelCollection>(),
	ProgramStateNames(this, {VINFO(n), VINFO(h), VINFO(meshRect), VINFO(layer_number_2d)}, {})
{
	pMesh = pMesh_;

	pSDemag = nullptr;

	//save ferromagnetic mesh dimensions now so we can tell later if they have changed (and uninitialize)
	n = pMesh->n;
	h = pMesh->h;
	meshRect = pMesh->meshRect;

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

void SDemag_Demag::Set_SDemag_Pointer(SDemag *pSDemag_)
{
	pSDemag = pSDemag_;

	UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

//Force the SDemag_Demag modules to calculate the layer_number_2d value, used for 2D layered convolution
	//this value is obtained from the minor id of this module held in *pMesh : the minor ids are guaranteed to be sequential and start from zero
	//thus if we've added enough of these SDemag_Demag modules, all layers will be covered
void SDemag_Demag::Set_2D_Layering(void)
{
	for (int idx = 0; idx < (*pMesh)().size(); idx++) {

		if ((*pMesh)[idx] == this) {

			INT2 moduleid = (*pMesh)().get_id_from_index(idx);

			layer_number_2d = moduleid.minor;
		}
	}
}

//initialize transfer object
BError SDemag_Demag::Initialize_Mesh_Transfer(void)
{
	BError error(__FUNCTION__);

	//convolution rectangle for this module
	Rect convolution_rect = pSDemag->get_convolution_rect(this);

	//common discretisation cellsize (may differ in thickness in 2D mode)
	DBL3 h_common = convolution_rect / pSDemag->n_common;

	//set correct size for transfer VEC
	if (!transfer.resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);

	//initialize transfer object

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

		if (!transfer.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//the Hdemag.size() check is needed : if in CUDA mode, this method will be called to initialize mesh transfer in transfer so the transfer info can be copied over to the gpu
		//In this case it's possible Hdemag does not have the correct memory allocated; if not in CUDA mode we first pass through Initialization method before calling this, in which case Hdemag will be sized correctly.
		if (pMesh->pSMesh->EvaluationSpeedup() && Hdemag.size() == transfer.size()) {

			//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
			if (!Hdemag.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
		}

		//transfer values from M - we need this to get number of non-empty cells
		transfer.transfer_in();

		non_empty_cells = transfer.get_nonempty_cells();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		if (!transfer.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
			{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
			MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//the Hdemag.size() check is needed : if in CUDA mode, this method will be called to initialize mesh transfer in transfer so the transfer info can be copied over to the gpu
		//In this case it's possible Hdemag does not have the correct memory allocated; if not in CUDA mode we first pass through Initialization method before calling this, in which case Hdemag will be sized correctly.
		if (pMesh->pSMesh->EvaluationSpeedup() && Hdemag.size() == transfer.size()) {

			//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
			if (!Hdemag.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
				{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
				MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
		}

		//transfer values from M - we need this to get number of non-empty cells
		//NOTE : do not use transfer_in_averaged here as for antiferromagnetic meshes in the ground state this will result in zero values everywhere, looking like there are no non-empty cells
		transfer.transfer_in();

		non_empty_cells = transfer.get_nonempty_cells();
	}

	//avoid division by zero
	if (!non_empty_cells) non_empty_cells = 1;

	return error;
}

BError SDemag_Demag::Initialize(void)
{
	BError error(CLASS_STR(SDemag_Demag));

	if (!pSDemag->IsInitialized()) error = pSDemag->Initialize();
	if (error) return error;
	
	//make sure to allocate memory for Hdemag if we need it
	if (pMesh->pSMesh->EvaluationSpeedup()) Hdemag.resize(pSDemag->get_convolution_rect(this) / pSDemag->n_common, pSDemag->get_convolution_rect(this));
	else Hdemag.clear();

	pSDemag->Hdemag_calculated = false;

	if (!initialized) {

		//calculate kernels for super-mesh convolution.
		if (!error) error = Calculate_Demag_Kernels(pSDemag->kernel_collection);
		if (error) return error(BERROR_OUTOFMEMORY_CRIT);

		//setup mesh transfer?

		//convolution rectangle for this module
		Rect convolution_rect = pSDemag->get_convolution_rect(this);

		//common discretisation cellsize (may differ in thickness in 2D mode)
		DBL3 h_common = convolution_rect / pSDemag->n_common;

		if (convolution_rect == meshRect && h_common == h) {

			//no transfer required
			do_transfer = false;
			transfer.clear();

			non_empty_cells = pMesh->M.get_nonempty_cells();
		}
		else {

			do_transfer = true;

			Initialize_Mesh_Transfer();
		}

		if (pSDemag->total_nonempty_volume) {

			energy_density_weight = non_empty_cells * h_common.dim() / pSDemag->total_nonempty_volume;
		}

		//avoid division by zero
		if (!non_empty_cells) non_empty_cells = 1;

		initialized = true;
	}

	return error;
}

BError SDemag_Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemag_Demag));

	//only need to uninitialize here if n, h or rectangle have changed
	if (n != pMesh->n || h != pMesh->h || meshRect != pMesh->meshRect) {

		//SDemag will see this module is uninitialized, which will force uninitialization of all other SDemag_Demag modules and resetting calculated kernels
		Uninitialize();

		n = pMesh->n;
		h = pMesh->h;
		meshRect = pMesh->meshRect;
	}

	//if memory needs to be allocated for Hdemag, it will be done through Initialize 
	Hdemag.clear();
	pSDemag->Hdemag_calculated = false;

	if (layer_number_2d >= 0) {

		//if we are in 2d layering mode, then must make sure we don't have extra SDemag_Demag we don't need : i.e. if layer_number_2d >= n.z, we must delete this module
		if (layer_number_2d >= n.z) {

			//call for this module to be deleted.
			pMesh->DelModule(this);

			//this module is no longer held in memory - must return immediately as there's nothing else we are allowed to do here.
			return BError();
		}
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError SDemag_Demag::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SDemag_Demag));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new SDemagCUDA_Demag(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

	Uninitialize();

#endif

	return error;
}

#endif
