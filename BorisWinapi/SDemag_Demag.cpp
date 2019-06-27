#include "stdafx.h"
#include "SDemag_Demag.h"

#ifdef MODULE_SDEMAG

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "SDemag.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SDemag_Demag::SDemag_Demag(Mesh *pMesh_) :
	Modules(),
	Convolution<DemagKernelCollection>(),
	ProgramStateNames(this, {VINFO(n), VINFO(h), VINFO(meshRect)}, {})
{
	pMesh = dynamic_cast<FMesh*>(pMesh_);

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
	if (!transfer.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

	//transfer values from M - we need this to get number of non-empty cells
	transfer.transfer_in();

	non_empty_cells = transfer.get_nonempty_cells();

	return error;
}

BError SDemag_Demag::Initialize(void)
{
	BError error(CLASS_STR(SDemag_Demag));

	if (!pSDemag->IsInitialized()) error = pSDemag->Initialize();
	if (error) return error;

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
	if (n != pMesh->n || h != pMesh->h || meshRect != pMesh->meshRect || cfgMessage == UPDATECONFIG_FORCEUPDATE) {

		//SDemag will see this module is uninitialized, which will force uninitialization of all other SDemag_Demag modules and resetting calculated kernels
		Uninitialize();

		n = pMesh->n;
		h = pMesh->h;
		meshRect = pMesh->meshRect;
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
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
		pModuleCUDA = new SDemagCUDA_Demag(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA), this);
		error = pModuleCUDA->Error_On_Create();
	}

	Uninitialize();

#endif

	return error;
}

#endif