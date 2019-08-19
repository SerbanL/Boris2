#include "stdafx.h"
#include "Demag.h"

#ifdef MODULE_DEMAG

#include "Mesh_Ferromagnetic.h"

#if COMPILECUDA == 1
#include "DemagCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Demag::Demag(Mesh *pMesh_) : 
	Modules(),
	Convolution<DemagKernel>(pMesh_->GetMeshSize(), pMesh_->GetMeshCellsize()),
	ProgramStateNames(this, {VINFO(demag_pbc_images)}, {})
{
	pMesh = dynamic_cast<FMesh*>(pMesh_);

	Uninitialize();

	error_on_create = Convolution_Error_on_Create();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

BError Demag::Initialize(void) 
{	
	BError error(CLASS_STR(Demag));

	if(!initialized) {
		
		error = Calculate_Demag_Kernels();

		if (!error) initialized = true;
	}

	return error;
}

BError Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag));

	//only need to uninitialize if n or h have changed
	if (!CheckDimensions(pMesh->n, pMesh->h, demag_pbc_images)) {
		
		Uninitialize();

		//Set convolution dimensions for embedded multiplication and required PBC conditions
		error = SetDimensions(pMesh->n, pMesh->h, true, demag_pbc_images);
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
	}
#endif

	return error;
}

//Set PBC
void Demag::Set_PBC(INT3 demag_pbc_images_)
{
	demag_pbc_images = demag_pbc_images_;

	//update will be needed if pbc settings have changed
	UpdateConfiguration();
}

BError Demag::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Demag));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new DemagCUDA(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA));
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Demag::UpdateField(void) 
{
	//convolute and get "energy" value
	energy = Convolute(pMesh->M, pMesh->Heff, false);

	//finish off energy value
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
	else energy = 0;

	return energy;
}

#endif


