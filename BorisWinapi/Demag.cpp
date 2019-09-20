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

Demag::~Demag()
{
	//when deleting the Demag module any pbc settings should no longer take effect in this mesh
	//thus must clear pbc flags in M

	pMesh->M.set_pbc(false, false, false);

	//same for the CUDA version if we are in cuda mode
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->M()->copyflags_from_cpuvec(pMesh->M);
	}
#endif
}

BError Demag::Initialize(void) 
{	
	BError error(CLASS_STR(Demag));

	if(!initialized) {
		
		error = Calculate_Demag_Kernels();

		if (!error) initialized = true;
	}

	//make sure to allocate memory for Hdemag if we need it
	if (pMesh->meshODE.EvaluationSpeedup()) Hdemag.resize(pMesh->h, pMesh->meshRect);
	else Hdemag.clear();

	Hdemag_calculated = false;
	
	return error;
}

BError Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(pMesh->n, pMesh->h, demag_pbc_images)) {
		
		Uninitialize();

		//Set convolution dimensions for embedded multiplication and required PBC conditions
		error = SetDimensions(pMesh->n, pMesh->h, true, demag_pbc_images);
	}

	//if memory needs to be allocated for Hdemag, it will be done through Initialize 
	Hdemag.clear();
	Hdemag_calculated = false;

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
	}
#endif

	return error;
}

//Set PBC
BError Demag::Set_PBC(INT3 demag_pbc_images_)
{
	BError error(__FUNCTION__);

	demag_pbc_images = demag_pbc_images_;

	//set pbc conditions in M : if any are zero then pbc is disabled in that dimension
	pMesh->M.set_pbc(demag_pbc_images.x, demag_pbc_images.y, demag_pbc_images.z);

	//same for the CUDA version if we are in cuda mode
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!pMesh->pMeshCUDA->M()->copyflags_from_cpuvec(pMesh->M)) return error(BERROR_GPUERROR_CRIT);
	}
#endif

	//update will be needed if pbc settings have changed
	error = UpdateConfiguration();

	return error;
}

BError Demag::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Demag));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new DemagCUDA(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA), this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Demag::UpdateField(void) 
{
	if (pMesh->meshODE.EvaluationSpeedup()) {

		//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

		int update_type = pMesh->meshODE.Check_Step_Update();

		if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

			//calculate field required

			if (update_type == EVALSPEEDUPSTEP_COMPUTE_AND_SAVE) {

				//calculate field and save it for next time : we'll need to use it (expecting update_type = EVALSPEEDUPSTEP_SKIP next time)

				//convolute and get "energy" value
				energy = Convolute(pMesh->M, Hdemag, true);

				Hdemag_calculated = true;

				//finish off energy value
				if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
				else energy = 0;
			}
			else {

				//calculate field but do not save it for next time : we'll need to recalculate it again (expecting update_type != EVALSPEEDUPSTEP_SKIP again next time : EVALSPEEDUPSTEP_COMPUTE_NO_SAVE or EVALSPEEDUPSTEP_COMPUTE_AND_SAVE)

				//convolute and get "energy" value
				energy = Convolute(pMesh->M, pMesh->Heff, false);

				//good practice to set this to false
				Hdemag_calculated = false;

				//finish off energy value
				if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
				else energy = 0;

				//return here to avoid adding Hdemag to Heff : we've already added demag field contribution
				return energy;
			}
		}

		//add contribution to Heff
		#pragma omp parallel for
		for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

			pMesh->Heff[idx] += Hdemag[idx];
		}
	}
	else {

		//don't use evaluation speedup, so no need to use Hdemag (this won't have memory allocated anyway)

		//convolute and get "energy" value
		energy = Convolute(pMesh->M, pMesh->Heff, false);

		//finish off energy value
		if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
		else energy = 0;
	}

	return energy;
}

#endif


