#include "stdafx.h"
#include "Atom_Demag.h"
#include "SuperMesh.h"

#if defined(MODULE_COMPILATION_DEMAG) && ATOMISTIC == 1

#include "Atom_Mesh.h"

#if COMPILECUDA == 1
#include "Atom_DemagCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Atom_Demag::Atom_Demag(Atom_Mesh *paMesh_) :
	Modules(),
	Convolution<Atom_Demag, DemagKernel>(paMesh_->n_dm, paMesh_->h_dm),
	ProgramStateNames(this, { VINFO(demag_pbc_images) }, {})
{
	paMesh = paMesh_;

	Uninitialize();

	error_on_create = Convolution_Error_on_Create();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_Demag::~Atom_Demag()
{
	//when deleting the Demag module any pbc settings should no longer take effect in this mesh
	//thus must clear pbc flags in M1

	paMesh->M1.set_pbc(false, false, false);

	//same for the CUDA version if we are in cuda mode
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		paMesh->paMeshCUDA->M1()->copyflags_from_cpuvec(paMesh->M1);
	}
#endif
}

//Initialize mesh transfer from atomistic mesh to micromagnetic mesh for demag field computation
BError Atom_Demag::Initialize_Mesh_Transfer(void)
{
	BError error(CLASS_STR(Atom_Demag));

	if (!M.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
	if (!Hdemag.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);

	if (!M.Initialize_MeshTransfer({ &paMesh->M1 }, {}, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
	if (!Hdemag.Initialize_MeshTransfer({}, { &paMesh->Heff1 }, MESHTRANSFERTYPE_ENLARGED)) return error(BERROR_OUTOFMEMORY_CRIT);

	return error;
}

BError Atom_Demag::Initialize(void)
{
	BError error(CLASS_STR(Atom_Demag));

	if (!initialized) {

		error = Calculate_Demag_Kernels();

		error = Initialize_Mesh_Transfer();

		if (!error) initialized = true;
	}

	//need to calculate non-empty cells here so we don't waste time during computations (M is a VEC, not a VEC_VC, which means non-empty cells need to be calculated on every call)
	M.transfer_in();
	non_empty_cells = M.get_nonempty_cells();

	Hdemag_calculated = false;

	return error;
}

BError Atom_Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Demag));

	Uninitialize();

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(paMesh->n_dm, paMesh->h_dm, demag_pbc_images) || cfgMessage == UPDATECONFIG_MESHCHANGE) {

		Uninitialize();

		//Set convolution dimensions for embedded multiplication and required PBC conditions
		error = SetDimensions(paMesh->n_dm, paMesh->h_dm, true, demag_pbc_images);

		Hdemag.clear();
		M.clear();
	}

	Hdemag_calculated = false;

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

//Set PBC
BError Atom_Demag::Set_PBC(INT3 demag_pbc_images_)
{
	BError error(__FUNCTION__);

	demag_pbc_images = demag_pbc_images_;

	paMesh->Set_Magnetic_PBC(demag_pbc_images);

	//update will be needed if pbc settings have changed
	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	return error;
}

BError Atom_Demag::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Demag));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_DemagCUDA(paMesh->paMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Demag::UpdateField(void)
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	if (paMesh->pSMesh->GetEvaluationSpeedup()) {

		//use evaluation speedup method (calculate Hdemag only when required, otherwise just update effective field with previously calculated Hdemag)

		int update_type = paMesh->pSMesh->Check_Step_Update();

		if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

			//calculate field required

			//transfer magnetic moments to magnetization mesh, converting from moment to magnetization in the process
			M.transfer_in();

			//convolute and get "energy" value
			energy = Convolute(M, Hdemag, true);

			Hdemag_calculated = true;

			//finish off energy value
			if (non_empty_cells) energy *= -MU0 / (2 * non_empty_cells);
			else energy = 0;
		}

		//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hdemag.transfer_out();
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//don't use evaluation speedup

		//transfer magnetic moments to magnetization mesh, converting from moment to magnetization in the process
		M.transfer_in();

		//convolute and get "energy" value
		energy = Convolute(M, Hdemag, true);

		//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hdemag.transfer_out();

		//finish off energy value
		if (non_empty_cells) energy *= -MU0 / (2 * non_empty_cells);
		else energy = 0;
	}

	return energy;
}

#endif


