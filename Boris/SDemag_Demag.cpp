#include "stdafx.h"
#include "SDemag_Demag.h"

#ifdef MODULE_COMPILATION_SDEMAG

#include "Mesh.h"
#include "SuperMesh.h"
#include "SDemag.h"
#include "MeshDefs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SDemag_Demag::SDemag_Demag(Mesh *pMesh_) :
	Modules(),
	Convolution<SDemag_Demag, DemagKernelCollection>(),
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

	//set correct size for transfer VEC (if using transfer)
	if (!transfer.resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);

	//initialize transfer object and any Hdemag objects used for evaluation speedup

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {
		
		if (!transfer.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
			{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
			MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//the Hdemag.size() checks are needed : if initializing in CUDA mode, this routine will be called so transfer objecty can calculate the mesh transfer. we don't need Hdemag ... from SDemag_Demag module in CUDA mode (but we do need them from the SDemagCUDA_Demag).
		if (pMesh->pSMesh->GetEvaluationSpeedup()) {

			if (pMesh->pSMesh->GetEvaluationSpeedup() >= 1 && Hdemag.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMesh->pSMesh->GetEvaluationSpeedup() >= 2 && Hdemag2.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag2.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMesh->pSMesh->GetEvaluationSpeedup() >= 3 && Hdemag3.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag3.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}
		}

		//transfer values from M - we need this to get number of non-empty cells
		//NOTE : do not use transfer_in_averaged here as for antiferromagnetic meshes in the ground state this will result in zero values everywhere, looking like there are no non-empty cells
		transfer.transfer_in();
		non_empty_cells = transfer.get_nonempty_cells();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		if (!transfer.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		if (pMesh->pSMesh->GetEvaluationSpeedup()) {

			if (pMesh->pSMesh->GetEvaluationSpeedup() >= 1 && Hdemag.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMesh->pSMesh->GetEvaluationSpeedup() >= 2 && Hdemag2.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag2.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMesh->pSMesh->GetEvaluationSpeedup() >= 3 && Hdemag3.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag3.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}
		}

		//transfer values from M - we need this to get number of non-empty cells
		transfer.transfer_in();
		non_empty_cells = transfer.get_nonempty_cells();
	}

	//avoid division by zero
	if (!non_empty_cells) non_empty_cells = 1;

	return error;
}

//allocate memory and initialize mesh transfer for module Heff and energy display data
BError SDemag_Demag::Initialize_Module_Display(void)
{
	BError error(CLASS_STR(SDemag_Demag));

	//Make sure display data has memory allocated (or freed) as required. Do this first, in case we need to setup a transfer.
	//NOTE : if working in 2d layering mode, then only enable module heff display for the first occurence of the MOD_SDEMAG_DEMAG, hence the layer_number_2d <= 0 check.
	if (layer_number_2d <= 0) {

		error = Update_Module_Display_VECs(
			pMesh->h, pMesh->meshRect, 
			(MOD_)pMesh->Get_ActualModule_Heff_Display() == MOD_SDEMAG_DEMAG || pMesh->IsOutputDataSet_withRect(DATA_E_DEMAG),
			(MOD_)pMesh->Get_ActualModule_Energy_Display() == MOD_SDEMAG_DEMAG || pMesh->IsOutputDataSet_withRect(DATA_E_DEMAG));
		if (error) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//set correct size for module Heff and energy transfers, if needed
	//Only setup module Heff transfer if memory already allocated for  Module_Heff
	//In case of 2d layering mode, then only the first occurence of the MOD_SDEMAG_DEMAG module has Module_Heff allocated, hence the use of MOD_SDEMAG_DEMAG below.
	//Note the indexing (MOD_SDEMAG_DEMAG) : this assumes minor id = 0, thus first occurence of the module.
	if ((*pMesh)(MOD_SDEMAG_DEMAG)->Get_Module_Heff().linear_size()) {

		//convolution rectangle for this module
		Rect convolution_rect = pSDemag->get_convolution_rect(this);

		//common discretisation cellsize (may differ in thickness in 2D mode)
		DBL3 h_common = convolution_rect / pSDemag->n_common;

		//only initialize the transfer objects if needed
		if (convolution_rect != meshRect || h_common != h) {

			if (!transfer_Module_Heff.resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);
			if (!transfer_Module_energy.resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);

			if (!transfer_Module_Heff.Initialize_MeshTransfer({}, { &(*pMesh)(MOD_SDEMAG_DEMAG)->Get_Module_Heff() }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			if (!transfer_Module_energy.Initialize_MeshTransfer({}, { &(*pMesh)(MOD_SDEMAG_DEMAG)->Get_Module_Energy() }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else {

			transfer_Module_Heff.clear();
			transfer_Module_energy.clear();
		}
	}

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

		selfDemagCoeff = DemagTFunc().SelfDemag_PBC(h_common, pSDemag->n_common, pSDemag->demag_pbc_images);

		//make sure to allocate memory for Hdemag if we need it
		if (pMesh->pSMesh->GetEvaluationSpeedup() >= 1) Hdemag.resize(h_common, convolution_rect);
		else Hdemag.clear();

		if (pMesh->pSMesh->GetEvaluationSpeedup() >= 2) Hdemag2.resize(h_common, convolution_rect);
		else Hdemag2.clear();

		if (pMesh->pSMesh->GetEvaluationSpeedup() >= 3) Hdemag3.resize(h_common, convolution_rect);
		else Hdemag3.clear();

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

	//initialize module display data if needed
	error = Initialize_Module_Display();
	if (error) initialized = false;

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
	Hdemag2.clear();
	Hdemag3.clear();

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

double SDemag_Demag::UpdateField(void)
{
	//demag field update done through the supermesh module.
	//here we need to zero the module display objects in case they are used : if we have to transfer data into them from the display transfer objects, this is done by adding.
	//cannot set output when transferring since we can have multiple transfer objects contributing to the display objects
	ZeroModuleVECs();

	//same goes for the total energy
	energy = 0.0;

	return 0.0;
}


#endif
