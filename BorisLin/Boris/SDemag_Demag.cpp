#include "stdafx.h"
#include "SDemag_Demag.h"

#ifdef MODULE_COMPILATION_SDEMAG

#include "SuperMesh.h"
#include "SDemag.h"
#include "MeshDefs.h"
#include "SimScheduleDefs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SDemag_Demag::SDemag_Demag(MeshBase *pMeshBase_) :
	Modules(),
	Convolution<SDemag_Demag, DemagKernelCollection>(),
	ProgramStateNames(this, {VINFO(n), VINFO(h), VINFO(meshRect), VINFO(layer_number_2d)}, {})
{
	pMeshBase = pMeshBase_;

	if (!pMeshBase->is_atomistic()) pMesh = dynamic_cast<Mesh*>(pMeshBase);
	else paMesh = dynamic_cast<Atom_Mesh*>(pMeshBase);

	//this will be set separately using Set_SDemag_Pointer
	pSDemag = nullptr;

	//save ferromagnetic mesh dimensions now so we can tell later if they have changed (and uninitialize)
	n = pMeshBase->n;
	h = pMeshBase->h;
	meshRect = pMeshBase->meshRect;

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMeshBase->cudaEnabled) {

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
	for (int idx = 0; idx < (*pMeshBase)().size(); idx++) {

		if ((*pMesh)[idx] == this) {

			INT2 moduleid = (*pMeshBase)().get_id_from_index(idx);

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

	if (pMeshBase->GetMeshType() == MESH_ANTIFERROMAGNETIC) {
		
		if (!transfer.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
			{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
			MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//the Hdemag.size() checks are needed : if initializing in CUDA mode, this routine will be called so transfer objecty can calculate the mesh transfer. we don't need Hdemag ... from SDemag_Demag module in CUDA mode (but we do need them from the SDemagCUDA_Demag).
		if (pMeshBase->pSMesh->GetEvaluationSpeedup()) {

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 1 && Hdemag.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 2 && Hdemag2.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag2.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 3 && Hdemag3.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag3.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 4 && Hdemag4.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag4.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 5 && Hdemag5.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag5.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
					{ &pMesh->M }, { &pMesh->M2 }, { &pMesh->Heff }, { &pMesh->Heff2 },
					MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 6 && Hdemag6.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag6.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(
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
	/////////////////////////////////////  FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMeshBase->GetMeshType() == MESH_FERROMAGNETIC) {

		if (!transfer.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		if (pMeshBase->pSMesh->GetEvaluationSpeedup()) {

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 1 && Hdemag.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 2 && Hdemag2.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag2.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 3 && Hdemag3.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag3.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 4 && Hdemag4.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag4.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 5 && Hdemag5.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag5.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 6 && Hdemag6.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag6.Initialize_MeshTransfer({ &pMesh->M }, { &pMesh->Heff }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			}
		}

		//transfer values from M - we need this to get number of non-empty cells
		transfer.transfer_in();
		non_empty_cells = transfer.get_nonempty_cells();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// ATOMISTIC MESH ///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		if (!transfer.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);

		if (pMeshBase->pSMesh->GetEvaluationSpeedup()) {

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 1 && Hdemag.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 2 && Hdemag2.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag2.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 3 && Hdemag3.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag3.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 4 && Hdemag4.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag4.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 5 && Hdemag5.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag5.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
			}

			if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 6 && Hdemag6.size() == transfer.size()) {

				//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
				if (!Hdemag6.Initialize_MeshTransfer({ &paMesh->M1 }, { &paMesh->Heff1 }, MESHTRANSFERTYPE_WDENSITY, MUB)) return error(BERROR_OUTOFMEMORY_CRIT);
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
			pMeshBase->h, pMeshBase->meshRect,
			(MOD_)pMeshBase->Get_ActualModule_Heff_Display() == MOD_SDEMAG_DEMAG || pMeshBase->IsOutputDataSet_withRect(DATA_E_DEMAG) || pMeshBase->IsStageSet(SS_MONTECARLO),
			(MOD_)pMeshBase->Get_ActualModule_Energy_Display() == MOD_SDEMAG_DEMAG || pMeshBase->IsOutputDataSet_withRect(DATA_E_DEMAG) || pMeshBase->IsStageSet(SS_MONTECARLO));
		if (error) return error(BERROR_OUTOFMEMORY_CRIT);

		//if a Monte Carlo stage is set then we need to compute fields
		if (pMeshBase->IsStageSet(SS_MONTECARLO)) pMeshBase->Set_Force_MonteCarlo_ComputeFields(true);
	}

	//set correct size for module Heff and energy transfers, if needed
	//Only setup module Heff transfer if memory already allocated for  Module_Heff
	//In case of 2d layering mode, then only the first occurence of the MOD_SDEMAG_DEMAG module has Module_Heff allocated, hence the use of MOD_SDEMAG_DEMAG below.
	//Note the indexing (MOD_SDEMAG_DEMAG) : this assumes minor id = 0, thus first occurence of the module.
	if ((*pMeshBase)(MOD_SDEMAG_DEMAG)->Get_Module_Heff().linear_size()) {

		//convolution rectangle for this module
		Rect convolution_rect = pSDemag->get_convolution_rect(this);

		//common discretisation cellsize (may differ in thickness in 2D mode)
		DBL3 h_common = convolution_rect / pSDemag->n_common;

		//only initialize the transfer objects if needed
		if (convolution_rect != meshRect || h_common != h) {

			if (!transfer_Module_Heff.resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);
			if (!transfer_Module_energy.resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);

			if (!transfer_Module_Heff.Initialize_MeshTransfer({}, { &(*pMeshBase)(MOD_SDEMAG_DEMAG)->Get_Module_Heff() }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
			if (!transfer_Module_energy.Initialize_MeshTransfer({}, { &(*pMeshBase)(MOD_SDEMAG_DEMAG)->Get_Module_Energy() }, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
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
		if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 1) Hdemag.resize(h_common, convolution_rect);
		else Hdemag.clear();

		if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 2) Hdemag2.resize(h_common, convolution_rect);
		else Hdemag2.clear();

		if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 3) Hdemag3.resize(h_common, convolution_rect);
		else Hdemag3.clear();

		if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 4) Hdemag4.resize(h_common, convolution_rect);
		else Hdemag4.clear();

		if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 5) Hdemag5.resize(h_common, convolution_rect);
		else Hdemag5.clear();

		if (pMeshBase->pSMesh->GetEvaluationSpeedup() >= 6) Hdemag6.resize(h_common, convolution_rect);
		else Hdemag6.clear();

		if (!pMeshBase->is_atomistic() && convolution_rect == meshRect && h_common == h) {

			//no transfer required (always force transfer for atomistic meshes)
			do_transfer = false;
			transfer.clear();
			
			non_empty_cells = pMeshBase->Get_NonEmpty_Magnetic_Cells();
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
	if (n != pMeshBase->n || h != pMeshBase->h || meshRect != pMeshBase->meshRect) {

		//SDemag will see this module is uninitialized, which will force uninitialization of all other SDemag_Demag modules and resetting calculated kernels
		Uninitialize();

		n = pMeshBase->n;
		h = pMeshBase->h;
		meshRect = pMeshBase->meshRect;
	}

	//...or if convolution settings have changed
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_DEMAG_CONVCHANGE)) Uninitialize();

	//if memory needs to be allocated for Hdemag, it will be done through Initialize 
	if (!initialized) {

		Hdemag.clear();
		Hdemag2.clear();
		Hdemag3.clear();
		Hdemag4.clear();
		Hdemag5.clear();
		Hdemag6.clear();
	}

	if (layer_number_2d >= 0) {

		//if we are in 2d layering mode, then must make sure we don't have extra SDemag_Demag we don't need : i.e. if layer_number_2d >= n.z, we must delete this module
		if (layer_number_2d >= n.z) {

			//call for this module to be deleted.
			pMeshBase->DelModule(this);

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

	if (pMeshBase->pMeshBaseCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new SDemagCUDA_Demag(pMeshBase->pMeshBaseCUDA, this);
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

//-------------------Energy methods

//FM mesh
double SDemag_Demag::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//Module_Heff needs to be calculated (done during a Monte Carlo simulation, where this method would be used)
	if (Module_Heff.linear_size()) {

		if (pMesh) {

			//do not divide by 2 as we are not double-counting here
			if (Mnew != DBL3()) return -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * (Mnew - pMesh->M[spin_index]);
			else return -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * pMesh->M[spin_index];
		}
		else if (paMesh) {

			if (Mnew != DBL3()) return -MUB_MU0 * Module_Heff[paMesh->M1.cellidx_to_position(spin_index)] * (Mnew - paMesh->M1[spin_index]);
			else return -MUB_MU0 * Module_Heff[paMesh->M1.cellidx_to_position(spin_index)] * paMesh->M1[spin_index];
		}
	}
	
	return 0.0;
}

//AFM mesh
DBL2 SDemag_Demag::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	//Module_Heff needs to be calculated (done during a Monte Carlo simulation, where this method would be used)
	if (pMesh && Module_Heff.linear_size() && Module_Heff2.linear_size()) {

		DBL3 M = (pMesh->M[spin_index] + pMesh->M2[spin_index]) / 2;
		DBL3 Mnew = (Mnew_A + Mnew_B) / 2;

		double energy_ = 0.0;

		//do not divide by 2 as we are not double-counting here
		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			energy_ = -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * (Mnew - M);
		}
		else {

			energy_ = -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * M;
		}

		return DBL2(energy_, energy_);
	}
	else return DBL2();
}

#endif
