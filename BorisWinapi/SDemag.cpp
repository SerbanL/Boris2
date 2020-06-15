#include "stdafx.h"
#include "SDemag.h"

#ifdef MODULE_SDEMAG

#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SDemag::SDemag(SuperMesh *pSMesh_) :
	Modules(),
	Convolution<DemagKernel>(),
	ProgramStateNames(this, { VINFO(use_multilayered_convolution), VINFO(n_common), VINFO(use_default_n), VINFO(force_2d_convolution), VINFO(demag_pbc_images) }, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

SDemag::~SDemag()
{
	//when deleting the SDemag module any pbc settings should no longer take effect in all meshes
	//thus must clear pbc flags in all M

	demag_pbc_images = INT3();
	Set_M_PBC();

	//RAII : SDemag_Demag module were created in the constructor, so delete them here in any remaining ferromagnetic meshes 
	//(some could have been deleted already if any ferromagnetic mesh was deleted in the mean-time)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			//Delete any SDemag_Demag modules in this mesh
			(*pSMesh)[idx]->DelModule(MOD_SDEMAG_DEMAG);
		}
	}
}

//------------------ Helpers for multi-layered convolution control

//when SDemag created, it needs to add one SDemag_Demag module to each (anti)ferromagnetic mesh (or multiple if the 2D layering option is enabled).
BError SDemag::Create_SDemag_Demag_Modules(void)
{
	BError error(CLASS_STR(SDemag));

	//identify all existing ferrommagnetic meshes (magnetic computation enabled)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {

			if (force_2d_convolution < 2) {

				//not 2d layered convolution, so there can only be 1 SDemag_Demag module in each layer

				//Make sure each (anti)ferromagnetic mesh has one SDemag_Demag module added
				if (!(*pSMesh)[idx]->IsModuleSet(MOD_SDEMAG_DEMAG)) {

					error = (*pSMesh)[idx]->AddModule(MOD_SDEMAG_DEMAG);
				}

				(*pSMesh)[idx]->CallModuleMethod(&SDemag_Demag::Set_SDemag_Pointer, this);
			}
			else {

				//2d layered convolution - in each mesh need Mesh::n.z SDemag_Demag modules - one for each layer.
				//here we must ensure we have at least n.z such modules - if there are more, they will be deleted through UpdateConfiguration method in the respective SDemag_Demag modules.
				
				//number of layers required
				int num_layers = (*pSMesh)[idx]->n.z;

				while ((*pSMesh)[idx]->GetNumModules(MOD_SDEMAG_DEMAG) < num_layers && !error) {

					//keep adding SDemag_Demag modules until we have enough to cover all layers
					error = (*pSMesh)[idx]->AddModule(MOD_SDEMAG_DEMAG, true);
				}

				//set SDemag pointer in all SDemag_Demag modules
				(*pSMesh)[idx]->CallAllModulesMethod(&SDemag_Demag::Set_SDemag_Pointer, this);

				//switch on 2D layering in all SDemag_Demag modules in this mesh, making sure each one knows exactly what layer it is
				(*pSMesh)[idx]->CallAllModulesMethod(&SDemag_Demag::Set_2D_Layering);
			}
		}
	}

	return error;
}

//delete all SDemag_Demag modules - these are only created by SDemag
void SDemag::Destroy_SDemag_Demag_Modules(void)
{
	//identify all existing ferrommagnetic meshes (magnetic computation enabled)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			//Delete any SDemag_Demag modules in this mesh
			(*pSMesh)[idx]->DelModule(MOD_SDEMAG_DEMAG);

			//Also delete any Demag modules as these should not be present when we have the SDemag module enabled
			(*pSMesh)[idx]->DelModule(MOD_DEMAG);
		}
	}

	FFT_Spaces_Input.clear();
	Rect_collection.clear();
	kernel_collection.clear();
	pSDemag_Demag.clear();
}

//make sure the pSDemag_Demag list is up to date : if any mismatches found return false
bool SDemag::Update_SDemag_Demag_List(void)
{
	vector<SDemag_Demag*> pSDemag_Demag_;

	//also make sure the FFT spaces and rectangles list is correct -> rebuild it
	FFT_Spaces_Input.clear();
	kernel_collection.clear();

	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {

			for (int layer_idx = 0; layer_idx < (*pSMesh)[idx]->GetNumModules(MOD_SDEMAG_DEMAG); layer_idx++) {

				SDemag_Demag *pSDemag_Demag_Module = dynamic_cast<SDemag_Demag*>((*pSMesh)[idx]->GetModule(MOD_SDEMAG_DEMAG, layer_idx));

				pSDemag_Demag_.push_back(pSDemag_Demag_Module);

				//build the fft spaces, rectangles list, and demag kernels
				FFT_Spaces_Input.push_back(pSDemag_Demag_.back()->Get_Input_Scratch_Space());
				kernel_collection.push_back(dynamic_cast<DemagKernelCollection*>(pSDemag_Demag_.back()));
			}
		}
	}

	//before transferring pSDemag_Demag_ to pSDemag_Demag, check no changes have occured in order to determine the return value
	if (pSDemag_Demag_.size() == pSDemag_Demag.size()) {

		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

			if (pSDemag_Demag[idx] != pSDemag_Demag_[idx]) {

				pSDemag_Demag = pSDemag_Demag_;
				return false;
			}
		}

		return true;
	}
	else {

		pSDemag_Demag = pSDemag_Demag_;
		return false;
	}
}

//set default value for n_common : largest value from all (anti)ferromagnetic meshes
void SDemag::set_default_n_common(void)
{
	SZ3 n_common_old = n_common;

	n_common = SZ3(1);

	bool meshes_2d = true;

	//are all the meshes 2d?
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {

			if ((*pSMesh)[idx]->n.z != 1) {

				meshes_2d = false;
				break;
			}
		}
	}

	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {

			if (n_common.x < (*pSMesh)[idx]->n.x) n_common.x = (*pSMesh)[idx]->n.x;
			if (n_common.y < (*pSMesh)[idx]->n.y) n_common.y = (*pSMesh)[idx]->n.y;
			if (n_common.z < (*pSMesh)[idx]->n.z) n_common.z = (*pSMesh)[idx]->n.z;
		}
	}

	if (meshes_2d || force_2d_convolution) {

		//all 2D meshes, or 2D layered meshes, forced or otherwise : common n.z must be 1, thus enabling exact computation for layers with arbitrary thicknesses.
		n_common.z = 1;
	}

	//uninitialize only if n_common has changed
	if (n_common_old != n_common) Uninitialize();
}

//adjust Rect_collection so the rectangles are matching for multi-layered convolution. n_common should be calculated already.
void SDemag::set_Rect_collection(void)
{
	//first get raw Rect_collection
	Rect_collection.clear();

	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {

			for (int layer_idx = 0; layer_idx < (*pSMesh)[idx]->GetNumModules(MOD_SDEMAG_DEMAG); layer_idx++) {

				if (force_2d_convolution == 2) {

					//if in 2d layering mode, then the rectangle set in the rect collection must be that of the respective layer, not the entire mesh
					Rect_collection.push_back((*pSMesh)[idx]->GetMeshRect().get_zlayer(layer_idx, (*pSMesh)[idx]->GetMeshCellsize().z));
				}
				else {

					Rect_collection.push_back((*pSMesh)[idx]->GetMeshRect());
				}
			}
		}
	}

	//it's possible there are no modules set for multilayered convolution (e.g. all set to be excluded)
	if (!Rect_collection.size()) return;

	//now adjust it
	DBL3 max_sizes;

	double min_x = Rect_collection[0].s.x;
	double min_y = Rect_collection[0].s.y;
	double min_z = Rect_collection[0].s.z;

	//adjust rectangles so they have matching dimensions for multi-layered convolution
	
	//first find maximum size in each dimension
	//also find smallest starting coordinates along each axis
	for (int idx = 0; idx < Rect_collection.size(); idx++) {

		if (max_sizes.x < Rect_collection[idx].size().x) max_sizes.x = Rect_collection[idx].size().x;
		if (max_sizes.y < Rect_collection[idx].size().y) max_sizes.y = Rect_collection[idx].size().y;
		if (max_sizes.z < Rect_collection[idx].size().z) max_sizes.z = Rect_collection[idx].size().z;

		if (Rect_collection[idx].s.x < min_x) min_x = Rect_collection[idx].s.x;
		if (Rect_collection[idx].s.y < min_y) min_y = Rect_collection[idx].s.y;
		if (Rect_collection[idx].s.z < min_z) min_z = Rect_collection[idx].s.z;
	}

	//now enlarge rectangles so they all have sizes max_sizes (except in 2D where they keep their thickness)
	//enlarge them by setting the starting points as close as possible to the smallest starting points found above, along each axis
	//ideally they all start at the same point, thus making multi-layered convolution most efficient
	for (int idx = 0; idx < Rect_collection.size(); idx++) {
		
		if (Rect_collection[idx].e.x - max_sizes.x < min_x) {

			Rect_collection[idx].s.x = min_x;
			Rect_collection[idx].e.x = Rect_collection[idx].s.x + max_sizes.x;
		}
		else Rect_collection[idx].s.x = Rect_collection[idx].e.x - max_sizes.x;

		if (Rect_collection[idx].e.y - max_sizes.y < min_y) {

			Rect_collection[idx].s.y = min_y;
			Rect_collection[idx].e.y = Rect_collection[idx].s.y + max_sizes.y;
		}
		else Rect_collection[idx].s.y = Rect_collection[idx].e.y - max_sizes.y;

		if (n_common.z != 1) {

			//3D convolution so also set the z sizes
			if (Rect_collection[idx].e.z - max_sizes.z < min_z) {

				Rect_collection[idx].s.z = min_z;
				Rect_collection[idx].e.z = Rect_collection[idx].s.z + max_sizes.z;
			}
			else Rect_collection[idx].s.z = Rect_collection[idx].e.z - max_sizes.z;
		}
	}
}

//get maximum cellsize for multi-layered convolution (use it to normalize dimensions)
double SDemag::get_maximum_cellsize(void)
{
	double h_max = 0.0;

	for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

		if (h_max < pSDemag_Demag[idx]->h.maxdim()) h_max = pSDemag_Demag[idx]->h.maxdim();
	}

	return h_max;
}


//get convolution rectangle for the given SDemag_Demag module (remember this might not be the rectangle of M in that mesh, but an adjusted rectangle to make the convolution work)
Rect SDemag::get_convolution_rect(SDemag_Demag* demag_demag)
{
	for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

		if (pSDemag_Demag[idx] == demag_demag) return Rect_collection[idx];
	}

	//empty rect : something not right
	return Rect();
}

//-------------------Setters

//change between demag calculation types : super-mesh (status = false) or multilayered (status = true)
BError SDemag::Set_Multilayered_Convolution(bool status)
{
	BError error(CLASS_STR(SDemag));

	use_multilayered_convolution = status;
	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	if (use_multilayered_convolution && use_default_n) set_default_n_common();

	UninitializeAll();

	return error;
}

//enable multi-layered convolution and force it to 2D for all layers
BError SDemag::Set_2D_Multilayered_Convolution(int status)
{
	BError error(CLASS_STR(SDemag));

	use_multilayered_convolution = true;

	force_2d_convolution = status;

	if (force_2d_convolution) n_common.z = 1;
	else if (use_default_n) set_default_n_common();

	//first clear all currently set SDemag_Demag modules - these will be created as required through the UpdateConfiguration() method below.
	Destroy_SDemag_Demag_Modules();

	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	UninitializeAll();

	return error;
}

//set n_common for multi-layered convolution
BError SDemag::Set_n_common(SZ3 n)
{
	BError error(CLASS_STR(SDemag));

	n_common = n;

	use_multilayered_convolution = true;
	use_default_n = false;

	if (n_common.z == 1) force_2d_convolution = 1;
	else force_2d_convolution = 0;

	//first clear all currently set SDemag_Demag modules - these will be created as required through the UpdateConfiguration() method below.
	//we need to do this since it's possible force_2d_convolution mode was changed e.g. from 2 to 1.
	Destroy_SDemag_Demag_Modules();

	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	UninitializeAll();

	return error;
}

//set status for use_default_n
BError SDemag::Set_Default_n_status(bool status)
{
	BError error(CLASS_STR(SDemag));

	use_multilayered_convolution = true;
	use_default_n = status;

	if (use_default_n) set_default_n_common();

	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	UninitializeAll();

	return error;
}

//Set PBC images for supermesh demag
BError SDemag::Set_PBC(INT3 demag_pbc_images_)
{
	BError error(__FUNCTION__);

	demag_pbc_images = demag_pbc_images_;

	error = Set_M_PBC();

	UninitializeAll();

	//update will be needed if pbc settings have changed
	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	return error;
}

//Set PBC settings for M in all meshes
BError SDemag::Set_M_PBC(void)
{
	BError error(__FUNCTION__);

	//set pbc conditions in all M : if any are zero then pbc is disabled in that dimension

	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		//Set PBC irrespective of demag exclusion setting
		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			(*pSMesh)[idx]->Set_Magnetic_PBC(demag_pbc_images);
		}
	}

	return error;
}

//-------------------Abstract base class method implementations

BError SDemag::Initialize(void)
{
	BError error(CLASS_STR(SDemag));
	
	if (!use_multilayered_convolution) {

		//make sure to allocate memory for Hdemag if we need it
		if (pSMesh->GetEvaluationSpeedup()) Hdemag.resize(pSMesh->h_fm, pSMesh->sMeshRect_fm);
		else Hdemag.clear();
	}

	Hdemag_calculated = false;

	//FFT Kernels are not so quick to calculate - if already initialized then we are guaranteed they are correct
	if (!initialized) {

		//calculate kernels for super-mesh convolution
		if (!use_multilayered_convolution) {

			error = Calculate_Demag_Kernels();
			if (error) return error;

			Initialize_Mesh_Transfer();
		}
		else {

			//in multi-layered convolution mode must make sure all convolution sizes are set correctly, and rect collections also set
			//SDemag_Demag modules are initialized before SDemag, so they must check if SDemag is not initialized, in which case must call this
			//This will happen in the first SDemag_Demag module to initialize, so after that everything is set correctly to calculate kernels

			//update common discretisation if needed
			if (use_default_n) set_default_n_common();

			//make sure Rect_collection is correct
			set_Rect_collection();
			
			double h_max = get_maximum_cellsize();

			for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

				//h_convolution may differ from h_common in 2D mode
				DBL3 h_convolution = Rect_collection[idx] / n_common;

				if (!pSDemag_Demag[idx]->CheckDimensions(n_common, h_convolution, demag_pbc_images)) {

					//set convolution dimensions using the common discretisation
					//kernel collection must be used without multiplcation embedding. Calling this also sets full sizes for S and S2 scratch spaces.
					error = pSDemag_Demag[idx]->SetDimensions(n_common, h_convolution, false, demag_pbc_images);
					if (error) return error;
				}

				//set all rect collections
				error = pSDemag_Demag[idx]->Set_Rect_Collection(Rect_collection, Rect_collection[idx], h_max);
				if (error) return error;
			}
			
			//now everything is set correctly, ready to calculate demag kernel collections
		}

		//initialized ok.
		initialized = true;
	}

	//calculate total_nonempty_volume from all meshes participating in convolution
	if (use_multilayered_convolution) {

		total_nonempty_volume = 0.0;

		for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

			if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {

				total_nonempty_volume += pSMesh->pMesh[idx]->Get_NonEmpty_Magnetic_Volume();
			}
		}
	}

	return error;
}

//initialize transfer object for supermesh convolution
BError SDemag::Initialize_Mesh_Transfer(void)
{
	BError error(CLASS_STR(SDemag));

	//now calculate data required for mesh transfers, as well as demag corrections
	vector< VEC<DBL3>* > pVal_from, pVal_from2;
	vector< VEC<DBL3>* > pVal_to, pVal_to2;

	antiferromagnetic_meshes_present = false;

	//identify all existing (anti)ferrommagnetic meshes (magnetic computation enabled)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->is_atomistic()) {

			pVal_from.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->M));
			pVal_to.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff));

			pVal_from2.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->M2));
			pVal_to2.push_back(&(dynamic_cast<Mesh*>((*pSMesh)[idx])->Heff2));

			if ((*pSMesh)[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) antiferromagnetic_meshes_present = true;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// ALL FERROMAGNETIC MESHES //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!antiferromagnetic_meshes_present) {

		//Initialize the mesh transfer object for convolution on the super-mesh

		//use built-in corrections based on interpolation
		if (!sm_Vals.Initialize_MeshTransfer(pVal_from, pVal_to, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//the Hdemag.size() check is needed : if in CUDA mode, this method will be called to initialize mesh transfer in sm_Vals so the transfer info can be copied over to the gpu
		//In this case it's possible Hdemag does not have the correct memory allocated; if not in CUDA mode we first pass through Initialization method before calling this, in which case Hdemag will be sized correctly.
		if (pSMesh->GetEvaluationSpeedup() && Hdemag.size() == sm_Vals.size()) {

			//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
			if (!Hdemag.Initialize_MeshTransfer(pVal_from, pVal_to, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
		}

		//transfer values from invidual M meshes to sm_Vals - we need this to get number of non-empty cells
		sm_Vals.transfer_in();

		non_empty_cells = sm_Vals.get_nonempty_cells();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// AT LEAST ONE ANTIFERROMAGNETIC MESH ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//Initialize the mesh transfer object for convolution on the super-mesh

		//use built-in corrections based on interpolation
		if (!sm_Vals.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//the Hdemag.size() check is needed : if in CUDA mode, this method will be called to initialize mesh transfer in sm_Vals so the transfer info can be copied over to the gpu
		//In this case it's possible Hdemag does not have the correct memory allocated; if not in CUDA mode we first pass through Initialization method before calling this, in which case Hdemag will be sized correctly.
		if (pSMesh->GetEvaluationSpeedup() && Hdemag.size() == sm_Vals.size()) {

			//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
			if (!Hdemag.Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
		}

		//transfer values from invidual M meshes to sm_Vals - we need this to get number of non-empty cells
		//NOTE : do not use transfer_in_averaged here as for antiferromagnetic meshes in the ground state this will result in zero values everywhere, looking like there are no non-empty cells
		sm_Vals.transfer_in();

		non_empty_cells = sm_Vals.get_nonempty_cells();
	}

	//avoid division by zero
	if (!non_empty_cells) non_empty_cells = 1;

	return error;
}

void SDemag::UninitializeAll(void)
{
	Uninitialize();

	//be careful when using UninitializeAll : pSDemag_Demag must be up to date
	if (use_multilayered_convolution) {

		//Must have called UpdateConfiguration before, which makes sure pSDemag_Demag vector is correct
		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

			pSDemag_Demag[idx]->Uninitialize();
		}
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		//Must have called UpdateConfiguration before - in turn it will have called the CUDA version of UpdateConfiguration, which makes sure the pSDemagCUDA_Demag vector is correct
		dynamic_cast<SDemagCUDA*>(pModuleCUDA)->UninitializeAll();
	}
#endif
}

BError SDemag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemag));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_DEMAG_CONVCHANGE, UPDATECONFIG_SMESH_CELLSIZE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hdemag.clear();
		Hdemag_calculated = false;

		if (!use_multilayered_convolution) {

			//don't need memory allocated for multi-layered convolution
			Destroy_SDemag_Demag_Modules();

			//for super-mesh convolution just need a single convolution and sm_Vals to be sized correctly

			//only need to uninitialize if n_fm or h_fm have changed
			if (!CheckDimensions(pSMesh->n_fm, pSMesh->h_fm, demag_pbc_images) || cfgMessage == UPDATECONFIG_DEMAG_CONVCHANGE) {

				Uninitialize();
				error = SetDimensions(pSMesh->n_fm, pSMesh->h_fm, true, demag_pbc_images);

				if (!sm_Vals.resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFMEMORY_CRIT);
			}
		}
		else {

			//don't need memory allocated for supermesh demag
			sm_Vals.clear();
			SetDimensions(SZ3(1), DBL3());

			if (use_default_n) set_default_n_common();

			//for multi-layered convolution need a convolution object for each mesh (make modules) - or multiple if the 2D layering option is enabled.

			//new (anti)ferromagnetic meshes could have been added
			Create_SDemag_Demag_Modules();

			//make sure the list of SDemag_Demag  modules is up to date. If it was not, must re-initialize.
			if (!Update_SDemag_Demag_List()) Uninitialize();

			//If SDemag or any SDemag_Demag modules are uninitialized, then Uninitialize all SDemag_Demag modules
			for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

				initialized &= pSDemag_Demag[idx]->IsInitialized();
			}
		}

		//if a new mesh has been added we must also set any possible pbc conditions for M
		error = Set_M_PBC();
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	//important this is at the end, so the CUDA version of UpdateConfiguration is executed before
	if (!initialized) UninitializeAll();

	return error;
}

BError SDemag::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SDemag));

#if COMPILECUDA == 1

	pModuleCUDA = new SDemagCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

	Uninitialize();

#endif

	return error;
}

double SDemag::UpdateField(void)
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// SUPERMESH CONVOLUTION ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!use_multilayered_convolution) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////// EVAL SPEEDUP - SUPERMESH ///////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (pSMesh->GetEvaluationSpeedup()) {

			//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

			int update_type = pSMesh->Check_Step_Update();

			if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

				//calculate field required

				if (update_type == EVALSPEEDUPSTEP_COMPUTE_AND_SAVE) {

					//calculate field and save it for next time : we'll need to use it (expecting update_type = EVALSPEEDUPSTEP_SKIP next time)

					if (!antiferromagnetic_meshes_present) {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals.transfer_in();
					}
					else {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals.transfer_in_averaged();
					}

					//convolute and get "energy" value
					energy = Convolute(sm_Vals, Hdemag, true);

					Hdemag_calculated = true;

					//finish off energy value
					energy *= -MU0 / (2 * non_empty_cells);
				}
				else {

					//calculate field but do not save it for next time : we'll need to recalculate it again (expecting update_type != EVALSPEEDUPSTEP_SKIP again next time : EVALSPEEDUPSTEP_COMPUTE_NO_SAVE or EVALSPEEDUPSTEP_COMPUTE_AND_SAVE)

					if (!antiferromagnetic_meshes_present) {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals.transfer_in();

						//convolute and get "energy" value
						energy = Convolute(sm_Vals, sm_Vals, true);

						//good practice to set this to false
						Hdemag_calculated = false;

						//finish off energy value
						energy *= -MU0 / (2 * non_empty_cells);

						//transfer to individual Heff meshes
						sm_Vals.transfer_out();
					}
					else {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals.transfer_in_averaged();

						//convolute and get "energy" value
						energy = Convolute(sm_Vals, sm_Vals, true);

						//good practice to set this to false
						Hdemag_calculated = false;

						//finish off energy value
						energy *= -MU0 / (2 * non_empty_cells);

						//transfer to individual Heff meshes
						sm_Vals.transfer_out_duplicated();
					}

					//return here to avoid adding Hdemag to Heff : we've already added demag field contribution
					return energy;
				}
			}

			if (!antiferromagnetic_meshes_present) {

				//transfer contribution to meshes Heff
				Hdemag.transfer_out();
			}
			else {

				//transfer contribution to meshes Heff
				Hdemag.transfer_out_duplicated();
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// NO SPEEDUP - SUPERMESH ///////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		else {

			//don't use evaluation speedup, so no need to use Hdemag (this won't have memory allocated anyway)

			if (!antiferromagnetic_meshes_present) {

				//transfer values from invidual M meshes to sm_Vals
				sm_Vals.transfer_in();

				//convolution with demag kernels, output overwrites in sm_Vals
				energy = Convolute(sm_Vals, sm_Vals, true);

				//finish off energy value
				energy *= -MU0 / (2 * non_empty_cells);

				//transfer to individual Heff meshes
				sm_Vals.transfer_out();
			}
			else {

				//transfer values from invidual M meshes to sm_Vals
				sm_Vals.transfer_in_averaged();

				//convolution with demag kernels, output overwrites in sm_Vals
				energy = Convolute(sm_Vals, sm_Vals, true);

				//finish off energy value
				energy *= -MU0 / (2 * non_empty_cells);

				//transfer to individual Heff meshes
				sm_Vals.transfer_out_duplicated();
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// MULTILAYERED CONVOLUTION //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		///////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////// EVAL SPEEDUP - MULTILAYERED ////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (pSMesh->GetEvaluationSpeedup()) {

			//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

			int update_type = pSMesh->Check_Step_Update();

			if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

				//calculate field required

				//Forward FFT for all ferromagnetic meshes
				for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemag_Demag[idx]->pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						if (pSDemag_Demag[idx]->do_transfer) {

							//transfer from M to common meshing
							pSDemag_Demag[idx]->transfer.transfer_in_averaged();

							//do forward FFT
							pSDemag_Demag[idx]->ForwardFFT(pSDemag_Demag[idx]->transfer);
						}
						else {

							pSDemag_Demag[idx]->ForwardFFT_AveragedInputs(pSDemag_Demag[idx]->pMesh->M, pSDemag_Demag[idx]->pMesh->M2);
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {

						if (pSDemag_Demag[idx]->do_transfer) {

							//transfer from M to common meshing
							pSDemag_Demag[idx]->transfer.transfer_in();

							//do forward FFT
							pSDemag_Demag[idx]->ForwardFFT(pSDemag_Demag[idx]->transfer);
						}
						else {

							pSDemag_Demag[idx]->ForwardFFT(pSDemag_Demag[idx]->pMesh->M);
						}
					}
				}
		
				//Kernel multiplications for multiple inputs. Reverse loop ordering improves cache use at both ends.
				for (int idx = pSDemag_Demag.size() - 1; idx >= 0; idx--) {

					pSDemag_Demag[idx]->KernelMultiplication_MultipleInputs(FFT_Spaces_Input);
				}

				if (update_type == EVALSPEEDUPSTEP_COMPUTE_AND_SAVE) {

					//calculate field and save it for next time : we'll need to use it (expecting update_type = EVALSPEEDUPSTEP_SKIP next time)

					energy = 0;

					//Inverse FFT
					for (int idx_mesh = 0; idx_mesh < pSDemag_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemag_Demag[idx_mesh]->pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemag_Demag[idx_mesh]->do_transfer) {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT(pSDemag_Demag[idx_mesh]->transfer, pSDemag_Demag[idx_mesh]->Hdemag, true) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;
							}
							else {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT_AveragedInputs(
									pSDemag_Demag[idx_mesh]->pMesh->M, pSDemag_Demag[idx_mesh]->pMesh->M2,
									pSDemag_Demag[idx_mesh]->Hdemag, true) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemag_Demag[idx_mesh]->do_transfer) {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT(pSDemag_Demag[idx_mesh]->transfer, pSDemag_Demag[idx_mesh]->Hdemag, true) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;
							}
							else {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT(pSDemag_Demag[idx_mesh]->pMesh->M, pSDemag_Demag[idx_mesh]->Hdemag, true) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;
							}
						}
					}

					//finish off energy value
					energy *= -MU0 / 2;

					Hdemag_calculated = true;
				}
				else {

					//calculate field but do not save it for next time : we'll need to recalculate it again (expecting update_type != EVALSPEEDUPSTEP_SKIP again next time : EVALSPEEDUPSTEP_COMPUTE_NO_SAVE or EVALSPEEDUPSTEP_COMPUTE_AND_SAVE)

					energy = 0;

					//Inverse FFT
					for (int idx_mesh = 0; idx_mesh < pSDemag_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemag_Demag[idx_mesh]->pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemag_Demag[idx_mesh]->do_transfer) {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT(pSDemag_Demag[idx_mesh]->transfer, pSDemag_Demag[idx_mesh]->transfer, true) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;

								//transfer to Heff in each mesh
								pSDemag_Demag[idx_mesh]->transfer.transfer_out_duplicated();
							}
							else {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT_AveragedInputs_DuplicatedOutputs(
									pSDemag_Demag[idx_mesh]->pMesh->M, pSDemag_Demag[idx_mesh]->pMesh->M2, 
									pSDemag_Demag[idx_mesh]->pMesh->Heff, pSDemag_Demag[idx_mesh]->pMesh->Heff2, false) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemag_Demag[idx_mesh]->do_transfer) {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT(pSDemag_Demag[idx_mesh]->transfer, pSDemag_Demag[idx_mesh]->transfer, true) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;

								//transfer to Heff in each mesh
								pSDemag_Demag[idx_mesh]->transfer.transfer_out();
							}
							else {

								//do inverse FFT and accumulate energy
								energy += (pSDemag_Demag[idx_mesh]->InverseFFT(pSDemag_Demag[idx_mesh]->pMesh->M, pSDemag_Demag[idx_mesh]->pMesh->Heff, false) / pSDemag_Demag[idx_mesh]->non_empty_cells) * pSDemag_Demag[idx_mesh]->energy_density_weight;
							}
						}
					}

					//finish off energy value
					energy *= -MU0 / 2;

					//good practice to set this to false
					Hdemag_calculated = false;

					//return here to avoid adding Hdemag to Heff : we've already added demag field contribution
					return energy;
				}
			}

			//add contribution to Heff in each mesh
			for (int idx_mesh = 0; idx_mesh < pSDemag_Demag.size(); idx_mesh++) {

				///////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				if (pSDemag_Demag[idx_mesh]->pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					if (pSDemag_Demag[idx_mesh]->do_transfer) {

						//transfer to Heff in each mesh
						pSDemag_Demag[idx_mesh]->Hdemag.transfer_out_duplicated();
					}
					else {

						//add contribution to Heff in each mesh
						#pragma omp parallel for
						for (int idx_point = 0; idx_point < pSDemag_Demag[idx_mesh]->Hdemag.linear_size(); idx_point++) {

							pSDemag_Demag[idx_mesh]->pMesh->Heff[idx_point] += pSDemag_Demag[idx_mesh]->Hdemag[idx_point];
							pSDemag_Demag[idx_mesh]->pMesh->Heff2[idx_point] += pSDemag_Demag[idx_mesh]->Hdemag[idx_point];
						}
					}
				}

				///////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				else {

					if (pSDemag_Demag[idx_mesh]->do_transfer) {

						//transfer to Heff in each mesh
						pSDemag_Demag[idx_mesh]->Hdemag.transfer_out();
					}
					else {

						//add contribution to Heff in each mesh
#pragma omp parallel for
						for (int idx_point = 0; idx_point < pSDemag_Demag[idx_mesh]->Hdemag.linear_size(); idx_point++) {

							pSDemag_Demag[idx_mesh]->pMesh->Heff[idx_point] += pSDemag_Demag[idx_mesh]->Hdemag[idx_point];
						}
					}
				}
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////// NO SPEEDUP - MULTILAYERED /////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		else {

			//don't use evaluation speedup, so no need to use Hdemag in SDemag_Demag modules (this won't have memory allocated anyway)
			
			//Forward FFT for all ferromagnetic meshes
			for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

				///////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				if (pSDemag_Demag[idx]->pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					if (pSDemag_Demag[idx]->do_transfer) {

						//transfer from M to common meshing
						pSDemag_Demag[idx]->transfer.transfer_in_averaged();

						//do forward FFT
						pSDemag_Demag[idx]->ForwardFFT(pSDemag_Demag[idx]->transfer);
					}
					else {

						pSDemag_Demag[idx]->ForwardFFT_AveragedInputs(pSDemag_Demag[idx]->pMesh->M, pSDemag_Demag[idx]->pMesh->M2);
					}
				}

				///////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				else {

					if (pSDemag_Demag[idx]->do_transfer) {

						//transfer from M to common meshing
						pSDemag_Demag[idx]->transfer.transfer_in();

						//do forward FFT
						pSDemag_Demag[idx]->ForwardFFT(pSDemag_Demag[idx]->transfer);
					}
					else {

						pSDemag_Demag[idx]->ForwardFFT(pSDemag_Demag[idx]->pMesh->M);
					}
				}
			}

			//Kernel multiplications for multiple inputs. Reverse loop ordering improves cache use at both ends.
			for (int idx = pSDemag_Demag.size() - 1; idx >= 0; idx--) {

				pSDemag_Demag[idx]->KernelMultiplication_MultipleInputs(FFT_Spaces_Input);
			}

			energy = 0;

			//Inverse FFT
			for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

				///////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				if (pSDemag_Demag[idx]->pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					if (pSDemag_Demag[idx]->do_transfer) {

						//do inverse FFT and accumulate energy
						energy += (pSDemag_Demag[idx]->InverseFFT(pSDemag_Demag[idx]->transfer, pSDemag_Demag[idx]->transfer, true) / pSDemag_Demag[idx]->non_empty_cells) * pSDemag_Demag[idx]->energy_density_weight;

						//transfer to Heff in each mesh
						pSDemag_Demag[idx]->transfer.transfer_out_duplicated();
					}
					else {

						//do inverse FFT and accumulate energy
						energy += (pSDemag_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
							pSDemag_Demag[idx]->pMesh->M, pSDemag_Demag[idx]->pMesh->M2,
							pSDemag_Demag[idx]->pMesh->Heff, pSDemag_Demag[idx]->pMesh->Heff2, false) / pSDemag_Demag[idx]->non_empty_cells) * pSDemag_Demag[idx]->energy_density_weight;
					}
				}

				///////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				else {

					if (pSDemag_Demag[idx]->do_transfer) {

						//do inverse FFT and accumulate energy
						energy += (pSDemag_Demag[idx]->InverseFFT(pSDemag_Demag[idx]->transfer, pSDemag_Demag[idx]->transfer, true) / pSDemag_Demag[idx]->non_empty_cells) * pSDemag_Demag[idx]->energy_density_weight;

						//transfer to Heff in each mesh
						pSDemag_Demag[idx]->transfer.transfer_out();
					}
					else {

						//do inverse FFT and accumulate energy
						energy += (pSDemag_Demag[idx]->InverseFFT(pSDemag_Demag[idx]->pMesh->M, pSDemag_Demag[idx]->pMesh->Heff, false) / pSDemag_Demag[idx]->non_empty_cells) * pSDemag_Demag[idx]->energy_density_weight;
					}
				}
			}

			//finish off energy value
			energy *= -MU0 / 2;
		}
	}

	return energy;
}

#endif