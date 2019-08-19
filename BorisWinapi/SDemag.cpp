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

	error_on_create = UpdateConfiguration();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

SDemag::~SDemag()
{
	//RAII : SDemag_Demag module were created in the constructor, so delete them here in any remaining ferromagnetic meshes 
	//(some could have been deleted already if any ferromagnetic mesh was deleted in the mean-time)
	Destroy_SDemag_Demag_Modules();
}

//------------------ Helpers for multi-layered convolution control

//when SDemag created, it needs to add one SDemag_Demag module to each ferromagnetic mesh
BError SDemag::Create_SDemag_Demag_Modules(void)
{
	BError error(CLASS_STR(SDemag));

	//identify all existing ferrommagnetic meshes (magnetic computation enabled)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			//Make sure each ferromagnetic mesh has one SDemag_Demag module added
			if (!(*pSMesh)[idx]->IsModuleSet(MOD_SDEMAG_DEMAG)) {

				error = (*pSMesh)[idx]->AddModule(MOD_SDEMAG_DEMAG);
			}

			(*pSMesh)[idx]->CallModuleMethod(&SDemag_Demag::Set_SDemag_Pointer, this);
		}
	}

	return error;
}

//when SDemag destroyed, it must destroy the SDemag_Demag module from each ferromagnetic mesh
void SDemag::Destroy_SDemag_Demag_Modules(void)
{
	//identify all existing ferrommagnetic meshes (magnetic computation enabled)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			(*pSMesh)[idx]->DelModule(MOD_SDEMAG_DEMAG);
		}
	}

	FFT_Spaces_Input.clear();
	FFT_Spaces_Output.clear();
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
	FFT_Spaces_Output.clear();
	kernel_collection.clear();

	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			SDemag_Demag *pSDemag_Demag_Module = dynamic_cast<SDemag_Demag*>((*pSMesh)[idx]->GetModule(MOD_SDEMAG_DEMAG));

			pSDemag_Demag_.push_back(pSDemag_Demag_Module);

			//build the fft spaces, rectangles list, and demag kernels
			FFT_Spaces_Input.push_back(pSDemag_Demag_.back()->Get_Input_Scratch_Space());
			FFT_Spaces_Output.push_back(pSDemag_Demag_.back()->Get_Output_Scratch_Space());
			Rect_collection.push_back((*pSMesh)[idx]->GetMeshRect());
			kernel_collection.push_back(dynamic_cast<DemagKernelCollection*>(pSDemag_Demag_.back()));
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

//set default value for n_common : largest value from all ferromagnetic meshes
void SDemag::set_default_n_common(void)
{
	SZ3 n_common_old = n_common;

	n_common = SZ3(1);

	bool layers_2d = true;

	//are all the layers 2d?
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			if ((*pSMesh)[idx]->n.z != 1) {

				layers_2d = false;
				break;
			}
		}
	}

	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			if (layers_2d || force_2d_convolution) {

				//2D layers
				if (n_common.x < (*pSMesh)[idx]->n.x) n_common.x = (*pSMesh)[idx]->n.x;
				if (n_common.y < (*pSMesh)[idx]->n.y) n_common.y = (*pSMesh)[idx]->n.y;
				n_common.z = 1;
			}
			else {

				if (n_common.x < (*pSMesh)[idx]->n.x) n_common.x = (*pSMesh)[idx]->n.x;
				if (n_common.y < (*pSMesh)[idx]->n.y) n_common.y = (*pSMesh)[idx]->n.y;
				if (n_common.z < (*pSMesh)[idx]->n.z) n_common.z = (*pSMesh)[idx]->n.z;
			}
		}
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

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			Rect_collection.push_back((*pSMesh)[idx]->GetMeshRect());
		}
	}

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
	error = UpdateConfiguration();

	if (use_multilayered_convolution && use_default_n) set_default_n_common();

	UninitializeAll();

	return error;
}

//enable multi-layered convolution and force it to 2D for all layers
BError SDemag::Set_2D_Multilayered_Convolution(bool status)
{
	BError error(CLASS_STR(SDemag));

	use_multilayered_convolution = true;
	force_2d_convolution = status;

	if (force_2d_convolution) n_common.z = 1;
	else if (use_default_n) set_default_n_common();

	error = UpdateConfiguration();

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

	if (n_common.z == 1) force_2d_convolution = true;
	else force_2d_convolution = false;

	error = UpdateConfiguration();

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

	error = UpdateConfiguration();

	UninitializeAll();

	return error;
}

//Set PBC images for supermesh demag
void SDemag::Set_PBC(INT3 demag_pbc_images_)
{
	demag_pbc_images = demag_pbc_images_;

	//update will be needed if pbc settings have changed
	UpdateConfiguration();
}

//-------------------Abstract base class method implementations

BError SDemag::Initialize(void)
{
	BError error(CLASS_STR(SDemag));

	/*
	//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE
	if (use_multilayered_convolution && initialized && !Kernels.size()) {

		//for multi-layered convolution, after all SDemag_Demag modules have initialized gather kernel collection sorted by kernel here. SDemag will already have been initialized.
		//Note it is important that Kernels vector was cleared on first initialization

		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

			//not yet
			if (!pSDemag_Demag[idx]->IsInitialized()) return error;
		}

		//at this point everything is initialized so we have everything we need to start sorting kernels

		//multiple inputs version : the space itself is the output

		auto add_kernel_entry = [&](int idx_in, int idx_out) -> void
		{
			//go through all existing entries in Kernels
			for (int idx_ker = 0; idx_ker < Kernels.size(); idx_ker++) {

				//if matching kernel found add new In, Out pair, then return
				if (Kernels[idx_ker].add_entry_if_kernel_matches(
					pSDemag_Demag[idx_out]->Get_Kernel(idx_in),
					pSDemag_Demag[idx_in]->Get_Input_Scratch_Space(),
					pSDemag_Demag[idx_out]->Get_Output_Scratch_Space(),
					pSDemag_Demag[idx_out]->is_inverse_shifted(idx_in))) {

					return;
				}
			}

			//no match found so add new kernel entry
			Kernels.push_back(
				KerTypeCollection(
					pSDemag_Demag[idx_out]->Get_Kernel(idx_in),
					pSDemag_Demag[idx_in]->Get_Input_Scratch_Space(),
					pSDemag_Demag[idx_out]->Get_Output_Scratch_Space(),
					pSDemag_Demag[idx_out]->is_inverse_shifted(idx_in)));
		};

		//make sure the first entries are the self demag kernels - these set the outputs, everything else add to outputs, so must be first
		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

			add_kernel_entry(idx, idx);
		}

		//now go through each kernel and add new entry in Kernels if no match found, else add new entry in KerTypeCollection if match found
		//not an efficient way to do this but it's simple; an efficient method is not needed since we don't have to deal with that many layers
		for (int idx_out = 0; idx_out < pSDemag_Demag.size(); idx_out++) {
			for (int idx_in = 0; idx_in < pSDemag_Demag.size(); idx_in++) {

				//skip diagonal elements (the self demag elements) as we've already done them
				if (idx_in == idx_out) continue;

				add_kernel_entry(idx_in, idx_out);
			}
		}

		return error;
	}
	*/

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

			//will collect Kernels at the end
			Kernels.clear();

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

	//make sure the energy density weights are correct
	if (use_multilayered_convolution) {

		double total_nonempty_volume = 0.0;

		for (int idx = 0; idx < (int)pSDemag_Demag.size(); idx++) {

			total_nonempty_volume += (double)pSDemag_Demag[idx]->pMesh->M.get_nonempty_cells() * pSDemag_Demag[idx]->pMesh->M.h.dim();
		}

		if (total_nonempty_volume) {

			for (int idx = 0; idx < (int)pSDemag_Demag.size(); idx++) {

				pSDemag_Demag[idx]->energy_density_weight = (double)pSDemag_Demag[idx]->pMesh->M.get_nonempty_cells() * pSDemag_Demag[idx]->pMesh->M.h.dim() / total_nonempty_volume;
			}
		}
	}

	return error;
}

//initialize transfer object
BError SDemag::Initialize_Mesh_Transfer(void)
{
	BError error(CLASS_STR(SDemag));

	//now calculate data required for mesh transfers, as well as demag corrections
	vector< VEC<DBL3>* > pVal_from;
	vector< VEC<DBL3>* > pVal_to;

	//identify all existing ferrommagnetic meshes (magnetic computation enabled)
	for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

		if ((*pSMesh)[idx]->MComputation_Enabled()) {

			pVal_from.push_back(&((*pSMesh)[idx]->M));
			pVal_to.push_back(&((*pSMesh)[idx]->Heff));
		}
	}

	//Initialize the mesh transfer object for convolution on the super-mesh

	//use built-in corrections based on interpolation
	if (!sm_Vals.Initialize_MeshTransfer(pVal_from, pVal_to, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

	//transfer values from invidual M meshes to sm_Vals - we need this to get number of non-empty cells
	sm_Vals.transfer_in();

	non_empty_cells = sm_Vals.get_nonempty_cells();

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
		reinterpret_cast<SDemagCUDA*>(pModuleCUDA)->UninitializeAll();
	}
#endif
}

BError SDemag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemag));

	if (!use_multilayered_convolution) {

		//don't need memory allocated for multi-layered convolution
		Destroy_SDemag_Demag_Modules();

		//for super-mesh convolution just need a single convolution and sm_Vals to be sized correctly

		//only need to uninitialize if n_fm or h_fm have changed
		if (!CheckDimensions(pSMesh->n_fm, pSMesh->h_fm, demag_pbc_images) || cfgMessage == UPDATECONFIG_FORCEUPDATE) {

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

		//for multi-layered convolution need a convolution object for each mesh (make modules)

		//new ferromagnetic meshes could have been added
		Create_SDemag_Demag_Modules();

		//make sure the list of SDemag_Demag  modules is up to date. If it was not, must re-initialize.
		if (!Update_SDemag_Demag_List()) Uninitialize();

		//If SDemag or any SDemag_Demag modules are uninitialized, then Uninitialize all SDemag_Demag modules
		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

			initialized &= pSDemag_Demag[idx]->IsInitialized();
		}
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
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
	if (!use_multilayered_convolution) {

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
		
		//Forward FFT for all ferromagnetic meshes
		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

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
		
		//Kernel multiplications for multiple inputs. Reverse loop ordering improves cache use at both ends.
		for (int idx = pSDemag_Demag.size() - 1; idx >= 0; idx--) {

			pSDemag_Demag[idx]->KernelMultiplication_MultipleInputs(FFT_Spaces_Input);
		}
		
		/*
		//Multiplication done by kernel type
		//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE
		if (n_common.z == 1) {

			for (int idx = 0; idx < Kernels.size(); idx++) {

				Kernels[idx].Kernel_Multiplication_2D();
			}
		}
		else {

			for (int idx = 0; idx < Kernels.size(); idx++) {

				Kernels[idx].Kernel_Multiplication_3D();
			}
		}
		*/

		energy = 0;

		//Inverse FFT
		for (int idx = 0; idx < pSDemag_Demag.size(); idx++) {

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

		//finish off energy value
		energy *= -MU0 / 2;
	}

	return energy;
}

#endif