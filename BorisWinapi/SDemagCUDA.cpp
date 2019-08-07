#include "stdafx.h"
#include "SDemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "ManagedMeshCUDA.h"

#include "SuperMesh.h"
#include "SDemag.h"

SDemagCUDA::SDemagCUDA(SuperMesh* pSMesh_, SDemag* pSDemag_) :
	ModulesCUDA(),
	ConvolutionCUDA<DemagKernelCUDA>()
{
	Uninitialize();

	pSMesh = pSMesh_;
	pSDemag = pSDemag_;

	error_on_create = UpdateConfiguration();
}

SDemagCUDA::~SDemagCUDA()
{
	//copy values back to cpu version
	if (Holder_Module_Available()) {

		if (!pSDemag->use_multilayered_convolution) {

			sm_Vals()->copy_to_cpuvec(pSDemag->sm_Vals);
		}

		//must force SDemag module to re-initialize as it's not properly initialized when CUDA module active even though its initialized flag is true.
		pSDemag->Uninitialize();
	}

	clear_vector(Kernels);
}

void SDemagCUDA::UninitializeAll(void)
{
	Uninitialize();

	//be careful when using UninitializeAll : pSDemagCUDA_Demag must be up to date
	if (pSDemag->use_multilayered_convolution) {

		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			pSDemagCUDA_Demag[idx]->Uninitialize();
		}
	}
}

BError SDemagCUDA::Initialize(void)
{
	BError error(CLASS_STR(SDemagCUDA));

	/*
	//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE
	if (pSDemag->use_multilayered_convolution && initialized && !Kernels.size()) {

		//for multi-layered convolution, after all SDemag_Demag modules have initialized gather kernel collection sorted by kernel here. SDemag will already have been initialized.
		//Note it is important that Kernels vector was cleared on first initialization

		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			//not yet
			if (!pSDemagCUDA_Demag[idx]->IsInitialized()) return error;
		}

		//at this point everything is initialized so we have everything we need to start sorting kernels

		//multiple inputs version : the space itself is the output

		auto add_kernel_entry = [&](int idx_in, int idx_out) -> void
		{
			//go through all existing entries in Kernels
			for (int idx_ker = 0; idx_ker < Kernels.size(); idx_ker++) {

				//if matching kernel found add new In, Out pair, then return
				if (Kernels[idx_ker]->add_entry_if_kernel_matches(
					pSDemagCUDA_Demag[idx_out]->Get_Kernel(idx_in),
					pSDemagCUDA_Demag[idx_in]->Get_Input_Scratch_Space_x(),
					pSDemagCUDA_Demag[idx_in]->Get_Input_Scratch_Space_y(),
					pSDemagCUDA_Demag[idx_in]->Get_Input_Scratch_Space_z(),
					pSDemagCUDA_Demag[idx_out]->Get_Output_Scratch_Space_x(),
					pSDemagCUDA_Demag[idx_out]->Get_Output_Scratch_Space_y(),
					pSDemagCUDA_Demag[idx_out]->Get_Output_Scratch_Space_z(),
					pSDemagCUDA_Demag[idx_out]->is_inverse_shifted(idx_in))) {

					return;
				}
			}
			
			//no match found so add new kernel entry
			Kernels.push_back(
				new KerTypeCollectionCUDA(
					pSDemagCUDA_Demag[idx_out]->Get_Kernel(idx_in),
					pSDemagCUDA_Demag[idx_in]->Get_Input_Scratch_Space_x(),
					pSDemagCUDA_Demag[idx_in]->Get_Input_Scratch_Space_y(),
					pSDemagCUDA_Demag[idx_in]->Get_Input_Scratch_Space_z(),
					pSDemagCUDA_Demag[idx_out]->Get_Output_Scratch_Space_x(),
					pSDemagCUDA_Demag[idx_out]->Get_Output_Scratch_Space_y(),
					pSDemagCUDA_Demag[idx_out]->Get_Output_Scratch_Space_z(),
					pSDemagCUDA_Demag[idx_out]->is_inverse_shifted(idx_in)));
		};

		//make sure the first entries are the self demag kernels - these set the outputs, everything else add to outputs, so must be first
		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			add_kernel_entry(idx, idx);
		}

		//now go through each kernel and add new entry in Kernels if no match found, else add new entry in KerTypeCollection if match found
		//not an efficient way to do this but it's simple; an efficient method is not needed since we don't have to deal with that many layers
		for (int idx_out = 0; idx_out < pSDemagCUDA_Demag.size(); idx_out++) {
			for (int idx_in = 0; idx_in < pSDemagCUDA_Demag.size(); idx_in++) {

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
		
		if (!pSDemag->use_multilayered_convolution) {

			error = Calculate_Demag_Kernels();
			if (error) return error;

			//array of pointers to input meshes (M) and oputput meshes (Heff) to transfer from and to
			cu_arr<cuVEC<cuReal3>> pVal_from;
			cu_arr<cuVEC<cuReal3>> pVal_to;

			//identify all existing ferrommagnetic meshes (magnetic computation enabled)
			for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

				if ((*pSMesh)[idx]->MComputation_Enabled()) {

					pVal_from.push_back((cuVEC<cuReal3>*&)(*pSMesh)[idx]->pMeshCUDA->M.get_managed_object());
					pVal_to.push_back((cuVEC<cuReal3>*&)(*pSMesh)[idx]->pMeshCUDA->Heff.get_managed_object());
				}
			}

			//Initialize the cpu mesh transfer object - note, SDemag::initialized is set true, but SDemag is not properly initialized - hence need the SDemagCUDA destructor to uninitialize SDemag
			error = pSDemag->Initialize_Mesh_Transfer();
			if (error) return error;

			//Now copy mesh transfer object to cuda version
			if (!sm_Vals()->copy_transfer_info(pVal_from, pVal_to, pSDemag->sm_Vals)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else {

			//in multi-layered convolution mode must make sure all convolution sizes are set correctly, and rect collections also set
			//SDemag_Demag modules are initialized before SDemag, so they must check if SDemag is not initialized, in which case must call this
			//This will happen in the first SDemag_Demag module to initialize, so after that everything is set correctly to calculate kernels

			//will collect Kernels at the end
			clear_vector(Kernels);

			//update common discretisation if needed
			if (pSDemag->use_default_n) pSDemag->set_default_n_common();

			//make sure Rect_collection is correct
			pSDemag->set_Rect_collection();

			//rectangle collection copy from SDemag
			Rect_collection.resize(pSDemag->Rect_collection.size());

			for (int idx = 0; idx < pSDemag->Rect_collection.size(); idx++) {

				Rect_collection[idx] = (cuRect)pSDemag->Rect_collection[idx];
			}

			cuReal h_max = (cuReal)pSDemag->get_maximum_cellsize();

			FFT_Spaces_x_Output.clear();
			FFT_Spaces_y_Output.clear();
			FFT_Spaces_z_Output.clear();

			for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

				//h_convolution may differ from h_common in 2D mode
				cuReal3 h_convolution = Rect_collection[idx] / (cuSZ3)pSDemag->n_common;

				if (!pSDemagCUDA_Demag[idx]->CheckDimensions((cuSZ3)pSDemag->n_common, h_convolution)) {

					//set convolution dimensions using the common discretisation
					//kernel collection must be used without multiplcation embedding. Calling this also sets full sizes for S and S2 scratch spaces.
					error = pSDemagCUDA_Demag[idx]->SetDimensions((cuSZ3)pSDemag->n_common, h_convolution, false);
					if (error) return error;
				}

				//set all rect collections
				error = pSDemagCUDA_Demag[idx]->Set_Rect_Collection(Rect_collection, Rect_collection[idx], h_max);
				if (error) return error;

				//now that output spaces have been allocated we can collect them in a cu_arr, ready to pass them to a __global__
				//This is important! FFT_Spaces_x_Output (etc.) is a cu_arr that stores the gpu memory locations of the FFT output arrays (cuS2)
				//These gpu memory locations can change when cuS2 is allocated, cleared, resized etc.
				//Thus only collect these after all memory has been allocated - above SetDimensions is the last call to change the cuS2 arrays memory locations
				FFT_Spaces_x_Output.push_back(pSDemagCUDA_Demag[idx]->Get_Output_Scratch_Space_x()->get_managed_array());
				FFT_Spaces_y_Output.push_back(pSDemagCUDA_Demag[idx]->Get_Output_Scratch_Space_y()->get_managed_array());
				FFT_Spaces_z_Output.push_back(pSDemagCUDA_Demag[idx]->Get_Output_Scratch_Space_z()->get_managed_array());
			}

			//now everything is set correctly, ready to calculate demag kernel collections
		}

		//initialized ok.
		initialized = true;
	}

	//make sure the energy density weights are correct
	if (pSDemag->use_multilayered_convolution) {

		double total_nonempty_volume = 0.0;

		for (int idx = 0; idx < (int)pSDemagCUDA_Demag.size(); idx++) {

			total_nonempty_volume += (double)pSDemag->pSDemag_Demag[idx]->pMesh->M.get_nonempty_cells() * pSDemag->pSDemag_Demag[idx]->pMesh->M.h.dim();
			pSDemagCUDA_Demag[idx]->energy_density_weight.from_cpu((cuReal)0.0);
		}

		if (total_nonempty_volume) {

			for (int idx = 0; idx < (int)pSDemagCUDA_Demag.size(); idx++) {

				double energy_density_weight = (double)pSDemag->pSDemag_Demag[idx]->pMesh->M.get_nonempty_cells() * pSDemag->pSDemag_Demag[idx]->pMesh->M.h.dim() / total_nonempty_volume;
				pSDemagCUDA_Demag[idx]->energy_density_weight.from_cpu((cuReal)energy_density_weight);
			}
		}
	}

	return error;
}

BError SDemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemagCUDA));
	
	if (!pSDemag->use_multilayered_convolution) {

		//only need to uninitialize if n or h have changed
		if (!CheckDimensions(pSMesh->n_fm, pSMesh->h_fm)) {

			Uninitialize();
			error = SetDimensions(pSMesh->n_fm, pSMesh->h_fm);

			if (!sm_Vals()->resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
	}
	else {

		//don't need memory allocated for supermesh demag
		SetDimensions(cuSZ3(1), cuReal3());
		sm_Vals()->clear();

		//rebuild kernel collection, and pSDemagCUDA_Demag
		kernel_collection.resize(pSDemag->kernel_collection.size());
		pSDemagCUDA_Demag.resize(pSDemag->pSDemag_Demag.size());
		FFT_Spaces_x_Input.resize(pSDemag->FFT_Spaces_Input.size());
		FFT_Spaces_y_Input.resize(pSDemag->FFT_Spaces_Input.size());
		FFT_Spaces_z_Input.resize(pSDemag->FFT_Spaces_Input.size());

		for (int idx = 0; idx < kernel_collection.size(); idx++) {

			kernel_collection[idx] = dynamic_cast<DemagKernelCollectionCUDA*>(pSDemag->pSDemag_Demag[idx]->pModuleCUDA);

			pSDemagCUDA_Demag[idx] = reinterpret_cast<SDemagCUDA_Demag*>(pSDemag->pSDemag_Demag[idx]->pModuleCUDA);

			FFT_Spaces_x_Input[idx] = pSDemagCUDA_Demag[idx]->Get_Input_Scratch_Space_x();
			FFT_Spaces_y_Input[idx] = pSDemagCUDA_Demag[idx]->Get_Input_Scratch_Space_y();
			FFT_Spaces_z_Input[idx] = pSDemagCUDA_Demag[idx]->Get_Input_Scratch_Space_z();
		}

		//mirror SDemag initialization flag
		initialized &= pSDemag->IsInitialized();

		//If SDemagCUDA or any SDemagCUDA_Demag modules are uninitialized, then Uninitialize all SDemagCUDA_Demag modules
		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			initialized &= pSDemagCUDA_Demag[idx]->IsInitialized();
		}

		if (!initialized) UninitializeAll();
	}
	
	return error;
}

void SDemagCUDA::UpdateField(void)
{
	if (!pSDemag->use_multilayered_convolution) {

		//transfer values from invidual M meshes to sm_Vals
		sm_Vals()->transfer_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());

		//only need energy after ode solver step finished
		if (pSMesh->CurrentTimeStepSolved()) ZeroEnergy();

		Convolute(sm_Vals, sm_Vals, energy, pSMesh->CurrentTimeStepSolved(), true);

		//transfer to individual Heff meshes
		sm_Vals()->transfer_out(pSDemag->sm_Vals.size_transfer_out());
	}
	else {

		//Forward FFT for all ferromagnetic meshes
		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			if (pSDemagCUDA_Demag[idx]->do_transfer) {

				//transfer from M to common meshing
				pSDemagCUDA_Demag[idx]->transfer()->transfer_in(pSDemag->pSDemag_Demag[idx]->transfer.linear_size(), pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_in());

				//do forward FFT
				pSDemagCUDA_Demag[idx]->ForwardFFT(pSDemagCUDA_Demag[idx]->transfer);
			}
			else {

				pSDemagCUDA_Demag[idx]->ForwardFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M);
			}
		}

		//Kernel multiplications for multiple inputs.
		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			pSDemagCUDA_Demag[idx]->KernelMultiplication_MultipleInputs(FFT_Spaces_x_Input, FFT_Spaces_y_Input, FFT_Spaces_z_Input);
		}
		
		/*
		Testing only : multiple inputs version is always quicker (up to twice quicker in tests).
		//Kernel multiplications for multiple outputs

		//first one sets the outputs
		pSDemagCUDA_Demag[0]->KernelMultiplication_MultipleOutputs_Set(FFT_Spaces_x_Output, FFT_Spaces_y_Output, FFT_Spaces_z_Output);

		//the others add into outputs
		for (int idx = 1; idx < pSDemagCUDA_Demag.size(); idx++) {

			pSDemagCUDA_Demag[idx]->KernelMultiplication_MultipleOutputs_Add(FFT_Spaces_x_Output, FFT_Spaces_y_Output, FFT_Spaces_z_Output);
		}
		*/

		/*
		//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE
		//Multiplication done by kernel type
		if (pSDemag->n_common.z == 1) {

			for (int idx = 0; idx < Kernels.size(); idx++) {

				Kernels[idx]->Kernel_Multiplication_2D(false);
			}
		}
		else {

			for (int idx = 0; idx < Kernels.size(); idx++) {

				Kernels[idx]->Kernel_Multiplication_3D();
			}
		}
		*/

		if (pSMesh->CurrentTimeStepSolved()) {

			//only need energy after ode solver step finished
			ZeroEnergy();

			//Inverse FFT
			for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

				if (pSDemagCUDA_Demag[idx]->do_transfer) {
					
					//do inverse FFT and accumulate energy
					pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);

					//transfer to Heff in each mesh
					pSDemagCUDA_Demag[idx]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
				}
				else {

					//do inverse FFT and accumulate energy. Add to Heff in each mesh
					pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, false);
				}
			}
		}
		else {

			//no energy contribution needed

			//Inverse FFT
			for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

				if (pSDemagCUDA_Demag[idx]->do_transfer) {

					//do inverse FFT
					pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, energy, false, true);

					//transfer to Heff in each mesh
					pSDemagCUDA_Demag[idx]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
				}
				else {

					//do inverse FFT. Add to Heff in each mesh
					pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, energy, false, false);
				}
			}
		}
	}
}

#endif

#endif