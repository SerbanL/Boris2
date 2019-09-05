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

			for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

				//h_convolution may differ from h_common in 2D mode
				cuReal3 h_convolution = Rect_collection[idx] / (cuSZ3)pSDemag->n_common;

				if (!pSDemagCUDA_Demag[idx]->CheckDimensions((cuSZ3)pSDemag->n_common, h_convolution, pSDemag->Get_PBC())) {

					//set convolution dimensions using the common discretisation
					//kernel collection must be used without multiplcation embedding. Calling this also sets full sizes for S and S2 scratch spaces.
					error = pSDemagCUDA_Demag[idx]->SetDimensions((cuSZ3)pSDemag->n_common, h_convolution, false, pSDemag->Get_PBC());
					if (error) return error;
				}

				//set all rect collections
				error = pSDemagCUDA_Demag[idx]->Set_Rect_Collection(Rect_collection, Rect_collection[idx], h_max);
				if (error) return error;
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
		if (!CheckDimensions(pSMesh->n_fm, pSMesh->h_fm, pSDemag->Get_PBC()) || cfgMessage == UPDATECONFIG_FORCEUPDATE) {

			Uninitialize();
			error = SetDimensions(pSMesh->n_fm, pSMesh->h_fm, true, pSDemag->Get_PBC());

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