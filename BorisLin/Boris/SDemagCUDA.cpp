#include "stdafx.h"
#include "SDemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "ManagedMeshCUDA.h"

#include "SuperMesh.h"
#include "SDemag.h"

SDemagCUDA::SDemagCUDA(SuperMesh* pSMesh_, SDemag* pSDemag_) :
	ModulesCUDA(),
	ConvolutionCUDA<SDemagCUDA, DemagKernelCUDA>()
{
	Uninitialize();

	pSMesh = pSMesh_;
	pSDemag = pSDemag_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
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
		
		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// SUPERMESH CONVOLUTION ////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (!pSDemag->use_multilayered_convolution) {

			error = Calculate_Demag_Kernels(true, pSMesh->Get_Kernel_Initialize_on_GPU());
			if (error) return error;

			//array of pointers to input meshes (M) and oputput meshes (Heff) to transfer from and to
			cu_arr<cuVEC<cuReal3>> pVal_from, pVal_from2;
			cu_arr<cuVEC<cuReal3>> pVal_to, pVal_to2;
			//atomistic meshes input / output
			cu_arr<cuVEC<cuReal3>> pVal_afrom, pVal_ato;

			//identify all existing ferrommagnetic meshes (magnetic computation enabled)
			for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

				if ((*pSMesh)[idx]->MComputation_Enabled()) {

					if (!(*pSMesh)[idx]->is_atomistic()) {

						//micromagnetic mesh

						pVal_from.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->M.get_managed_object());
						pVal_to.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff.get_managed_object());

						pVal_from2.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->M2.get_managed_object());
						pVal_to2.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff2.get_managed_object());
					}
					else {

						//atomistic mesh
						
						pVal_afrom.push_back((cuVEC<cuReal3>*&)dynamic_cast<Atom_Mesh*>((*pSMesh)[idx])->paMeshCUDA->M1.get_managed_object());
						pVal_ato.push_back((cuVEC<cuReal3>*&)dynamic_cast<Atom_Mesh*>((*pSMesh)[idx])->paMeshCUDA->Heff1.get_managed_object());
					}
				}
			}

			//Initialize the cpu mesh transfer object - note, SDemag::initialized is set true, but SDemag is not properly initialized - hence need the SDemagCUDA destructor to uninitialize SDemag
			error = pSDemag->Initialize_Mesh_Transfer();
			if (error) return error;

			if (pVal_from.size()) {

				///////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////// ALL FERROMAGNETIC MESHES //////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				if (!pSDemag->antiferromagnetic_meshes_present) {

					//Now copy mesh transfer object to cuda version
					if (!sm_Vals()->copy_transfer_info(pVal_from, pVal_to, pSDemag->sm_Vals.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				///////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////// AT LEAST ONE ANTIFERROMAGNETIC MESH ///////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				else {

					//Now copy mesh transfer object to cuda version
					if (!sm_Vals()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag->sm_Vals.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////// ATOMISTIC MESHES //////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			if (pVal_afrom.size()) {

				//Now copy mesh transfer object to cuda version
				if (!sm_Vals()->copy_transfer2_info(pVal_afrom, pVal_ato, pSDemag->sm_Vals.get_transfer2())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////// MULTILAYERED CONVOLUTION //////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

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

			cuBReal h_max = (cuBReal)pSDemag->get_maximum_cellsize();

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
				error = pSDemagCUDA_Demag[idx]->Set_Rect_Collection(Rect_collection, Rect_collection[idx], h_max, idx);
				if (error) return error;
			}

			//now everything is set correctly, ready to calculate demag kernel collections
		}

		//initialized ok.
		initialized = true;

		//mirror SDemag initialized flag
		pSDemag->initialized = true;
	}

	//calculate total_nonempty_volume from all meshes participating in convolution
	if (pSDemag->use_multilayered_convolution) {

		total_nonempty_volume = 0.0;

		for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

			if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->Get_Demag_Exclusion()) {
				
				total_nonempty_volume += pSMesh->pMesh[idx]->Get_NonEmpty_Magnetic_Volume();
			}
		}
	}

	num_Hdemag_saved = 0;

	return error;
}

BError SDemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemagCUDA));
	
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_DEMAG_CONVCHANGE, UPDATECONFIG_SMESH_CELLSIZE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {

		if (!pSDemag->use_multilayered_convolution) {

			//only need to uninitialize if n or h have changed
			if (!CheckDimensions(pSMesh->n_fm, pSMesh->h_fm, pSDemag->Get_PBC()) || cfgMessage == UPDATECONFIG_DEMAG_CONVCHANGE || cfgMessage == UPDATECONFIG_MESHCHANGE) {

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

				pSDemagCUDA_Demag[idx] = dynamic_cast<SDemagCUDA_Demag*>(pSDemag->pSDemag_Demag[idx]->pModuleCUDA);

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
	}
	
	return error;
}

void SDemagCUDA::UpdateField(void)
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// SUPERMESH CONVOLUTION ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!pSDemag->use_multilayered_convolution) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// NO SPEEDUP - SUPERMESH ///////////////////////////////////
		// eval speedup not used for supermesh convolution - possible, but not worth implementing it //
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (!pSDemag->antiferromagnetic_meshes_present) {

			//transfer values from invidual M meshes to sm_Vals
			if (pSDemag->sm_Vals.size_transfer_in()) sm_Vals()->transfer_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());
			//transfer from atomistic mesh (if any) - clear input only if there was no transfer from micromagnetic meshes, else add in
			if (pSDemag->sm_Vals.size_transfer2_in()) sm_Vals()->transfer2_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer2_in(), pSDemag->sm_Vals.size_transfer_in() == 0);

			//only need energy after ode solver step finished
			if (pSMesh->CurrentTimeStepSolved()) ZeroEnergy();

			Convolute(sm_Vals, sm_Vals, energy, pSMesh->CurrentTimeStepSolved(), true);

			//transfer to individual Heff meshes (micromagnetic and atomistc meshes)
			if (pSDemag->sm_Vals.size_transfer_out()) sm_Vals()->transfer_out(pSDemag->sm_Vals.size_transfer_out());
			if (pSDemag->sm_Vals.size_transfer2_out()) sm_Vals()->transfer2_out(pSDemag->sm_Vals.size_transfer2_out());
		}
		else {

			//transfer values from invidual M meshes to sm_Vals
			if (pSDemag->sm_Vals.size_transfer_in()) sm_Vals()->transfer_in_averaged(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());
			//transfer from atomistic mesh (if any) - clear input only if there was no transfer from micromagnetic meshes, else add in
			if (pSDemag->sm_Vals.size_transfer2_in()) sm_Vals()->transfer2_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer2_in(), pSDemag->sm_Vals.size_transfer_in() == 0);

			//only need energy after ode solver step finished
			if (pSMesh->CurrentTimeStepSolved()) ZeroEnergy();

			Convolute(sm_Vals, sm_Vals, energy, pSMesh->CurrentTimeStepSolved(), true);

			//transfer to individual Heff meshes (micromagnetic and atomistc meshes)
			if (pSDemag->sm_Vals.size_transfer_out()) sm_Vals()->transfer_out_duplicated(pSDemag->sm_Vals.size_transfer_out());
			if (pSDemag->sm_Vals.size_transfer2_out()) sm_Vals()->transfer2_out(pSDemag->sm_Vals.size_transfer2_out());
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// MULTILAYERED CONVOLUTION //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//Forward FFT for all ferromagnetic meshes
		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			///////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			if (pSDemagCUDA_Demag[idx]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				if (pSDemagCUDA_Demag[idx]->do_transfer) {

					//transfer from M to common meshing
					pSDemagCUDA_Demag[idx]->transfer()->transfer_in_averaged(pSDemag->pSDemag_Demag[idx]->transfer.linear_size(), pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_in());

					//do forward FFT
					pSDemagCUDA_Demag[idx]->ForwardFFT(pSDemagCUDA_Demag[idx]->transfer);
				}
				else {

					pSDemagCUDA_Demag[idx]->ForwardFFT_AveragedInputs(pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2);
				}
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			else {

				if (pSDemagCUDA_Demag[idx]->do_transfer) {

					//transfer from M to common meshing
					pSDemagCUDA_Demag[idx]->transfer()->transfer_in(pSDemag->pSDemag_Demag[idx]->transfer.linear_size(), pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_in());

					//do forward FFT
					pSDemagCUDA_Demag[idx]->ForwardFFT(pSDemagCUDA_Demag[idx]->transfer);
				}
				else {

					//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh
					pSDemagCUDA_Demag[idx]->ForwardFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M);
				}
			}
		}

		//Kernel multiplications for multiple inputs.
		for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

			pSDemagCUDA_Demag[idx]->KernelMultiplication_MultipleInputs(FFT_Spaces_x_Input, FFT_Spaces_y_Input, FFT_Spaces_z_Input);
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////// NO SPEEDUP - MULTILAYERED /////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (!pSMesh->GetEvaluationSpeedup() || (num_Hdemag_saved < pSMesh->GetEvaluationSpeedup() && !pSMesh->Check_Step_Update())) {

			//only need energy after ode solver step finished
			if (pSMesh->CurrentTimeStepSolved()) {

				double energy_cpu = 0.0;

				//Inverse FFT
				for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemagCUDA_Demag[idx]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						if (pSDemagCUDA_Demag[idx]->do_transfer) {

							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx]->InverseFFT(
									pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, true, true, &pSDemagCUDA_Demag[idx]->transfer_Module_Heff, &pSDemagCUDA_Demag[idx]->transfer_Module_energy);

								pSDemagCUDA_Demag[idx]->transfer_Module_Heff()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.size_transfer_out());
								pSDemagCUDA_Demag[idx]->transfer_Module_energy()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_energy.size_transfer_out());
							}
							else {

								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, true, true);
							}

							//transfer to Heff in each mesh
							pSDemagCUDA_Demag[idx]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
						}
						else {

							if (pSDemag->pSDemag_Demag[idx]->Module_Heff.linear_size()) {

								//do inverse FFT and accumulate energy. Add to Heff in each mesh
								pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2,
									pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx]->energy, true, false, &pSDemagCUDA_Demag[idx]->Module_Heff, &pSDemagCUDA_Demag[idx]->Module_energy);
							}
							else {

								//do inverse FFT and accumulate energy. Add to Heff in each mesh
								pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2,
									pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx]->energy, true, false);
							}
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {

						if (pSDemagCUDA_Demag[idx]->do_transfer) {

							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx]->InverseFFT(
									pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, true, true, &pSDemagCUDA_Demag[idx]->transfer_Module_Heff, &pSDemagCUDA_Demag[idx]->transfer_Module_energy);

								pSDemagCUDA_Demag[idx]->transfer_Module_Heff()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.size_transfer_out());
								pSDemagCUDA_Demag[idx]->transfer_Module_energy()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_energy.size_transfer_out());
							}
							else {

								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, true, true);
							}

							//transfer to Heff in each mesh
							pSDemagCUDA_Demag[idx]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
						}
						else {

							//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

							//do inverse FFT and accumulate energy. Add to Heff in each mesh
							if (pSDemag->pSDemag_Demag[idx]->Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx]->InverseFFT(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->energy, true, false, &pSDemagCUDA_Demag[idx]->Module_Heff, &pSDemagCUDA_Demag[idx]->Module_energy);
							}
							else {

								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->energy, true, false);
							}
						}
					}

					//build total energy
					energy_cpu += pSDemagCUDA_Demag[idx]->energy.to_cpu() * pSDemag->pSDemag_Demag[idx]->energy_density_weight;
				}

				energy.from_cpu(energy_cpu);
			}

			//no energy contribution needed
			else {

				//Inverse FFT
				for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemagCUDA_Demag[idx]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						if (pSDemagCUDA_Demag[idx]->do_transfer) {

							//do inverse FFT
							if (pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx]->InverseFFT(
									pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, false, true, &pSDemagCUDA_Demag[idx]->transfer_Module_Heff, &pSDemagCUDA_Demag[idx]->transfer_Module_energy);

								pSDemagCUDA_Demag[idx]->transfer_Module_Heff()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.size_transfer_out());
								pSDemagCUDA_Demag[idx]->transfer_Module_energy()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_energy.size_transfer_out());
							}
							else {

								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, false, true);
							}

							//transfer to Heff in each mesh
							pSDemagCUDA_Demag[idx]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
						}
						else {

							if (pSDemag->pSDemag_Demag[idx]->Module_Heff.linear_size()) {

								//do inverse FFT. Add to Heff in each mesh
								pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2,
									pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx]->energy, false, false, &pSDemagCUDA_Demag[idx]->Module_Heff, &pSDemagCUDA_Demag[idx]->Module_energy);
							}
							else {

								//do inverse FFT. Add to Heff in each mesh
								pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2,
									pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx]->energy, false, false);
							}
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {

						if (pSDemagCUDA_Demag[idx]->do_transfer) {

							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx]->InverseFFT(
									pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, false, true, &pSDemagCUDA_Demag[idx]->transfer_Module_Heff, &pSDemagCUDA_Demag[idx]->transfer_Module_energy);

								pSDemagCUDA_Demag[idx]->transfer_Module_Heff()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_Heff.size_transfer_out());
								pSDemagCUDA_Demag[idx]->transfer_Module_energy()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer_Module_energy.size_transfer_out());
							}
							else {

								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->energy, false, true);
							}

							//transfer to Heff in each mesh
							pSDemagCUDA_Demag[idx]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
						}
						else {

							//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

							//do inverse FFT and accumulate energy. Add to Heff in each mesh
							if (pSDemag->pSDemag_Demag[idx]->Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx]->InverseFFT(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->energy, false, false, &pSDemagCUDA_Demag[idx]->Module_Heff, &pSDemagCUDA_Demag[idx]->Module_energy);
							}
							else {

								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->energy, false, false);
							}
						}
					}
				}
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////// EVAL SPEEDUP - MULTILAYERED ////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		else {
		
			//update if required by ODE solver or if we don't have enough previous evaluations saved to extrapolate
			if (pSMesh->Check_Step_Update() || num_Hdemag_saved < pSMesh->GetEvaluationSpeedup()) {
				
				double energy_cpu = 0.0;

				for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

					cu_obj<cuVEC<cuReal3>>* pHdemag;

					if (num_Hdemag_saved < pSMesh->GetEvaluationSpeedup()) {

						//don't have enough evaluations, so save next one
						switch (num_Hdemag_saved)
						{
						case 0:
							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
							if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag1 = pSMesh->Get_EvalStep_Time();
							break;
						case 1:
							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag2;
							if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag2 = pSMesh->Get_EvalStep_Time();
							break;
						case 2:
							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag3;
							if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag3 = pSMesh->Get_EvalStep_Time();
							break;
						case 3:
							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag4;
							if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag4 = pSMesh->Get_EvalStep_Time();
							break;
						case 4:
							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag5;
							if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag5 = pSMesh->Get_EvalStep_Time();
							break;
						case 5:
							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag6;
							if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag6 = pSMesh->Get_EvalStep_Time();
							break;
						}

						if (idx_mesh == pSDemagCUDA_Demag.size() - 1) num_Hdemag_saved++;
					}
					else {

						//have enough evaluations saved, so just cycle between them now

						//QUINTIC
						if (pSMesh->GetEvaluationSpeedup() == 6) {

							//1, 2, 3, 4, 5, 6 -> next is 1
							if (time_demag6 > time_demag5 && time_demag5 > time_demag4 && time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag1 = pSMesh->Get_EvalStep_Time();
							}
							//2, 3, 4, 5, 6, 1 -> next is 2
							else if (time_demag1 > time_demag2) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag2;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag2 = pSMesh->Get_EvalStep_Time();
							}
							//3, 4, 5, 6, 1, 2 -> next is 3
							else if (time_demag2 > time_demag3) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag3;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag3 = pSMesh->Get_EvalStep_Time();
							}
							//4, 5, 6, 1, 2, 3 -> next is 4
							else if (time_demag3 > time_demag4) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag4;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag4 = pSMesh->Get_EvalStep_Time();
							}
							//5, 6, 1, 2, 3, 4 -> next is 5
							else if (time_demag4 > time_demag5) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag5;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag5 = pSMesh->Get_EvalStep_Time();
							}
							else {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag6;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag6 = pSMesh->Get_EvalStep_Time();
							}
						}
						//QUARTIC
						else if (pSMesh->GetEvaluationSpeedup() == 5) {

							//1, 2, 3, 4, 5 -> next is 1
							if (time_demag5 > time_demag4 && time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag1 = pSMesh->Get_EvalStep_Time();
							}
							//2, 3, 4, 5, 1 -> next is 2
							else if (time_demag1 > time_demag2) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag2;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag2 = pSMesh->Get_EvalStep_Time();
							}
							//3, 4, 5, 1, 2 -> next is 3
							else if (time_demag2 > time_demag3) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag3;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag3 = pSMesh->Get_EvalStep_Time();
							}
							//4, 5, 1, 2, 3 -> next is 4
							else if (time_demag3 > time_demag4) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag4;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag4 = pSMesh->Get_EvalStep_Time();
							}
							else {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag5;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag5 = pSMesh->Get_EvalStep_Time();
							}
						}
						//CUBIC
						else if (pSMesh->GetEvaluationSpeedup() == 4) {

							//1, 2, 3, 4 -> next is 1
							if (time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag1 = pSMesh->Get_EvalStep_Time();
							}
							//2, 3, 4, 1 -> next is 2
							else if (time_demag1 > time_demag2) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag2;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag2 = pSMesh->Get_EvalStep_Time();
							}
							//3, 4, 1, 2 -> next is 3
							else if (time_demag2 > time_demag3) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag3;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag3 = pSMesh->Get_EvalStep_Time();
							}
							else {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag4;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag4 = pSMesh->Get_EvalStep_Time();
							}
						}
						//QUADRATIC
						else if (pSMesh->GetEvaluationSpeedup() == 3) {

							//1, 2, 3 -> next is 1
							if (time_demag3 > time_demag2 && time_demag2 > time_demag1) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag1 = pSMesh->Get_EvalStep_Time();
							}
							//2, 3, 1 -> next is 2
							else if (time_demag3 > time_demag2 && time_demag1 > time_demag2) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag2;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag2 = pSMesh->Get_EvalStep_Time();
							}
							//3, 1, 2 -> next is 3, leading to 1, 2, 3 again
							else {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag3;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag3 = pSMesh->Get_EvalStep_Time();
							}
						}
						//LINEAR
						else if (pSMesh->GetEvaluationSpeedup() == 2) {

							//1, 2 -> next is 1
							if (time_demag2 > time_demag1) {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag1 = pSMesh->Get_EvalStep_Time();
							}
							//2, 1 -> next is 2, leading to 1, 2 again
							else {

								pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag2;
								if (idx_mesh == pSDemagCUDA_Demag.size() - 1) time_demag2 = pSMesh->Get_EvalStep_Time();
							}
						}
						//STEP
						else {

							pHdemag = &pSDemagCUDA_Demag[idx_mesh]->Hdemag;
						}
					}
					
					//Inverse FFT
					
					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {
						
						if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx_mesh]->transfer_Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT(
									pSDemagCUDA_Demag[idx_mesh]->transfer, *pHdemag, pSDemagCUDA_Demag[idx_mesh]->energy, true, true, &pSDemagCUDA_Demag[idx_mesh]->transfer_Module_Heff, &pSDemagCUDA_Demag[idx_mesh]->transfer_Module_energy);

								pSDemagCUDA_Demag[idx_mesh]->transfer_Module_Heff()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer_Module_Heff.size_transfer_out());
								pSDemagCUDA_Demag[idx_mesh]->transfer_Module_energy()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer_Module_energy.size_transfer_out());
							}
							else {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT(pSDemagCUDA_Demag[idx_mesh]->transfer, *pHdemag, pSDemagCUDA_Demag[idx_mesh]->energy, true, true);
							}

							//transfer to Heff in each mesh.
							(*pHdemag)()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							
							//remove self demag contribution
							SDemag_EvalSpeedup_SubSelf(pSDemag->n_common.dim(), *pHdemag, pSDemagCUDA_Demag[idx_mesh]->transfer, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
						}
						else {

							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx_mesh]->Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT_AveragedInputs(
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2,
									*pHdemag,
									pSDemagCUDA_Demag[idx_mesh]->energy, true, true,
									&pSDemagCUDA_Demag[idx_mesh]->Module_Heff, &pSDemagCUDA_Demag[idx_mesh]->Module_energy);
							}
							else {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT_AveragedInputs(
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2,
									*pHdemag,
									pSDemagCUDA_Demag[idx_mesh]->energy, true, true);
							}

							//add contribution to Heff and Heff2 then remove self demag contribution
							SDemag_EvalSpeedup_AddField_SubSelf(pSDemag->n_common.dim(),
								pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2, 
								*pHdemag, 
								pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {
						
						if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {
							
							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx_mesh]->transfer_Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT(
									pSDemagCUDA_Demag[idx_mesh]->transfer, *pHdemag, pSDemagCUDA_Demag[idx_mesh]->energy, true, true, &pSDemagCUDA_Demag[idx_mesh]->transfer_Module_Heff, &pSDemagCUDA_Demag[idx_mesh]->transfer_Module_energy);

								pSDemagCUDA_Demag[idx_mesh]->transfer_Module_Heff()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer_Module_Heff.size_transfer_out());
								pSDemagCUDA_Demag[idx_mesh]->transfer_Module_energy()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer_Module_energy.size_transfer_out());
							}
							else {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT(pSDemagCUDA_Demag[idx_mesh]->transfer, *pHdemag, pSDemagCUDA_Demag[idx_mesh]->energy, true, true);
							}

							//transfer to Heff in each mesh.
							(*pHdemag)()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());

							//remove self demag contribution
							SDemag_EvalSpeedup_SubSelf(pSDemag->n_common.dim(), *pHdemag, pSDemagCUDA_Demag[idx_mesh]->transfer, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
						}
						else {
							
							//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

							//do inverse FFT and accumulate energy
							if (pSDemag->pSDemag_Demag[idx_mesh]->Module_Heff.linear_size()) {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT(
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, *pHdemag, pSDemagCUDA_Demag[idx_mesh]->energy, true, true,
									&pSDemagCUDA_Demag[idx_mesh]->Module_Heff, &pSDemagCUDA_Demag[idx_mesh]->Module_energy);
							}
							else {

								pSDemagCUDA_Demag[idx_mesh]->InverseFFT(pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, *pHdemag, pSDemagCUDA_Demag[idx_mesh]->energy, true, true);
							}
							
							//add contribution to Heff then remove self demag contribution
							SDemag_EvalSpeedup_AddField_SubSelf(pSDemag->n_common.dim(),
								pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
								*pHdemag,
								pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
						}
					}
					
					//build total energy
					energy_cpu += pSDemagCUDA_Demag[idx_mesh]->energy.to_cpu() * pSDemag->pSDemag_Demag[idx_mesh]->energy_density_weight;
				}

				energy.from_cpu(energy_cpu);
			}

			else {
			
				//not required to update, and we have enough previous evaluations: use previous Hdemag saves to extrapolate for current evaluation
				
				double a1 = 1.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
				double time = pSMesh->Get_EvalStep_Time();

				//QUINTIC
				if (pSMesh->GetEvaluationSpeedup() == 6) {

					a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5) * (time_demag1 - time_demag6));
					a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5) * (time_demag2 - time_demag6));
					a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5) * (time_demag3 - time_demag6));
					a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) * (time - time_demag6) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5) * (time_demag4 - time_demag6));
					a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag6) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4) * (time_demag5 - time_demag6));
					a6 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag6 - time_demag1) * (time_demag6 - time_demag2) * (time_demag6 - time_demag3) * (time_demag6 - time_demag4) * (time_demag6 - time_demag5));

					for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////
						
						if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5, pSDemagCUDA_Demag[idx_mesh]->Hdemag6,
									a1, a2, a3, a4, a5, a6,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5, pSDemagCUDA_Demag[idx_mesh]->Hdemag6,
									a1, a2, a3, a4, a5, a6,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
						
						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5, pSDemagCUDA_Demag[idx_mesh]->Hdemag6,
									a1, a2, a3, a4, a5, a6,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5, pSDemagCUDA_Demag[idx_mesh]->Hdemag6,
									a1, a2, a3, a4, a5, a6,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
					}
				}
				//QUARTIC
				else if (pSMesh->GetEvaluationSpeedup() == 5) {

					a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5));
					a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5));
					a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5));
					a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5));
					a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4));

					for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5,
									a1, a2, a3, a4, a5,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5,
									a1, a2, a3, a4, a5,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5,
									a1, a2, a3, a4, a5,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4, pSDemagCUDA_Demag[idx_mesh]->Hdemag5,
									a1, a2, a3, a4, a5,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
					}
				}
				//CUBIC
				else if (pSMesh->GetEvaluationSpeedup() == 4) {

					a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4));
					a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4));
					a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4));
					a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3));

					for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4,
									a1, a2, a3, a4,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4,
									a1, a2, a3, a4,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4,
									a1, a2, a3, a4,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3, pSDemagCUDA_Demag[idx_mesh]->Hdemag4,
									a1, a2, a3, a4,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
					}
				}
				//QUADRATIC
				else if (pSMesh->GetEvaluationSpeedup() == 3) {

					if (time_demag2 != time_demag1 && time_demag2 != time_demag3 && time_demag1 != time_demag3) {

						a1 = (time - time_demag2) * (time - time_demag3) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3));
						a2 = (time - time_demag1) * (time - time_demag3) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3));
						a3 = (time - time_demag1) * (time - time_demag2) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2));
					}

					for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3,
									a1, a2, a3,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3,
									a1, a2, a3,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3,
									a1, a2, a3,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2, pSDemagCUDA_Demag[idx_mesh]->Hdemag3,
									a1, a2, a3,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
					}
				}
				//LINEAR
				else if (pSMesh->GetEvaluationSpeedup() == 2) {

					if (time_demag2 != time_demag1) {

						a1 = (time - time_demag2) / (time_demag1 - time_demag2);
						a2 = (time - time_demag1) / (time_demag2 - time_demag1);
					}

					for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2,
									a1, a2,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2,
									a1, a2,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
						
						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2,
									a1, a2,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag, pSDemagCUDA_Demag[idx_mesh]->Hdemag2,
									a1, a2,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
					}
				}
				//STEP
				else {

					for (int idx_mesh = 0; idx_mesh < pSDemagCUDA_Demag.size(); idx_mesh++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx_mesh]->pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff2,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M2, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
						
						///////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////// OTHER MAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx_mesh]->do_transfer) {

								//construct demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_SetExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->transfer,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag,
									pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx_mesh]->transfer()->transfer_out(pSDemag->pSDemag_Demag[idx_mesh]->transfer.size_transfer_out());
							}
							else {

								//transfer is forced for atomistic meshes, so if no transfer required, this must mean a micromagnetic mesh

								//add demag field approximation, including self demag contribution
								SDemag_EvalSpeedup_AddExtrapField_AddSelf(pSDemag->n_common.dim(),
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->Heff,
									pSDemagCUDA_Demag[idx_mesh]->Hdemag,
									pSDemagCUDA_Demag[idx_mesh]->pMeshCUDA->M, pSDemagCUDA_Demag[idx_mesh]->selfDemagCoeff);
							}
						}
					}
				}
			}
		}
	}
}

#endif

#endif