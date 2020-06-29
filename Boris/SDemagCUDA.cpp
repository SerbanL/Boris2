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

	if (!pSDemag->use_multilayered_convolution) {

		//make sure to allocate memory for Hdemag if we need it
		if (pSMesh->GetEvaluationSpeedup()) Hdemag()->resize(pSMesh->h_fm, pSMesh->sMeshRect_fm);
		else Hdemag()->clear();
	}

	Hdemag_calculated = false;

	//FFT Kernels are not so quick to calculate - if already initialized then we are guaranteed they are correct
	if (!initialized) {
		
		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// SUPERMESH CONVOLUTION ////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (!pSDemag->use_multilayered_convolution) {

			error = Calculate_Demag_Kernels();
			if (error) return error;

			//array of pointers to input meshes (M) and oputput meshes (Heff) to transfer from and to
			cu_arr<cuVEC<cuReal3>> pVal_from, pVal_from2;
			cu_arr<cuVEC<cuReal3>> pVal_to, pVal_to2;

			//identify all existing ferrommagnetic meshes (magnetic computation enabled)
			for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

				if ((*pSMesh)[idx]->MComputation_Enabled() && !(*pSMesh)[idx]->is_atomistic()) {

					pVal_from.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->M.get_managed_object());
					pVal_to.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff.get_managed_object());

					pVal_from2.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->M2.get_managed_object());
					pVal_to2.push_back((cuVEC<cuReal3>*&)dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->Heff2.get_managed_object());
				}
			}

			//Initialize the cpu mesh transfer object - note, SDemag::initialized is set true, but SDemag is not properly initialized - hence need the SDemagCUDA destructor to uninitialize SDemag
			error = pSDemag->Initialize_Mesh_Transfer();
			if (error) return error;

			if (!pSDemag->antiferromagnetic_meshes_present) {

				//Now copy mesh transfer object to cuda version
				if (!sm_Vals()->copy_transfer_info(pVal_from, pVal_to, pSDemag->sm_Vals)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pSMesh->GetEvaluationSpeedup()) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info(pVal_from, pVal_to, pSDemag->sm_Vals)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}
			else {

				//Now copy mesh transfer object to cuda version
				if (!sm_Vals()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag->sm_Vals)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pSMesh->GetEvaluationSpeedup()) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag->sm_Vals)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
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
				error = pSDemagCUDA_Demag[idx]->Set_Rect_Collection(Rect_collection, Rect_collection[idx], h_max);
				if (error) return error;
			}

			//now everything is set correctly, ready to calculate demag kernel collections
		}

		//initialized ok.
		initialized = true;
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

	return error;
}

BError SDemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemagCUDA));
	
	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_DEMAG_CONVCHANGE, UPDATECONFIG_SMESH_CELLSIZE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hdemag()->clear();
		Hdemag_calculated = false;

		if (!pSDemag->use_multilayered_convolution) {

			//only need to uninitialize if n or h have changed
			if (!CheckDimensions(pSMesh->n_fm, pSMesh->h_fm, pSDemag->Get_PBC())) {

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
		////////////////////////////////// EVAL SPEEDUP - SUPERMESH ///////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		if (pSMesh->GetEvaluationSpeedup()) {

			//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

			int update_type = pSMesh->Check_Step_Update();

			if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

				//calculate field required

				if (update_type == EVALSPEEDUPSTEP_COMPUTE_AND_SAVE) {

					//calculate field and save it for next time : we'll need to use it (expecting update_type = EVALSPEEDUPSTEP_SKIP next time)

					if (!pSDemag->antiferromagnetic_meshes_present) {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals()->transfer_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());				
					}
					else {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals()->transfer_in_averaged(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());
					}

					ZeroEnergy();
					Convolute(sm_Vals, Hdemag, energy, true, true);

					Hdemag_calculated = true;
				}
				else {

					//calculate field but do not save it for next time : we'll need to recalculate it again (expecting update_type != EVALSPEEDUPSTEP_SKIP again next time : EVALSPEEDUPSTEP_COMPUTE_NO_SAVE or EVALSPEEDUPSTEP_COMPUTE_AND_SAVE)

					if (!pSDemag->antiferromagnetic_meshes_present) {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals()->transfer_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());

						ZeroEnergy();
						Convolute(sm_Vals, sm_Vals, energy, true, true);

						//good practice to set this to false
						Hdemag_calculated = false;

						//transfer to individual Heff meshes
						sm_Vals()->transfer_out(pSDemag->sm_Vals.size_transfer_out());
					}
					else {

						//transfer values from invidual M meshes to sm_Vals
						sm_Vals()->transfer_in_averaged(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());

						ZeroEnergy();
						Convolute(sm_Vals, sm_Vals, energy, true, true);

						//good practice to set this to false
						Hdemag_calculated = false;

						//transfer to individual Heff meshes
						sm_Vals()->transfer_out_duplicated(pSDemag->sm_Vals.size_transfer_out());
					}

					//return here to avoid adding Hdemag to Heff : we've already added demag field contribution
					return;
					
				}
			}

			if (!pSDemag->antiferromagnetic_meshes_present) {

				//transfer contribution to meshes Heff
				Hdemag()->transfer_out(pSDemag->sm_Vals.size_transfer_out());
			}
			else {

				//transfer contribution to meshes Heff
				Hdemag()->transfer_out_duplicated(pSDemag->sm_Vals.size_transfer_out());
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// NO SPEEDUP - SUPERMESH ///////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		else {

			//don't use evaluation speedup, so no need to use Hdemag (this won't have memory allocated anyway)

			if (!pSDemag->antiferromagnetic_meshes_present) {

				//transfer values from invidual M meshes to sm_Vals
				sm_Vals()->transfer_in(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());

				//only need energy after ode solver step finished
				if (pSMesh->CurrentTimeStepSolved()) ZeroEnergy();

				Convolute(sm_Vals, sm_Vals, energy, pSMesh->CurrentTimeStepSolved(), true);

				//transfer to individual Heff meshes
				sm_Vals()->transfer_out(pSDemag->sm_Vals.size_transfer_out());
			}
			else {

				//transfer values from invidual M meshes to sm_Vals
				sm_Vals()->transfer_in_averaged(pSDemag->sm_Vals.linear_size(), pSDemag->sm_Vals.size_transfer_in());

				//only need energy after ode solver step finished
				if (pSMesh->CurrentTimeStepSolved()) ZeroEnergy();

				Convolute(sm_Vals, sm_Vals, energy, pSMesh->CurrentTimeStepSolved(), true);

				//transfer to individual Heff meshes
				sm_Vals()->transfer_out_duplicated(pSDemag->sm_Vals.size_transfer_out());
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
				for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

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
					////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {

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
				}

				//Kernel multiplications for multiple inputs.
				for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

					pSDemagCUDA_Demag[idx]->KernelMultiplication_MultipleInputs(FFT_Spaces_x_Input, FFT_Spaces_y_Input, FFT_Spaces_z_Input);
				}

				if (update_type == EVALSPEEDUPSTEP_COMPUTE_AND_SAVE) {

					//calculate field and save it for next time : we'll need to use it (expecting update_type = EVALSPEEDUPSTEP_SKIP next time)

					ZeroEnergy();

					//Inverse FFT
					for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx]->do_transfer) {

								//do inverse FFT and accumulate energy
								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->Hdemag, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);
							}
							else {

								//do inverse FFT and accumulate energy.
								pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2, 
									pSDemagCUDA_Demag[idx]->Hdemag, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

							if (pSDemagCUDA_Demag[idx]->do_transfer) {

								//do inverse FFT and accumulate energy
								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->Hdemag, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);
							}
							else {

								//do inverse FFT and accumulate energy.
								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->Hdemag, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);
							}
						}
					}

					Hdemag_calculated = true;
				}
				else {

					//calculate field but do not save it for next time : we'll need to recalculate it again (expecting update_type != EVALSPEEDUPSTEP_SKIP again next time : EVALSPEEDUPSTEP_COMPUTE_NO_SAVE or EVALSPEEDUPSTEP_COMPUTE_AND_SAVE)

					ZeroEnergy();

					//Inverse FFT
					for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

						///////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

							if (pSDemagCUDA_Demag[idx]->do_transfer) {

								//do inverse FFT and accumulate energy
								pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);

								//transfer to Heff in each mesh
								pSDemagCUDA_Demag[idx]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
							}
							else {

								//do inverse FFT and accumulate energy. Add to Heff in each mesh
								pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
									pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2, 
									pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2, 
									energy, pSDemagCUDA_Demag[idx]->energy_density_weight, false);
							}
						}

						///////////////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////

						else {

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
					
					//good practice to set this to false
					Hdemag_calculated = false;

					//return here to avoid adding Hdemag to Heff : we've already added demag field contribution
					return;
				}
			}

			//add contribution to Heff in each mesh
			for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

				///////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					if (pSDemagCUDA_Demag[idx]->do_transfer) {

						//transfer to Heff in each mesh
						pSDemagCUDA_Demag[idx]->Hdemag()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx]->Hdemag.size_transfer_out());
					}
					else {

						//add contribution to Heff
						pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff()->add_values(pSDemagCUDA_Demag[idx]->pMeshCUDA->n.dim(), pSDemagCUDA_Demag[idx]->Hdemag);
						pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2()->add_values(pSDemagCUDA_Demag[idx]->pMeshCUDA->n.dim(), pSDemagCUDA_Demag[idx]->Hdemag);
					}
				}

				///////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				else {

					if (pSDemagCUDA_Demag[idx]->do_transfer) {

						//transfer to Heff in each mesh
						pSDemagCUDA_Demag[idx]->Hdemag()->transfer_out(pSDemag->pSDemag_Demag[idx]->Hdemag.size_transfer_out());
					}
					else {

						//add contribution to Heff
						pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff()->add_values(pSDemagCUDA_Demag[idx]->pMeshCUDA->n.dim(), pSDemagCUDA_Demag[idx]->Hdemag);
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
			for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

				///////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

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
				////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////

				else {

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

					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						if (pSDemagCUDA_Demag[idx]->do_transfer) {

							//do inverse FFT and accumulate energy
							pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, energy, pSDemagCUDA_Demag[idx]->energy_density_weight, true);

							//transfer to Heff in each mesh
							pSDemagCUDA_Demag[idx]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
						}
						else {

							//do inverse FFT and accumulate energy. Add to Heff in each mesh
							pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
								pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2, 
								pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2, 
								energy, pSDemagCUDA_Demag[idx]->energy_density_weight, false);
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////
					////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {

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
			}
			else {

				//no energy contribution needed

				//Inverse FFT
				for (int idx = 0; idx < pSDemagCUDA_Demag.size(); idx++) {

					///////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					if (pSDemagCUDA_Demag[idx]->pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						if (pSDemagCUDA_Demag[idx]->do_transfer) {

							//do inverse FFT
							pSDemagCUDA_Demag[idx]->InverseFFT(pSDemagCUDA_Demag[idx]->transfer, pSDemagCUDA_Demag[idx]->transfer, energy, false, true);

							//transfer to Heff in each mesh
							pSDemagCUDA_Demag[idx]->transfer()->transfer_out_duplicated(pSDemag->pSDemag_Demag[idx]->transfer.size_transfer_out());
						}
						else {

							//do inverse FFT. Add to Heff in each mesh
							pSDemagCUDA_Demag[idx]->InverseFFT_AveragedInputs_DuplicatedOutputs(
								pSDemagCUDA_Demag[idx]->pMeshCUDA->M, pSDemagCUDA_Demag[idx]->pMeshCUDA->M2, 
								pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff, pSDemagCUDA_Demag[idx]->pMeshCUDA->Heff2, 
								energy, false, false);
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////
					////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////

					else {

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
			}
		}	
	}
}

#endif

#endif