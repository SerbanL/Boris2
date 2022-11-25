#include "stdafx.h"
#include "SDemagCUDA_Demag.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "SDemag.h"
#include "SDemag_Demag.h"
#include "MeshDefs.h"
#include "SuperMesh.h"

SDemagCUDA_Demag::SDemagCUDA_Demag(MeshBaseCUDA* pMeshBaseCUDA_, SDemag_Demag *pSDemag_Demag_) :
	ModulesCUDA(),
	ConvolutionCUDA<SDemagCUDA_Demag, DemagKernelCollectionCUDA>()
{
	Uninitialize();

	pMeshBaseCUDA = pMeshBaseCUDA_;

	if (!pMeshBaseCUDA->is_atomistic()) pMeshCUDA = dynamic_cast<MeshCUDA*>(pMeshBaseCUDA);
	else paMeshCUDA = dynamic_cast<Atom_MeshCUDA*>(pMeshBaseCUDA);

	pSDemag_Demag = pSDemag_Demag_;
}

SDemagCUDA_Demag::~SDemagCUDA_Demag() 
{
}

BError SDemagCUDA_Demag::Initialize(void)
{
	BError error(CLASS_STR(SDemagCUDA_Demag));

	//no energy density contribution here
	ZeroEnergy();

	//pointer to cpu SDemag object
	SDemag* pSDemag = pSDemag_Demag->pSDemag;

	//pointer to gpu SDemagCUDA object
	SDemagCUDA* pSDemagCUDA = dynamic_cast<SDemagCUDA*>(pSDemag->pModuleCUDA);

	if (!pSDemagCUDA->IsInitialized()) error = pSDemagCUDA->Initialize();
	if (error) return error;

	//convolution rectangle for this module
	cuRect convolution_rect = (cuRect)pSDemag->get_convolution_rect(pSDemag_Demag);

	//common discretisation cellsize (may differ in thickness in 2D mode)
	cuReal3 h_common = convolution_rect / (cuSZ3)pSDemag->n_common;

	if (!initialized) {
		
		//calculate kernels for super-mesh convolution.
		if (!error) error = Calculate_Demag_Kernels(pSDemagCUDA->kernel_collection, pSDemag_Demag->pSDemag->pSMesh->Get_Kernel_Initialize_on_GPU());
		if (error) return error(BERROR_OUTOFMEMORY_CRIT);

		//setup mesh transfer?

		selfDemagCoeff.from_cpu(DemagTFunc().SelfDemag_PBC(h_common, (cuSZ3)pSDemag->n_common, pSDemag->demag_pbc_images));

		//make sure to allocate memory for Hdemag if we need it
		if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 1) Hdemag()->resize(h_common, convolution_rect);
		else Hdemag()->clear();

		if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 2) Hdemag2()->resize(h_common, convolution_rect);
		else Hdemag2()->clear();

		if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 3) Hdemag3()->resize(h_common, convolution_rect);
		else Hdemag3()->clear();

		if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 4) Hdemag4()->resize(h_common, convolution_rect);
		else Hdemag4()->clear();

		if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 5) Hdemag5()->resize(h_common, convolution_rect);
		else Hdemag5()->clear();

		if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 6) Hdemag6()->resize(h_common, convolution_rect);
		else Hdemag6()->clear();

		int non_empty_cells = 0;

		if (!pMeshBaseCUDA->is_atomistic() && convolution_rect == (cuRect)pSDemag_Demag->meshRect && h_common == (cuReal3)pSDemag_Demag->h) {

			//no transfer required (always force transfer for atomistic meshes)
			do_transfer = false;
			transfer()->clear();
			non_empty_cells = pMeshCUDA->M()->get_nonempty_cells_cpu();
		}
		else {

			do_transfer = true;

			//set correct size for transfer cuVEC
			if (!transfer()->resize(h_common, convolution_rect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

			//now calculate in cpu memory
			pSDemag_Demag->Initialize_Mesh_Transfer();

			//Now copy mesh transfer object to cuda version
			cu_arr<cuVEC<cuReal3>> pVal_from, pVal_from2;
			cu_arr<cuVEC<cuReal3>> pVal_to, pVal_to2;

			cu_arr<cuVEC<cuReal3>> pVal_afrom, pVal_ato;

			if (!pMeshBaseCUDA->is_atomistic()) {

				pVal_from.push_back((cuVEC<cuReal3>*&)pMeshCUDA->M.get_managed_object());
				pVal_to.push_back((cuVEC<cuReal3>*&)pMeshCUDA->Heff.get_managed_object());
				pVal_from2.push_back((cuVEC<cuReal3>*&)pMeshCUDA->M2.get_managed_object());
				pVal_to2.push_back((cuVEC<cuReal3>*&)pMeshCUDA->Heff2.get_managed_object());
			}
			else {

				pVal_afrom.push_back((cuVEC<cuReal3>*&)paMeshCUDA->M1.get_managed_object());
				pVal_ato.push_back((cuVEC<cuReal3>*&)paMeshCUDA->Heff1.get_managed_object());
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			if (pMeshBaseCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

				if (!transfer()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 1) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 2) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag2()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 3) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag3()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 4) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag4()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 5) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag5()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 6) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag6()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			else if (pMeshBaseCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				if (!transfer()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 1) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 2) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag2()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 3) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag3()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 4) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag4()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 5) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag5()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 6) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag6()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////// ATOMISTIC MESH ///////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			else {

				if (!transfer()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 1) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 2) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag2()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 3) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag3()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 4) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag4()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 5) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag5()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}

				if (pMeshBaseCUDA->GetEvaluationSpeedup() >= 6) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag6()->copy_transfer_info(pVal_afrom, pVal_ato, pSDemag_Demag->transfer.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}

			non_empty_cells = pSDemag_Demag->transfer.get_nonempty_cells();
		}
		
		//setup energy density weights
		if (pSDemagCUDA->total_nonempty_volume) {

			energy_density_weight.from_cpu((cuBReal)(non_empty_cells * h_common.dim() / pSDemagCUDA->total_nonempty_volume));
		}
		
		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required : mirror cpu module version
	error = pSDemag_Demag->Initialize_Module_Display();
	if (error) { initialized = false; return error(BERROR_OUTOFMEMORY_CRIT); }

	if (pSDemag_Demag->Get_Module_Heff().linear_size()) {

		error = Update_Module_Display_VECs((cuReal3)pMeshBaseCUDA->h, (cuRect)pMeshBaseCUDA->meshRect, true, true);
		if (error) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		set_SDemag_DemagCUDA_pointers();
	}

	//allocate memory and copy transfer info for module display transfer objects if needed
	if (pSDemag_Demag->transfer_Module_Heff.linear_size()) {

		if (!transfer_Module_Heff()->resize(h_common, convolution_rect)) { initialized = false; return error(BERROR_OUTOFGPUMEMORY_CRIT); }

		cu_arr<cuVEC<cuReal3>> pVal_from;
		cu_arr<cuVEC<cuReal3>> pVal_to;

		if (!pMeshBaseCUDA->is_atomistic()) {

			pVal_to.push_back((cuVEC<cuReal3>*&)(*(pSDemag_Demag->pMesh))(MOD_SDEMAG_DEMAG)->Get_Module_HeffCUDA().get_managed_object());
		}
		else {

			pVal_to.push_back((cuVEC<cuReal3>*&)(*(pSDemag_Demag->paMesh))(MOD_SDEMAG_DEMAG)->Get_Module_HeffCUDA().get_managed_object());
		}

		if (!transfer_Module_Heff()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer_Module_Heff.get_transfer())) { initialized = false; return error(BERROR_OUTOFGPUMEMORY_CRIT); }
	}

	if (pSDemag_Demag->transfer_Module_energy.linear_size()) {

		if (!transfer_Module_energy()->resize(h_common, convolution_rect)) { initialized = false; return error(BERROR_OUTOFGPUMEMORY_CRIT); }

		cu_arr<cuVEC<cuBReal>> pVal_from;
		cu_arr<cuVEC<cuBReal>> pVal_to;

		if (!pMeshBaseCUDA->is_atomistic()) {

			pVal_to.push_back((cuVEC<cuBReal>*&)(*(pSDemag_Demag->pMesh))(MOD_SDEMAG_DEMAG)->Get_Module_EnergyCUDA().get_managed_object());
		}
		else {

			pVal_to.push_back((cuVEC<cuBReal>*&)(*(pSDemag_Demag->paMesh))(MOD_SDEMAG_DEMAG)->Get_Module_EnergyCUDA().get_managed_object());
		}

		if (!transfer_Module_energy()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer_Module_energy.get_transfer())) { initialized = false; return error(BERROR_OUTOFGPUMEMORY_CRIT); }
	}

	return error;
}

BError SDemagCUDA_Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemagCUDA_Demag));

	//just mirror the initialization flag in the cpu version module
	if (!pSDemag_Demag->IsInitialized()) Uninitialize();

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_DEMAG_CONVCHANGE)) Uninitialize();

	//if memory needs to be allocated for Hdemag, it will be done through Initialize 
	if (!initialized) {

		Hdemag()->clear();
		Hdemag2()->clear();
		Hdemag3()->clear();
		Hdemag4()->clear();
		Hdemag5()->clear();
		Hdemag6()->clear();
	}

	return error;
}

void SDemagCUDA_Demag::UpdateField(void)
{
	//demag field update done through the supermesh module.
	//here we need to zero the module display objects in case they are used : if we have to transfer data into them from the display transfer objects, this is done by adding.
	//cannot set output when transferring since we can have multiple transfer objects contributing to the display objects
	ZeroModuleVECs();

	//same goes for the total energy
	ZeroEnergy();
}

#endif

#endif