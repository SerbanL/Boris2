#include "stdafx.h"
#include "SDemagCUDA_Demag.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "MeshCUDA.h"
#include "SDemag.h"
#include "SDemag_Demag.h"
#include "MeshDefs.h"

SDemagCUDA_Demag::SDemagCUDA_Demag(MeshCUDA* pMeshCUDA_, SDemag_Demag *pSDemag_Demag_) :
	ModulesCUDA(),
	ConvolutionCUDA<SDemagCUDA_Demag, DemagKernelCollectionCUDA>()
{
	Uninitialize();

	pMeshCUDA = pMeshCUDA_;

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

	//make sure to allocate memory for Hdemag if we need it
	if (pMeshCUDA->GetEvaluationSpeedup()) Hdemag()->resize((cuRect)pSDemag->get_convolution_rect(pSDemag_Demag) / (cuSZ3)pSDemag->n_common, (cuRect)pSDemag->get_convolution_rect(pSDemag_Demag));
	else Hdemag()->clear();

	if (!initialized) {

		//calculate kernels for super-mesh convolution.
		if (!error) error = Calculate_Demag_Kernels(pSDemagCUDA->kernel_collection);
		if (error) return error(BERROR_OUTOFMEMORY_CRIT);

		//setup mesh transfer?

		//convolution rectangle for this module
		cuRect convolution_rect = (cuRect)pSDemag->get_convolution_rect(pSDemag_Demag);

		//common discretisation cellsize (may differ in thickness in 2D mode)
		cuReal3 h_common = convolution_rect / (cuSZ3)pSDemag->n_common;

		int non_empty_cells = 0;

		if (convolution_rect == (cuRect)pSDemag_Demag->meshRect && h_common == (cuReal3)pSDemag_Demag->h) {

			//no transfer required
			do_transfer = false;
			transfer()->clear();

			non_empty_cells = pMeshCUDA->M()->get_nonempty_cells_cpu();
		}
		else {

			do_transfer = true;

			//set correct size for transfer VEC

			//first calculate in cpu memory
			pSDemag_Demag->Initialize_Mesh_Transfer();

			//set correct size for transfer cuVEC
			if (!transfer()->resize(h_common, convolution_rect)) error_on_create(BERROR_OUTOFMEMORY_CRIT);

			//Now copy mesh transfer object to cuda version
			cu_arr<cuVEC<cuReal3>> pVal_from, pVal_from2;
			cu_arr<cuVEC<cuReal3>> pVal_to, pVal_to2;

			pVal_from.push_back((cuVEC<cuReal3>*&)pMeshCUDA->M.get_managed_object());
			pVal_to.push_back((cuVEC<cuReal3>*&)pMeshCUDA->Heff.get_managed_object());
			pVal_from2.push_back((cuVEC<cuReal3>*&)pMeshCUDA->M2.get_managed_object());
			pVal_to2.push_back((cuVEC<cuReal3>*&)pMeshCUDA->Heff2.get_managed_object());

			///////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

				if (!transfer()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pMeshCUDA->GetEvaluationSpeedup()) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info(pVal_from, pVal_to, pSDemag_Demag->transfer)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

			else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				if (!transfer()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				if (pMeshCUDA->GetEvaluationSpeedup()) {

					//initialize mesh transfer for Hdemag as well if we are using evaluation speedup
					if (!Hdemag()->copy_transfer_info_averagedinputs_duplicatedoutputs(pVal_from, pVal_from2, pVal_to, pVal_to2, pSDemag_Demag->transfer)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
				}
			}

			non_empty_cells = pSDemag_Demag->transfer.get_nonempty_cells();
		}
		
		if (pSDemagCUDA->total_nonempty_volume) {

			energy_density_weight.from_cpu((cuBReal)(non_empty_cells * h_common.dim() / pSDemagCUDA->total_nonempty_volume));
		}
		
		
		initialized = true;
	}

	return error;
}

BError SDemagCUDA_Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SDemagCUDA_Demag));

	//just mirror the initialization flag in the cpu version module
	if (!pSDemag_Demag->IsInitialized()) Uninitialize();

	//if memory needs to be allocated for Hdemag, it will be done through Initialize 
	Hdemag()->clear();

	return error;
}

#endif

#endif