#include "stdafx.h"
#include "DemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DEMAG

#include "MeshCUDA.h"
#include "Mesh.h"
#include "Demag.h"

DemagCUDA::DemagCUDA(MeshCUDA* pMeshCUDA_, Demag *pDemag_) :
	ModulesCUDA(), 
	ConvolutionCUDA<DemagCUDA, DemagKernelCUDA>(pMeshCUDA_->n, pMeshCUDA_->h)
{
	Uninitialize();

	pMeshCUDA = pMeshCUDA_;

	pDemag = pDemag_;

	error_on_create = Convolution_Error_on_Create();
}

DemagCUDA::~DemagCUDA() {}

BError DemagCUDA::Initialize(void)
{
	BError error(CLASS_STR(DemagCUDA));

	if (!initialized) {

		error = Calculate_Demag_Kernels();

		if (!error) initialized = true;
	}

	//make sure to allocate memory for Hdemag if we need it
	if (pMeshCUDA->GetEvaluationSpeedup()) {

		Hdemag()->resize(pMeshCUDA->h, pMeshCUDA->meshRect);
	}
	else {

		Hdemag()->clear();
	}

	Hdemag_calculated = false;

	return error;
}

BError DemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DemagCUDA));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(pMeshCUDA->n, pMeshCUDA->h, pDemag->Get_PBC())) {

		Uninitialize();

		//Set convolution dimensions and required PBC conditions
		error = SetDimensions(pMeshCUDA->n, pMeshCUDA->h, true, pDemag->Get_PBC());
	}

	//if memory needs to be allocated for Hdemag, it will be done through Initialize 
	Hdemag()->clear();
	Hdemag_calculated = false;

	return error;
}

void DemagCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetEvaluationSpeedup()) {

		//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

		int update_type = pMeshCUDA->Check_Step_Update();

		if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

			//calculate field required

			if (update_type == EVALSPEEDUPSTEP_COMPUTE_AND_SAVE) {

				//calculate field and save it for next time : we'll need to use it (expecting update_type = EVALSPEEDUPSTEP_SKIP next time)

				//convolute and get energy value
				ZeroEnergy();

				if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

					Convolute(pMeshCUDA->M, Hdemag, energy, true, true);
				}
				else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					Convolute_AveragedInputs(pMeshCUDA->M, pMeshCUDA->M2, Hdemag, energy, true, true);
				}

				Hdemag_calculated = true;
			}
			else {

				//calculate field but do not save it for next time : we'll need to recalculate it again (expecting update_type != EVALSPEEDUPSTEP_SKIP again next time : EVALSPEEDUPSTEP_COMPUTE_NO_SAVE or EVALSPEEDUPSTEP_COMPUTE_AND_SAVE)

				//convolute and get energy value
				ZeroEnergy();

				if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

					Convolute(pMeshCUDA->M, pMeshCUDA->Heff, energy, true, false);
				}
				else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					Convolute_AveragedInputs_DuplicatedOutputs(pMeshCUDA->M, pMeshCUDA->M2, pMeshCUDA->Heff, pMeshCUDA->Heff2, energy, true, true);
				}

				//good practice to set this to false
				Hdemag_calculated = false;

				//return here to avoid adding Hdemag to Heff : we've already added demag field contribution
				return;
			}
		}

		//add contribution to Heff
		pMeshCUDA->Heff()->add_values(pMeshCUDA->n.dim(), Hdemag);

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) pMeshCUDA->Heff2()->add_values(pMeshCUDA->n.dim(), Hdemag);
	}
	else {

		if (pMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

		//Convolute_AveragedInputs_DuplicatedOutputs

		if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

			Convolute(pMeshCUDA->M, pMeshCUDA->Heff, energy, pMeshCUDA->CurrentTimeStepSolved(), false);
		}

		else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			Convolute_AveragedInputs_DuplicatedOutputs(pMeshCUDA->M, pMeshCUDA->M2, pMeshCUDA->Heff, pMeshCUDA->Heff2, energy, pMeshCUDA->CurrentTimeStepSolved(), false);
		}
	}
}

#endif

#endif