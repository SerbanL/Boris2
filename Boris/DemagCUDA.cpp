#include "stdafx.h"
#include "DemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DEMAG

#include "SimScheduleDefs.h"

#include "MeshCUDA.h"
#include "Mesh.h"
#include "Demag.h"
#include "DataDefs.h"
#include "SuperMesh.h"

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

		error = Calculate_Demag_Kernels(true, pDemag->pMesh->pSMesh->Get_Kernel_Initialize_on_GPU());

		selfDemagCoeff.from_cpu(DemagTFunc().SelfDemag_PBC(pMeshCUDA->h, pMeshCUDA->n, pDemag->Get_PBC()));

		if (!error) initialized = true;
	}

	//make sure to allocate memory for Hdemag if we need it
	if (pMeshCUDA->GetEvaluationSpeedup() >= 3) { if (!Hdemag3()->resize(pMeshCUDA->h, pMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
	else Hdemag3()->clear();

	if (pMeshCUDA->GetEvaluationSpeedup() >= 2) { if (!Hdemag2()->resize(pMeshCUDA->h, pMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
	else Hdemag2()->clear();

	if (pMeshCUDA->GetEvaluationSpeedup() >= 1) { if (!Hdemag()->resize(pMeshCUDA->h, pMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
	else Hdemag()->clear();

	num_Hdemag_saved = 0;

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_DEMAG || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG) || pMeshCUDA->IsStageSet(SS_MONTECARLO),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_DEMAG || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG) || pMeshCUDA->IsStageSet(SS_MONTECARLO));
	if (error) initialized = false;

	if (initialized) set_DemagCUDA_pointers();

	//if a Monte Carlo stage is set then we need to compute fields
	if (pMeshCUDA->IsStageSet(SS_MONTECARLO)) pMeshCUDA->Set_Force_MonteCarlo_ComputeFields(true);

	return error;
}

BError DemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DemagCUDA));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(pMeshCUDA->n, pMeshCUDA->h, pDemag->Get_PBC()) || cfgMessage == UPDATECONFIG_DEMAG_CONVCHANGE) {

		Uninitialize();

		//Set convolution dimensions and required PBC conditions
		error = SetDimensions(pMeshCUDA->n, pMeshCUDA->h, true, pDemag->Get_PBC());

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hdemag()->clear();
		Hdemag2()->clear();
		Hdemag3()->clear();
	}

	num_Hdemag_saved = 0;

	return error;
}

void DemagCUDA::UpdateField(void)
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!pMeshCUDA->GetEvaluationSpeedup() || (num_Hdemag_saved < pMeshCUDA->GetEvaluationSpeedup() && !pMeshCUDA->Check_Step_Update())) {

		//don't use evaluation speedup, so no need to use Hdemag (this won't have memory allocated anyway) - or else we are using speedup but don't yet have enough previous evaluations at steps where we should be extrapolating

		if (pMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

		if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {
			
			if (Module_Heff_size) Convolute(pMeshCUDA->M, pMeshCUDA->Heff, energy, pMeshCUDA->CurrentTimeStepSolved(), false, &Module_Heff, &Module_energy);
			else Convolute(pMeshCUDA->M, pMeshCUDA->Heff, energy, pMeshCUDA->CurrentTimeStepSolved(), false);
		}

		else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			if (Module_Heff_size) Convolute_AveragedInputs_DuplicatedOutputs(pMeshCUDA->M, pMeshCUDA->M2, pMeshCUDA->Heff, pMeshCUDA->Heff2, energy, pMeshCUDA->CurrentTimeStepSolved(), false, &Module_Heff, &Module_energy);
			else Convolute_AveragedInputs_DuplicatedOutputs(pMeshCUDA->M, pMeshCUDA->M2, pMeshCUDA->Heff, pMeshCUDA->Heff2, energy, pMeshCUDA->CurrentTimeStepSolved(), false);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

		//update if required by ODE solver or if we don't have enough previous evaluations saved to extrapolate
		if (pMeshCUDA->Check_Step_Update() || num_Hdemag_saved < pMeshCUDA->GetEvaluationSpeedup()) {

			cu_obj<cuVEC<cuReal3>>* pHdemag;

			if (num_Hdemag_saved < pMeshCUDA->GetEvaluationSpeedup()) {

				//don't have enough evaluations, so save next one
				switch (num_Hdemag_saved)
				{
				case 0:
					pHdemag = &Hdemag;
					time_demag1 = pMeshCUDA->Get_EvalStep_Time();
					break;
				case 1:
					pHdemag = &Hdemag2;
					time_demag2 = pMeshCUDA->Get_EvalStep_Time();
					break;
				case 2:
					pHdemag = &Hdemag3;
					time_demag3 = pMeshCUDA->Get_EvalStep_Time();
					break;
				}

				num_Hdemag_saved++;
			}
			else {

				//have enough evaluations saved, so just cycle between them now

				//QUADRATIC
				if (pMeshCUDA->GetEvaluationSpeedup() == 3) {

					//1, 2, 3 -> next is 1
					if (time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMeshCUDA->Get_EvalStep_Time();
					}
					//2, 3, 1 -> next is 2
					else if (time_demag3 > time_demag2 && time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = pMeshCUDA->Get_EvalStep_Time();
					}
					//3, 1, 2 -> next is 3, leading to 1, 2, 3 again
					else {

						pHdemag = &Hdemag3;
						time_demag3 = pMeshCUDA->Get_EvalStep_Time();
					}
				}
				//LINEAR
				else if (pMeshCUDA->GetEvaluationSpeedup() == 2) {

					//1, 2 -> next is 1
					if (time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMeshCUDA->Get_EvalStep_Time();
					}
					//2, 1 -> next is 2, leading to 1, 2 again
					else {

						pHdemag = &Hdemag2;
						time_demag2 = pMeshCUDA->Get_EvalStep_Time();
					}
				}
				//STEP
				else {

					pHdemag = &Hdemag;
				}
			}

			//do evaluation
			ZeroEnergy();

			if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

				if (Module_Heff_size) Convolute(pMeshCUDA->M, *pHdemag, energy, true, true, &Module_Heff, &Module_energy);
				else Convolute(pMeshCUDA->M, *pHdemag, energy, true, true);
			}

			else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				if (Module_Heff_size) Convolute_AveragedInputs(pMeshCUDA->M, pMeshCUDA->M2, *pHdemag, energy, true, true, &Module_Heff, &Module_energy);
				else Convolute_AveragedInputs(pMeshCUDA->M, pMeshCUDA->M2, *pHdemag, energy, true, true);
			}

			if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				//add contribution to Heff and Heff2 and subtract self demag from *pHDemag
				Demag_EvalSpeedup_AddField_SubSelf(pMeshCUDA->Heff, pMeshCUDA->Heff2, *pHdemag, pMeshCUDA->M, pMeshCUDA->M2);
			}
			else {

				//add contribution to Heff and subtract self demag from *pHDemag
				Demag_EvalSpeedup_AddField_SubSelf(pMeshCUDA->Heff, *pHdemag, pMeshCUDA->M);
			}
		}
		else {

			//not required to update, and we have enough previous evaluations: use previous Hdemag saves to extrapolate for current evaluation

			cuBReal a1 = 1.0, a2 = 0.0, a3 = 0.0;
			cuBReal time = pMeshCUDA->Get_EvalStep_Time();

			//QUADRATIC
			if (pMeshCUDA->GetEvaluationSpeedup() == 3) {

				if (time_demag2 != time_demag1 && time_demag2 != time_demag3 && time_demag1 != time_demag3) {

					a1 = (time - time_demag2) * (time - time_demag3) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3));
					a2 = (time - time_demag1) * (time - time_demag3) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3));
					a3 = (time - time_demag1) * (time - time_demag2) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2));
				}

				if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2, together with self demag
					Demag_EvalSpeedup_AddExtrapField_AddSelf(pMeshCUDA->Heff, pMeshCUDA->Heff2, a1, a2, a3, pMeshCUDA->M, pMeshCUDA->M2);

				}
				else {

					//add contribution to Heff, together with self demag
					Demag_EvalSpeedup_AddExtrapField_AddSelf(pMeshCUDA->Heff, a1, a2, a3, pMeshCUDA->M);
				}
			}
			//LINEAR
			else if (pMeshCUDA->GetEvaluationSpeedup() == 2) {

				if (time_demag2 != time_demag1) {

					a1 = (time - time_demag2) / (time_demag1 - time_demag2);
					a2 = (time - time_demag1) / (time_demag2 - time_demag1);
				}

				if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2, together with self demag
					Demag_EvalSpeedup_AddExtrapField_AddSelf(pMeshCUDA->Heff, pMeshCUDA->Heff2, a1, a2, pMeshCUDA->M, pMeshCUDA->M2);
				}
				else {

					//add contribution to Heff, together with self demag
					Demag_EvalSpeedup_AddExtrapField_AddSelf(pMeshCUDA->Heff, a1, a2, pMeshCUDA->M);
				}
			}
			//STEP
			else {

				if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2, together with self demag
					Demag_EvalSpeedup_AddExtrapField_AddSelf(pMeshCUDA->Heff, pMeshCUDA->Heff2, pMeshCUDA->M, pMeshCUDA->M2);
				}
				else {

					//add contribution to Heff, together with self demag
					Demag_EvalSpeedup_AddExtrapField_AddSelf(pMeshCUDA->Heff, pMeshCUDA->M);
				}
			}
		}
	}
}

#endif

#endif