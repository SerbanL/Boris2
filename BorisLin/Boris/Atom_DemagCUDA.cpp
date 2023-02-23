#include "stdafx.h"
#include "Atom_DemagCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_DEMAG) && ATOMISTIC == 1

#include "SimScheduleDefs.h"

#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "Atom_Demag.h"
#include "DataDefs.h"
#include "SuperMesh.h"

Atom_DemagCUDA::Atom_DemagCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_Demag *paDemag_) :
	ModulesCUDA(),
	ConvolutionCUDA<Atom_DemagCUDA, DemagKernelCUDA>(paMeshCUDA_->n_dm, paMeshCUDA_->h_dm)
{
	Uninitialize();

	paMeshCUDA = paMeshCUDA_;

	paDemag = paDemag_;

	error_on_create = Convolution_Error_on_Create();
}

Atom_DemagCUDA::~Atom_DemagCUDA() {}

BError Atom_DemagCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_DemagCUDA));

	if (!initialized) {

		error = Calculate_Demag_Kernels(true, paDemag->paMesh->pSMesh->Get_Kernel_Initialize_on_GPU());

		selfDemagCoeff.from_cpu(DemagTFunc().SelfDemag_PBC(paMeshCUDA->h_dm, paMeshCUDA->n_dm, paDemag->Get_PBC()));
		
		if (!M()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!Hd()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		//make sure to allocate memory for Hdemag if we need it
		if (paMeshCUDA->GetEvaluationSpeedup() >= 6) { if (!Hdemag6()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
		else Hdemag6()->clear();

		if (paMeshCUDA->GetEvaluationSpeedup() >= 5) { if (!Hdemag5()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
		else Hdemag5()->clear();

		if (paMeshCUDA->GetEvaluationSpeedup() >= 4) { if (!Hdemag4()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
		else Hdemag4()->clear();

		if (paMeshCUDA->GetEvaluationSpeedup() >= 3) { if (!Hdemag3()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
		else Hdemag3()->clear();

		if (paMeshCUDA->GetEvaluationSpeedup() >= 2) { if (!Hdemag2()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
		else Hdemag2()->clear();

		if (paMeshCUDA->GetEvaluationSpeedup() >= 1) { if (!Hdemag()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT); }
		else Hdemag()->clear();

		//initialize mesh transfer on CPU
		error = paDemag->Initialize_Mesh_Transfer();

		if (!error) {

			cu_arr<cuVEC<cuReal3>> pVal_from_M, pVal_to_M;
			cu_arr<cuVEC<cuReal3>> pVal_from_H, pVal_to_H;

			pVal_from_M.push_back((cuVEC<cuReal3>*&)paMeshCUDA->M1.get_managed_object());
			pVal_to_H.push_back((cuVEC<cuReal3>*&)paMeshCUDA->Heff1.get_managed_object());

			//initialize mesh transfer for M and Hd in GPU memory as well
			if (!M()->copy_transfer_info(pVal_from_M, pVal_to_M, paDemag->M.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Hd()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hd.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

			if (paMeshCUDA->GetEvaluationSpeedup() >= 1) if (!Hdemag()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 2) if (!Hdemag2()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag2.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 3) if (!Hdemag3()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag3.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 4) if (!Hdemag4()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag4.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 5) if (!Hdemag5()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag5.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 6) if (!Hdemag6()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag6.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}

		if (!error) initialized = true;
	}

	num_Hdemag_saved = 0;
	
	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h_dm, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_DEMAG || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG) || paMeshCUDA->IsStageSet(SS_MONTECARLO),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_DEMAG || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG) || paMeshCUDA->IsStageSet(SS_MONTECARLO));
	if (error)	initialized = false;

	if (initialized) set_Atom_DemagCUDA_pointers();

	//if a Monte Carlo stage is set then we need to compute fields
	if (paMeshCUDA->IsStageSet(SS_MONTECARLO)) paMeshCUDA->Set_Force_MonteCarlo_ComputeFields(true);

	return error;
}

BError Atom_DemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DemagCUDA));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, paDemag->Get_PBC()) || cfgMessage == UPDATECONFIG_MESHCHANGE || cfgMessage == UPDATECONFIG_DEMAG_CONVCHANGE) {

		Uninitialize();

		//Set convolution dimensions and required PBC conditions
		error = SetDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, true, paDemag->Get_PBC());

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hd()->clear();
		M()->clear();

		Hdemag()->clear();
		Hdemag2()->clear();
		Hdemag3()->clear();
		Hdemag4()->clear();
		Hdemag5()->clear();
		Hdemag6()->clear();
	}

	num_Hdemag_saved = 0;

	return error;
}

void Atom_DemagCUDA::UpdateField(void)
{
	//transfer magnetic moments to magnetization mesh, converting from moment to magnetization in the process
	M()->transfer_in(paDemag->M.linear_size(), paDemag->M.size_transfer_in());

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!paMeshCUDA->GetEvaluationSpeedup() || (num_Hdemag_saved < paMeshCUDA->GetEvaluationSpeedup() && !paMeshCUDA->Check_Step_Update())) {

		if (paMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

		//convolute and get energy density value
		if (Module_Heff_size) Convolute(M, Hd, energy, paMeshCUDA->CurrentTimeStepSolved(), true, &Module_Heff, &Module_energy);
		else Convolute(M, Hd, energy, paMeshCUDA->CurrentTimeStepSolved(), true);

		//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hd()->transfer_out(paDemag->Hd.size_transfer_out());
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

		//update if required by ODE solver or if we don't have enough previous evaluations saved to extrapolate
		if (paMeshCUDA->Check_Step_Update() || num_Hdemag_saved < paMeshCUDA->GetEvaluationSpeedup()) {

			cu_obj<cuVEC<cuReal3>>* pHdemag;

			if (num_Hdemag_saved < paMeshCUDA->GetEvaluationSpeedup()) {

				//don't have enough evaluations, so save next one
				switch (num_Hdemag_saved)
				{
				case 0:
					pHdemag = &Hdemag;
					time_demag1 = paMeshCUDA->Get_EvalStep_Time();
					break;
				case 1:
					pHdemag = &Hdemag2;
					time_demag2 = paMeshCUDA->Get_EvalStep_Time();
					break;
				case 2:
					pHdemag = &Hdemag3;
					time_demag3 = paMeshCUDA->Get_EvalStep_Time();
					break;
				case 3:
					pHdemag = &Hdemag4;
					time_demag4 = paMeshCUDA->Get_EvalStep_Time();
					break;
				case 4:
					pHdemag = &Hdemag5;
					time_demag5 = paMeshCUDA->Get_EvalStep_Time();
					break;
				case 5:
					pHdemag = &Hdemag6;
					time_demag6 = paMeshCUDA->Get_EvalStep_Time();
					break;
				}

				num_Hdemag_saved++;
			}
			else {

				//have enough evaluations saved, so just cycle between them now

				//QUINTIC
				if (paMeshCUDA->GetEvaluationSpeedup() == 6) {

					//1, 2, 3, 4, 5, 6 -> next is 1
					if (time_demag6 > time_demag5 && time_demag5 > time_demag4 && time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMeshCUDA->Get_EvalStep_Time();
					}
					//2, 3, 4, 5, 6, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = paMeshCUDA->Get_EvalStep_Time();
					}
					//3, 4, 5, 6, 1, 2 -> next is 3
					else if (time_demag2 > time_demag3) {

						pHdemag = &Hdemag3;
						time_demag3 = paMeshCUDA->Get_EvalStep_Time();
					}
					//4, 5, 6, 1, 2, 3 -> next is 4
					else if (time_demag3 > time_demag4) {

						pHdemag = &Hdemag4;
						time_demag4 = paMeshCUDA->Get_EvalStep_Time();
					}
					//5, 6, 1, 2, 3, 4 -> next is 5
					else if (time_demag4 > time_demag5) {

						pHdemag = &Hdemag5;
						time_demag5 = paMeshCUDA->Get_EvalStep_Time();
					}
					else {

						pHdemag = &Hdemag6;
						time_demag6 = paMeshCUDA->Get_EvalStep_Time();
					}
				}
				//QUARTIC
				else if (paMeshCUDA->GetEvaluationSpeedup() == 5) {

					//1, 2, 3, 4, 5 -> next is 1
					if (time_demag5 > time_demag4 && time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMeshCUDA->Get_EvalStep_Time();
					}
					//2, 3, 4, 5, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = paMeshCUDA->Get_EvalStep_Time();
					}
					//3, 4, 5, 1, 2 -> next is 3
					else if (time_demag2 > time_demag3) {

						pHdemag = &Hdemag3;
						time_demag3 = paMeshCUDA->Get_EvalStep_Time();
					}
					//4, 5, 1, 2, 3 -> next is 4
					else if (time_demag3 > time_demag4) {

						pHdemag = &Hdemag4;
						time_demag4 = paMeshCUDA->Get_EvalStep_Time();
					}
					else {

						pHdemag = &Hdemag5;
						time_demag5 = paMeshCUDA->Get_EvalStep_Time();
					}
				}
				//CUBIC
				else if (paMeshCUDA->GetEvaluationSpeedup() == 4) {

					//1, 2, 3, 4 -> next is 1
					if (time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMeshCUDA->Get_EvalStep_Time();
					}
					//2, 3, 4, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = paMeshCUDA->Get_EvalStep_Time();
					}
					//3, 4, 1, 2 -> next is 3
					else if (time_demag2 > time_demag3) {

						pHdemag = &Hdemag3;
						time_demag3 = paMeshCUDA->Get_EvalStep_Time();
					}
					else {

						pHdemag = &Hdemag4;
						time_demag4 = paMeshCUDA->Get_EvalStep_Time();
					}
				}
				//QUADRATIC
				else if (paMeshCUDA->GetEvaluationSpeedup() == 3) {

					//1, 2, 3 -> next is 1
					if (time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMeshCUDA->Get_EvalStep_Time();
					}
					//2, 3, 1 -> next is 2
					else if (time_demag3 > time_demag2 && time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = paMeshCUDA->Get_EvalStep_Time();
					}
					//3, 1, 2 -> next is 3, leading to 1, 2, 3 again
					else {

						pHdemag = &Hdemag3;
						time_demag3 = paMeshCUDA->Get_EvalStep_Time();
					}
				}
				//LINEAR
				else if (paMeshCUDA->GetEvaluationSpeedup() == 2) {

					//1, 2 -> next is 1
					if (time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMeshCUDA->Get_EvalStep_Time();
					}
					//2, 1 -> next is 2, leading to 1, 2 again
					else {

						pHdemag = &Hdemag2;
						time_demag2 = paMeshCUDA->Get_EvalStep_Time();
					}
				}
				//STEP
				else {

					pHdemag = &Hdemag;
				}
			}

			//do evaluation
			ZeroEnergy();
			
			if (Module_Heff_size) Convolute(M, *pHdemag, energy, true, true, &Module_Heff, &Module_energy);
			else Convolute(M, *pHdemag, energy, true, true);

			//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			(*pHdemag)()->transfer_out(paDemag->Hd.size_transfer_out());

			//subtract self demag from *pHDemag
			Atom_Demag_EvalSpeedup_SubSelf(*pHdemag);
		}
		else {

			//not required to update, and we have enough previous evaluations: use previous Hdemag saves to extrapolate for current evaluation

			cuBReal a1 = 1.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
			cuBReal time = paMeshCUDA->Get_EvalStep_Time();

			//QUINTIC
			if (paMeshCUDA->GetEvaluationSpeedup() == 6) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5) * (time_demag1 - time_demag6));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5) * (time_demag2 - time_demag6));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5) * (time_demag3 - time_demag6));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) * (time - time_demag6) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5) * (time_demag4 - time_demag6));
				a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag6) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4) * (time_demag5 - time_demag6));
				a6 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag6 - time_demag1) * (time_demag6 - time_demag2) * (time_demag6 - time_demag3) * (time_demag6 - time_demag4) * (time_demag6 - time_demag5));

				//construct effective field approximation
				Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3, a4, a5, a6);

				//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
				Hd()->transfer_out(paDemag->Hd.size_transfer_out());
			}
			//QUARTIC
			else if (paMeshCUDA->GetEvaluationSpeedup() == 5) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5));
				a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4));

				//construct effective field approximation
				Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3, a4, a5);

				//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
				Hd()->transfer_out(paDemag->Hd.size_transfer_out());
			}
			//CUBIC
			else if (paMeshCUDA->GetEvaluationSpeedup() == 4) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3));

				//construct effective field approximation
				Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3, a4);

				//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
				Hd()->transfer_out(paDemag->Hd.size_transfer_out());
			}
			//QUADRATIC
			else if (paMeshCUDA->GetEvaluationSpeedup() == 3) {

				if (time_demag2 != time_demag1 && time_demag2 != time_demag3 && time_demag1 != time_demag3) {

					a1 = (time - time_demag2) * (time - time_demag3) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3));
					a2 = (time - time_demag1) * (time - time_demag3) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3));
					a3 = (time - time_demag1) * (time - time_demag2) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2));
				}

				//construct effective field approximation
				Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3);

				//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
				Hd()->transfer_out(paDemag->Hd.size_transfer_out());
			}
			//LINEAR
			else if (paMeshCUDA->GetEvaluationSpeedup() == 2) {

				if (time_demag2 != time_demag1) {

					a1 = (time - time_demag2) / (time_demag1 - time_demag2);
					a2 = (time - time_demag1) / (time_demag2 - time_demag1);
				}

				//construct effective field approximation
				Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2);

				//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
				Hd()->transfer_out(paDemag->Hd.size_transfer_out());
			}
			//STEP
			else {

				//construct effective field approximation
				Atom_Demag_EvalSpeedup_SetExtrapField_AddSelf(Hd);

				//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
				Hd()->transfer_out(paDemag->Hd.size_transfer_out());
			}
		}
	}
}

#endif

#endif
