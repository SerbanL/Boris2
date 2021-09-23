#include "stdafx.h"
#include "Atom_DipoleDipoleCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "SimScheduleDefs.h"

#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "Atom_DipoleDipole.h"
#include "DataDefs.h"

Atom_DipoleDipoleCUDA::Atom_DipoleDipoleCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_DipoleDipole *paDipoleDipole_) :
	ModulesCUDA(),
	ConvolutionCUDA<Atom_DipoleDipoleCUDA, DipoleDipoleKernelCUDA>(paMeshCUDA_->n_dm, paMeshCUDA_->h_dm)
{
	Uninitialize();

	paMeshCUDA = paMeshCUDA_;

	paDipoleDipole = paDipoleDipole_;

	error_on_create = Convolution_Error_on_Create();
}

Atom_DipoleDipoleCUDA::~Atom_DipoleDipoleCUDA() {}

BError Atom_DipoleDipoleCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_DipoleDipoleCUDA));

	if (!initialized) {

		//initialize mesh transfer on CPU
		error = paDipoleDipole->Initialize_Mesh_Transfer();
		using_macrocell = paDipoleDipole->using_macrocell;

		if (using_macrocell) {

			//M used only if using macrocell

			if (!M()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Hd()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

			cu_arr<cuVEC<cuReal3>> pVal_from_M, pVal_to_M;
			pVal_from_M.push_back((cuVEC<cuReal3>*&)paMeshCUDA->M1.get_managed_object());
			if (!M()->copy_transfer_info(pVal_from_M, pVal_to_M, paDipoleDipole->M)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

			cu_arr<cuVEC<cuReal3>> pVal_from_H, pVal_to_H;
			pVal_to_H.push_back((cuVEC<cuReal3>*&)paMeshCUDA->Heff1.get_managed_object());
			if (!Hd()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hd)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else {

			M()->clear();
			Hd()->clear();
		}

		if (paMeshCUDA->GetEvaluationSpeedup()) {
			
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

			cu_arr<cuVEC<cuReal3>> pVal_from_H, pVal_to_H;
			pVal_to_H.push_back((cuVEC<cuReal3>*&)paMeshCUDA->Heff1.get_managed_object());

			if (paMeshCUDA->GetEvaluationSpeedup() >= 1) if (!Hdemag()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hdemag)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 2) if (!Hdemag2()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hdemag2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 3) if (!Hdemag3()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hdemag3)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 4) if (!Hdemag4()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hdemag4)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 5) if (!Hdemag5()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hdemag5)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (paMeshCUDA->GetEvaluationSpeedup() >= 6) if (!Hdemag6()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hdemag6)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else {

			Hdemag()->clear();
			Hdemag2()->clear();
			Hdemag3()->clear();
			Hdemag4()->clear();
			Hdemag5()->clear();
			Hdemag6()->clear();
		}

		error = Calculate_DipoleDipole_Kernels(using_macrocell);

		if (using_macrocell) selfDemagCoeff.from_cpu(paDipoleDipole->selfDemagCoeff);

		if (!error) initialized = true;
	}

	num_Hdemag_saved = 0;

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(using_macrocell ? (cuReal3)paMeshCUDA->h_dm : (cuReal3)paMeshCUDA->h), (cuRect)paMeshCUDA->meshRect,
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ATOM_DIPOLEDIPOLE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG) || paMeshCUDA->IsStageSet(SS_MONTECARLO),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ATOM_DIPOLEDIPOLE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG) || paMeshCUDA->IsStageSet(SS_MONTECARLO));
	if (error)	initialized = false;

	if (initialized) set_Atom_DipoleDipoleCUDA_pointers();

	//if a Monte Carlo stage is set then we need to compute fields
	if (paMeshCUDA->IsStageSet(SS_MONTECARLO)) paMeshCUDA->Set_Force_MonteCarlo_ComputeFields(true);

	return error;
}

BError Atom_DipoleDipoleCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DipoleDipoleCUDA));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, paDipoleDipole->Get_PBC()) || cfgMessage == UPDATECONFIG_MESHCHANGE || cfgMessage == UPDATECONFIG_DEMAG_CONVCHANGE) {

		Uninitialize();

		//Set convolution dimensions and required PBC conditions
		error = SetDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, true, paDipoleDipole->Get_PBC());

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

void Atom_DipoleDipoleCUDA::UpdateField(void)
{
	//transfer magnetic moments to macrocell mesh
	if (using_macrocell) M()->transfer_in(paDipoleDipole->M.linear_size(), paDipoleDipole->M.size_transfer_in());

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!paMeshCUDA->GetEvaluationSpeedup() || (num_Hdemag_saved < paMeshCUDA->GetEvaluationSpeedup() && !paMeshCUDA->Check_Step_Update())) {

		//don't use evaluation speedup

		if (using_macrocell) {

			//convolute; energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
			if (Module_Heff_size) Convolute(M, Hd, energy, false, true, &Module_Heff, &Module_energy);
			else Convolute(M, Hd, energy, false, true);

			//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
		}
		else {

			if (paMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

			//not using macrocell so get moments directly from mesh

			//convolute and get energy density value
			if (Module_Heff_size) Convolute(paMeshCUDA->M1, paMeshCUDA->Heff1, energy, paMeshCUDA->CurrentTimeStepSolved(), false, &Module_Heff, &Module_energy);
			else Convolute(paMeshCUDA->M1, paMeshCUDA->Heff1, energy, paMeshCUDA->CurrentTimeStepSolved(), false);

			if (paMeshCUDA->CurrentTimeStepSolved()) Energy_to_EnergyDensity(paMeshCUDA->Heff1);
		}
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

			if (using_macrocell) {

				//convolute; energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
				if (Module_Heff_size) Convolute(M, *pHdemag, energy, false, true, &Module_Heff, &Module_energy);
				else Convolute(M, *pHdemag, energy, false, true);
			}
			else {

				//not using macrocell so get moments directly from mesh

				//convolute and get energy density value
				if (Module_Heff_size) Convolute(paMeshCUDA->M1, *pHdemag, energy, true, true, &Module_Heff, &Module_energy);
				else Convolute(paMeshCUDA->M1, *pHdemag, energy, true, true);

				//finish off energy value
				Energy_to_EnergyDensity(*pHdemag);
			}

			//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			(*pHdemag)()->transfer_out(paDipoleDipole->Hdemag.size_transfer_out());

			//subtract self contribution from *pHDemag
			if (using_macrocell) Atom_DipoleDipole_EvalSpeedup_SubSelf(*pHdemag);
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
				if (using_macrocell) {

					Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3, a4, a5, a6);

					//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
					Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
				}
				else {

					Atom_DipoleDipole_EvalSpeedup_AddExtrapField(paMeshCUDA->Heff1, a1, a2, a3, a4, a5, a6);
				}
			}
			//QUARTIC
			else if (paMeshCUDA->GetEvaluationSpeedup() == 5) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5));
				a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4));

				//construct effective field approximation
				if (using_macrocell) {

					Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3, a4, a5);

					//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
					Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
				}
				else {

					Atom_DipoleDipole_EvalSpeedup_AddExtrapField(paMeshCUDA->Heff1, a1, a2, a3, a4, a5);
				}
			}
			//CUBIC
			else if (paMeshCUDA->GetEvaluationSpeedup() == 4) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3));

				//construct effective field approximation
				if (using_macrocell) {

					Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3, a4);

					//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
					Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
				}
				else {

					Atom_DipoleDipole_EvalSpeedup_AddExtrapField(paMeshCUDA->Heff1, a1, a2, a3, a4);
				}
			}
			//QUADRATIC
			else if (paMeshCUDA->GetEvaluationSpeedup() == 3) {

				a1 = (time - time_demag2) * (time - time_demag3) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3));
				a2 = (time - time_demag1) * (time - time_demag3) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3));
				a3 = (time - time_demag1) * (time - time_demag2) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2));

				//construct effective field approximation
				if (using_macrocell) {

					Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2, a3);

					//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
					Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
				}
				else {

					Atom_DipoleDipole_EvalSpeedup_AddExtrapField(paMeshCUDA->Heff1, a1, a2, a3);
				}
			}
			//LINEAR
			else if (paMeshCUDA->GetEvaluationSpeedup() == 2) {

				a1 = (time - time_demag2) / (time_demag1 - time_demag2);
				a2 = (time - time_demag1) / (time_demag2 - time_demag1);

				//construct effective field approximation
				if (using_macrocell) {

					Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(Hd, a1, a2);

					//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
					Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
				}
				else {

					Atom_DipoleDipole_EvalSpeedup_AddExtrapField(paMeshCUDA->Heff1, a1, a2);
				}
			}
			//STEP
			else {

				//construct effective field approximation
				if (using_macrocell) {

					Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(Hd);

					//transfer field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
					Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
				}
				else {

					Atom_DipoleDipole_EvalSpeedup_AddExtrapField(paMeshCUDA->Heff1);
				}
			}
		}
	}
}

#endif

#endif
