#include "stdafx.h"
#include "Atom_DipoleDipoleCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "Atom_DipoleDipole.h"

Atom_DipoleDipoleCUDA::Atom_DipoleDipoleCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_DipoleDipole *paDipoleDipole_) :
	ModulesCUDA(),
	ConvolutionCUDA<DipoleDipoleKernelCUDA>(paMeshCUDA_->n_dm, paMeshCUDA_->h_dm)
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
			
			cu_arr<cuVEC<cuReal3>> pVal_from_M, pVal_to_M;
			pVal_from_M.push_back((cuVEC<cuReal3>*&)paMeshCUDA->M1.get_managed_object());
			if (!M()->copy_transfer_info(pVal_from_M, pVal_to_M, paDipoleDipole->M)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else M()->clear();

		if (using_macrocell || paMeshCUDA->GetEvaluationSpeedup()) {

			//Hd used if using macrocell, or if using eval speedup

			if (!Hd()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			
			cu_arr<cuVEC<cuReal3>> pVal_from_H, pVal_to_H;
			pVal_to_H.push_back((cuVEC<cuReal3>*&)paMeshCUDA->Heff1.get_managed_object());
			if (!Hd()->copy_transfer_info(pVal_from_H, pVal_to_H, paDipoleDipole->Hd)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else Hd()->clear();

		error = Calculate_DipoleDipole_Kernels(using_macrocell);

		if (!error) initialized = true;
	}

	Hd_calculated = false;

	return error;
}

BError Atom_DipoleDipoleCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DipoleDipoleCUDA));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, paDipoleDipole->Get_PBC()) || cfgMessage == UPDATECONFIG_MESHCHANGE) {

		Uninitialize();

		//Set convolution dimensions and required PBC conditions
		error = SetDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, true, paDipoleDipole->Get_PBC());

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hd()->clear();
		M()->clear();
	}

	Hd_calculated = false;

	return error;
}

void Atom_DipoleDipoleCUDA::UpdateField(void)
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (paMeshCUDA->GetEvaluationSpeedup()) {

		//use evaluation speedup method (calculate Hdemag only when required, otherwise just update effective field with previously calculated Hdemag)

		int update_type = paMeshCUDA->Check_Step_Update();

		if (update_type != EVALSPEEDUPSTEP_SKIP || !Hd_calculated) {

			//calculate field required

			if (using_macrocell) {

				//transfer magnetic moments to macrocell mesh
				M()->transfer_in(paDipoleDipole->M.linear_size(), paDipoleDipole->M.size_transfer_in());

				//convolute; energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
				Convolute(M, Hd, energy, false, true);
			}
			else {

				//not using macrocell so get moments directly from mesh

				if (paMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

				//convolute and get energy density value
				Convolute(paMeshCUDA->M1, Hd, energy, paMeshCUDA->CurrentTimeStepSolved(), true);

				//finish off energy value
				if (paMeshCUDA->CurrentTimeStepSolved()) Energy_to_EnergyDensity(Hd);
			}

			Hd_calculated = true;
		}

		//transfer dipole-dipole field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//don't use evaluation speedup

		if (using_macrocell) {

			//transfer magnetic moments to magnetisation mesh, converting from moment to magnetisation in the process
			M()->transfer_in(paDipoleDipole->M.linear_size(), paDipoleDipole->M.size_transfer_in());

			//convolute; energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
			Convolute(M, Hd, energy, false , true);

			//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			Hd()->transfer_out(paDipoleDipole->Hd.size_transfer_out());
		}
		else {

			if (paMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

			//not using macrocell so get moments directly from mesh

			//convolute and get energy density value
			Convolute(paMeshCUDA->M1, paMeshCUDA->Heff1, energy, paMeshCUDA->CurrentTimeStepSolved(), false);

			if (paMeshCUDA->CurrentTimeStepSolved()) Energy_to_EnergyDensity(paMeshCUDA->Heff1);
		}
	}
}

#endif

#endif
