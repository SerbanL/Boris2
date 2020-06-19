#include "stdafx.h"
#include "Atom_DemagCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_DEMAG) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "Atom_Demag.h"

Atom_DemagCUDA::Atom_DemagCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_Demag *paDemag_) :
	ModulesCUDA(),
	ConvolutionCUDA<DemagKernelCUDA>(paMeshCUDA_->n_dm, paMeshCUDA_->h_dm)
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

		error = Calculate_Demag_Kernels();

		if (!M()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!Hdemag()->resize(paMeshCUDA->h_dm, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		//initialize mesh transfer on CPU
		error = paDemag->Initialize_Mesh_Transfer();

		if (!error) {

			cu_arr<cuVEC<cuReal3>> pVal_from_M, pVal_to_M;
			cu_arr<cuVEC<cuReal3>> pVal_from_H, pVal_to_H;

			pVal_from_M.push_back((cuVEC<cuReal3>*&)paMeshCUDA->M1.get_managed_object());
			pVal_to_H.push_back((cuVEC<cuReal3>*&)paMeshCUDA->Heff1.get_managed_object());

			//initialize mesh transfer for M and Hdemag in GPU memory as well
			if (!M()->copy_transfer_info(pVal_from_M, pVal_to_M, paDemag->M)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Hdemag()->copy_transfer_info(pVal_from_H, pVal_to_H, paDemag->Hdemag)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}

		if (!error) initialized = true;
	}

	Hdemag_calculated = false;

	return error;
}

BError Atom_DemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DemagCUDA));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, paDemag->Get_PBC()) || cfgMessage == UPDATECONFIG_MESHCHANGE) {

		Uninitialize();

		//Set convolution dimensions and required PBC conditions
		error = SetDimensions(paMeshCUDA->n_dm, paMeshCUDA->h_dm, true, paDemag->Get_PBC());

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hdemag()->clear();
		M()->clear();
	}

	Hdemag_calculated = false;

	return error;
}

void Atom_DemagCUDA::UpdateField(void)
{
	if (paMeshCUDA->GetEvaluationSpeedup()) {

		//use evaluation speedup method (calculate Hdemag only when required, otherwise just update effective field with previously calculated Hdemag)

		int update_type = paMeshCUDA->Check_Step_Update();

		if (update_type != EVALSPEEDUPSTEP_SKIP || !Hdemag_calculated) {

			//calculate field required

			//transfer magnetic moments to magnetisation mesh, converting from moment to magnetisation in the process
			M()->transfer_in(paDemag->M.linear_size(), paDemag->M.size_transfer_in());

			//convolute and get energy density value
			ZeroEnergy();
			Convolute(M, Hdemag, energy, paMeshCUDA->CurrentTimeStepSolved(), true);

			Hdemag_calculated = true;
		}

		//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hdemag()->transfer_out(paDemag->Hdemag.size_transfer_out());
	}
	else {

		if (paMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

		//transfer magnetic moments to magnetisation mesh, converting from moment to magnetisation in the process
		M()->transfer_in(paDemag->M.linear_size(), paDemag->M.size_transfer_in());

		//convolute and get energy density value
		Convolute(M, Hdemag, energy, paMeshCUDA->CurrentTimeStepSolved(), true);

		//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hdemag()->transfer_out(paDemag->Hdemag.size_transfer_out());
	}
}

#endif

#endif
