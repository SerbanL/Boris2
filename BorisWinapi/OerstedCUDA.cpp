#include "stdafx.h"
#include "OerstedCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_OERSTED

#include "SuperMesh.h"
#include "Oersted.h"

OerstedCUDA::OerstedCUDA(SuperMesh* pSMesh_, Oersted* pOersted_) :
	ModulesCUDA(),
	ConvolutionCUDA<OerstedKernelCUDA>(pSMesh_->n_e, pSMesh_->h_e)
{
	Uninitialize();

	pSMesh = pSMesh_;
	pOersted = pOersted_;

	error_on_create = Convolution_Error_on_Create();

	//set from cpu version of sm_vals
	if (!sm_Vals()->set_from_cpuvec(pOersted->sm_Vals)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
}

OerstedCUDA::~OerstedCUDA()
{
	//copy values back to cpu version
	if (Holder_Module_Available()) {

		sm_Vals()->copy_to_cpuvec(pOersted->sm_Vals);
	}
}

BError OerstedCUDA::Initialize(void)
{
	BError error(CLASS_STR(OerstedCUDA));

	//not counting this to the total energy density for now
	ZeroEnergy();

	//FFT Kernels are not so quick to calculate - if already initialized then we are guaranteed they are correct
	if (!initialized) {

		vector< VEC<DBL3>* > pVal_from_cpu;
		vector< VEC<DBL3>* > pVal_to_cpu;

		//array of pointers to input meshes (M) and oputput meshes (Heff) to transfer from and to
		cu_arr<cuVEC<cuReal3>> pVal_from;
		cu_arr<cuVEC<cuReal3>> pVal_to;

		//identify all existing ferrommagnetic meshes (magnetic computation enabled)
		for (int idx = 0; idx < (int)pSMesh->pMesh.size(); idx++) {

			if ((*pSMesh)[idx]->EComputation_Enabled()) {

				pVal_from_cpu.push_back(&((*pSMesh)[idx]->Jc));
				pVal_from.push_back((cuVEC<cuReal3>*&)(*pSMesh)[idx]->pMeshCUDA->Jc.get_managed_object());
			}

			if ((*pSMesh)[idx]->MComputation_Enabled()) {

				pVal_to_cpu.push_back(&((*pSMesh)[idx]->Heff));
				pVal_to.push_back((cuVEC<cuReal3>*&)(*pSMesh)[idx]->pMeshCUDA->Heff.get_managed_object());
			}
		}

		//Initialize the mesh transfer object.
		if (!pOersted->sm_Vals.Initialize_MeshTransfer(pVal_from_cpu, pVal_to_cpu, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		//Now copy mesh transfer object to cuda version
		if (!sm_Vals()->copy_transfer_info(pVal_from, pVal_to, pOersted->sm_Vals)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		error = Calculate_Oersted_Kernels();

		if (!error) initialized = true;
	}

	return error;
}

BError OerstedCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(OerstedCUDA));

	//only need to uninitialize if n or h have changed
	if (!CheckDimensions(pSMesh->n_e, pSMesh->h_e)) {

		Uninitialize();
		error = SetDimensions(pSMesh->n_e, pSMesh->h_e);

		if (!sm_Vals()->resize(pSMesh->h_e, pSMesh->sMeshRect_e)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	return error;
}

void OerstedCUDA::UpdateField(void)
{
	//only recalculate Oersted field if there was a significant change in current density (judged based on transport solver iterations)
	if (pSMesh->CallModuleMethod(&STransport::Transport_Recalculated) || !pOersted->oefield_computed) {

		//transfer values from invidual Jc meshes to sm_Vals
		sm_Vals()->transfer_in(pOersted->sm_Vals.linear_size(), pOersted->sm_Vals.size_transfer_in());

		//only need energy after ode solver step finished
		if (pSMesh->CurrentTimeStepSolved()) ZeroEnergy();

		Convolute(sm_Vals, sm_Vals, energy, false, true);

		//transfer to individual Heff meshes
		sm_Vals()->transfer_out(pOersted->sm_Vals.size_transfer_out());

		pOersted->oefield_computed = true;
	}
	else {

		//transfer to individual Heff meshes
		sm_Vals()->transfer_out(pOersted->sm_Vals.size_transfer_out());
	}
}

#endif

#endif