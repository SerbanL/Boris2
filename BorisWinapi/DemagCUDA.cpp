#include "stdafx.h"
#include "DemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_DEMAG

#include "Mesh_FerromagneticCUDA.h"
#include "Mesh_Ferromagnetic.h"

DemagCUDA::DemagCUDA(FMeshCUDA* pMeshCUDA_) : 
	ModulesCUDA(), 
	ConvolutionCUDA<DemagKernelCUDA>(pMeshCUDA_->n, pMeshCUDA_->h)
{
	Uninitialize();

	pMeshCUDA = pMeshCUDA_;

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

	return error;
}

BError DemagCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DemagCUDA));

	//only need to uninitialize if n or h have changed
	if (!CheckDimensions(pMeshCUDA->n, pMeshCUDA->h)) {

		Uninitialize();
		error = SetDimensions(pMeshCUDA->n, pMeshCUDA->h);
	}

	return error;
}

void DemagCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) ZeroEnergy();

	Convolute(pMeshCUDA->M, pMeshCUDA->Heff, energy, pMeshCUDA->CurrentTimeStepSolved(), false);
}

#endif

#endif