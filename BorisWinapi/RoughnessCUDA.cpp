#include "stdafx.h"
#include "RoughnessCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ROUGHNESS

#include "Mesh_FerromagneticCUDA.h"
#include "Roughness.h"

RoughnessCUDA::RoughnessCUDA(FMeshCUDA* pMeshCUDA_, Roughness* pRough_) :
	ModulesCUDA()
{
	Uninitialize();

	pMeshCUDA = pMeshCUDA_;
	pRough = pRough_;

	error_on_create = UpdateConfiguration();
}

RoughnessCUDA::~RoughnessCUDA() {}

BError RoughnessCUDA::Initialize(void)
{
	BError error(CLASS_STR(RoughnessCUDA));

	if (!pRough->initialized) {

		//Calculate Fmul_rough and Fomul_rough on cpu...
		pRough->Initialize();
	}

	if (!initialized) {

		//...then transfer to GPU
		Fmul_rough()->set_from_cpuvec(pRough->Fmul_rough);
		Fomul_rough()->set_from_cpuvec(pRough->Fomul_rough);
	}

	initialized = true;

	return error;
}

BError RoughnessCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(RoughnessCUDA));

	Uninitialize();

	return error;
}

#endif

#endif