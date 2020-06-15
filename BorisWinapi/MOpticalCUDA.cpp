#include "stdafx.h"
#include "MOpticalCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_MOPTICAL

#include "MOptical.h"
#include "MeshCUDA.h"

MOpticalCUDA::MOpticalCUDA(MeshCUDA* pMeshCUDA_) :
	ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

MOpticalCUDA::~MOpticalCUDA()
{}

BError MOpticalCUDA::Initialize(void)
{
	BError error(CLASS_STR(MOpticalCUDA));

	initialized = true;

	return error;
}

BError MOpticalCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MOpticalCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();
		Initialize();
	}

	return error;
}

void MOpticalCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
}

#endif

#endif