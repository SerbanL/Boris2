#include "stdafx.h"
#include "DMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_DMEXCHANGE

#include "Mesh_FerromagneticCUDA.h"

DMExchangeCUDA::DMExchangeCUDA(FMeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

DMExchangeCUDA::~DMExchangeCUDA()
{}

BError DMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(DMExchangeCUDA));

	initialized = true;

	return error;
}

BError DMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DMExchangeCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif

