#include "stdafx.h"
#include "iDMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_IDMEXCHANGE

#include "Mesh_FerromagneticCUDA.h"

iDMExchangeCUDA::iDMExchangeCUDA(FMeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

iDMExchangeCUDA::~iDMExchangeCUDA()
{}

BError iDMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(iDMExchangeCUDA));

	initialized = true;

	return error;
}

BError iDMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(iDMExchangeCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif

