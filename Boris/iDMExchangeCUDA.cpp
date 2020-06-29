#include "stdafx.h"
#include "iDMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_IDMEXCHANGE

#include "MeshCUDA.h"
#include "iDMExchange.h"

iDMExchangeCUDA::iDMExchangeCUDA(MeshCUDA* pMeshCUDA_, iDMExchange* piDMExchange)
	:
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, piDMExchange)
{
	pMeshCUDA = pMeshCUDA_;
}

iDMExchangeCUDA::~iDMExchangeCUDA()
{}

BError iDMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(iDMExchangeCUDA));

	error = ExchangeBaseCUDA::Initialize();

	if (!error) initialized = true;

	return error;
}

BError iDMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(iDMExchangeCUDA));

	Uninitialize();

	return error;
}

#endif

#endif

