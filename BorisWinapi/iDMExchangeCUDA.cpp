#include "stdafx.h"
#include "iDMExchangeCUDA.h"
#include "iDMExchange.h"

#if COMPILECUDA == 1

#ifdef MODULE_IDMEXCHANGE

#include "Mesh_FerromagneticCUDA.h"

iDMExchangeCUDA::iDMExchangeCUDA(FMeshCUDA* pMeshCUDA_, iDMExchange* piDMExchange)
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

