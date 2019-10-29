#include "stdafx.h"
#include "DMExchangeCUDA.h"
#include "DMExchange.h"

#if COMPILECUDA == 1

#ifdef MODULE_DMEXCHANGE

#include "MeshCUDA.h"

DMExchangeCUDA::DMExchangeCUDA(MeshCUDA* pMeshCUDA_, DMExchange* pDMExchange)
	: 
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, pDMExchange)
{
	pMeshCUDA = pMeshCUDA_;
}

DMExchangeCUDA::~DMExchangeCUDA()
{}

BError DMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(DMExchangeCUDA));

	error = ExchangeBaseCUDA::Initialize();

	if (!error) initialized = true;

	return error;
}

BError DMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DMExchangeCUDA));

	Uninitialize();

	return error;
}

#endif

#endif

