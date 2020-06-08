#include "stdafx.h"
#include "ExchangeCUDA.h"
#include "Exchange.h"

#ifdef MODULE_EXCHANGE

#include "MeshCUDA.h"

#if COMPILECUDA == 1

Exch_6ngbr_NeuCUDA::Exch_6ngbr_NeuCUDA(MeshCUDA* pMeshCUDA_, Exch_6ngbr_Neu* pExch_6ngbr_Neu) : 
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, pExch_6ngbr_Neu)
{
	pMeshCUDA = pMeshCUDA_;
}

Exch_6ngbr_NeuCUDA::~Exch_6ngbr_NeuCUDA()
{}

BError Exch_6ngbr_NeuCUDA::Initialize(void)
{
	BError error(CLASS_STR(Exch_6ngbr_NeuCUDA));

	error = ExchangeBaseCUDA::Initialize();

	if (!error) initialized = true;

	return error;
}

BError Exch_6ngbr_NeuCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Exch_6ngbr_NeuCUDA));

	Uninitialize();

	return error;
}

#endif

#endif

