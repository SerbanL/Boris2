#include "stdafx.h"
#include "ExchangeCUDA.h"

#ifdef MODULE_EXCHANGE

#include "Mesh_FerromagneticCUDA.h"

#if COMPILECUDA == 1

Exch_6ngbr_NeuCUDA::Exch_6ngbr_NeuCUDA(FMeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Exch_6ngbr_NeuCUDA::~Exch_6ngbr_NeuCUDA()
{}

BError Exch_6ngbr_NeuCUDA::Initialize(void)
{
	BError error(CLASS_STR(Exch_6ngbr_NeuCUDA));

	initialized = true;

	return error;
}

BError Exch_6ngbr_NeuCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Exch_6ngbr_NeuCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif

