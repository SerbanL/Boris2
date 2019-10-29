#include "stdafx.h"
#include "Demag_NCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_DEMAG_N

#include "MeshCUDA.h"

Demag_NCUDA::Demag_NCUDA(MeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Demag_NCUDA::~Demag_NCUDA()
{}

BError Demag_NCUDA::Initialize(void)
{
	BError error(CLASS_STR(Demag_NCUDA));

	initialized = true;

	return error;
}

BError Demag_NCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag_NCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif
