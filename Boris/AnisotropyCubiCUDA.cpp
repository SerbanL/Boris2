#include "stdafx.h"
#include "AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANICUBI

#include "MeshCUDA.h"

//--------------- CUBIC

Anisotropy_CubicCUDA::Anisotropy_CubicCUDA(MeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Anisotropy_CubicCUDA::~Anisotropy_CubicCUDA()
{}

BError Anisotropy_CubicCUDA::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_CubicCUDA));

	initialized = true;

	return error;
}

BError Anisotropy_CubicCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_CubicCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif