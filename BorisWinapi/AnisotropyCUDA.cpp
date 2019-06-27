#include "stdafx.h"
#include "AnisotropyCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ANIUNI

#include "Mesh_FerromagneticCUDA.h"

//--------------- UNIAXIAL

Anisotropy_UniaxialCUDA::Anisotropy_UniaxialCUDA(FMeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Anisotropy_UniaxialCUDA::~Anisotropy_UniaxialCUDA()
{}

BError Anisotropy_UniaxialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_UniaxialCUDA));

	initialized = true;

	return error;
}

BError Anisotropy_UniaxialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_UniaxialCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif