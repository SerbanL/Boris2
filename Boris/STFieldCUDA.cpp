#include "stdafx.h"
#include "STFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STFIELD

#include "STField.h"
#include "MeshCUDA.h"

STFieldCUDA::STFieldCUDA(MeshCUDA* pMeshCUDA_, STField* pHolderModule)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

STFieldCUDA::~STFieldCUDA()
{}

BError STFieldCUDA::Initialize(void)
{
	BError error(CLASS_STR(STFieldCUDA));

	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError STFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SOTFieldCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif
