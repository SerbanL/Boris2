#include "stdafx.h"
#include "SOTFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SOTFIELD

#include "SOTField.h"
#include "Mesh_FerromagneticCUDA.h"

SOTFieldCUDA::SOTFieldCUDA(FMeshCUDA* pMeshCUDA_, SOTField* pHolderModule)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

SOTFieldCUDA::~SOTFieldCUDA()
{}

BError SOTFieldCUDA::Initialize(void)
{
	BError error(CLASS_STR(SOTFieldCUDA));

	initialized = true;

	return error;
}

BError SOTFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SOTFieldCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif
