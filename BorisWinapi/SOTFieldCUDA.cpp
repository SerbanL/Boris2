#include "stdafx.h"
#include "SOTFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SOTFIELD

#include "SOTField.h"
#include "MeshCUDA.h"

SOTFieldCUDA::SOTFieldCUDA(MeshCUDA* pMeshCUDA_, SOTField* pHolderModule)
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

	//no energy density contribution here
	ZeroEnergy();

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
