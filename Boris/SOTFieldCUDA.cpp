#include "stdafx.h"
#include "SOTFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SOTFIELD

#include "SOTField.h"
#include "MeshCUDA.h"
#include "MeshDefs.h"

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

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs((cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, (MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_SOTFIELD, (MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_SOTFIELD, pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError SOTFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SOTFieldCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
