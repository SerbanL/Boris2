#include "stdafx.h"
#include "STFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STFIELD

#include "STField.h"
#include "MeshCUDA.h"
#include "MeshDefs.h"

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

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs((cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, (MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_STFIELD, (MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_STFIELD, pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError STFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SOTFieldCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
