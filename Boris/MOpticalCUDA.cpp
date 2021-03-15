#include "stdafx.h"
#include "MOpticalCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MOPTICAL

#include "MOptical.h"
#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

MOpticalCUDA::MOpticalCUDA(MeshCUDA* pMeshCUDA_) :
	ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

MOpticalCUDA::~MOpticalCUDA()
{}

BError MOpticalCUDA::Initialize(void)
{
	BError error(CLASS_STR(MOpticalCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_MOPTICAL || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_MOPTICAL),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_MOPTICAL || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_MOPTICAL),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError MOpticalCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MOpticalCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();
	}

	return error;
}

void MOpticalCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
}

#endif

#endif