#include "stdafx.h"
#include "RoughnessCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ROUGHNESS

#include "MeshCUDA.h"
#include "Roughness.h"
#include "MeshDefs.h"
#include "DataDefs.h"

RoughnessCUDA::RoughnessCUDA(MeshCUDA* pMeshCUDA_, Roughness* pRough_) :
	ModulesCUDA()
{
	Uninitialize();

	pMeshCUDA = pMeshCUDA_;
	pRough = pRough_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

RoughnessCUDA::~RoughnessCUDA() {}

BError RoughnessCUDA::Initialize(void)
{
	BError error(CLASS_STR(RoughnessCUDA));

	if (!pRough->initialized) {

		//Calculate Fmul_rough and Fomul_rough on cpu...
		pRough->Initialize();
	}

	if (!initialized) {

		//...then transfer to GPU
		Fmul_rough()->set_from_cpuvec(pRough->Fmul_rough);
		Fomul_rough()->set_from_cpuvec(pRough->Fomul_rough);
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ROUGHNESS || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ROUGH),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ROUGHNESS || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ROUGH),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	if (initialized) set_RoughnessCUDA_pointers();

	return error;
}

BError RoughnessCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(RoughnessCUDA));

	Uninitialize();

	return error;
}

#endif

#endif