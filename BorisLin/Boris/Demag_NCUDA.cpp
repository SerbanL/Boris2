#include "stdafx.h"
#include "Demag_NCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DEMAG_N

#include "MeshCUDA.h"
#include "DataDefs.h"

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

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_DEMAG_N || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_DEMAG_N || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG));
	if (!error)	initialized = true;

	return error;
}

BError Demag_NCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag_NCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
