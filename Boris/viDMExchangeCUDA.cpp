#include "stdafx.h"
#include "viDMExchangeCUDA.h"
#include "DataDefs.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_VIDMEXCHANGE

#include "MeshCUDA.h"
#include "viDMExchange.h"
#include "MeshDefs.h"

viDMExchangeCUDA::viDMExchangeCUDA(MeshCUDA* pMeshCUDA_, viDMExchange* pviDMExchange)
	:
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, pviDMExchange)
{
	pMeshCUDA = pMeshCUDA_;
}

viDMExchangeCUDA::~viDMExchangeCUDA()
{}

BError viDMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(viDMExchangeCUDA));

	error = ExchangeBaseCUDA::Initialize();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect,
		(MOD_)pMeshCUDA->Get_ActualModule_Heff_Display() == MOD_VIDMEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)pMeshCUDA->Get_ActualModule_Energy_Display() == MOD_VIDMEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError viDMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(viDMExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 viDMExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif