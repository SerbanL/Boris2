#include "stdafx.h"
#include "DMExchangeCUDA.h"
#include "DataDefs.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DMEXCHANGE

#include "MeshCUDA.h"
#include "DMExchange.h"
#include "MeshDefs.h"

DMExchangeCUDA::DMExchangeCUDA(MeshCUDA* pMeshCUDA_, DMExchange* pDMExchange)
	: 
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, pDMExchange)
{
	pMeshCUDA = pMeshCUDA_;
}

DMExchangeCUDA::~DMExchangeCUDA()
{}

BError DMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(DMExchangeCUDA));

	error = ExchangeBaseCUDA::Initialize();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_ActualModule_Heff_Display() == MOD_DMEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)pMeshCUDA->Get_ActualModule_Energy_Display() == MOD_DMEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError DMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DMExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 DMExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif

