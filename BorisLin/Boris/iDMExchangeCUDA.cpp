#include "stdafx.h"
#include "iDMExchangeCUDA.h"
#include "DataDefs.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_IDMEXCHANGE

#include "MeshCUDA.h"
#include "iDMExchange.h"
#include "MeshDefs.h"

iDMExchangeCUDA::iDMExchangeCUDA(MeshCUDA* pMeshCUDA_, iDMExchange* piDMExchange)
	:
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, piDMExchange)
{
	pMeshCUDA = pMeshCUDA_;
}

iDMExchangeCUDA::~iDMExchangeCUDA()
{}

BError iDMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(iDMExchangeCUDA));

	error = ExchangeBaseCUDA::Initialize();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_ActualModule_Heff_Display() == MOD_IDMEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)pMeshCUDA->Get_ActualModule_Energy_Display() == MOD_IDMEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError iDMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(iDMExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 iDMExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif

