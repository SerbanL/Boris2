#include "stdafx.h"
#include "ExchangeCUDA.h"
#include "DataDefs.h"

#ifdef MODULE_COMPILATION_EXCHANGE

#include "MeshCUDA.h"
#include "Exchange.h"
#include "MeshDefs.h"

#if COMPILECUDA == 1

Exch_6ngbr_NeuCUDA::Exch_6ngbr_NeuCUDA(MeshCUDA* pMeshCUDA_, Exch_6ngbr_Neu* pExch_6ngbr_Neu) : 
	ModulesCUDA(),
	ExchangeBaseCUDA(pMeshCUDA_, pExch_6ngbr_Neu)
{
	pMeshCUDA = pMeshCUDA_;
}

Exch_6ngbr_NeuCUDA::~Exch_6ngbr_NeuCUDA()
{}

BError Exch_6ngbr_NeuCUDA::Initialize(void)
{
	BError error(CLASS_STR(Exch_6ngbr_NeuCUDA));

	error = ExchangeBaseCUDA::Initialize();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_EXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_EXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError Exch_6ngbr_NeuCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Exch_6ngbr_NeuCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Exch_6ngbr_NeuCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif

