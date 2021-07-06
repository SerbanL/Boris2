#include "stdafx.h"
#include "SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "SurfExchange.h"
#include "Mesh_Ferromagnetic.h"
#include "Mesh_FerromagneticCUDA.h"
#include "DataDefs.h"

SurfExchangeCUDA::SurfExchangeCUDA(MeshCUDA* pMeshCUDA_, SurfExchange* pSurfExch_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
	pSurfExch = pSurfExch_;
}

SurfExchangeCUDA::~SurfExchangeCUDA()
{}

BError SurfExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(SurfExchangeCUDA));

	//clear cu_arrs then rebuild them from information in SurfExchange module
	pMeshFM_Bot.clear();
	pMeshFM_Top.clear();
	pMeshAFM_Bot.clear();
	pMeshAFM_Top.clear();

	//make sure information in SurfExchange module is up to date
	error = pSurfExch->Initialize();

	if (!error) {

		for (int idx = 0; idx < pSurfExch->pMesh_Bot.size(); idx++) {

			if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

				pMeshFM_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				pMeshAFM_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		for (int idx = 0; idx < pSurfExch->pMesh_Top.size(); idx++) {

			if (pSurfExch->pMesh_Top[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

				pMeshFM_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Top[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				pMeshAFM_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		//copy number of coupled cells to gpu memory : SurfExchange has been initialized so the count is correct
		coupled_cells.from_cpu(pSurfExch->coupled_cells);

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_SURFEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_SURFEXCH) || pMeshCUDA->IsOutputDataSet(DATA_T_SURFEXCH),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_SURFEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_SURFEXCH));
	if (error) initialized = false;

	if (initialized) set_SurfExchangeCUDA_pointers();

	return error;
}

BError SurfExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SurfExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 SurfExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif

