#include "stdafx.h"
#include "SurfExchangeCUDA_AFM.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "SurfExchange_AFM.h"
#include "Mesh_AntiFerromagnetic.h"
#include "Mesh_AntiFerromagneticCUDA.h"

SurfExchangeCUDA_AFM::SurfExchangeCUDA_AFM(MeshCUDA* pMeshCUDA_, SurfExchange_AFM* pSurfExch_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
	pSurfExch = pSurfExch_;
}

SurfExchangeCUDA_AFM::~SurfExchangeCUDA_AFM()
{}

BError SurfExchangeCUDA_AFM::Initialize(void)
{
	BError error(CLASS_STR(SurfExchangeCUDA_AFM));

	//clear cu_arrs then rebuild them from information in SurfExchange module
	pMeshFM_Bot.clear();
	pMeshFM_Top.clear();
	pMeshAFM_Bot.clear();
	pMeshAFM_Top.clear();
	pMeshDia_Bot.clear();
	pMeshDia_Top.clear();

	//make sure information in SurfExchange module is up to date
	error = pSurfExch->Initialize();

	if (!error) {

		for (int idx = 0; idx < pSurfExch->pMesh_Bot.size(); idx++) {

			if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_DIAMAGNETIC) {

				pMeshDia_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

				pMeshFM_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				pMeshAFM_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		for (int idx = 0; idx < pSurfExch->pMesh_Top.size(); idx++) {

			if (pSurfExch->pMesh_Top[idx]->GetMeshType() == MESH_DIAMAGNETIC) {

				pMeshDia_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Top[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

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

	return error;
}

BError SurfExchangeCUDA_AFM::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SurfExchangeCUDA_AFM));

	Uninitialize();

	return error;
}

#endif

#endif

