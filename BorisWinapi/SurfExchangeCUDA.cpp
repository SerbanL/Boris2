#include "stdafx.h"
#include "SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SURFEXCHANGE

#include "SurfExchange.h"
#include "Mesh_Ferromagnetic.h"
#include "Mesh_FerromagneticCUDA.h"

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
	pMesh_Bot.clear();
	pMesh_Top.clear();
	pMeshDia_Bot.clear();
	pMeshDia_Top.clear();

	//make sure information in SurfExchange module is up to date
	error = pSurfExch->Initialize();

	if (!error) {

		for (int idx = 0; idx < pSurfExch->pMesh_Bot.size(); idx++) {

			if (pSurfExch->pMesh_Bot[idx]->GetMeshType() != MESH_DIAMAGNETIC) {

				pMesh_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else {

				pMeshDia_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		for (int idx = 0; idx < pSurfExch->pMesh_Top.size(); idx++) {

			if (pSurfExch->pMesh_Top[idx]->GetMeshType() != MESH_DIAMAGNETIC) {

				pMesh_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else {

				pMeshDia_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		//copy number of coupled cells to gpu memory : SurfExchange has been initialized so the count is correct
		coupled_cells.from_cpu(pSurfExch->coupled_cells);

		initialized = true;
	}

	return error;
}

BError SurfExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SurfExchangeCUDA));

	Uninitialize();

	return error;
}

#endif

#endif

