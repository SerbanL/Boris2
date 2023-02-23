#include "stdafx.h"
#include "SurfExchangeCUDA_AFM.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "SurfExchange_AFM.h"

#include "Mesh_AntiFerromagnetic.h"
#include "Mesh_AntiFerromagneticCUDA.h"

#include "Atom_Mesh.h"
#include "Atom_MeshCUDA.h"

#include "DataDefs.h"

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
	pMeshAtom_Bot.clear();
	pMeshAtom_Top.clear();

	//make sure information in SurfExchange module is up to date
	error = pSurfExch->Initialize();

	if (!error) {

		//micromagnetic meshes BOTTOM
		for (int idx = 0; idx < pSurfExch->pMesh_Bot.size(); idx++) {

			if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

				pMeshFM_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Bot[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				pMeshAFM_Bot.push_back(pSurfExch->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		//atomistic meshes BOTTOM
		for (int idx = 0; idx < pSurfExch->paMesh_Bot.size(); idx++) {

			pMeshAtom_Bot.push_back(pSurfExch->paMesh_Bot[idx]->paMeshCUDA->cuaMesh.get_managed_object());
		}

		//micromagnetic meshes TOP
		for (int idx = 0; idx < pSurfExch->pMesh_Top.size(); idx++) {

			if (pSurfExch->pMesh_Top[idx]->GetMeshType() == MESH_FERROMAGNETIC) {

				pMeshFM_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
			else if (pSurfExch->pMesh_Top[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				pMeshAFM_Top.push_back(pSurfExch->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}
		}

		//atomistic meshes TOP
		for (int idx = 0; idx < pSurfExch->paMesh_Top.size(); idx++) {

			pMeshAtom_Top.push_back(pSurfExch->paMesh_Top[idx]->paMeshCUDA->cuaMesh.get_managed_object());
		}

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect,
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_SURFEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_SURFEXCH) || pMeshCUDA->IsOutputDataSet(DATA_T_SURFEXCH),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_SURFEXCHANGE || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_SURFEXCH));
	if (error) initialized = false;

	if (initialized) set_SurfExchangeCUDA_AFM_pointers();

	return error;
}

BError SurfExchangeCUDA_AFM::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SurfExchangeCUDA_AFM));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 SurfExchangeCUDA_AFM::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif

