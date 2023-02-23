#include "stdafx.h"
#include "STFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STFIELD

#include "STField.h"
#include "MeshDefs.h"

#include "Mesh_Ferromagnetic.h"
#include "Mesh_FerromagneticCUDA.h"

#include "Atom_Mesh.h"
#include "Atom_MeshCUDA.h"

STFieldCUDA::STFieldCUDA(MeshCUDA* pMeshCUDA_, STField* pSTField_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
	pSTField = pSTField_;
}

STFieldCUDA::~STFieldCUDA()
{}

BError STFieldCUDA::Initialize(void)
{
	BError error(CLASS_STR(STFieldCUDA));

	//clear cu_arrs then rebuild them from information in STField module
	pMeshFM_Bot.clear();
	pMeshFM_Top.clear();
	pMeshAtom_Bot.clear();
	pMeshAtom_Top.clear();

	if (pSTField->pMesh->STp.get0() == DBL3()) {

		fixed_polarization = false;

		//make sure information in STField module is up to date
		error = pSTField->Initialize();

		if (!error) {

			for (int idx = 0; idx < pSTField->pMesh_Bot.size(); idx++) {

				pMeshFM_Bot.push_back(pSTField->pMesh_Bot[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}

			for (int idx = 0; idx < pSTField->pMesh_Top.size(); idx++) {

				pMeshFM_Top.push_back(pSTField->pMesh_Top[idx]->pMeshCUDA->cuMesh.get_managed_object());
			}

			for (int idx = 0; idx < pSTField->paMesh_Bot.size(); idx++) {

				pMeshAtom_Bot.push_back(pSTField->paMesh_Bot[idx]->paMeshCUDA->cuaMesh.get_managed_object());
			}

			for (int idx = 0; idx < pSTField->paMesh_Top.size(); idx++) {

				pMeshAtom_Top.push_back(pSTField->paMesh_Top[idx]->paMeshCUDA->cuaMesh.get_managed_object());
			}

			initialized = true;
		}
	}
	else fixed_polarization = true;

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_STFIELD, 
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_STFIELD, 
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError STFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(STFieldCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
