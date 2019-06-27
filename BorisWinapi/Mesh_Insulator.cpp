#include "stdafx.h"
#include "Mesh_Insulator.h"
#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

InsulatorMesh::InsulatorMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_INSULATOR, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n_t), VINFO(h_t), VINFO(Temp), VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(base_temperature), VINFO(thermCond), VINFO(density), VINFO(shc)
		},
		{
			//Modules Implementations
			IINFO(Heat)
		})
{}

InsulatorMesh::InsulatorMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Mesh(MESH_INSULATOR, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n_t), VINFO(h_t), VINFO(Temp), VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(base_temperature), VINFO(thermCond), VINFO(density), VINFO(shc)
		},
		{
			//Modules Implementations
			IINFO(Heat)
		})
{
	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_NONE;

	meshRect = meshRect_;

	h_t = h_;

	error_on_create = UpdateConfiguration();

	//--------------------------

	//If cuda is enabled we also need to make the cuda mesh version
	if (cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

void InsulatorMesh::RepairObjectState(void)
{

}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
BError InsulatorMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(string(CLASS_STR(InsulatorMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (!error) error = pMeshCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	//change mesh dimensions in all currently set effective field modules
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) {

			error = pMod[idx]->UpdateConfiguration(cfgMessage);
		}
	}

	//update material parameters spatial dependence as cellsize and rectangle could have changed
	if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

	return error;
}

BError InsulatorMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(InsulatorMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		//delete MeshCUDA object and null (just in case)
		if (pMeshCUDA) delete pMeshCUDA;
		pMeshCUDA = nullptr;

		//then make MeshCUDA object, copying over currently held cpu data
		pMeshCUDA = new InsulatorMeshCUDA(this);
		error = pMeshCUDA->Error_On_Create();
		if (!error) error = pMeshCUDA->cuMesh()->set_pointers(pMeshCUDA);
	}
	else {

		//delete MeshCUDA object and null
		if (pMeshCUDA) delete pMeshCUDA;
		pMeshCUDA = nullptr;
	}

	//--------------------------------------------

	//SwitchCUDA state for all active modules in this mesh
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) error = pMod[idx]->SwitchCUDAState(cudaState);
	}

#endif

	return error;
}

//----------------------------------- VARIOUS GET/SET METHODS