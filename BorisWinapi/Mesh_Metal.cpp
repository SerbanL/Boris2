#include "stdafx.h"
#include "Mesh_Metal.h"
#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

MetalMesh::MetalMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_METAL, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(elecCond), VINFO(De), VINFO(SHA), VINFO(iSHA), VINFO(l_sf), VINFO(Gi), VINFO(Gmix), VINFO(base_temperature), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Transport), IINFO(Heat)
		})
{}

MetalMesh::MetalMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Mesh(MESH_METAL, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(elecCond), VINFO(De), VINFO(SHA), VINFO(iSHA), VINFO(l_sf), VINFO(Gi), VINFO(Gmix), VINFO(base_temperature), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Transport), IINFO(Heat)
		})
{
	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_NONE;

	meshRect = meshRect_;
	
	h_e = h_;
	h_t = h_;

	error_on_create = UpdateConfiguration();

	//--------------------------

	//If cuda is enabled we also need to make the cuda mesh version
	if (cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

void MetalMesh::RepairObjectState(void)
{

}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
BError MetalMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(string(CLASS_STR(MetalMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

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

BError MetalMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(MetalMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		//delete MeshCUDA object and null (just in case)
		if (pMeshCUDA) delete pMeshCUDA;
		pMeshCUDA = nullptr;

		//then make MeshCUDA object, copying over currently held cpu data
		pMeshCUDA = new MetalMeshCUDA(this);
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