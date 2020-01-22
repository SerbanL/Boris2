#include "stdafx.h"
#include "Mesh_Insulator.h"
#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

InsulatorMesh::InsulatorMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_INSULATOR, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n_t), VINFO(h_t), VINFO(Temp), VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(base_temperature), VINFO(T_equation), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
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
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n_t), VINFO(h_t), VINFO(Temp), VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(base_temperature), VINFO(T_equation), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
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

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

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

	///////////////////////////////////////////////////////
	//Mesh specific configuration
	///////////////////////////////////////////////////////

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		//update material parameters spatial dependence as cellsize and rectangle could have changed
		if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

		//update any text equations used in mesh parameters (dependence on mesh dimensions possible)
		if (!error) update_meshparam_equations();
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (!error) error = pMeshCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	///////////////////////////////////////////////////////
	//Update configuration for mesh modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) {

			error = pMod[idx]->UpdateConfiguration(cfgMessage);
		}
	}

	return error;
}

void InsulatorMesh::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	///////////////////////////////////////////////////////
	//Update configuration in this mesh
	///////////////////////////////////////////////////////

	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		//Update text equations other than those used in mesh parameters
		UpdateTEquationUserConstants();

		//update any text equations used in mesh parameters
		update_meshparam_equations();
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (T_equation.is_set()) T_equation.clear();
	}

	///////////////////////////////////////////////////////
	//Update configuration for mesh modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx]) {

			pMod[idx]->UpdateConfiguration_Values(cfgMessage);
		}
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		pMeshCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError InsulatorMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(InsulatorMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshCUDA = new InsulatorMeshCUDA(this);
			error = pMeshCUDA->Error_On_Create();
			if (!error) error = pMeshCUDA->cuMesh()->set_pointers(pMeshCUDA);
		}
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