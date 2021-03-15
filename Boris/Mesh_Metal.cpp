#include "stdafx.h"
#include "Mesh_Metal.h"

#ifdef MESH_COMPILATION_METAL

#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

MetalMesh::MetalMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_METAL, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), 
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(Module_Heff_Display), VINFO(Module_Energy_Display),
			VINFO(meshRect), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), 
			VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(Temp_l), 
			VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(elecCond), VINFO(De), VINFO(n_density), VINFO(SHA), VINFO(iSHA), VINFO(l_sf), VINFO(Gi), VINFO(Gmix), 
			VINFO(base_temperature), VINFO(T_equation), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q)
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
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId),
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(Module_Heff_Display), VINFO(Module_Energy_Display),
			VINFO(meshRect), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m),
			VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(Temp_l),
			VINFO(pMod),
			//Members in this derived class

			//Material Parameters
			VINFO(elecCond), VINFO(De), VINFO(n_density), VINFO(SHA), VINFO(iSHA), VINFO(l_sf), VINFO(Gi), VINFO(Gmix), 
			VINFO(base_temperature), VINFO(T_equation), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q)
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
	h_m = h_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

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

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError MetalMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(std::string(CLASS_STR(MetalMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	///////////////////////////////////////////////////////
	//Mesh specific configuration
	///////////////////////////////////////////////////////

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		//update material parameters spatial dependence as cellsize and rectangle could have changed
		if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

		//update any text equations used in mesh parameters (dependence on mesh dimensions possible)
		if (!error) update_all_meshparam_equations();
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

void MetalMesh::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	///////////////////////////////////////////////////////
	//Update configuration in this mesh
	///////////////////////////////////////////////////////

	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		//Update text equations other than those used in mesh parameters
		UpdateTEquationUserConstants();

		//update any text equations used in mesh parameters
		update_all_meshparam_equations();
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

BError MetalMesh::SwitchCUDAState(bool cudaState)
{
	BError error(std::string(CLASS_STR(MetalMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshBaseCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshBaseCUDA = new MetalMeshCUDA(this);
			pMeshCUDA = dynamic_cast<MeshCUDA*>(pMeshBaseCUDA);

			error = pMeshBaseCUDA->Error_On_Create();
			if (!error) error = pMeshCUDA->cuMesh()->set_pointers(pMeshCUDA);
		}
	}
	else {

		//delete MeshCUDA object and null
		if (pMeshBaseCUDA) delete pMeshBaseCUDA;
		pMeshBaseCUDA = nullptr;
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

#endif