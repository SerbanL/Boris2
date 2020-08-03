#include "stdafx.h"
#include "Mesh_Diamagnetic.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

DiaMesh::DiaMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_DIAMAGNETIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId),
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar),
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), 
			VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(Temp_l), VINFO(u_disp), VINFO(strain_diag), VINFO(strain_odiag), VINFO(pMod), VINFO(exclude_from_multiconvdemag),
			//Members in this derived class
			//Material Parameters
			VINFO(neta_dia), VINFO(susrel), VINFO(cHA), VINFO(cHmo),
			VINFO(elecCond), VINFO(De), VINFO(n_density), VINFO(SHA), VINFO(iSHA), VINFO(l_sf), VINFO(Gi), VINFO(Gmix),
			VINFO(base_temperature), VINFO(T_equation), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Demag), IINFO(SDemag_Demag), IINFO(SurfExchange),
			IINFO(Zeeman), IINFO(MOptical),
			IINFO(Transport), IINFO(Heat)
		}),
	meshODE(this)
{
	atomic_moment = 0.0;
	T_Curie = 0.0;
	T_Curie_material = 0.0;
}

DiaMesh::DiaMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Mesh(MESH_DIAMAGNETIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId),
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar),
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), 
			VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(Temp_l), VINFO(u_disp), VINFO(strain_diag), VINFO(strain_odiag), VINFO(pMod), VINFO(exclude_from_multiconvdemag),
			//Members in this derived class
			//Material Parameters
			VINFO(neta_dia), VINFO(susrel), VINFO(cHA), VINFO(cHmo),
			VINFO(elecCond), VINFO(De), VINFO(n_density), VINFO(SHA), VINFO(iSHA), VINFO(l_sf), VINFO(Gi), VINFO(Gmix),
			VINFO(base_temperature), VINFO(T_equation), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Demag), IINFO(SDemag_Demag), IINFO(SurfExchange),
			IINFO(Zeeman), IINFO(MOptical),
			IINFO(Transport), IINFO(Heat)
		}),
	meshODE(this)
{
	atomic_moment = 0.0;
	T_Curie = 0.0;
	T_Curie_material = 0.0;

	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_MAGNETIZATION;

	meshRect = meshRect_;

	h = h_;
	h_e = h_;
	h_t = h_;
	h_m = h_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//default modules configuration
	if (!error_on_create) error_on_create = AddModule(MOD_ZEEMAN);

	//--------------------------

	//If cuda is enabled we also need to make the cuda mesh version
	if (cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

DiaMesh::~DiaMesh()
{
}

void DiaMesh::RepairObjectState(void)
{
	//at this point Heff is empty and must not be since MComputation_Enabled will report wrong
	Heff.assign(h, meshRect, DBL3(0, 0, 0));

	//Make sure Zeeman module is always the first one in the list : Zeeman module sets Heff (if Zeeman module disabled then PrepareIteration clears Heff)
	if (IsModuleSet(MOD_ZEEMAN)) {

		int idxZeeman = pMod.get_index_from_ID(MOD_ZEEMAN);
		if (idxZeeman != 0) pMod.move(idxZeeman);
	}
}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError DiaMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(string(CLASS_STR(DiaMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	///////////////////////////////////////////////////////
	//Mesh specific configuration
	///////////////////////////////////////////////////////

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		//resize arrays held in Mesh - if changing mesh to new size then use resize : this maps values from old mesh to new mesh size (keeping magnitudes). Otherwise assign magnetization along x in whole mesh.
		if (M.linear_size()) {

			if (!M.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else if (!M.assign(h, meshRect, DBL3(0, 0, 0))) return error(BERROR_OUTOFMEMORY_CRIT);

		if (!Heff.assign(h, meshRect, DBL3(0, 0, 0))) return error(BERROR_OUTOFMEMORY_CRIT);

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

	//------------------------

	///////////////////////////////////////////////////////
	//Update configuration for mesh modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) {

			error = pMod[idx]->UpdateConfiguration(cfgMessage);
		}
	}

	///////////////////////////////////////////////////////
	//Update configuration for mesh ode solver
	///////////////////////////////////////////////////////

	if (!error) error = meshODE.UpdateConfiguration(cfgMessage);

	return error;
}

void DiaMesh::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
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

	///////////////////////////////////////////////////////
	//Update configuration for mesh ode solver
	///////////////////////////////////////////////////////

	meshODE.UpdateConfiguration_Values(cfgMessage);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		pMeshCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError DiaMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(DiaMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshBaseCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshBaseCUDA = new DiaMeshCUDA(this);
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

	//SwitchCUDA state for differential equation
	if (!error) error = meshODE.SwitchCUDAState(cudaState);

	//--------------------------------------------

	//SwitchCUDA state for all active modules in this mesh
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) error = pMod[idx]->SwitchCUDAState(cudaState);
	}

#endif

	return error;
}

#endif