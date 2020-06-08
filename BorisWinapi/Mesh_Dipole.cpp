#include "stdafx.h"
#include "Mesh_Dipole.h"

#ifdef MESH_COMPILATION_DIPOLE

#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

DipoleMesh::DipoleMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_DIPOLE, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), 
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), 
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class
			//Material Parameters
			VINFO(Ms), VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(De), VINFO(n_density), VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix), 
			VINFO(base_temperature), VINFO(T_equation), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
		},
		{
			IINFO(Transport), IINFO(Heat)
		})
{
	M.assign(SZ3(1), DBL3(-Ms, 0, 0));

	recalculateStrayField = true;
}

DipoleMesh::DipoleMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Mesh(MESH_DIPOLE, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), 
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), 
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class
			//Material Parameters
			VINFO(Ms), VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(De), VINFO(n_density), VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix), 
			VINFO(base_temperature), VINFO(T_equation), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
		},
		{
			IINFO(Transport), IINFO(Heat)
		})
{
	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_NONE;

	meshRect = meshRect_;
	n = SZ3(1);
	
	h = meshRect.size();
	h_e = h_;
	h_t = h_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//--------------------------

	//If cuda is enabled we also need to make the cuda mesh version
	if (cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

void DipoleMesh::RepairObjectState(void)
{
	//calculate scaling function (Curie Weiss law)
	if (T_Curie > 0) {

		pCurieWeiss->Initialize_CurieWeiss(0.0, T_Curie);
	}
	else {

		pCurieWeiss->Initialize_CurieWeiss(0.0, 1.0);
		pLongRelSus->Initialize_LongitudinalRelSusceptibility(pCurieWeiss->get_data(), atomic_moment, 1.0);
	}
}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError DipoleMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(string(CLASS_STR(DipoleMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	n = SZ3(1);
	h = meshRect.size();

	///////////////////////////////////////////////////////
	//Mesh specific configuration
	///////////////////////////////////////////////////////

	if (M.linear_size()) M.resize(h, meshRect);
	else M.assign(h, meshRect, DBL3(-Ms, 0, 0));

	//update material parameters spatial dependence as cellsize and rectangle could have changed
	if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

	//update any text equations used in mesh parameters (dependence on mesh dimensions possible)
	if (!error) update_all_meshparam_equations();

	//set dipole value from Ms - this could have changed
	Reset_Mdipole();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (!error) error = pMeshCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	///////////////////////////////////////////////////////
	//Update configuration for mesh modules
	///////////////////////////////////////////////////////

	//update configuration in all currently set modules in this mesh
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) {

			error = pMod[idx]->UpdateConfiguration(cfgMessage);
		}
	}

	return error;
}

void DipoleMesh::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
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

BError DipoleMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(DipoleMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshCUDA = new DipoleMeshCUDA(this);
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

	//--------------------------------------------

#endif

	return error;
}

//----------------------------------- VARIOUS SET METHODS

//set magnitude for Mdipole (also setting recalculateStrayField flag)
void DipoleMesh::Reset_Mdipole(void)
{
	//set dipole value from Ms - this could have changed
	if (!Temp.linear_size()) M.renormalize((double)Ms);
	else {

		double Temperature = Temp.average_nonempty_omp();
		M.renormalize((double)Ms.get(Temperature));
	}

	recalculateStrayField = true;

#if COMPILECUDA == 1
	if (pMeshCUDA) reinterpret_cast<DipoleMeshCUDA*>(pMeshCUDA)->Reset_Mdipole();
#endif
}

void DipoleMesh::SetMagAngle(double polar, double azim, Rect rectangle)
{
#if COMPILECUDA
	if (pMeshCUDA) {

		reinterpret_cast<DipoleMeshCUDA*>(pMeshCUDA)->SetMagAngle(polar, azim);
		return;
	}
#endif

	M[0] = Polar_to_Cartesian(DBL3(Ms, polar, azim)); 
	recalculateStrayField = true;
}

void DipoleMesh::SetCurieTemperature(double Tc, bool set_default_dependences)
{
	//Curie temperature is a constant in mesh parameter equations, so update them
	if (Tc != T_Curie) update_all_meshparam_equations();

	if (Tc > 0) {

		T_Curie = Tc;

		if (set_default_dependences) {

			//set default temperature dependences
			Ms.set_t_scaling_equation(string("me(T/Tc)"), userConstants, T_Curie, base_temperature);

			//make sure to also update them - this method can be called during a simulation, e.g. if field changes.
			Ms.update(base_temperature);
		}
	}
	else {

		//turn it off
		T_Curie = 0.0;

		if (set_default_dependences) {

			//reset temperature dependencies for affected parameters
			Ms.clear_t_scaling();
		}
	}

	Reset_Mdipole();

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		//T_Curie changed : sync with cuda version
		pMeshCUDA->T_Curie.from_cpu((cuBReal)T_Curie);
	}
#endif
}

//----------------------------------- VARIOUS GET METHODS

bool DipoleMesh::Check_recalculateStrayField(void)
{
	//recalculate stray field if flag is set (Mdipole has changed) or non-uniform temperature is enabled and Ms has a temperature dependence (in this case always recalculate)

	if (Temp.linear_size() && Ms.is_tdep()) {

		double Temperature = Temp.average_nonempty_omp();
		M.renormalize((double)Ms.get(Temperature));

		return true;
	}

	return recalculateStrayField;
}

#endif