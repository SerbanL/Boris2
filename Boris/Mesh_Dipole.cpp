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
			VINFO(dipole_velocity), VINFO(dipole_shift_debt), VINFO(dipole_shift_clip), VINFO(dipole_last_time),
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
			VINFO(dipole_velocity), VINFO(dipole_shift_debt), VINFO(dipole_shift_clip), VINFO(dipole_last_time),
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
	BError error(std::string(CLASS_STR(DipoleMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	n = SZ3(1);
	h = meshRect.size();
	dipole_last_time = pSMesh->GetTime();

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

	//reset time for dipole shifting algorithm
	dipole_last_time = pSMesh->GetTime();

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

//called at the start of each iteration
void DipoleMesh::PrepareNewIteration(void)
{
	Dipole_Shifting_Algorithm();

	if (recalculateStrayField && !strayField_recalculated) strayField_recalculated = true;
	else if (strayField_recalculated) { strayField_recalculated = false; recalculateStrayField = false; }
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
	BError error(std::string(CLASS_STR(DipoleMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshBaseCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshBaseCUDA = new DipoleMeshCUDA(this);
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
	if (pMeshCUDA) dynamic_cast<DipoleMeshCUDA*>(pMeshCUDA)->Reset_Mdipole();
#endif
}

void DipoleMesh::SetMagAngle(double polar, double azim, Rect rectangle)
{
#if COMPILECUDA
	if (pMeshCUDA) {

		dynamic_cast<DipoleMeshCUDA*>(pMeshCUDA)->SetMagAngle(polar, azim);
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
			Ms.set_t_scaling_equation(std::string("me(T/Tc)"), userConstants, T_Curie, base_temperature);

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

//shift dipole mesh rectangle by given amount
void DipoleMesh::Shift_Dipole(DBL3 shift)
{
	meshRect += shift;

	M.rect += shift;
	if (V.linear_size()) V.rect += shift;
	if (S.linear_size()) S.rect += shift;
	if (elC.linear_size()) elC.rect += shift;
	if (Temp.linear_size()) Temp.rect += shift;

#if COMPILECUDA == 1
	if (pMeshCUDA) dynamic_cast<DipoleMeshCUDA*>(pMeshCUDA)->Shift_Dipole();
#endif

	recalculateStrayField = true;
	strayField_recalculated = false;
}

//at the start of each iteration see if we need to implement a moving dipole
void DipoleMesh::Dipole_Shifting_Algorithm(void)
{
	if (dipole_velocity != DBL3() && pSMesh->CurrentTimeStepSolved()) {

		//current time so we can calculate required shift
		double dipole_current_time = pSMesh->GetTime();

		//if current time less than stored previous time then something is wrong (e.g. ode was reset - reset shifting debt as well)
		if (dipole_current_time < dipole_last_time) {

			dipole_last_time = dipole_current_time;
			dipole_shift_debt = DBL3();
		}

		//add to total amount of shiting which hasn't yet been executed (the shift debt)
		dipole_shift_debt += (dipole_current_time - dipole_last_time) * dipole_velocity;

		//clip the shift to execute if required
		DBL3 shift = DBL3(
			dipole_shift_clip.x > 0.0 ? 0.0 : dipole_shift_debt.x,
			dipole_shift_clip.y > 0.0 ? 0.0 : dipole_shift_debt.y,
			dipole_shift_clip.z > 0.0 ? 0.0 : dipole_shift_debt.z);

		if (dipole_shift_clip.x > 0.0 && fabs(dipole_shift_debt.x) > dipole_shift_clip.x)
			shift.x = floor(fabs(dipole_shift_debt.x) / dipole_shift_clip.x) * dipole_shift_clip.x * get_sign(dipole_shift_debt.x);

		if (dipole_shift_clip.y > 0.0 && fabs(dipole_shift_debt.y) > dipole_shift_clip.y)
			shift.y = floor(fabs(dipole_shift_debt.y) / dipole_shift_clip.y) * dipole_shift_clip.y * get_sign(dipole_shift_debt.y);

		if (dipole_shift_clip.z > 0.0 && fabs(dipole_shift_debt.z) > dipole_shift_clip.z)
			shift.z = floor(fabs(dipole_shift_debt.z) / dipole_shift_clip.z) * dipole_shift_clip.z * get_sign(dipole_shift_debt.z);

		//execute shift if needed
		if (shift != DBL3()) {

			Shift_Dipole(shift);
			dipole_shift_debt -= shift;
		}

		dipole_last_time = dipole_current_time;
	}
}

//----------------------------------- VARIOUS GET METHODS

bool DipoleMesh::CheckRecalculateStrayField(void)
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