#include "stdafx.h"
#include "Mesh_Dipole.h"
#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

DipoleMesh::DipoleMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_DIPOLE, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class
			//Material Parameters
			VINFO(Ms), VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(De), VINFO(n_density), VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix), VINFO(base_temperature), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
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
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class
			//Material Parameters
			VINFO(Ms), VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(De), VINFO(n_density), VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix), VINFO(base_temperature), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
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

//----------------------------------- IMPORTANT CONTROL METHODS

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
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

BError DipoleMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(DipoleMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		//delete MeshCUDA object and null (just in case)
		if (pMeshCUDA) delete pMeshCUDA;
		pMeshCUDA = nullptr;

		//then make MeshCUDA object, copying over currently held cpu data
		pMeshCUDA = new DipoleMeshCUDA(this);
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

void DipoleMesh::SetMagnetisationAngle(double polar, double azim, Rect rectangle)
{
#if COMPILECUDA
	if (pMeshCUDA) {

		reinterpret_cast<DipoleMeshCUDA*>(pMeshCUDA)->SetMagnetisationAngle(polar, azim);
		return;
	}
#endif

	M[0] = Polar_to_Cartesian(DBL3(Ms, polar, azim)); 
	recalculateStrayField = true;
}

void DipoleMesh::SetCurieTemperature(double Tc)
{
	if (Tc >= 1) {

		T_Curie = Tc;

		std::vector<double> t_scaling_me;

		t_scaling_me.assign(MAX_TEMPERATURE + 1, 1.0);

		int chunk = (int)floor_epsilon(MAX_TEMPERATURE / OmpThreads);

		#pragma omp parallel for
		for (int thread = 0; thread < OmpThreads; thread++) {

			double me = 1.0;

			for (int t_value = 1 + thread * chunk; t_value <= (thread != OmpThreads - 1 ? (thread + 1) * chunk : (int)MAX_TEMPERATURE); t_value++) {

				double c = 3 * Tc / t_value;

				//This is F(x) = B(cx + d) - x
				std::function<double(double)> F = [&](double x) -> double {

					double arg = c * x;
					return ((1 + exp(-2 * arg)) / (1 - exp(-2 * arg))) - 1 / arg - x;
				};

				//This is F'(x) = B'(cx + d) * c - 1
				std::function<double(double)> Fdiff = [&](double x) -> double {

					double arg = c * x;
					double coth = (1 + exp(-2 * arg)) / (1 - exp(-2 * arg));
					return (((1 - coth * coth) + 1 / (arg*arg)) * c - 1);
				};

				//solve for me at temperature = t_value : Ms0 scaling
				//next starting value for the root finder is me - this makes the NewtonRaphson algorithm extremely efficient (typically a single iteration is required to keep within 1e-6 accuracy!)
				me = Root_NewtonRaphson(F, Fdiff, me, 1e-6);

				//Ms scaling
				t_scaling_me[t_value] = me;
			}
		}

		Ms.set_precalculated_scaling_array(t_scaling_me);
		Ms.update(base_temperature);
	}
	else {

		//turn it off

		//reset temperature dependencies for affected parameters
		Ms.clear_t_scaling();
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