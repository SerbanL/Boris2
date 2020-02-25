#include "stdafx.h"
#include "Demag_N.h"

#ifdef MODULE_DEMAG_N

#include "Mesh.h"

#if COMPILECUDA == 1
#include "Demag_NCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Demag_N::Demag_N(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Demag_N::~Demag_N()
{
}

BError Demag_N::Initialize(void)
{
	initialized = true;

	return BError(CLASS_STR(Demag_N));
}

BError Demag_N::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag_N));

	Uninitialize();

	Initialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Demag_N::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Demag_N));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Demag_NCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Demag_N::UpdateField(void)
{
	double energy = 0;

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction (+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
				DBL2 Nxy = pMesh->Nxy;

				DBL3 Heff_value = DBL3(-Nxy.x * pMesh->M[idx].x, -Nxy.y * pMesh->M[idx].y, -(1 - Nxy.x - Nxy.y) * pMesh->M[idx].z);

				pMesh->Heff[idx] += Heff_value;

				energy += pMesh->M[idx] * Heff_value;
			}
		}
	}

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction (+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
				DBL2 Nxy = pMesh->Nxy;

				DBL3 Mval = (pMesh->M[idx] + pMesh->M2[idx]) / 2;

				DBL3 Heff_value = DBL3(-Nxy.x * Mval.x, -Nxy.y * Mval.y, -(1 - Nxy.x - Nxy.y) * Mval.z);

				pMesh->Heff[idx] += Heff_value;
				pMesh->Heff2[idx] += Heff_value;

				energy += Mval * Heff_value;
			}
		}
	}

	if (pMesh->M.get_nonempty_cells()) this->energy = -energy * MU0 / (2 * pMesh->M.get_nonempty_cells());
	else this->energy = 0;

	return this->energy;
}

#endif