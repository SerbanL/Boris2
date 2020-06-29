#include "stdafx.h"
#include "Atom_Demag_N.h"

#if defined(MODULE_COMPILATION_DEMAG_N) && ATOMISTIC == 1

#include "Atom_Mesh.h"

#if COMPILECUDA == 1
#include "Atom_Demag_NCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Atom_Demag_N::Atom_Demag_N(Atom_Mesh *paMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	paMesh = paMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_Demag_N::~Atom_Demag_N()
{
}

BError Atom_Demag_N::Initialize(void)
{
	initialized = true;

	return BError(CLASS_STR(Atom_Demag_N));
}

BError Atom_Demag_N::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Demag_N));

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

BError Atom_Demag_N::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Demag_N));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_Demag_NCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Demag_N::UpdateField(void)
{
	double energy = 0;

	//used to convert moment to magnetisation in each atomistic unit cell
	double conversion = MUB / paMesh->M1.h.dim();

#pragma omp parallel for reduction (+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
			DBL2 Nxy = paMesh->Nxy;

			DBL3 Heff_value = DBL3(-Nxy.x * paMesh->M1[idx].x, -Nxy.y * paMesh->M1[idx].y, -(1 - Nxy.x - Nxy.y) * paMesh->M1[idx].z) * conversion;

			paMesh->Heff1[idx] += Heff_value;

			//energy density contribution (to be scaled)
			energy += paMesh->M1[idx] * Heff_value;
		}
	}

	if (paMesh->M1.get_nonempty_cells()) this->energy = -energy * MUB_MU0 / (2 * paMesh->M1.get_nonempty_cells());
	else this->energy = 0;

	return this->energy;
}

#endif