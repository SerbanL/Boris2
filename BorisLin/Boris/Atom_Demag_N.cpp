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
	BError error(CLASS_STR(Atom_Demag_N));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_DEMAG_N || paMesh->IsOutputDataSet_withRect(DATA_E_DEMAG), 
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_DEMAG_N || paMesh->IsOutputDataSet_withRect(DATA_E_DEMAG));
	if (!error)	initialized = true;

	return error;
}

BError Atom_Demag_N::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Demag_N));

	Uninitialize();

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

	//used to convert moment to magnetization in each atomistic unit cell
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

			if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
			if (Module_energy.linear_size()) Module_energy[idx] = -MUB_MU0 * paMesh->M1[idx] * Heff_value / (2 * paMesh->M1.h.dim());
		}
	}

	if (paMesh->M1.get_nonempty_cells()) this->energy = -energy * MU0 * conversion / (2 * paMesh->M1.get_nonempty_cells());
	else this->energy = 0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_Demag_N::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
		DBL2 Nxy = paMesh->Nxy;

		double Nz = (1 - Nxy.x - Nxy.y);

		double r = MUB / paMesh->h.dim();
		DBL3 S = paMesh->M1[spin_index];

		if (Mnew != DBL3()) return (MUB_MU0 / 2) * r * (Nxy.x * (Mnew.x*Mnew.x - S.x*S.x) + Nxy.y * (Mnew.y*Mnew.y - S.y*S.y) + Nz * (Mnew.z*Mnew.z - S.z*S.z));
		else return (MUB_MU0 / 2) * r * (Nxy.x * S.x*S.x + Nxy.y * S.y*S.y + Nz * S.z*S.z);
	}
	else return 0.0;
}

#endif