#include "stdafx.h"
#include "Demag_N.h"

#ifdef MODULE_COMPILATION_DEMAG_N

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
	BError error(CLASS_STR(Demag_N));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_DEMAG_N || pMesh->IsOutputDataSet_withRect(DATA_E_DEMAG),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_DEMAG_N || pMesh->IsOutputDataSet_withRect(DATA_E_DEMAG));
	if (!error)	initialized = true;

	return error;
}

BError Demag_N::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag_N));

	Uninitialize();

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

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * (pMesh->M[idx] * Heff_value) / 2;
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

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * (Mval * Heff_value) / 2;
			}
		}
	}

	if (pMesh->M.get_nonempty_cells()) this->energy = -energy * MU0 / (2 * pMesh->M.get_nonempty_cells());
	else this->energy = 0;

	return this->energy;
}

//-------------------Energy methods

//FM mesh
double Demag_N::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
		DBL2 Nxy = pMesh->Nxy;
		double Nz = (1 - Nxy.x - Nxy.y);

		if (Mnew != DBL3()) {

			return (MU0 / 2) * pMesh->h.dim() * (
				(Mnew * DBL3(Nxy.x * Mnew.x, Nxy.y * Mnew.y, Nz * Mnew.z)) -
				(pMesh->M[spin_index] * DBL3(Nxy.x * pMesh->M[spin_index].x, Nxy.y * pMesh->M[spin_index].y, Nz * pMesh->M[spin_index].z)));
		}
		else return (MU0 / 2) * pMesh->h.dim() * (pMesh->M[spin_index] * DBL3(Nxy.x * pMesh->M[spin_index].x, Nxy.y * pMesh->M[spin_index].y, Nz * pMesh->M[spin_index].z));
	}
	else return 0.0;
}

//AFM mesh
DBL2 Demag_N::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
		DBL2 Nxy = pMesh->Nxy;
		double Nz = (1 - Nxy.x - Nxy.y);

		DBL3 M = (pMesh->M[spin_index] + pMesh->M2[spin_index]) / 2;
		DBL3 Mnew = (Mnew_A + Mnew_B) / 2;

		double energy_ = 0.0;

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			energy_ = (MU0 / 2) * pMesh->h.dim() * (
				(Mnew * DBL3(Nxy.x * Mnew.x, Nxy.y * Mnew.y, Nz * Mnew.z)) -
				(M * DBL3(Nxy.x * M.x, Nxy.y * M.y, Nz * M.z)));
		}
		else energy_ = (MU0 / 2) * pMesh->h.dim() * (M * DBL3(Nxy.x * M.x, Nxy.y * M.y, Nz * M.z));

		return DBL2(energy_, energy_);
	}
	else return DBL2();
}

#endif