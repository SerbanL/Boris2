#include "stdafx.h"
#include "Atom_Exchange.h"

#if defined(MODULE_COMPILATION_EXCHANGE) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_ExchangeCUDA.h"
#endif

Atom_Exchange::Atom_Exchange(Atom_Mesh *paMesh_) :
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

BError Atom_Exchange::Initialize(void)
{
	BError error(CLASS_STR(Atom_Exchange));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_EXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_EXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_EXCH));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_Exchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Exchange));

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

BError Atom_Exchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Exchange));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_ExchangeCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Exchange::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double J = paMesh->J;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->J, J);
		
			//update effective field with the Heisenberg exchange field
			DBL3 Heff_value = (J / (MUB_MU0*mu_s)) * paMesh->M1.ngbr_dirsum(idx);

			paMesh->Heff1[idx] += Heff_value;

			//update energy E = -mu_s * Bex. Will finish off at the end with prefactors.
			energy += paMesh->M1[idx] * Heff_value;

			if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
			if (Module_energy.linear_size()) Module_energy[idx] = -MUB_MU0 * paMesh->M1[idx] * Heff_value / (2  * paMesh->M1.h.dim());
		}
	}

	//convert to energy density and return. Divide by two since in the Hamiltonian the sum is performed only once for every pair of spins, but if you use the M.H expression each sum appears twice.
	//Also note, this energy density is not the same as the micromagnetic one, due to different zero-energy points.
	//To obtain the micromagnetic energy density you also have to subtract the energy density obtained at saturation from the Heisenberg Hamiltonian.
	//Could be done here easily, but decided not to (add J * number of neigbhors to energy).
	if (non_empty_volume) this->energy = -MUB_MU0 * energy / (2*non_empty_volume);
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_Exchange::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double J = paMesh->J;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->J, J);

		//local spin energy given -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
		return -J * ((Mnew.normalized() - paMesh->M1[spin_index].normalized()) * paMesh->M1.ngbr_dirsum(spin_index));
	}
	else return 0.0;
}

double Atom_Exchange::Get_Energy(int spin_index)
{
	if (paMesh->M1.is_not_empty(spin_index)) {

		double J = paMesh->J;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->J, J);

		return -J * paMesh->M1[spin_index].normalized() * paMesh->M1.ngbr_dirsum(spin_index);
	}
	else return 0.0;
}

//-------------------Torque methods

DBL3 Atom_Exchange::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<Atom_ExchangeCUDA*>(pModuleCUDA)->GetTorque(avRect);
#endif

	return CalculateTorque(paMesh->M1, avRect);
}

#endif