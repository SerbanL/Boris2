#include "stdafx.h"
#include "Atom_Anisotropy.h"

#if defined(MODULE_COMPILATION_ANIUNI) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_AnisotropyCUDA.h"
#endif

Atom_Anisotropy_Uniaxial::Atom_Anisotropy_Uniaxial(Atom_Mesh *paMesh_) :
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

BError Atom_Anisotropy_Uniaxial::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_Uniaxial));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_ANIUNI || paMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_ANIUNI || paMesh->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_Anisotropy_Uniaxial::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_Uniaxial));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Atom_Anisotropy_Uniaxial::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_Uniaxial));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_Anisotropy_UniaxialCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Anisotropy_Uniaxial::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double K1 = paMesh->K1;
			double K2 = paMesh->K2;
			DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->K1, K1, paMesh->K2, K2, paMesh->mcanis_ea1, mcanis_ea1);

			//calculate m.ea dot product
			double dotprod = (paMesh->M1[idx] * mcanis_ea1) / mu_s;

			//update effective field with the anisotropy field
			DBL3 Heff_value = (2 / (MUB_MU0*mu_s)) * dotprod * (K1 + 2 * K2 * (1 - dotprod * dotprod)) * mcanis_ea1;

			paMesh->Heff1[idx] += Heff_value;

			energy += (K1 + K2 * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod);

			if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
			if (Module_energy.linear_size()) Module_energy[idx] = (K1 + K2 * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod) / paMesh->M1.h.dim();
		}
	}
	
	//convert to energy density and return
	if (non_empty_volume) this->energy = energy / non_empty_volume;
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_Anisotropy_Uniaxial::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double K1 = paMesh->K1;
		double K2 = paMesh->K2;
		DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->K1, K1, paMesh->K2, K2, paMesh->mcanis_ea1, mcanis_ea1);

		//calculate m.ea dot product
		double dotprod = paMesh->M1[spin_index].normalized() * mcanis_ea1;
		double dpsq = dotprod * dotprod;

		if (Mnew != DBL3()) {

			double dotprod_new = Mnew.normalized() * mcanis_ea1;
			double dpsq_new = dotprod_new * dotprod_new;

			//Hamiltonian contribution as -Ku * (S * ea)^2, where S is the local spin direction
			return -K1 * (dpsq_new - dpsq) - K2 * (dpsq_new * dpsq_new - dpsq * dpsq);
		}
		else return -K1 * dpsq - K2 * dpsq * dpsq;
	}
	else return 0.0;
}

//-------------------Torque methods

DBL3 Atom_Anisotropy_Uniaxial::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(paMesh->M1, avRect);
}

#endif