#include "stdafx.h"
#include "Atom_AnisotropyCubi.h"

#if defined(MODULE_COMPILATION_ANICUBI) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_AnisotropyCubiCUDA.h"
#endif

Atom_Anisotropy_Cubic::Atom_Anisotropy_Cubic(Atom_Mesh *paMesh_) :
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

BError Atom_Anisotropy_Cubic::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_Cubic));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_ANICUBI || paMesh->IsOutputDataSet_withRect(DATA_E_ANIS), 
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_ANICUBI || paMesh->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_Anisotropy_Cubic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_Cubic));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Atom_Anisotropy_Cubic::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_Cubic));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_Anisotropy_CubiCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Anisotropy_Cubic::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double K1 = paMesh->K1;
			double K2 = paMesh->K2;
			DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
			DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
			DBL3 mcanis_ea3 = paMesh->mcanis_ea3;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->K1, K1, paMesh->K2, K2, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2, paMesh->mcanis_ea3, mcanis_ea3);

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			double d1 = (paMesh->M1[idx] * mcanis_ea1) / mu_s;
			double d2 = (paMesh->M1[idx] * mcanis_ea2) / mu_s;
			double d3 = (paMesh->M1[idx] * mcanis_ea3) / mu_s;

			//terms for K1 contribution
			double a1 = d1 * (d2*d2 + d3*d3);
			double a2 = d2 * (d1*d1 + d3*d3);
			double a3 = d3 * (d1*d1 + d2*d2);

			//terms for K2 contribution
			double d123 = d1*d2*d3;

			double b1 = d123 * d2*d3;
			double b2 = d123 * d1*d3;
			double b3 = d123 * d1*d2;

			//update effective field with the anisotropy field
			DBL3 Heff_value = DBL3(
				(-2 * K1 / (MUB_MU0*mu_s)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
				+ (-2 * K2 / (MUB_MU0*mu_s)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),
				(-2 * K1 / (MUB_MU0*mu_s)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
				+ (-2 * K2 / (MUB_MU0*mu_s)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),
				(-2 * K1 / (MUB_MU0*mu_s)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
				+ (-2 * K2 / (MUB_MU0*mu_s)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
			);

			paMesh->Heff1[idx] += Heff_value;

			//update energy
			energy += K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123;

			if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
			if (Module_energy.linear_size()) Module_energy[idx] = (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123) / paMesh->M1.h.dim();
		}
	}

	//convert to energy density and return
	if (non_empty_volume) this->energy = energy / non_empty_volume;
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_Anisotropy_Cubic::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double K1 = paMesh->K1;
		double K2 = paMesh->K2;
		DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = paMesh->mcanis_ea3;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->K1, K1, paMesh->K2, K2, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2, paMesh->mcanis_ea3, mcanis_ea3);

		DBL3 S = paMesh->M1[spin_index].normalized();
		DBL3 S_new = Mnew.normalized();

		//calculate m.ea1, m.ea2 and m.ea3 dot products
		double d1 = S * mcanis_ea1;
		double d2 = S * mcanis_ea2;
		double d3 = S * mcanis_ea3;

		if (Mnew != DBL3()) {

			double d1_new = S_new * mcanis_ea1;
			double d2_new = S_new * mcanis_ea2;
			double d3_new = S_new * mcanis_ea3;

			//Hamiltonian contribution as K * (Sx^2*Sy^2 + Sx^2*Sz^2 + Sy^2*Sz^2), where S is the local spin direction (for easy axes coinciding with the xyz system)
			//This is equivalent to the form -K/2 * (Sx^4 + Sy^4 + Sz^4) - energy zero point differs but that's immaterial.
			//Also note the correct signs here for given easy axes (need to be careful, some publications have this wrong).
			return K1 * (d1_new*d1_new*d2_new*d2_new + d1_new*d1_new*d3_new*d3_new + d2_new*d2_new*d3_new*d3_new - d1*d1*d2*d2 - d1*d1*d3*d3 - d2*d2*d3*d3)
				+ K2 * (d1_new*d2_new*d3_new*d1_new*d2_new*d3_new - d1*d2*d3*d1*d2*d3);
		}
		else return K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d1*d2*d3*d1*d2*d3;
	}
	else return 0.0;
}

//-------------------Torque methods

DBL3 Atom_Anisotropy_Cubic::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(paMesh->M1, avRect);
}

#endif