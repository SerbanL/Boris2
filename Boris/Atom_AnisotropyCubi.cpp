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

	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_Anisotropy_Cubic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_Cubic));

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
			double K = paMesh->K;
			DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
			DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->K, K, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			DBL3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			double d1 = (paMesh->M1[idx] * mcanis_ea1) / mu_s;
			double d2 = (paMesh->M1[idx] * mcanis_ea2) / mu_s;
			double d3 = (paMesh->M1[idx] * mcanis_ea3) / mu_s;

			//terms for K contribution
			double a1 = d1 * (d2*d2 + d3*d3);
			double a2 = d2 * (d1*d1 + d3*d3);
			double a3 = d3 * (d1*d1 + d2*d2);

			//update effective field with the anisotropy field
			DBL3 Heff_value = DBL3(
				(-2 * K / (MUB_MU0*mu_s)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3),
				(-2 * K / (MUB_MU0*mu_s)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3),
				(-2 * K / (MUB_MU0*mu_s)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
			);

			paMesh->Heff1[idx] += Heff_value;

			//update energy
			energy += K * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3);
		}
	}

	//convert to energy density and return
	if (non_empty_volume) this->energy = energy / non_empty_volume;
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy density methods

double Atom_Anisotropy_Cubic::GetEnergyDensity(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return dynamic_cast<Atom_Anisotropy_CubiCUDA*>(pModuleCUDA)->GetEnergyDensity(avRect);
#endif

	double energy = 0;

	int num_points = 0;

#pragma omp parallel for reduction(+:energy, num_points)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		//only average over values in given rectangle
		if (!avRect.contains(paMesh->M1.cellidx_to_position(idx))) continue;

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double K = paMesh->K;
			DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
			DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->K, K, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			DBL3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			double d1 = (paMesh->M1[idx] * mcanis_ea1) / mu_s;
			double d2 = (paMesh->M1[idx] * mcanis_ea2) / mu_s;
			double d3 = (paMesh->M1[idx] * mcanis_ea3) / mu_s;

			//update energy
			energy += K * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3);
			num_points++;
		}
	}

	//convert to energy density and return
	if (num_points) energy = energy / (num_points * paMesh->M1.h.dim());
	else energy = 0.0;

	return energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_Anisotropy_Cubic::Get_Atomistic_Energy(int spin_index)
{
	//For CUDA there are separate device functions using by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double K = paMesh->K;
		DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->K, K, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2);

		//vector product of ea1 and ea2 : the third orthogonal axis
		DBL3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

		DBL3 S = paMesh->M1[spin_index].normalized();

		//calculate m.ea1, m.ea2 and m.ea3 dot products
		double d1 = S * mcanis_ea1;
		double d2 = S * mcanis_ea2;
		double d3 = S * mcanis_ea3;

		//Hamiltonian contribution as K * (Sx^2*Sy^2 + Sx^2*Sz^2 + Sy^2*Sz^2), where S is the local spin direction (for easy axes coinciding with the xyz system)
		//This is equivalent to the form -K/2 * (Sx^4 + Sy^4 + Sz^4) - energy zero point differs but that's immaterial.
		//Also note the correct signs here for given easy axes (need to be careful, some publications have this wrong).
		return K * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3);
	}
	else return 0.0;
}

#endif