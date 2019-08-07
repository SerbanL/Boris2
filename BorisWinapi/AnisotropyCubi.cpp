#include "stdafx.h"
#include "AnisotropyCubi.h"

#ifdef MODULE_ANICUBI

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "AnisotropyCubiCUDA.h"
#endif

//Cubic anisotropy
//
//(E/V) = K1 * (a^2 * b^2 + b^2 * g^2 + g^2 * a^2) + K2 * a^2 * b^2 * g^2   (J/m^3)
//
//where a = sin(theta) * cos(phi)
//      b = sin(theta) * sin(phi)
//	    g = cos(theta)
//
// Theta angle between z axis and magnetization
// Phi xy-plane angle of magnetization projection on xy-plane
//

Anisotropy_Cubic::Anisotropy_Cubic(Mesh *pMesh) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	this->pMesh = dynamic_cast<FMesh*>(pMesh);

	error_on_create = UpdateConfiguration();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

BError Anisotropy_Cubic::Initialize(void)
{
	initialized = true;

	return BError(CLASS_STR(Anisotropy_Cubic));
}

BError Anisotropy_Cubic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_Cubic));

	Uninitialize();

	Initialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
	}
#endif

	return error;
}

BError Anisotropy_Cubic::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Anisotropy_Cubic));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Anisotropy_CubicCUDA(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA));
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Anisotropy_Cubic::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			double K1 = pMesh->K1;
			double K2 = pMesh->K2;
			DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
			DBL3 mcanis_ea2 = pMesh->mcanis_ea2;

			pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			DBL3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			double d1 = (pMesh->M[idx] * mcanis_ea1) / Ms;
			double d2 = (pMesh->M[idx] * mcanis_ea2) / Ms;
			double d3 = (pMesh->M[idx] * mcanis_ea3) / Ms;

			//terms for K1 contribution
			double a1 = d1 * (d2*d2 + d3 * d3);
			double a2 = d2 * (d1*d1 + d3 * d3);
			double a3 = d3 * (d1*d1 + d2 * d2);

			//terms for K2 contribution
			double d123 = d1 * d2*d3;

			double b1 = d123 * d2*d3;
			double b2 = d123 * d1*d3;
			double b3 = d123 * d1*d2;

			//update effective field with the anisotropy field
			DBL3 Heff_value = DBL3(
				(-2 * K1 / (MU0*Ms)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
				+ (-2 * K2 / (MU0*Ms)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),

				(-2 * K1 / (MU0*Ms)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
				+ (-2 * K2 / (MU0*Ms)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),

				(-2 * K1 / (MU0*Ms)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
				+ (-2 * K2 / (MU0*Ms)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
			);

			pMesh->Heff[idx] += Heff_value;

			//update energy (E/V)
			energy += K1 * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3) + K2 * d123*d123;
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy /= pMesh->M.get_nonempty_cells();
	else energy = 0;

	this->energy = energy;

	return energy;
}

#endif