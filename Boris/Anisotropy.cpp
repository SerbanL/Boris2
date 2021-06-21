#include "stdafx.h"
#include "Anisotropy.h"

#ifdef MODULE_COMPILATION_ANIUNI

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "AnisotropyCUDA.h"
#endif

//
//
// Let theta be the angle between magnetization and z axis
// Let phi be the xy-plane angle between the projection of magnetizaiton on the xy-plane and the x axis
//
// Spherical polar coordinates:
//
// alpha = cos(phi) * sin (theta)
// beta = sin(phi) * sin(theta)
// gamma = cos(theta)
//
// Hk = - (1 / MU0) dE/dM = - (1 / MU0 Ms) dE/dm  ( dE/dM = x^ dE/dMx + y^ dE/dMy + z^ dE/dMz )

//Uniaxial anisotropy (simple type):
//
//(E/V) = K1 * sin^2(theta) (J/m^3)
//
//i.e. theta can be taken as the angle between magnetization and easy axis. Thus in this case we define K1 and an easy axis direction
//
//Thus (E/V) = K1 * [ 1 - (M.ea/Ms)^2]  where ea is the easy axis direction. Theta is then the angle between M and the easy axis.
//
//Thus can calculate the anisotropy field Hk (which is added into the Heff field) as :
//
//Hk = 2K1/(MU0*Ms) * ea * (ea.m)
//
//Thus best to specify K1 with 4 components : K1, eax, eay, eaz
//
//For uniaxial anisotropy with 4th order term : (E/V) = K1 * sin^2(theta) + K2 * sin^4(theta) = K1 * (1 - (m.ea)^2) + K2 * (1 - (m.ea)^2)^2 (J/m^3)
//
//As above, calculate the anisotropy field as :
//
//Hk = 2/(MU0*Ms) [ K1 * ea * (ea.m) + 2 * K2 * ea * (ea.m) - 2 * K2 * ea * (ea.m)^3 ]
//

Anisotropy_Uniaxial::Anisotropy_Uniaxial(Mesh *pMesh_) : 
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

BError Anisotropy_Uniaxial::Initialize(void) 
{ 	
	BError error(CLASS_STR(Anisotropy_Uniaxial));
	
	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_ANIUNI || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_ANIUNI || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}
	
BError Anisotropy_Uniaxial::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_Uniaxial));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Anisotropy_Uniaxial::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Anisotropy_Uniaxial));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Anisotropy_UniaxialCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Anisotropy_Uniaxial::UpdateField(void)
{
	double energy = 0;

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double K1 = pMesh->K1;
				double K2 = pMesh->K2;
				DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->mcanis_ea1, mcanis_ea1);

				//calculate m.ea dot product
				double dotprod = (pMesh->M[idx] * mcanis_ea1) / Ms;

				//update effective field with the anisotropy field
				DBL3 Heff_value = (2 / (MU0*Ms)) * dotprod * (K1 + 2 * K2 * (1 - dotprod * dotprod)) * mcanis_ea1;

				pMesh->Heff[idx] += Heff_value;

				//update energy (E/V) = K1 * sin^2(theta) + K2 * sin^4(theta) = K1 * [ 1 - dotprod*dotprod ] + K2 * [1 - dotprod * dotprod]^2
				energy += (K1 + K2 * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod);

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_energy.linear_size()) Module_energy[idx] = (K1 + K2 * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod);
			}
		}
	}

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				DBL2 K1_AFM = pMesh->K1_AFM;
				DBL2 K2_AFM = pMesh->K2_AFM;
				DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->K1_AFM, K1_AFM, pMesh->K2_AFM, K2_AFM, pMesh->mcanis_ea1, mcanis_ea1);

				//calculate m.ea dot product
				double dotprod = (pMesh->M[idx] * mcanis_ea1) / Ms_AFM.i;
				double dotprod2 = (pMesh->M2[idx] * mcanis_ea1) / Ms_AFM.j;

				//update effective field with the anisotropy field
				DBL3 Heff_value = (2 / (MU0*Ms_AFM.i)) * dotprod * (K1_AFM.i + 2 * K2_AFM.i * (1 - dotprod * dotprod)) * mcanis_ea1;
				DBL3 Heff_value2 = (2 / (MU0*Ms_AFM.j)) * dotprod2 * (K1_AFM.j + 2 * K2_AFM.j * (1 - dotprod2 * dotprod2)) * mcanis_ea1;

				pMesh->Heff[idx] += Heff_value;
				pMesh->Heff2[idx] += Heff_value2;

				//update energy (E/V) = K1 * sin^2(theta) + K2 * sin^4(theta) = K1 * [ 1 - dotprod*dotprod ] + K2 * [1 - dotprod * dotprod]^2
				energy += ((K1_AFM.i + K2_AFM.i * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod) + (K1_AFM.j + K2_AFM.j * (1 - dotprod2 * dotprod2)) * (1 - dotprod2 * dotprod2)) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Heff_value2;
				if (Module_energy.linear_size()) Module_energy[idx] = (K1_AFM.i + K2_AFM.i * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod);
				if (Module_energy2.linear_size()) Module_energy2[idx] = (K1_AFM.j + K2_AFM.j * (1 - dotprod2 * dotprod2)) * (1 - dotprod2 * dotprod2);
			}
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy /= pMesh->M.get_nonempty_cells();
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//-------------------Energy methods

double Anisotropy_Uniaxial::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double K1 = pMesh->K1;
		double K2 = pMesh->K2;
		double Ms = pMesh->Ms;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->mcanis_ea1, mcanis_ea1);

		//calculate m.ea dot product
		double dotprod = (pMesh->M[spin_index] * mcanis_ea1) / Ms;
		double dotprod_new = Mnew * mcanis_ea1 / Ms;
		double dpsq = dotprod * dotprod;
		double dpsq_new = dotprod_new * dotprod_new;

		return pMesh->h.dim() * ((K1 + K2 * (1 - dpsq_new)) * (1 - dpsq_new) - (K1 + K2 * (1 - dpsq)) * (1 - dpsq));
	}
	else return 0.0;
}

double Anisotropy_Uniaxial::Get_Energy(int spin_index)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double K1 = pMesh->K1;
		double K2 = pMesh->K2;
		double Ms = pMesh->Ms;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->mcanis_ea1, mcanis_ea1);

		//calculate m.ea dot product
		double dotprod = (pMesh->M[spin_index] * mcanis_ea1) / Ms;
		double dpsq = dotprod * dotprod;

		return pMesh->h.dim() * (K1 + K2 * (1 - dpsq)) * (1 - dpsq);
	}
	else return 0.0;
}

//-------------------Torque methods

DBL3 Anisotropy_Uniaxial::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<Anisotropy_UniaxialCUDA*>(pModuleCUDA)->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif