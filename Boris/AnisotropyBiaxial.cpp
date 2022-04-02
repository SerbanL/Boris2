#include "stdafx.h"
#include "AnisotropyBiaxial.h"

#ifdef MODULE_COMPILATION_ANIBI

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "AnisotropyBiaxialCUDA.h"
#endif

//This module implements combination of uniaxial (2nd order) and biaxial anisotropies (fourth order) as:

//E/V = K1 * (1 - (m.e1)^2) + K2 * ((m.e2)^2 * (m.e3)^2)
//
//K1, K2, e1, e2, e3 are the usual user-controllable mesh parameters
//

Anisotropy_Biaxial::Anisotropy_Biaxial(Mesh *pMesh_) :
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

BError Anisotropy_Biaxial::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_Biaxial));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_ANIBI || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_ANIBI || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError Anisotropy_Biaxial::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_Biaxial));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Anisotropy_Biaxial::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Anisotropy_Biaxial));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Anisotropy_BiaxialCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Anisotropy_Biaxial::UpdateField(void)
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
				DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
				DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

				//calculate m.ea1 dot product (uniaxial contribution)
				double u1 = (pMesh->M[idx] * mcanis_ea1) / Ms;

				//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
				double b1 = (pMesh->M[idx] * mcanis_ea2) / Ms;
				double b2 = (pMesh->M[idx] * mcanis_ea3) / Ms;

				//update effective field with the anisotropy field
				DBL3 Heff_value = (2 / (MU0*Ms)) * (K1 * u1 * mcanis_ea1 - K2 * (b1*b2*b2 * mcanis_ea2 + b1*b1*b2 * mcanis_ea3));

				pMesh->Heff[idx] += Heff_value;

				energy += K1 * (1 - u1*u1) + K2 * b1*b1*b2*b2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_energy.linear_size()) Module_energy[idx] = K1 * (1 - u1*u1) + K2 * b1*b1*b2*b2;
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
				DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
				DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->K1_AFM, K1_AFM, pMesh->K2_AFM, K2_AFM, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

				//calculate m.ea1 dot product (uniaxial contribution)
				double u1_A = (pMesh->M[idx] * mcanis_ea1) / Ms_AFM.i;
				double u1_B = (pMesh->M2[idx] * mcanis_ea1) / Ms_AFM.j;

				//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
				double b1_A = (pMesh->M[idx] * mcanis_ea2) / Ms_AFM.i;
				double b1_B = (pMesh->M2[idx] * mcanis_ea2) / Ms_AFM.j;
				
				double b2_A = (pMesh->M[idx] * mcanis_ea3) / Ms_AFM.i;
				double b2_B = (pMesh->M2[idx] * mcanis_ea3) / Ms_AFM.j;

				//update effective field with the anisotropy field
				DBL3 Heff_value = (2 / (MU0*Ms_AFM.i)) * (K1_AFM.i * u1_A * mcanis_ea1 - K2_AFM.i * (b1_A*b2_A*b2_A * mcanis_ea2 + b1_A*b1_A*b2_A * mcanis_ea3));
				DBL3 Heff_value2 = (2 / (MU0*Ms_AFM.j)) * (K1_AFM.j * u1_B * mcanis_ea1 - K2_AFM.j * (b1_B*b2_B*b2_B * mcanis_ea2 + b1_B*b1_B*b2_B * mcanis_ea3));

				pMesh->Heff[idx] += Heff_value;
				pMesh->Heff2[idx] += Heff_value2;

				energy += (K1_AFM.i * (1 - u1_A*u1_A) + K2_AFM.i * b1_A*b1_A*b2_A*b2_A + K1_AFM.j * (1 - u1_B*u1_B) + K2_AFM.j * b1_B*b1_B*b2_B*b2_B) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Heff_value2;
				if (Module_energy.linear_size()) Module_energy[idx] = K1_AFM.i * (1 - u1_A*u1_A) + K2_AFM.i * b1_A*b1_A*b2_A*b2_A;
				if (Module_energy2.linear_size()) Module_energy2[idx] = K1_AFM.j * (1 - u1_B*u1_B) + K2_AFM.j * b1_B*b1_B*b2_B*b2_B;
			}
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy /= pMesh->M.get_nonempty_cells();
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//-------------------Energy methods

//FM Mesh
double Anisotropy_Biaxial::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double Ms = pMesh->Ms;
		double K1 = pMesh->K1;
		double K2 = pMesh->K2;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		//calculate m.ea1 dot product (uniaxial contribution)
		double u1 = pMesh->M[spin_index] * mcanis_ea1 / Ms;

		//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
		double b1 = pMesh->M[spin_index] * mcanis_ea2 / Ms;
		double b2 = pMesh->M[spin_index] * mcanis_ea3 / Ms;

		if (Mnew != DBL3()) {

			//calculate m.ea1 dot product (uniaxial contribution)
			double u1new = Mnew * mcanis_ea1 / Ms;

			//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
			double b1new = Mnew * mcanis_ea2 / Ms;
			double b2new = Mnew * mcanis_ea3 / Ms;

			return pMesh->h.dim() * ((K1 * (1 - u1new * u1new) + K2 * b1new*b1new*b2new*b2new) - (K1 * (1 - u1 * u1) + K2 * b1*b1*b2*b2));
		}
		else return pMesh->h.dim() * (K1 * (1 - u1 * u1) + K2 * b1*b1*b2*b2);
	}
	else return 0.0;
}

//AFM mesh
DBL2 Anisotropy_Biaxial::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL2 K1_AFM = pMesh->K1_AFM;
		DBL2 K2_AFM = pMesh->K2_AFM;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms_AFM, Ms_AFM, pMesh->K1_AFM, K1_AFM, pMesh->K2_AFM, K2_AFM, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		//calculate m.ea1 dot product (uniaxial contribution)
		double u1_A = (pMesh->M[spin_index] * mcanis_ea1) / Ms_AFM.i;
		double u1_B = (pMesh->M2[spin_index] * mcanis_ea1) / Ms_AFM.j;

		//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
		double b1_A = (pMesh->M[spin_index] * mcanis_ea2) / Ms_AFM.i;
		double b1_B = (pMesh->M2[spin_index] * mcanis_ea2) / Ms_AFM.j;

		double b2_A = (pMesh->M[spin_index] * mcanis_ea3) / Ms_AFM.i;
		double b2_B = (pMesh->M2[spin_index] * mcanis_ea3) / Ms_AFM.j;

		double energyA = K1_AFM.i * (1 - u1_A * u1_A) + K2_AFM.i * b1_A*b1_A*b2_A*b2_A;
		double energyB = K1_AFM.j * (1 - u1_B * u1_B) + K2_AFM.j * b1_B*b1_B*b2_B*b2_B;

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			//calculate m.ea1 dot product (uniaxial contribution)
			double u1_Anew = (Mnew_A * mcanis_ea1) / Ms_AFM.i;
			double u1_Bnew = (Mnew_B * mcanis_ea1) / Ms_AFM.j;

			//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
			double b1_Anew = (Mnew_A * mcanis_ea2) / Ms_AFM.i;
			double b1_Bnew = (Mnew_B * mcanis_ea2) / Ms_AFM.j;

			double b2_Anew = (Mnew_A * mcanis_ea3) / Ms_AFM.i;
			double b2_Bnew = (Mnew_B * mcanis_ea3) / Ms_AFM.j;

			double energyAnew = K1_AFM.i * (1 - u1_Anew * u1_Anew) + K2_AFM.i * b1_Anew*b1_Anew*b2_Anew*b2_Anew;
			double energyBnew = K1_AFM.j * (1 - u1_Bnew * u1_Bnew) + K2_AFM.j * b1_Bnew*b1_Bnew*b2_Bnew*b2_Bnew;

			return pMesh->h.dim() * DBL2(energyAnew - energyA, energyBnew - energyB);
		}
		else return pMesh->h.dim() * DBL2(energyA, energyB);
	}
	else return DBL2();
}

//-------------------Torque methods

DBL3 Anisotropy_Biaxial::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif