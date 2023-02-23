#include "stdafx.h"
#include "AnisotropyCubi.h"

#ifdef MODULE_COMPILATION_ANICUBI

#include "Mesh.h"
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
	this->pMesh = pMesh;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

BError Anisotropy_Cubic::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_Cubic));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_ANICUBI || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_ANICUBI || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS), 
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError Anisotropy_Cubic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_Cubic));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
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
		pModuleCUDA = new Anisotropy_CubicCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Anisotropy_Cubic::UpdateField(void)
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

				//calculate m.ea1, m.ea2 and m.ea3 dot products
				double d1 = (pMesh->M[idx] * mcanis_ea1) / Ms;
				double d2 = (pMesh->M[idx] * mcanis_ea2) / Ms;
				double d3 = (pMesh->M[idx] * mcanis_ea3) / Ms;

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
					(-2 * K1 / (MU0*Ms)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
					+ (-2 * K2 / (MU0*Ms)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),

					(-2 * K1 / (MU0*Ms)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
					+ (-2 * K2 / (MU0*Ms)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),

					(-2 * K1 / (MU0*Ms)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
					+ (-2 * K2 / (MU0*Ms)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
				);

				pMesh->Heff[idx] += Heff_value;

				//update energy (E/V)
				energy += K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_energy.linear_size()) Module_energy[idx] = K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123;
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

				//calculate m.ea1, m.ea2 and m.ea3 dot products
				double d1 = (pMesh->M[idx] * mcanis_ea1) / Ms_AFM.i;
				double d2 = (pMesh->M[idx] * mcanis_ea2) / Ms_AFM.i;
				double d3 = (pMesh->M[idx] * mcanis_ea3) / Ms_AFM.i;

				//terms for K1 contribution
				double a1 = d1 * (d2 * d2 + d3 * d3);
				double a2 = d2 * (d1 * d1 + d3 * d3);
				double a3 = d3 * (d1 * d1 + d2 * d2);

				//terms for K2 contribution
				double d123 = d1 * d2*d3;

				double b1 = d123 * d2*d3;
				double b2 = d123 * d1*d3;
				double b3 = d123 * d1*d2;

				//same thing for sub-lattice B

				double d1B = (pMesh->M2[idx] * mcanis_ea1) / Ms_AFM.j;
				double d2B = (pMesh->M2[idx] * mcanis_ea2) / Ms_AFM.j;
				double d3B = (pMesh->M2[idx] * mcanis_ea3) / Ms_AFM.j;

				double a1B = d1B * (d2B * d2B + d3B * d3B);
				double a2B = d2B * (d1B * d1B + d3B * d3B);
				double a3B = d3B * (d1B * d1B + d2B * d2B);

				double d123B = d1B * d2B*d3B;

				double b1B = d123B * d2B*d3B;
				double b2B = d123B * d1B*d3B;
				double b3B = d123B * d1B*d2B;

				//update effective field with the anisotropy field
				DBL3 Heff_value = DBL3(
					(-2 * K1_AFM.i / (MU0*Ms_AFM.i)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
					+ (-2 * K2_AFM.i / (MU0*Ms_AFM.i)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),

					(-2 * K1_AFM.i / (MU0*Ms_AFM.i)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
					+ (-2 * K2_AFM.i / (MU0*Ms_AFM.i)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),

					(-2 * K1_AFM.i / (MU0*Ms_AFM.i)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
					+ (-2 * K2_AFM.i / (MU0*Ms_AFM.i)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
				);

				pMesh->Heff[idx] += Heff_value;

				//same thing for sub-lattice B

				DBL3 Heff_value2 = DBL3(
					(-2 * K1_AFM.j / (MU0*Ms_AFM.j)) * (mcanis_ea1.i * a1B + mcanis_ea2.i * a2B + mcanis_ea3.i * a3B)
					+ (-2 * K2_AFM.j / (MU0*Ms_AFM.j)) * (mcanis_ea1.i * b1B + mcanis_ea2.i * b2B + mcanis_ea3.i * b3B),

					(-2 * K1_AFM.j / (MU0*Ms_AFM.j)) * (mcanis_ea1.j * a1B + mcanis_ea2.j * a2B + mcanis_ea3.j * a3B)
					+ (-2 * K2_AFM.j / (MU0*Ms_AFM.j)) * (mcanis_ea1.j * b1B + mcanis_ea2.j * b2B + mcanis_ea3.j * b3B),

					(-2 * K1_AFM.j / (MU0*Ms_AFM.j)) * (mcanis_ea1.k * a1B + mcanis_ea2.k * a2B + mcanis_ea3.k * a3B)
					+ (-2 * K2_AFM.j / (MU0*Ms_AFM.j)) * (mcanis_ea1.k * b1B + mcanis_ea2.k * b2B + mcanis_ea3.k * b3B)
				);

				pMesh->Heff2[idx] += Heff_value2;

				//update energy (E/V)
				energy += (K1_AFM.i * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2_AFM.i * d123*d123 + K1_AFM.j * (d1B*d1B*d2B*d2B + d1B*d1B*d3B*d3B + d2B*d2B*d3B*d3B) + K2_AFM.j * d123B*d123B) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Heff_value2;
				if (Module_energy.linear_size()) Module_energy[idx] = K1_AFM.i * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3) + K2_AFM.i * d123*d123;
				if (Module_energy2.linear_size()) Module_energy2[idx] = K2_AFM.i * d123*d123 + K1_AFM.j * (d1B*d1B*d2B*d2B + d1B * d1B*d3B*d3B + d2B * d2B*d3B*d3B) + K2_AFM.j * d123B*d123B;
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
double Anisotropy_Cubic::Get_EnergyChange(int spin_index, DBL3 Mnew)
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

		DBL3 S = pMesh->M[spin_index] / Ms;
		DBL3 S_new = Mnew / Ms;

		//calculate m.ea1, m.ea2 and m.ea3 dot products
		double d1 = S * mcanis_ea1;
		double d2 = S * mcanis_ea2;
		double d3 = S * mcanis_ea3;
		double d123 = d1 * d2 * d3;

		if (Mnew != DBL3()) {

			double d1_new = S_new * mcanis_ea1;
			double d2_new = S_new * mcanis_ea2;
			double d3_new = S_new * mcanis_ea3;
			double d123_new = d1_new * d2_new * d3_new;

			return pMesh->h.dim() * (
				(K1 * (d1_new*d1_new*d2_new*d2_new + d1_new*d1_new*d3_new*d3_new + d2_new*d2_new*d3_new*d3_new) + K2 * d123_new*d123_new) 
				- (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123));
		}
		else return pMesh->h.dim() * (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123);
	}
	else return 0.0;
}

//AFM mesh
DBL2 Anisotropy_Cubic::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL2 K1_AFM = pMesh->K1_AFM;
		DBL2 K2_AFM = pMesh->K2_AFM;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms_AFM, Ms_AFM, pMesh->K1_AFM, K1_AFM, pMesh->K2_AFM, K2_AFM, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		//calculate m.ea1, m.ea2 and m.ea3 dot products
		double d1 = (pMesh->M[spin_index] * mcanis_ea1) / Ms_AFM.i;
		double d2 = (pMesh->M[spin_index] * mcanis_ea2) / Ms_AFM.i;
		double d3 = (pMesh->M[spin_index] * mcanis_ea3) / Ms_AFM.i;
		double d123 = d1*d2*d3;

		//same thing for sub-lattice B

		double d1B = (pMesh->M2[spin_index] * mcanis_ea1) / Ms_AFM.j;
		double d2B = (pMesh->M2[spin_index] * mcanis_ea2) / Ms_AFM.j;
		double d3B = (pMesh->M2[spin_index] * mcanis_ea3) / Ms_AFM.j;
		double d123B = d1B*d2B*d3B;

		double energyA = K1_AFM.i * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2_AFM.i * d123*d123;
		double energyB = K1_AFM.j * (d1B*d1B*d2B*d2B + d1B*d1B*d3B*d3B + d2B*d2B*d3B*d3B) + K2_AFM.j * d123B*d123B;

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			double d1new = (pMesh->M[spin_index] * mcanis_ea1) / Ms_AFM.i;
			double d2new = (pMesh->M[spin_index] * mcanis_ea2) / Ms_AFM.i;
			double d3new = (pMesh->M[spin_index] * mcanis_ea3) / Ms_AFM.i;
			double d123new = d1new*d2new*d3new;

			//same thing for sub-lattice B

			double d1Bnew = (pMesh->M2[spin_index] * mcanis_ea1) / Ms_AFM.j;
			double d2Bnew = (pMesh->M2[spin_index] * mcanis_ea2) / Ms_AFM.j;
			double d3Bnew = (pMesh->M2[spin_index] * mcanis_ea3) / Ms_AFM.j;
			double d123Bnew = d1Bnew*d2Bnew*d3Bnew;

			double energyAnew = K1_AFM.i * (d1new*d1new*d2new*d2new + d1new*d1new*d3new*d3new + d2new*d2new*d3new*d3new) + K2_AFM.i * d123new*d123new;
			double energyBnew = K1_AFM.j * (d1Bnew*d1Bnew*d2Bnew*d2Bnew + d1Bnew*d1Bnew*d3Bnew*d3Bnew + d2Bnew*d2Bnew*d3Bnew*d3Bnew) + K2_AFM.j * d123Bnew*d123Bnew;

			return pMesh->h.dim() * DBL2(energyAnew - energyA, energyBnew - energyB);
		}
		else return pMesh->h.dim() * DBL2(energyA, energyB);
	}
	else return DBL2();
}

//-------------------Torque methods

DBL3 Anisotropy_Cubic::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif