#include "stdafx.h"
#include "AnisotropyTensorial.h"

#ifdef MODULE_COMPILATION_ANITENS

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "AnisotropyTensorialCUDA.h"
#endif

Anisotropy_Tensorial::Anisotropy_Tensorial(Mesh *pMesh_) :
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

BError Anisotropy_Tensorial::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_Tensorial));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_ANITENS || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_ANITENS || pMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError Anisotropy_Tensorial::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_Tensorial));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_PARAMCHANGED)) {

		Uninitialize();
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Anisotropy_Tensorial::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Anisotropy_Tensorial));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Anisotropy_TensorialCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Anisotropy_Tensorial::UpdateField(void)
{
	double energy = 0;

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double K1 = pMesh->K1;
				double K2 = pMesh->K2;
				double K3 = pMesh->K3;
				DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
				DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
				DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->K3, K3, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

				//calculate dot products
				double a = (pMesh->M[idx] * mcanis_ea1) / Ms;
				double b = (pMesh->M[idx] * mcanis_ea2) / Ms;
				double c = (pMesh->M[idx] * mcanis_ea3) / Ms;

				DBL3 Heff_value;
				double energy_ = 0.0;

				for (int tidx = 0; tidx < pMesh->Kt.size(); tidx++) {

					//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
					//(-d / mu0Ms) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

					double ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
					double ap = 0.0, bp = 0.0, cp = 0.0;
					if (pMesh->Kt[tidx].j > 0) { ap1 = pow(a, pMesh->Kt[tidx].j - 1); ap = ap1 * a; } else ap = pow(a, pMesh->Kt[tidx].j);
					if (pMesh->Kt[tidx].k > 0) { bp1 = pow(b, pMesh->Kt[tidx].k - 1); bp = bp1 * b; } else bp = pow(b, pMesh->Kt[tidx].k);
					if (pMesh->Kt[tidx].l > 0) { cp1 = pow(c, pMesh->Kt[tidx].l - 1); cp = cp1 * c; } else cp = pow(c, pMesh->Kt[tidx].l);

					double coeff;
					int order = pMesh->Kt[tidx].j + pMesh->Kt[tidx].k + pMesh->Kt[tidx].l;
					if (order == 2) coeff = -K1 * pMesh->Kt[tidx].i / (MU0*Ms);
					else if (order == 4) coeff = -K2 * pMesh->Kt[tidx].i / (MU0*Ms);
					else if (order == 6) coeff = -K3 * pMesh->Kt[tidx].i / (MU0*Ms);
					else coeff = -pMesh->Kt[tidx].i / (MU0*Ms);

					Heff_value += coeff * (pMesh->Kt[tidx].j * ap1*bp*cp * mcanis_ea1 + pMesh->Kt[tidx].k * ap*bp1*cp * mcanis_ea2 + pMesh->Kt[tidx].l * ap*bp*cp1 * mcanis_ea3);

					energy_ += -coeff * MU0*Ms * ap*bp*cp;
				}
				
				pMesh->Heff[idx] += Heff_value;

				energy += energy_;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_energy.linear_size()) Module_energy[idx] = energy_;
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
				DBL2 K3_AFM = pMesh->K3_AFM;
				DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
				DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
				DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->K1_AFM, K1_AFM, pMesh->K2_AFM, K2_AFM, pMesh->K3_AFM, K3_AFM, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

				//calculate dot products
				double a = (pMesh->M[idx] * mcanis_ea1) / Ms_AFM.i;
				double b = (pMesh->M[idx] * mcanis_ea2) / Ms_AFM.i;
				double c = (pMesh->M[idx] * mcanis_ea3) / Ms_AFM.i;

				double a2 = (pMesh->M2[idx] * mcanis_ea1) / Ms_AFM.j;
				double b2 = (pMesh->M2[idx] * mcanis_ea2) / Ms_AFM.j;
				double c2 = (pMesh->M2[idx] * mcanis_ea3) / Ms_AFM.j;

				DBL3 Heff_value, Heff_value2;
				double energy_ = 0.0, energy2_ = 0.0;

				for (int tidx = 0; tidx < pMesh->Kt.size(); tidx++) {

					//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
					//(-d / mu0Ms) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

					double ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
					double ap = 0.0, bp = 0.0, cp = 0.0;
					if (pMesh->Kt[tidx].j > 0) { ap1 = pow(a, pMesh->Kt[tidx].j - 1); ap = ap1 * a; }
					else ap = pow(a, pMesh->Kt[tidx].j);
					if (pMesh->Kt[tidx].k > 0) { bp1 = pow(b, pMesh->Kt[tidx].k - 1); bp = bp1 * b; }
					else bp = pow(b, pMesh->Kt[tidx].k);
					if (pMesh->Kt[tidx].l > 0) { cp1 = pow(c, pMesh->Kt[tidx].l - 1); cp = cp1 * c; }
					else cp = pow(c, pMesh->Kt[tidx].l);

					double coeff;
					int order = pMesh->Kt[tidx].j + pMesh->Kt[tidx].k + pMesh->Kt[tidx].l;
					if (order == 2) coeff = -K1_AFM.i * pMesh->Kt[tidx].i / (MU0*Ms_AFM.i);
					else if (order == 4) coeff = -K2_AFM.i * pMesh->Kt[tidx].i / (MU0*Ms_AFM.i);
					else if (order == 6) coeff = -K3_AFM.i * pMesh->Kt[tidx].i / (MU0*Ms_AFM.i);
					else coeff = -pMesh->Kt[tidx].i / (MU0*Ms_AFM.i);

					Heff_value += coeff * (pMesh->Kt[tidx].j * ap1*bp*cp * mcanis_ea1 + pMesh->Kt[tidx].k * ap*bp1*cp * mcanis_ea2 + pMesh->Kt[tidx].l * ap*bp*cp1 * mcanis_ea3);

					energy_ += -coeff * MU0*Ms_AFM.i * ap*bp*cp;
				}

				for (int tidx = 0; tidx < pMesh->Kt2.size(); tidx++) {

					//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
					//(-d / mu0Ms) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

					double ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
					double ap = 0.0, bp = 0.0, cp = 0.0;
					if (pMesh->Kt2[tidx].j > 0) { ap1 = pow(a2, pMesh->Kt2[tidx].j - 1); ap = ap1 * a2; }
					else ap = pow(a2, pMesh->Kt2[tidx].j);
					if (pMesh->Kt2[tidx].k > 0) { bp1 = pow(b2, pMesh->Kt2[tidx].k - 1); bp = bp1 * b2; }
					else bp = pow(b2, pMesh->Kt2[tidx].k);
					if (pMesh->Kt2[tidx].l > 0) { cp1 = pow(c2, pMesh->Kt2[tidx].l - 1); cp = cp1 * c2; }
					else cp = pow(c2, pMesh->Kt2[tidx].l);

					double coeff;
					int order = pMesh->Kt2[tidx].j + pMesh->Kt2[tidx].k + pMesh->Kt2[tidx].l;
					if (order == 2) coeff = -K1_AFM.j * pMesh->Kt2[tidx].i / (MU0*Ms_AFM.j);
					else if (order == 4) coeff = -K2_AFM.j * pMesh->Kt2[tidx].i / (MU0*Ms_AFM.j);
					else if (order == 6) coeff = -K3_AFM.j * pMesh->Kt2[tidx].i / (MU0*Ms_AFM.j);
					else coeff = -pMesh->Kt2[tidx].i / (MU0*Ms_AFM.j);

					Heff_value2 += coeff * (pMesh->Kt2[tidx].j * ap1*bp*cp * mcanis_ea1 + pMesh->Kt2[tidx].k * ap*bp1*cp * mcanis_ea2 + pMesh->Kt2[tidx].l * ap*bp*cp1 * mcanis_ea3);

					energy2_ += -coeff * MU0*Ms_AFM.j * ap*bp*cp;
				}

				pMesh->Heff[idx] += Heff_value;
				pMesh->Heff2[idx] += Heff_value2;

				energy += (energy_ + energy2_) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Heff_value2;
				if (Module_energy.linear_size()) Module_energy[idx] = energy_;
				if (Module_energy2.linear_size()) Module_energy2[idx] = energy2_;
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
double Anisotropy_Tensorial::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double Ms = pMesh->Ms;
		double K1 = pMesh->K1;
		double K2 = pMesh->K2;
		double K3 = pMesh->K3;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms, Ms, pMesh->K1, K1, pMesh->K2, K2, pMesh->K3, K3, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		auto Get_Energy = [&](double a, double b, double c) -> double {

			double energy_ = 0.0;

			for (int tidx = 0; tidx < pMesh->Kt.size(); tidx++) {

				double coeff;
				int order = pMesh->Kt[tidx].j + pMesh->Kt[tidx].k + pMesh->Kt[tidx].l;
				if (order == 2) coeff = K1 * pMesh->Kt[tidx].i;
				else if (order == 4) coeff = K2 * pMesh->Kt[tidx].i;
				else if (order == 6) coeff = K3 * pMesh->Kt[tidx].i;
				else coeff = pMesh->Kt[tidx].i;

				energy_ += coeff * pow(a, pMesh->Kt[tidx].j)*pow(b, pMesh->Kt[tidx].k)*pow(c, pMesh->Kt[tidx].l);
			}

			return energy_;
		};

		//calculate dot products
		double a = pMesh->M[spin_index] * mcanis_ea1 / Ms;
		double b = pMesh->M[spin_index] * mcanis_ea2 / Ms;
		double c = pMesh->M[spin_index] * mcanis_ea3 / Ms;

		double energy_ = Get_Energy(a, b, c);

		if (Mnew != DBL3()) {

			double anew = Mnew * mcanis_ea1 / Ms;
			double bnew = Mnew * mcanis_ea2 / Ms;
			double cnew = Mnew * mcanis_ea3 / Ms;

			double energynew_ = Get_Energy(anew, bnew, cnew);

			return pMesh->h.dim() * (energynew_ - energy_);
		}
		else return pMesh->h.dim() * energy_;
	}
	else return 0.0;
}

//AFM mesh
DBL2 Anisotropy_Tensorial::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL2 K1_AFM = pMesh->K1_AFM;
		DBL2 K2_AFM = pMesh->K2_AFM;
		DBL2 K3_AFM = pMesh->K3_AFM;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms_AFM, Ms_AFM, pMesh->K1_AFM, K1_AFM, pMesh->K2_AFM, K2_AFM, pMesh->K3_AFM, K3_AFM, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		auto Get_Energy = [&](double a, double b, double c, double a2, double b2, double c2) -> DBL2 {

			double energyA = 0.0, energyB = 0.0;

			for (int tidx = 0; tidx < pMesh->Kt.size(); tidx++) {

				double coeff;
				int order = pMesh->Kt[tidx].j + pMesh->Kt[tidx].k + pMesh->Kt[tidx].l;
				if (order == 2) coeff = K1_AFM.i * pMesh->Kt[tidx].i;
				else if (order == 4) coeff = K2_AFM.i * pMesh->Kt[tidx].i;
				else if (order == 6) coeff = K3_AFM.i * pMesh->Kt[tidx].i;
				else coeff = pMesh->Kt[tidx].i;

				energyA += coeff * pow(a, pMesh->Kt[tidx].j)*pow(b, pMesh->Kt[tidx].k)*pow(c, pMesh->Kt[tidx].l);
			}

			for (int tidx = 0; tidx < pMesh->Kt2.size(); tidx++) {

				double coeff;
				int order = pMesh->Kt2[tidx].j + pMesh->Kt2[tidx].k + pMesh->Kt2[tidx].l;
				if (order == 2) coeff = K1_AFM.j * pMesh->Kt2[tidx].i;
				else if (order == 4) coeff = K2_AFM.j * pMesh->Kt2[tidx].i;
				else if (order == 6) coeff = K3_AFM.j * pMesh->Kt2[tidx].i;
				else coeff = pMesh->Kt2[tidx].i;

				energyB += coeff * pow(a2, pMesh->Kt2[tidx].j)*pow(b2, pMesh->Kt2[tidx].k)*pow(c2, pMesh->Kt2[tidx].l);
			}

			return DBL2(energyA, energyB);
		};

		//calculate dot products
		double a = (pMesh->M[spin_index] * mcanis_ea1) / Ms_AFM.i;
		double b = (pMesh->M[spin_index] * mcanis_ea2) / Ms_AFM.i;
		double c = (pMesh->M[spin_index] * mcanis_ea3) / Ms_AFM.i;

		double a2 = (pMesh->M2[spin_index] * mcanis_ea1) / Ms_AFM.j;
		double b2 = (pMesh->M2[spin_index] * mcanis_ea2) / Ms_AFM.j;
		double c2 = (pMesh->M2[spin_index] * mcanis_ea3) / Ms_AFM.j;

		DBL2 energy_ = Get_Energy(a, b, c, a2, b2, c2);

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			double anew = Mnew_A * mcanis_ea1 / Ms_AFM.j;
			double bnew = Mnew_A * mcanis_ea2 / Ms_AFM.j;
			double cnew = Mnew_A * mcanis_ea3 / Ms_AFM.j;

			double a2new = Mnew_B * mcanis_ea1 / Ms_AFM.j;
			double b2new = Mnew_B * mcanis_ea2 / Ms_AFM.j;
			double c2new = Mnew_B * mcanis_ea3 / Ms_AFM.j;

			DBL2 energynew_ = Get_Energy(anew, bnew, cnew, a2new, b2new, c2new);

			return pMesh->h.dim() * (energynew_ - energy_);
		}
		else return pMesh->h.dim() * energy_;
	}
	else return DBL2();
}

//-------------------Torque methods

DBL3 Anisotropy_Tensorial::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif