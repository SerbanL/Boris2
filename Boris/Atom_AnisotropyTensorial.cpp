#include "stdafx.h"
#include "Atom_AnisotropyTensorial.h"

#if defined(MODULE_COMPILATION_ANITENS) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_AnisotropyTensorialCUDA.h"
#endif

Atom_Anisotropy_Tensorial::Atom_Anisotropy_Tensorial(Atom_Mesh *paMesh_) :
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

BError Atom_Anisotropy_Tensorial::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_Tensorial));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect,
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_ANITENS || paMesh->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_ANITENS || paMesh->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_Anisotropy_Tensorial::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_Tensorial));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Atom_Anisotropy_Tensorial::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_Tensorial));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_Anisotropy_TensorialCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_Anisotropy_Tensorial::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double K1 = paMesh->K1;
			double K2 = paMesh->K2;
			double K3 = paMesh->K3;
			DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
			DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
			DBL3 mcanis_ea3 = paMesh->mcanis_ea3;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->K1, K1, paMesh->K2, K2, paMesh->K3, K3, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2, paMesh->mcanis_ea3, mcanis_ea3);

			//calculate dot products
			double a = (paMesh->M1[idx] * mcanis_ea1) / mu_s;
			double b = (paMesh->M1[idx] * mcanis_ea2) / mu_s;
			double c = (paMesh->M1[idx] * mcanis_ea3) / mu_s;

			DBL3 Heff_value;
			double energy_ = 0.0;

			for (int tidx = 0; tidx < paMesh->Kt.size(); tidx++) {

				//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
				//(-d / mu0 mu_s) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

				double ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
				double ap = 0.0, bp = 0.0, cp = 0.0;
				if (paMesh->Kt[tidx].j > 0) { ap1 = pow(a, paMesh->Kt[tidx].j - 1); ap = ap1 * a; }
				else ap = pow(a, paMesh->Kt[tidx].j);
				if (paMesh->Kt[tidx].k > 0) { bp1 = pow(b, paMesh->Kt[tidx].k - 1); bp = bp1 * b; }
				else bp = pow(b, paMesh->Kt[tidx].k);
				if (paMesh->Kt[tidx].l > 0) { cp1 = pow(c, paMesh->Kt[tidx].l - 1); cp = cp1 * c; }
				else cp = pow(c, paMesh->Kt[tidx].l);

				double coeff;
				int order = paMesh->Kt[tidx].j + paMesh->Kt[tidx].k + paMesh->Kt[tidx].l;
				if (order == 2) coeff = -K1 * paMesh->Kt[tidx].i / (MUB_MU0*mu_s);
				else if (order == 4) coeff = -K2 * paMesh->Kt[tidx].i / (MUB_MU0*mu_s);
				else if (order == 6) coeff = -K3 * paMesh->Kt[tidx].i / (MUB_MU0*mu_s);
				else coeff = -paMesh->Kt[tidx].i / (MUB_MU0*mu_s);

				Heff_value += coeff * (paMesh->Kt[tidx].j * ap1*bp*cp * mcanis_ea1 + paMesh->Kt[tidx].k * ap*bp1*cp * mcanis_ea2 + paMesh->Kt[tidx].l * ap*bp*cp1 * mcanis_ea3);

				energy_ += -coeff * MUB_MU0*mu_s * ap*bp*cp;
			}

			paMesh->Heff1[idx] += Heff_value;

			energy += energy_;

			if (Module_Heff.linear_size()) Module_Heff[idx] = Heff_value;
			if (Module_energy.linear_size()) Module_energy[idx] = energy_ / paMesh->M1.h.dim();
		}
	}

	//convert to energy density and return
	if (non_empty_volume) this->energy = energy / non_empty_volume;
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_Anisotropy_Tensorial::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double K1 = paMesh->K1;
		double K2 = paMesh->K2;
		double K3 = paMesh->K3;
		DBL3 mcanis_ea1 = paMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = paMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = paMesh->mcanis_ea3;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->K1, K1, paMesh->K2, K2, paMesh->K3, K3, paMesh->mcanis_ea1, mcanis_ea1, paMesh->mcanis_ea2, mcanis_ea2, paMesh->mcanis_ea3, mcanis_ea3);
		
		auto Get_Energy = [&](double a, double b, double c) -> double {

			double energy_ = 0.0;

			for (int tidx = 0; tidx < paMesh->Kt.size(); tidx++) {

				double coeff;
				int order = paMesh->Kt[tidx].j + paMesh->Kt[tidx].k + paMesh->Kt[tidx].l;
				if (order == 2) coeff = K1 * paMesh->Kt[tidx].i;
				else if (order == 4) coeff = K2 * paMesh->Kt[tidx].i;
				else if (order == 6) coeff = K3 * paMesh->Kt[tidx].i;
				else coeff = paMesh->Kt[tidx].i;

				energy_ += coeff * pow(a, paMesh->Kt[tidx].j)*pow(b, paMesh->Kt[tidx].k)*pow(c, paMesh->Kt[tidx].l);
			}

			return energy_;
		};

		//calculate dot products
		double a = paMesh->M1[spin_index].normalized() * mcanis_ea1;
		double b = paMesh->M1[spin_index].normalized() * mcanis_ea2;
		double c = paMesh->M1[spin_index].normalized() * mcanis_ea3;

		double energy_ = Get_Energy(a, b, c);

		if (Mnew != DBL3()) {

			double anew = Mnew.normalized() * mcanis_ea1;
			double bnew = Mnew.normalized() * mcanis_ea2;
			double cnew = Mnew.normalized() * mcanis_ea3;

			double energynew_ = Get_Energy(anew, bnew, cnew);

			return energynew_ - energy_;
		}
		else return energy_;
	}
	else return 0.0;
}

//-------------------Torque methods

DBL3 Atom_Anisotropy_Tensorial::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(paMesh->M1, avRect);
}

#endif