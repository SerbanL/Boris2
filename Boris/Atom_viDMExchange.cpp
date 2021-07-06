#include "stdafx.h"
#include "Atom_viDMExchange.h"

#if defined(MODULE_COMPILATION_VIDMEXCHANGE) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_viDMExchangeCUDA.h"
#endif

Atom_viDMExchange::Atom_viDMExchange(Atom_Mesh *paMesh_) :
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

BError Atom_viDMExchange::Initialize(void)
{
	BError error(CLASS_STR(Atom_viDMExchange));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect,
		(MOD_)paMesh->Get_ActualModule_Heff_Display() == MOD_VIDMEXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)paMesh->Get_ActualModule_Heff_Display() == MOD_VIDMEXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_EXCH));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_viDMExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_viDMExchange));

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

BError Atom_viDMExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_viDMExchange));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_viDMExchangeCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_viDMExchange::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double J = paMesh->J;
			double D = paMesh->D;
			DBL3 D_dir = paMesh->D_dir;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->J, J, paMesh->D, D, paMesh->D_dir, D_dir);

			//update effective field with the Heisenberg and iDMI exchange field
			DBL3 Hexch_A = J * paMesh->M1.ngbr_dirsum(idx) / (MUB_MU0*mu_s);

			DBL3 hexch_D_x = paMesh->M1.xanisotropic_ngbr_dirsum(idx);
			DBL3 hexch_D_y = paMesh->M1.yanisotropic_ngbr_dirsum(idx);
			DBL3 hexch_D_z = paMesh->M1.zanisotropic_ngbr_dirsum(idx);

			DBL3 Hexch_D = D * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z) / (MUB_MU0*mu_s);

			paMesh->Heff1[idx] += (Hexch_A + Hexch_D);

			//update energy E = -mu_s * Bex
			energy += paMesh->M1[idx] * (Hexch_A + Hexch_D);

			//spatial dependence display of effective field and energy density
			if (Module_Heff.linear_size() && Module_energy.linear_size()) {

				if ((MOD_)paMesh->Get_Module_Heff_Display() == MOD_EXCHANGE) {

					//total : direct and DMI
					Module_Heff[idx] = Hexch_A + Hexch_D;
					Module_energy[idx] = -MUB_MU0 * (paMesh->M1[idx] * (Hexch_A + Hexch_D)) / (2 * paMesh->M1.h.dim());
				}
				else {

					//just DMI
					Module_Heff[idx] = Hexch_D;
					Module_energy[idx] = -MUB_MU0 * (paMesh->M1[idx] * Hexch_D) / (2 * paMesh->M1.h.dim());
				}
			}
		}
	}

	//convert to energy density and return. Divide by two since in the Hamiltonian the sum is performed only once for every pair of spins, but if you use the M.H expression each sum appears twice.
	//Also note, this energy density is not the same as the micromagnetic one, due to different zero-energy points.
	if (non_empty_volume) this->energy = -MUB_MU0 * energy / (2 * non_empty_volume);
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_viDMExchange::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double J = paMesh->J;
		double D = paMesh->D;
		DBL3 D_dir = paMesh->D_dir;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->J, J, paMesh->D, D, paMesh->D_dir, D_dir);

		DBL3 hexch_D_x = paMesh->M1.xanisotropic_ngbr_dirsum(spin_index);
		DBL3 hexch_D_y = paMesh->M1.yanisotropic_ngbr_dirsum(spin_index);
		DBL3 hexch_D_z = paMesh->M1.zanisotropic_ngbr_dirsum(spin_index);
		DBL3 hexch_D = D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z;

		if (Mnew != DBL3()) return (Mnew.normalized() - paMesh->M1[spin_index].normalized()) * (-J * paMesh->M1.ngbr_dirsum(spin_index) - D * hexch_D);
		else return paMesh->M1[spin_index].normalized() * (-J * paMesh->M1.ngbr_dirsum(spin_index) - D * hexch_D);
	}
	else return 0.0;
}

//-------------------Torque methods

DBL3 Atom_viDMExchange::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<Atom_viDMExchangeCUDA*>(pModuleCUDA)->GetTorque(avRect);
#endif

	return CalculateTorque(paMesh->M1, avRect);
}

#endif