#include "stdafx.h"
#include "Atom_iDMExchange.h"

#if defined(MODULE_COMPILATION_IDMEXCHANGE) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_iDMExchangeCUDA.h"
#endif

Atom_iDMExchange::Atom_iDMExchange(Atom_Mesh *paMesh_) :
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

BError Atom_iDMExchange::Initialize(void)
{
	BError error(CLASS_STR(Atom_iDMExchange));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_ActualModule_Heff_Display() == MOD_IDMEXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)paMesh->Get_ActualModule_Heff_Display() == MOD_IDMEXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_EXCH));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_iDMExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_iDMExchange));

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

BError Atom_iDMExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_iDMExchange));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_iDMExchangeCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_iDMExchange::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			double mu_s = paMesh->mu_s;
			double J = paMesh->J;
			double D = paMesh->D;
			paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->J, J, paMesh->D, D);

			//update effective field with the Heisenberg and iDMI exchange field
			DBL3 Hexch_A = J * paMesh->M1.ngbr_dirsum(idx) / (MUB_MU0*mu_s);
			DBL3 Hexch_D = D * paMesh->M1.zanisotropic_ngbr_dirsum(idx) / (MUB_MU0*mu_s);

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
	if (non_empty_volume) this->energy = -MUB_MU0 * energy / (2*non_empty_volume);
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_iDMExchange::Get_Atomistic_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double J = paMesh->J;
		double D = paMesh->D;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->J, J, paMesh->D, D);

		//local spin energy given:
		//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
		//2) iDM exchange: D * Sum_over_neighbors_j((rij x z) . (Si x Sj))

		//Now zanisotropic_ngbr_dirsum returns (rij x z) x Sj, and Si . ((rij x z) x Sj) = -Si. (Sj x (rij x z)) = -(rij x z) . (Si x Sj)
		//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.zanisotropic_ngbr_dirsum(spin_index));

		return (Mnew.normalized() - paMesh->M1[spin_index].normalized()) * (-J * paMesh->M1.ngbr_dirsum(spin_index) - D * paMesh->M1.zanisotropic_ngbr_dirsum(spin_index));
	}
	else return 0.0;
}

#endif