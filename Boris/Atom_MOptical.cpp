#include "stdafx.h"
#include "Atom_MOptical.h"

#if defined(MODULE_COMPILATION_MOPTICAL) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_MOpticalCUDA.h"
#endif

Atom_MOptical::Atom_MOptical(Atom_Mesh *paMesh_) :
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

Atom_MOptical::~Atom_MOptical()
{
}

BError Atom_MOptical::Initialize(void)
{
	BError error(CLASS_STR(Atom_MOptical));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_MOPTICAL || paMesh->IsOutputDataSet_withRect(DATA_E_MOPTICAL),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_MOPTICAL || paMesh->IsOutputDataSet_withRect(DATA_E_MOPTICAL));
	if (!error)	initialized = true;

	non_empty_volume = paMesh->Get_NonEmpty_Magnetic_Volume();

	return error;
}

BError Atom_MOptical::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MOptical));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

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

void Atom_MOptical::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError Atom_MOptical::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_MOptical));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_MOpticalCUDA(paMesh->paMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_MOptical::UpdateField(void)
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		double cHmo = paMesh->cHmo;
		paMesh->update_parameters_mcoarse(idx, paMesh->cHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		paMesh->Heff1[idx] += DBL3(0, 0, cHmo);

		energy += -MUB * paMesh->M1[idx] * MU0 * DBL3(0, 0, cHmo);

		if (Module_Heff.linear_size()) Module_Heff[idx] = DBL3(0, 0, cHmo);
		if (Module_energy.linear_size()) Module_energy[idx] = -MUB * paMesh->M1[idx] * MU0 * DBL3(0, 0, cHmo) / paMesh->M1.h.dim();
	}

	//convert to energy density and return
	if (non_empty_volume) this->energy = energy / non_empty_volume;
	else this->energy = 0.0;

	return this->energy;
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_MOptical::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (paMesh->M1.is_not_empty(spin_index)) {

		double cHmo = paMesh->cHmo;
		paMesh->update_parameters_mcoarse(spin_index, paMesh->cHmo, cHmo);

		if (Mnew != DBL3()) return -MUB * (Mnew - paMesh->M1[spin_index]) * MU0 * DBL3(0, 0, cHmo);
		else return -MUB * paMesh->M1[spin_index] * MU0 * DBL3(0, 0, cHmo);
	}
	else return 0.0;
}

#endif