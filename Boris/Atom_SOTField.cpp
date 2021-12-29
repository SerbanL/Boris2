#include "stdafx.h"
#include "Atom_SOTField.h"

#if defined(MODULE_COMPILATION_SOTFIELD) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_SOTFieldCUDA.h"
#endif

Atom_SOTField::Atom_SOTField(Atom_Mesh* paMesh_) :
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

Atom_SOTField::~Atom_SOTField()
{
}

BError Atom_SOTField::Initialize(void)
{
	BError error(CLASS_STR(Atom_SOTField));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_SOTFIELD, 
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_SOTFIELD);
	if (!error)	initialized = true;

	return error;
}

BError Atom_SOTField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_SOTField));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Atom_SOTField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_SOTField));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_SOTFieldCUDA(paMesh->paMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_SOTField::UpdateField(void)
{
	if (!paMesh->E.linear_size()) return 0.0;

	double conv = paMesh->M1.h.dim() / MUB;

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		double grel = paMesh->grel;
		double mu_s = paMesh->mu_s;
		double SHA = paMesh->SHA;
		double flSOT = paMesh->flSOT;
		DBL3 STp = paMesh->STp;
		paMesh->update_parameters_mcoarse(idx, paMesh->grel, grel, paMesh->mu_s, mu_s, paMesh->SHA, SHA, paMesh->flSOT, flSOT, paMesh->STp, STp);

		if (IsNZ(grel)) {

			double a_const = -conv * (SHA * MUB_E / (GAMMA * grel)) / (mu_s * mu_s * paMesh->GetMeshDimensions().z);

			int idx_E = paMesh->E.position_to_cellidx(paMesh->M1.cellidx_to_position(idx));
			DBL3 p_vec = (STp ^ paMesh->E[idx_E]) * paMesh->elC[idx_E];

			DBL3 SOTField = a_const * ((paMesh->M1[idx] ^ p_vec) + flSOT * mu_s * p_vec);
			paMesh->Heff1[idx] += SOTField;

			if (Module_Heff.linear_size()) Module_Heff[idx] = SOTField;
		}
	}

	//don't count this as a contribution to the total energy density
	return 0.0;
}

#endif