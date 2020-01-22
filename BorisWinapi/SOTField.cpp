#include "stdafx.h"
#include "SOTField.h"

#ifdef MODULE_SOTFIELD

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "SOTFieldCUDA.h"
#endif

SOTField::SOTField(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = dynamic_cast<FMesh*>(pMesh_);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

SOTField::~SOTField()
{
}

BError SOTField::Initialize(void)
{
	BError error(CLASS_STR(SOTField));

	initialized = true;

	return error;
}

BError SOTField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SOTField));

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

BError SOTField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SOTField));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new SOTFieldCUDA(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA), this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double SOTField::UpdateField(void)
{
	if (!pMesh->E.linear_size()) return 0.0;

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		double grel = pMesh->grel;
		double Ms = pMesh->Ms;
		double SHA = pMesh->SHA;
		double flSOT = pMesh->flSOT;
		pMesh->update_parameters_mcoarse(idx, pMesh->grel, grel, pMesh->Ms, Ms, pMesh->SHA, SHA, pMesh->flSOT, flSOT);

		if (IsNZ(grel)) {

			double a_const = -(SHA * MUB_E / (GAMMA * grel)) / (Ms * Ms * pMesh->GetMeshDimensions().z);

			int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(idx));
			DBL3 p_vec = (DBL3(0, 0, 1) ^ pMesh->E[idx_E]) * pMesh->elC[idx_E];

			pMesh->Heff[idx] += a_const * ((pMesh->M[idx] ^ p_vec) + flSOT * Ms * p_vec);
		}
	}

	//don't count this as a contribution to the total energy density
	return 0.0;
}

#endif