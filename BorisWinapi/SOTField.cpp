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

	error_on_create = UpdateConfiguration();

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

		if (!error) error = pModuleCUDA->UpdateConfiguration();
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

void SOTField::UpdateField(void)
{
	if (!pMesh->Jc.linear_size()) return;

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		double Ms = pMesh->Ms;
		double SHA = pMesh->SHA;
		double flSOT = pMesh->flSOT;
		pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->SHA, SHA, pMesh->flSOT, flSOT);

		double a_const = -(HBAR_E / MU0) * SHA / (2 * Ms * Ms * pMesh->GetMeshDimensions().z);

		int idx_Jc = pMesh->Jc.position_to_cellidx(pMesh->M.cellidx_to_position(idx));
		DBL3 p_vec = DBL3(0, 0, 1) ^ pMesh->Jc[idx_Jc];

		pMesh->Heff[idx] += a_const * ((pMesh->M[idx] ^ p_vec) + flSOT * Ms * p_vec);
	}
}

#endif