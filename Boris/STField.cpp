#include "stdafx.h"
#include "STField.h"

#ifdef MODULE_COMPILATION_STFIELD

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "STFieldCUDA.h"
#endif

STField::STField(Mesh *pMesh_) :
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

STField::~STField()
{
}

BError STField::Initialize(void)
{
	BError error(CLASS_STR(STField));

	initialized = true;

	return error;
}

BError STField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(STField));

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

BError STField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(STField));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new STFieldCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double STField::UpdateField(void)
{
	if (!pMesh->E.linear_size()) return 0.0;
		
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		double grel = pMesh->grel;
		double Ms = pMesh->Ms;
		double flSOT = pMesh->flSOT;
		DBL2 STq = pMesh->STq;
		DBL2 STa = pMesh->STa;
		DBL3 STp = pMesh->STp;
		pMesh->update_parameters_mcoarse(idx, pMesh->grel, grel, pMesh->Ms, Ms, pMesh->flSOT, flSOT, pMesh->STq, STq, pMesh->STa, STa, pMesh->STp, STp);

		if (IsNZ(grel)) {

			double dotprod = (pMesh->M[idx] * STp) / Ms;
			double neta = STq.i / (STa.i + STq.j * dotprod) + STq.j / (STa.i - STq.j * dotprod);

			int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(idx));
			//z component of J
			double Jc = pMesh->E[idx_E].z * pMesh->elC[idx_E];

			double a_const = -(neta * MUB_E * Jc / (GAMMA * grel)) / (Ms * Ms * pMesh->GetMeshDimensions().z);

			pMesh->Heff[idx] += a_const * ((pMesh->M[idx] ^ STp) + flSOT * Ms * STp);
		}
	}

	//don't count this as a contribution to the total energy density
	return 0.0;
}

#endif