#include "stdafx.h"
#include "MOptical.h"

#ifdef MODULE_COMPILATION_MOPTICAL

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "MOpticalCUDA.h"
#endif

MOptical::MOptical(Mesh *pMesh_) :
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

MOptical::~MOptical()
{
}

BError MOptical::Initialize(void)
{
	BError error(CLASS_STR(MOptical));

	initialized = true;

	return error;
}

BError MOptical::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MOptical));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();
		Initialize();
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void MOptical::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError MOptical::MakeCUDAModule(void)
{
	BError error(CLASS_STR(MOptical));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new MOpticalCUDA(pMesh->pMeshCUDA);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double MOptical::UpdateField(void)
{
	if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double cHmo = pMesh->cHmo;
			pMesh->update_parameters_mcoarse(idx, pMesh->cHmo, cHmo);

			//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
			pMesh->Heff[idx] += DBL3(0, 0, cHmo);
			pMesh->Heff2[idx] += DBL3(0, 0, cHmo);
		}
	}

	else {

		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double cHmo = pMesh->cHmo;
			pMesh->update_parameters_mcoarse(idx, pMesh->cHmo, cHmo);

			//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
			pMesh->Heff[idx] += DBL3(0, 0, cHmo);
		}
	}

	//no energy density returned
	return 0.0;
}

#endif