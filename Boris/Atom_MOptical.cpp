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

	initialized = true;

	return error;
}

BError Atom_MOptical::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
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
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		double cHmo = paMesh->cHmo;
		paMesh->update_parameters_mcoarse(idx, paMesh->cHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		paMesh->Heff1[idx] += DBL3(0, 0, cHmo);
	}

	//no energy density returned
	return 0.0;
}

#endif