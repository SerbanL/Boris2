#include "stdafx.h"
#include "AnisotropyTensorialCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANITENS

#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Anisotropy_TensorialCUDA::Anisotropy_TensorialCUDA(MeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Anisotropy_TensorialCUDA::~Anisotropy_TensorialCUDA()
{}

BError Anisotropy_TensorialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_TensorialCUDA));

	if (!initialized) {

		if (!pMeshCUDA->Kt()->assign(cuSZ3(pMeshCUDA->get_tensorial_anisotropy().size(), 1, 1), cuReal4())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!pMeshCUDA->Kt2()->assign(cuSZ3(pMeshCUDA->get_tensorial_anisotropy2().size(), 1, 1), cuReal4())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		if (!pMeshCUDA->Kt()->copy_from_vector(pMeshCUDA->get_tensorial_anisotropy())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!pMeshCUDA->Kt2()->copy_from_vector(pMeshCUDA->get_tensorial_anisotropy2())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect,
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ANITENS || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ANITENS || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (error)	initialized = false;

	return error;
}

BError Anisotropy_TensorialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_TensorialCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_PARAMCHANGED)) {

		Uninitialize();
	}

	return error;
}

#endif

#endif