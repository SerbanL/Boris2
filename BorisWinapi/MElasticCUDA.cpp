#include "stdafx.h"
#include "MElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_MELASTIC

#include "MElastic.h"
#include "MeshCUDA.h"
#include "Mesh.h"

MElasticCUDA::MElasticCUDA(Mesh* pMesh_, MElastic* pMElastic_) :
	ModulesCUDA()
{
	pMesh = pMesh_;
	pMeshCUDA = pMesh->pMeshCUDA;
	pMElastic = pMElastic_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

MElasticCUDA::~MElasticCUDA()
{
	//Transport managed quantities must be transferred back to their cpu versions - values only as sizes must match
	if (Holder_Module_Available()) {

		pMeshCUDA->u_disp()->copy_to_cpuvec(pMesh->u_disp);
		pMeshCUDA->strain_diag()->copy_to_cpuvec(pMesh->strain_diag);
		pMeshCUDA->strain_odiag()->copy_to_cpuvec(pMesh->strain_odiag);
	}

	//clear mesh quantities as they're not used any more
	pMeshCUDA->u_disp()->clear();
	pMeshCUDA->strain_diag()->clear();
	pMeshCUDA->strain_odiag()->clear();
}

BError MElasticCUDA::Initialize(void)
{
	BError error(CLASS_STR(MElasticCUDA));

	initialized = true;

	return error;
}

BError MElasticCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MElasticCUDA));

	//make sure correct memory is assigned for mechanical quantities

	bool success = true;

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		if (pMeshCUDA->u_disp()->size_cpu().dim()) {

			success &= pMeshCUDA->u_disp()->resize(pMeshCUDA->h_m, pMeshCUDA->meshRect);
		}
		else {

			success &= pMeshCUDA->u_disp()->set_from_cpuvec(pMesh->u_disp);
		}

		//strain tensor - set empty cells using information in u_disp
		if (pMeshCUDA->u_disp()->size_cpu().dim()) {

			success &= pMeshCUDA->strain_diag()->resize(pMeshCUDA->h_m, pMeshCUDA->meshRect, (cuVEC_VC<cuReal3>&)pMeshCUDA->strain_diag);
			success &= pMeshCUDA->strain_odiag()->resize(pMeshCUDA->h_m, pMeshCUDA->meshRect, (cuVEC_VC<cuReal3>&)pMeshCUDA->strain_odiag);
		}
		else {

			success &= pMeshCUDA->strain_diag()->assign(pMeshCUDA->h_m, pMeshCUDA->meshRect, cuReal3(), (cuVEC_VC<cuReal3>&)pMeshCUDA->strain_diag);
			success &= pMeshCUDA->strain_odiag()->assign(pMeshCUDA->h_m, pMeshCUDA->meshRect, cuReal3(), (cuVEC_VC<cuReal3>&)pMeshCUDA->strain_odiag);
		}
	}

	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

void MElasticCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

	}
}

//copy all required mechanical VECs from their cpu versions
BError MElasticCUDA::copy_VECs_to_GPU(void)
{
	BError error(CLASS_STR(MElasticCUDA));

	bool success = true;

	success &= pMeshCUDA->u_disp()->set_from_cpuvec(pMesh->u_disp);
	success &= pMeshCUDA->strain_diag()->set_from_cpuvec(pMesh->strain_diag);
	success &= pMeshCUDA->strain_odiag()->set_from_cpuvec(pMesh->strain_odiag);

	if (!success) return error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif