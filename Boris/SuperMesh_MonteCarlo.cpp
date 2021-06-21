#include "stdafx.h"
#include "SuperMesh.h"

//switch to serial (true) or parallel (false) in given mesh - all if meshName is the supermesh handle
BError SuperMesh::Set_MonteCarlo_Serial(bool status, std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	//serial mode only possible with cuda off
	if (cudaEnabled && status == true) return error(BERROR_INCORRECTCONFIG);

	if (meshName == superMeshHandle) {

		//all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->Set_MonteCarlo_Serial(status);
		}
	}
	else {

		//named mesh only
		pMesh[meshName]->Set_MonteCarlo_Serial(status);
	}

	return error;
}


//switch to constrained Monnte-Carlo (true) or classical (false) in given mesh - all if meshName is the supermesh handle; if constrained, then use cmc_n direction.
BError SuperMesh::Set_MonteCarlo_Constrained(bool status, DBL3 cmc_n, std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all atomistic meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (status && !cmc_n.IsNull()) pMesh[idx]->Set_MonteCarlo_Constrained(cmc_n);
			else pMesh[idx]->Set_MonteCarlo_Constrained(DBL3());
		}
	}
	else {

		//named mesh only
		if (status && !cmc_n.IsNull()) pMesh[meshName]->Set_MonteCarlo_Constrained(cmc_n);
		else pMesh[meshName]->Set_MonteCarlo_Constrained(DBL3());
	}

	return error;
}

//Disable/enable MC iteration in named mesh
BError SuperMesh::Set_MonteCarlo_Disabled(bool status, std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->Set_MonteCarlo_Disabled(status);

	return error;
}