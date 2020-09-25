#include "stdafx.h"
#include "SuperMesh.h"

//switch to serial (true) or parallel (false) in given mesh - all if meshName is the supermesh handle
BError SuperMesh::Set_MonteCarlo_Serial(bool status, string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	//serial mode only possible with cuda off
	if (cudaEnabled && status == true) return error(BERROR_INCORRECTCONFIG);

	if (meshName == superMeshHandle) {

		//all atomistic meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (pMesh[idx]->is_atomistic()) {

				dynamic_cast<Atom_Mesh*>(pMesh[idx])->Set_MonteCarlo_Serial(status);
			}
		}
	}
	else {

		//named mesh only
		if(!pMesh[meshName]->is_atomistic()) return error(BERROR_INCORRECTCONFIG);

		dynamic_cast<Atom_Mesh*>(pMesh[meshName])->Set_MonteCarlo_Serial(status);
	}

	return error;
}


//switch to constrained Monnte-Carlo (true) or classical (false) in given mesh - all if meshName is the supermesh handle; if constrained, then use cmc_n direction.
BError SuperMesh::Set_MonteCarlo_Constrained(bool status, DBL3 cmc_n, string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//all atomistic meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (pMesh[idx]->is_atomistic()) {

				if (status && !cmc_n.IsNull()) dynamic_cast<Atom_Mesh*>(pMesh[idx])->Set_MonteCarlo_Constrained(cmc_n);
				else dynamic_cast<Atom_Mesh*>(pMesh[idx])->Set_MonteCarlo_Constrained(DBL3());
			}
		}
	}
	else {

		//named mesh only
		if (!pMesh[meshName]->is_atomistic()) return error(BERROR_INCORRECTCONFIG);

		if (status && !cmc_n.IsNull()) dynamic_cast<Atom_Mesh*>(pMesh[meshName])->Set_MonteCarlo_Constrained(cmc_n);
		else dynamic_cast<Atom_Mesh*>(pMesh[meshName])->Set_MonteCarlo_Constrained(DBL3());
	}

	return error;
}
