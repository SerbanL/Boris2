#include "stdafx.h"
#include "SuperMesh.h"

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

vector<PhysQ> SuperMesh::FetchOnScreenPhysicalQuantity(double detail_level)
{
	vector<PhysQ> physQ;

#if COMPILECUDA == 1
	if (pSMeshCUDA) {

		//if super-mesh display quantities are set with CUDA enabled then get them from PSMeshCUDA
		physQ = pSMeshCUDA->FetchOnScreenPhysicalQuantity(detail_level);
		if (physQ.size()) return physQ;
	}
#endif

	//get anything displayed on super-mesh
	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_SM_DEMAG:

		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			physQ.push_back(PhysQ(&(reinterpret_cast<SDemag*>(pSMod(MODS_SDEMAG))->GetDemagField()), displayedPhysicalQuantity).set_focus(true, superMeshHandle));
			return physQ;
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		if (IsSuperMeshModuleSet(MODS_OERSTED)) {

			physQ.push_back(PhysQ(&(reinterpret_cast<Oersted*>(pSMod(MODS_OERSTED))->GetOerstedField()), displayedPhysicalQuantity).set_focus(true, superMeshHandle));
			return physQ;
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		if (IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			physQ.push_back(PhysQ(&(reinterpret_cast<StrayField*>(pSMod(MODS_STRAYFIELD))->GetStrayField()), displayedPhysicalQuantity).set_focus(true, superMeshHandle));
			return physQ;
		}
		break;
	}

	//get anything displayed in individual meshes
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		physQ.push_back(pMesh[idx]->FetchOnScreenPhysicalQuantity(detail_level).set_focus(pMesh.get_key_from_index(idx) == activeMeshName, pMesh.get_key_from_index(idx)));
	}

	return physQ;
}

BError SuperMesh::SetDisplayedPhysicalQuantity(string meshName, int displayedPhysicalQuantity_)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle) {

		MESH_ meshType = pMesh[meshName]->GetMeshType();

		if (displayedPhysicalQuantity_ >= MESHDISPLAY_NONE && vector_contains(meshAllowedDisplay(meshType), (MESHDISPLAY_)displayedPhysicalQuantity_)) {

			pMesh[meshName]->SetDisplayedPhysicalQuantity(displayedPhysicalQuantity_);
		}
	}
	else {

		if (displayedPhysicalQuantity_ >= MESHDISPLAY_NONE && vector_contains(meshAllowedDisplay(MESH_SUPERMESH), (MESHDISPLAY_)displayedPhysicalQuantity_)) {

			displayedPhysicalQuantity = displayedPhysicalQuantity_;
		}
	}

	return error;
}