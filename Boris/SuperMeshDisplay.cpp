#include "stdafx.h"
#include "SuperMesh.h"

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

std::vector<PhysQ> SuperMesh::FetchOnScreenPhysicalQuantity(double detail_level)
{
	std::vector<PhysQ> physQ;

	bool cudaSupermesh = false;

#if COMPILECUDA == 1

	//Commands are executed on newly spawned threads, so if cuda is on and we are not using device 0 (default device) we must switch to required device, otherwise 0 will be used
	if (cudaEnabled && cudaDeviceSelect != 0) {

		int device = 0;
		cudaGetDevice(&device);
		if (device != cudaDeviceSelect) cudaSetDevice(cudaDeviceSelect);
	}

	if (pSMeshCUDA) {

		//if super-mesh display quantities are set with CUDA enabled then get them from PSMeshCUDA
		physQ = pSMeshCUDA->FetchOnScreenPhysicalQuantity(detail_level);
		cudaSupermesh = true;
	}
#endif

	//get anything displayed on super-mesh
	if (!cudaSupermesh) {

		switch (displayedPhysicalQuantity) {

		case MESHDISPLAY_NONE:
			break;

		case MESHDISPLAY_SM_DEMAG:

			if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

				physQ.push_back(PhysQ(&(dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->GetDemagField()), displayedPhysicalQuantity, (VEC3REP_)vec3rep).set_focus(true, superMeshHandle));
			}
			break;

		case MESHDISPLAY_SM_OERSTED:

			if (IsSuperMeshModuleSet(MODS_OERSTED)) {

				physQ.push_back(PhysQ(&(dynamic_cast<Oersted*>(pSMod(MODS_OERSTED))->GetOerstedField()), displayedPhysicalQuantity, (VEC3REP_)vec3rep).set_focus(true, superMeshHandle));
			}
			break;

		case MESHDISPLAY_SM_STRAYH:

			if (IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

				physQ.push_back(PhysQ(&(dynamic_cast<StrayField*>(pSMod(MODS_STRAYFIELD))->GetStrayField()), displayedPhysicalQuantity, (VEC3REP_)vec3rep).set_focus(true, superMeshHandle));
			}
			break;
		}
	}

	bool supermeshDisplay = physQ.size();
	//allow dual display with supermesh display as the foreground - need to check if any individual meshes have a background display enabled.
	if (supermeshDisplay) physQ.back().set_display_props(displayTransparency.i, displayThresholds, displayThresholdTrigger);

	//get anything displayed in individual meshes
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if (pMesh[idx]->IsDisplayBackgroundEnabled()) {

			//get foreground and set required transparency and thresholds : not applicable if supermesh quantity is in the foreground
			if (!supermeshDisplay) physQ.push_back(pMesh[idx]->FetchOnScreenPhysicalQuantity(detail_level).set_focus(pMesh.get_key_from_index(idx) == activeMeshName, pMesh.get_key_from_index(idx)).set_display_props(displayTransparency.i, displayThresholds, displayThresholdTrigger));

			//get background and set required transparency (thresholds do not apply)
			physQ.push_back(pMesh[idx]->FetchOnScreenPhysicalQuantity(detail_level, true).set_focus(false, pMesh.get_key_from_index(idx)).set_transparency(displayTransparency.j));
		}
		else if (!supermeshDisplay) {

			//get quantity and set thresholds (transparency does not apply) : not applicable if supermesh quantity is in the foreground
			physQ.push_back(pMesh[idx]->FetchOnScreenPhysicalQuantity(detail_level).set_focus(pMesh.get_key_from_index(idx) == activeMeshName, pMesh.get_key_from_index(idx)).set_thresholds(displayThresholds, displayThresholdTrigger));
		}
	}

	return physQ;
}

//save the quantity currently displayed on screen in an ovf2 file using the specified format
BError SuperMesh::SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType)
{
#if COMPILECUDA == 1
	if (pSMeshCUDA) { return pSMeshCUDA->SaveOnScreenPhysicalQuantity(fileName, ovf2_dataType); }
#endif

	BError error(__FUNCTION__);

	OVF2 ovf2;

	//get anything displayed on super-mesh
	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		//get anything displayed in active mesh
		error = active_mesh()->SaveOnScreenPhysicalQuantity(fileName, ovf2_dataType);
		break;

	case MESHDISPLAY_SM_DEMAG:

		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->GetDemagField(), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		if (IsSuperMeshModuleSet(MODS_OERSTED)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Oersted*>(pSMod(MODS_OERSTED))->GetOerstedField(), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		if (IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<StrayField*>(pSMod(MODS_STRAYFIELD))->GetStrayField(), ovf2_dataType);
		}
		break;
	}

	return error;
}

//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
void SuperMesh::PrepareDisplayedMeshValue(void)
{
	switch (displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
	{
		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->PrepareDisplayedMeshValue();
		}
	}
	break;

	default:

		////////////////
		//use a quantity displayed on the supermesh
		////////////////

#if COMPILECUDA == 1
		if (pSMeshCUDA) {

			pSMeshCUDA->PrepareDisplayedMeshValue();
		}
#endif
		break;
	}
}

//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
//Return an Any as the displayed quantity could be either a scalar or a vector.
Any SuperMesh::GetDisplayedMeshValue(DBL3 abs_pos)
{
#if COMPILECUDA == 1
	if (pSMeshCUDA) {

		//if super-mesh display quantities are set with CUDA enabled then get value from pSMeshCUDA
		return pSMeshCUDA->GetDisplayedMeshValue(abs_pos);
	}
#endif
	
	//get anything displayed on super-mesh
	switch (displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
	{
		//find which mesh holds abs_pos, if any, and return value displayed for that mesh
		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (pMesh[idx]->meshRect.contains(abs_pos)) return pMesh[idx]->GetDisplayedMeshValue(abs_pos);
		}
	}
		break;

		////////////////
		//use a quantity displayed on the supermesh
		////////////////

	case MESHDISPLAY_SM_DEMAG:

		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->GetDemagField()[abs_pos - sMeshRect_fm.s];
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		if (IsSuperMeshModuleSet(MODS_OERSTED)) {

			return dynamic_cast<Oersted*>(pSMod(MODS_OERSTED))->GetOerstedField()[abs_pos - sMeshRect_e.s];
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		if (IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			return dynamic_cast<StrayField*>(pSMod(MODS_STRAYFIELD))->GetStrayField()[abs_pos - sMeshRect_fm.s];
		}
		break;
	}

	return Any();
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any  SuperMesh::GetAverageDisplayedMeshValue(Rect rel_rect)
{
#if COMPILECUDA == 1
	if (pSMeshCUDA) {

		//if super-mesh display quantities are set with CUDA enabled then get value from pSMeshCUDA
		return pSMeshCUDA->GetAverageDisplayedMeshValue(rel_rect);
	}
#endif

	//get anything displayed on super-mesh
	switch (displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
	{
		//get average value from currently focused mesh instead
		return active_mesh()->GetAverageDisplayedMeshValue(rel_rect);
	}
	break;

	////////////////
	//use a quantity displayed on the supermesh
	////////////////

	case MESHDISPLAY_SM_DEMAG:

		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->GetDemagField().average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		if (IsSuperMeshModuleSet(MODS_OERSTED)) {

			return dynamic_cast<Oersted*>(pSMod(MODS_OERSTED))->GetOerstedField().average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		if (IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			return dynamic_cast<StrayField*>(pSMod(MODS_STRAYFIELD))->GetStrayField().average_nonempty_omp(rel_rect);
		}
		break;
	}

	return Any();
}

BError SuperMesh::SetDisplayedPhysicalQuantity(std::string meshName, int displayedPhysicalQuantity_)
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

BError SuperMesh::SetDisplayedBackgroundPhysicalQuantity(std::string meshName, int displayedBackgroundPhysicalQuantity_)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	MESH_ meshType = pMesh[meshName]->GetMeshType();

	if (displayedBackgroundPhysicalQuantity_ >= MESHDISPLAY_NONE && vector_contains(meshAllowedDisplay(meshType), (MESHDISPLAY_)displayedBackgroundPhysicalQuantity_)) {

		pMesh[meshName]->SetDisplayedBackgroundPhysicalQuantity(displayedBackgroundPhysicalQuantity_);
	}

	return error;
}

//Get/Set vectorial quantity representation options in named mesh (which could be the supermesh)
BError SuperMesh::SetVEC3Rep(std::string meshName, int vec3rep_)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle) {

		 pMesh[meshName]->SetVEC3Rep(vec3rep_);
	}
	else {

		vec3rep = vec3rep_;
	}

	return error;
}

int SuperMesh::GetVEC3Rep(std::string meshName)
{
	if (!contains(meshName) && meshName != superMeshHandle) return vec3rep;

	if (meshName != superMeshHandle) {

		return pMesh[meshName]->GetVEC3Rep();
	}
	else {

		return vec3rep;
	}
}