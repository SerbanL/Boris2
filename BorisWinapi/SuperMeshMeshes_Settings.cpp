#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- MESH HANDLING - SETTINGS

BError SuperMesh::SetMagAngle(string meshName, double polar, double azim)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->SetMagAngle(polar, azim);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetMagAngle(polar, azim);
		}
	}

	return error;
}

BError SuperMesh::SetMagAngle_Rect(string meshName, double polar, double azim, Rect rectangle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	Rect mesh_relrect = Rect(pMesh[meshName]->GetMeshDimensions());
	if (!mesh_relrect.contains(rectangle)) return error(BERROR_PARAMOUTOFBOUNDS);

	pMesh[meshName]->SetMagAngle(polar, azim, rectangle);

	return error;
}

//Invert magnetisation direction in given mesh (must be magnetic)
BError SuperMesh::SetInvertedMag(string meshName, bool x, bool y, bool z)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetInvertedMag(x, y, z);

	return error;
}

//Mirror magnetisation in given axis (literal x, y, or z) in given mesh (must be magnetic)
BError SuperMesh::SetMirroredMag(string meshName, string axis)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMirroredMag(axis);

	return error;
}

//Set random magentisation distribution in given mesh (must be magnetic)
BError SuperMesh::SetRandomMag(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetRandomMag();

	return error;
}

BError SuperMesh::SetMagDomainWall(string meshName, string longitudinal, string transverse, double width, double position)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	auto set_component_from_literals = [](string literal) -> int {

		if (!literal.length()) return 0;

		vector<string> literals = { "-z", "-y", "-x", "", "x", "y", "z" };

		for (int idx = 0; idx < (int)literals.size(); idx++)
			if (literal == literals[idx]) return (idx - 3);

		return 0;
	};

	int longi = set_component_from_literals(longitudinal), trans = set_component_from_literals(transverse);

	if (!longi || !trans) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMagDomainWall(longi, trans, width, position);

	return error;
}

//Set Neel skyrmion in given mesh with chirality, diameter and relative x-y plane position
BError SuperMesh::SetSkyrmion(string meshName, int orientation, int chirality, double diameter, DBL2 position)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	DBL3 mesh_dim = pMesh[meshName]->GetMeshDimensions();

	Rect skyrmion_relrect = Rect(DBL3(position.x - diameter / 2, position.y - diameter / 2, 0.0), DBL3(position.x + diameter / 2, position.y + diameter / 2, mesh_dim.z));

	Rect mesh_relrect = pMesh[meshName]->GetMeshRect() - pMesh[meshName]->GetOrigin();

	//the skyrmion rectangle must be contained in the mesh
	if (!mesh_relrect.contains(skyrmion_relrect)) return error(BERROR_PARAMOUTOFBOUNDS);

	pMesh[meshName]->SetSkyrmion(orientation, chirality, skyrmion_relrect);

	return error;
}

//Set Bloch skyrmion in given mesh with chirality, diameter and relative x-y plane position
BError SuperMesh::SetSkyrmionBloch(string meshName, int orientation, int chirality, double diameter, DBL2 position)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	DBL3 mesh_dim = pMesh[meshName]->GetMeshDimensions();

	Rect skyrmion_relrect = Rect(DBL3(position.x - diameter / 2, position.y - diameter / 2, 0.0), DBL3(position.x + diameter / 2, position.y + diameter / 2, mesh_dim.z));

	Rect mesh_relrect = pMesh[meshName]->GetMeshRect() - pMesh[meshName]->GetOrigin();

	//the skyrmion rectangle must be contained in the mesh
	if (!mesh_relrect.contains(skyrmion_relrect)) return error(BERROR_PARAMOUTOFBOUNDS);

	pMesh[meshName]->SetSkyrmionBloch(orientation, chirality, skyrmion_relrect);

	return error;
}

BError SuperMesh::SetField(string meshName, DBL3 field_cartesian)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->CallModuleMethod(&ZeemanBase::SetField, field_cartesian);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->CallModuleMethod(&ZeemanBase::SetField, field_cartesian);
		}
	}

	return error;
}

//Set uniform applied mechanical stress
BError SuperMesh::SetUniformStress(string meshName, DBL3 stress_cartesian)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->CallModuleMethod(&MElastic::SetUniformStress, stress_cartesian);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->CallModuleMethod(&MElastic::SetUniformStress, stress_cartesian);
		}
	}

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - transverse wall
BError SuperMesh::PrepareMovingMesh(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set dipoles left and right, size 4 times the mesh length
	string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](string& name) {

		string adjusted_name = name;

		int num = 1;
		while (contains(adjusted_name)) {

			adjusted_name = name + ToString(num++);
		}

		name = adjusted_name;
	};

	adjust_name(left_dipole);
	adjust_name(right_dipole);

	Rect meshRect = pMesh[meshName]->GetMeshRect();
	DBL3 sizes = meshRect.size();

	Rect meshRect_left = Rect(meshRect.s - (sizes & DBL3(1, 0, 0)) * 4, meshRect.s + (sizes & DBL3(0, 1, 1)));
	Rect meshRect_right = Rect(meshRect.e - (sizes & DBL3(0, 1, 1)), meshRect.e + (sizes & DBL3(1, 0, 0)) * 4);

	//check for intersection with any other meshes : abort if any found. Do not use IsNZ !!! volume is well below the epsilon so will not work correctly.
	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshRect().intersection_volume(meshRect_left) || pMesh[idx]->GetMeshRect().intersection_volume(meshRect_right))
			return error(BERROR_INCORRECTACTION);
	}

	if (!error) error = AddMesh(left_dipole, MESH_DIPOLE, meshRect_left);
	if (!error) error = SetMagAngle(left_dipole, 90.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagAngle(right_dipole, 90.0, 180.0);

	//2. Set moving mesh trigger, antisymmetric type and threshold
	SetMoveMeshTrigger(true, meshName);
	SetMoveMeshAntisymmetric(true);
	SetMoveMeshThreshold(MOVEMESH_ANTISYMMETRIC_THRESHOLD);

	//3. Set transverse domain wall structure along the x direction. Set the ends along x direction.
	SetMagDomainWall(meshName, "x", "y", pMesh[meshName]->GetMeshDimensions().x, 0.0);

	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	pMesh[meshName]->SetMagAngle(90, 0, Rect(DBL3(0), DBL3(MOVEMESH_ENDRATIO * dim.x, dim.y, dim.z)));
	pMesh[meshName]->SetMagAngle(90, 180, Rect(DBL3((1 - MOVEMESH_ENDRATIO) * dim.x, 0, 0), dim));

	//4. add stray field module
	if (!error) error = AddModule(superMeshHandle, MODS_STRAYFIELD);

	//5. set scalemeshrects to true
	scale_rects = true;

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - bloch wall
BError SuperMesh::PrepareMovingMesh_Bloch(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set moving mesh trigger, antisymmetric type and threshold
	SetMoveMeshTrigger(true, meshName);
	SetMoveMeshAntisymmetric(true);
	SetMoveMeshThreshold(MOVEMESH_ANTISYMMETRIC_THRESHOLD);

	//2. Set Bloch domain wall structure along the x direction. Set the ends along x direction.
	SetMagDomainWall(meshName, "z", "y", pMesh[meshName]->GetMeshDimensions().x, 0.0);

	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	pMesh[meshName]->SetMagAngle(0, 0, Rect(DBL3(0), DBL3(MOVEMESH_ENDRATIO * dim.x, dim.y, dim.z)));
	pMesh[meshName]->SetMagAngle(180, 0, Rect(DBL3((1 - MOVEMESH_ENDRATIO) * dim.x, 0, 0), dim));

	//3. Set dipoles left and right, size 4 times the mesh length
	string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](string& name) {

		string adjusted_name = name;

		int num = 1;
		while (contains(adjusted_name)) {

			adjusted_name = name + ToString(num++);
		}

		name = adjusted_name;
	};

	adjust_name(left_dipole);
	adjust_name(right_dipole);

	Rect meshRect = pMesh[meshName]->GetMeshRect();
	DBL3 sizes = meshRect.size();

	Rect meshRect_left = Rect(meshRect.s - (sizes & DBL3(1, 0, 0)) * 4, meshRect.s + (sizes & DBL3(0, 1, 1)));
	Rect meshRect_right = Rect(meshRect.e - (sizes & DBL3(0, 1, 1)), meshRect.e + (sizes & DBL3(1, 0, 0)) * 4);

	//check for intersection with any other meshes : abort if any found. Do not use IsNZ !!! volume is well below the epsilon so will not work correctly.
	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshRect().intersection_volume(meshRect_left) || pMesh[idx]->GetMeshRect().intersection_volume(meshRect_right))
			return error(BERROR_INCORRECTACTION);
	}

	if (!error) error = AddMesh(left_dipole, MESH_DIPOLE, meshRect_left);
	if (!error) error = SetMagAngle(left_dipole, 0.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagAngle(right_dipole, 180.0, 0.0);

	//4. add stray field module
	if (!error) error = AddModule(superMeshHandle, MODS_STRAYFIELD);

	//5. set scalemeshrects to true
	scale_rects = true;

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - neel wall
BError SuperMesh::PrepareMovingMesh_Neel(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set dipoles left and right, size 4 times the mesh length
	string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](string& name) {

		string adjusted_name = name;

		int num = 1;
		while (contains(adjusted_name)) {

			adjusted_name = name + ToString(num++);
		}

		name = adjusted_name;
	};

	adjust_name(left_dipole);
	adjust_name(right_dipole);

	Rect meshRect = pMesh[meshName]->GetMeshRect();
	DBL3 sizes = meshRect.size();

	Rect meshRect_left = Rect(meshRect.s - (sizes & DBL3(1, 0, 0)) * 4, meshRect.s + (sizes & DBL3(0, 1, 1)));
	Rect meshRect_right = Rect(meshRect.e - (sizes & DBL3(0, 1, 1)), meshRect.e + (sizes & DBL3(1, 0, 0)) * 4);

	//check for intersection with any other meshes : abort if any found. Do not use IsNZ !!! volume is well below the epsilon so will not work correctly.
	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshRect().intersection_volume(meshRect_left) || pMesh[idx]->GetMeshRect().intersection_volume(meshRect_right))
			return error(BERROR_INCORRECTACTION);
	}

	if (!error) error = AddMesh(left_dipole, MESH_DIPOLE, meshRect_left);
	if (!error) error = SetMagAngle(left_dipole, 0.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagAngle(right_dipole, 180.0, 0.0);

	//2. Set moving mesh trigger, antisymmetric type and threshold
	SetMoveMeshTrigger(true, meshName);
	SetMoveMeshAntisymmetric(true);
	SetMoveMeshThreshold(MOVEMESH_ANTISYMMETRIC_THRESHOLD);

	//3. Set Bloch domain wall structure along the x direction. Set the ends along x direction.
	SetMagDomainWall(meshName, "z", "x", pMesh[meshName]->GetMeshDimensions().x, 0.0);

	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	pMesh[meshName]->SetMagAngle(0, 0, Rect(DBL3(0), DBL3(MOVEMESH_ENDRATIO * dim.x, dim.y, dim.z)));
	pMesh[meshName]->SetMagAngle(180, 0, Rect(DBL3((1 - MOVEMESH_ENDRATIO) * dim.x, 0, 0), dim));

	//4. add stray field module
	if (!error) error = AddModule(superMeshHandle, MODS_STRAYFIELD);

	//5. set scalemeshrects to true
	scale_rects = true;

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - skyrmion
BError SuperMesh::PrepareMovingMesh_Skyrmion(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set dipoles left and right, size 4 times the mesh length
	string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](string& name) {

		string adjusted_name = name;

		int num = 1;
		while (contains(adjusted_name)) {

			adjusted_name = name + ToString(num++);
		}

		name = adjusted_name;
	};

	adjust_name(left_dipole);
	adjust_name(right_dipole);

	Rect meshRect = pMesh[meshName]->GetMeshRect();
	DBL3 sizes = meshRect.size();

	Rect meshRect_left = Rect(meshRect.s - (sizes & DBL3(1, 0, 0)) * 4, meshRect.s + (sizes & DBL3(0, 1, 1)));
	Rect meshRect_right = Rect(meshRect.e - (sizes & DBL3(0, 1, 1)), meshRect.e + (sizes & DBL3(1, 0, 0)) * 4);

	//check for intersection with any other meshes : abort if any found. Do not use IsNZ !!! volume is well below the epsilon so will not work correctly.
	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshRect().intersection_volume(meshRect_left) || pMesh[idx]->GetMeshRect().intersection_volume(meshRect_right))
			return error(BERROR_INCORRECTACTION);
	}

	if (!error) error = AddMesh(left_dipole, MESH_DIPOLE, meshRect_left);
	if (!error) error = SetMagAngle(left_dipole, 0.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagAngle(right_dipole, 0.0, 0.0);

	//2. Set moving mesh trigger, symmetric type and threshold
	SetMoveMeshTrigger(true, meshName);
	SetMoveMeshAntisymmetric(true);
	SetMoveMeshThreshold(MOVEMESH_SYMMETRIC_THRESHOLD);

	//3. Set skyrmion in the centre amd up direction everywhere else
	pMesh[meshName]->SetMagAngle(0, 0);

	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	Rect sky_rect = Rect(DBL3(dim.x / 2 - dim.y / 4, dim.y / 4, 0), DBL3(dim.x / 2 + dim.y / 4, dim.y * 3 / 4, dim.z));
	pMesh[meshName]->SetSkyrmion(-1, 1, sky_rect);

	//4. add stray field module
	if (!error) error = AddModule(superMeshHandle, MODS_STRAYFIELD);

	//5. set scalemeshrects to true
	scale_rects = true;

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//Clear moving mesh settings made by a prepare method
void SuperMesh::ClearMovingMesh(void)
{
	//1. Delete all dipole meshes which contain dip_left or dip_right in the name
	string left_dipole = "dip_left", right_dipole = "dip_right";

	vector<string> meshes_to_delete;

	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshType() == MESH_DIPOLE) {

			string meshName = pMesh.get_key_from_index(idx);

			if (meshName.find(left_dipole) != std::string::npos || meshName.find(right_dipole) != std::string::npos) {

				meshes_to_delete.push_back(meshName);
			}
		}
	}

	for (int idx = 0; idx < meshes_to_delete.size(); idx++) {

		DelMesh(meshes_to_delete[idx]);
	}

	//2. Disable moving mesh trigger
	SetMoveMeshTrigger(false);

	//3. N/A

	//4. Delete Stray Field module
	DelModule(superMeshHandle, MODS_STRAYFIELD);

	//5. set scalemeshrects to false
	scale_rects = false;

	UpdateConfiguration(UPDATECONFIG_GENERIC);
}

//set ferromagnetic mesh roughness refinement if Roughness module enabled in given mesh
BError SuperMesh::SetMeshRoughnessRefinement(string meshName, INT3 refine)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	if (!pMesh[meshName]->IsModuleSet(MOD_ROUGHNESS)) return error(BERROR_INCORRECTMODCONFIG);

	pMesh[meshName]->CallModuleMethod(&Roughness::set_refine, refine);

	return error;
}

//Set periodic boundary conditions for magnetization
//possible flags: x, y, z
BError SuperMesh::Set_PBC(string meshName, string flag, int images)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle) {

		if (!pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

		//pbc setting for individual mesh demag module

		if (pMesh[meshName]->Is_Demag_Enabled()) {

			INT3 pbc_images = pMesh[meshName]->CallModuleMethod(&DemagBase::Get_PBC);

			if (flag == "x") pbc_images.x = images;
			else if (flag == "y") pbc_images.y = images;
			else if (flag == "z") pbc_images.z = images;
			else return error(BERROR_INCORRECTVALUE);

			pMesh[meshName]->CallModuleMethod(&DemagBase::Set_PBC, pbc_images);
		}
		else return error(BERROR_INCORRECTNAME);
	}
	else {

		//pbc setting for supermesh demag module
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			INT3 pbc_images = dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->Get_PBC();

			if (flag == "x") pbc_images.x = images;
			else if (flag == "y") pbc_images.y = images;
			else if (flag == "z") pbc_images.z = images;
			else return error(BERROR_INCORRECTVALUE);

			error = dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->Set_PBC(pbc_images);
		}
		else return error(BERROR_INCORRECTNAME);
	}

	return error;
}

//set exchange coupling to neighboring meshes - an exchange-type module (i.e. inherit from ExchangeBase) must be enabled in the named mesh
BError SuperMesh::Set_ExchangeCoupledMeshes(bool status, string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMeshExchangeCoupling(status);

	return error;
}

//Set/Get multilayered demag exclusion : will need to call UpdateConfiguration when the flag is changed, so the correct SDemag_Demag modules and related settings are set from the SDemag module.
BError SuperMesh::Set_Demag_Exclusion(bool exclude_from_multiconvdemag, string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->Set_Demag_Exclusion(exclude_from_multiconvdemag);

	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	return error;
}

//set link_stochastic flag in named mesh, or all meshes if supermesh handle given
BError SuperMesh::SetLinkStochastic(bool link_stochastic, string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle) {

		if (!pMesh[meshName]->MComputation_Enabled() || pMesh[meshName]->is_atomistic()) return error(BERROR_INCORRECTNAME);

		//only applicable for micromagnetic meshes; for atomistic meshes stochasticity is applied at the individual atomic moment level
		error = dynamic_cast<Mesh*>(pMesh[meshName])->SetLinkStochastic(link_stochastic);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (pMesh[idx]->MComputation_Enabled() && !pMesh[meshName]->is_atomistic()) error = dynamic_cast<Mesh*>(pMesh[idx])->SetLinkStochastic(link_stochastic);
		}
	}

	return error;
}