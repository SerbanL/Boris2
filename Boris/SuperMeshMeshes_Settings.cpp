#include "stdafx.h"
#include "SuperMesh.h"
#include "SimSchedule.h"

//--------------------------------------------------------- MESH HANDLING - SETTINGS

BError SuperMesh::SetMagAngle(std::string meshName, double polar, double azim)
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

	//Need to renormalize magnetization
	error = UpdateConfiguration(UPDATECONFIG_PARAMVALUECHANGED_MLENGTH);

	return error;
}

BError SuperMesh::SetMagAngle_Rect(std::string meshName, double polar, double azim, Rect rectangle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	Rect mesh_relrect = Rect(pMesh[meshName]->GetMeshDimensions());
	if (!mesh_relrect.contains(rectangle)) return error(BERROR_PARAMOUTOFBOUNDS);

	pMesh[meshName]->SetMagAngle(polar, azim, rectangle);

	//Need to renormalize magnetization
	error = UpdateConfiguration(UPDATECONFIG_PARAMVALUECHANGED_MLENGTH);

	return error;
}

//set magnetization angle in given mesh using given shape
BError SuperMesh::SetMagAngle_Shape(std::string meshName, double polar, double azim, std::vector<MeshShape> shapes)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->SetMagAngle_Shape(polar, azim, shapes);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetMagAngle_Shape(polar, azim, shapes);
		}
	}

	//Need to renormalize magnetization
	error = UpdateConfiguration(UPDATECONFIG_PARAMVALUECHANGED_MLENGTH);

	return error;
}

//Set magnetization angle in solid object only containing given relative position uniformly using polar coordinates
BError SuperMesh::SetMagAngle_Object(std::string meshName, double polar, double azim, DBL3 position)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->SetMagAngle_Object(polar, azim, position);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetMagAngle_Object(polar, azim, position);
		}
	}

	//Need to renormalize magnetization
	error = UpdateConfiguration(UPDATECONFIG_PARAMVALUECHANGED_MLENGTH);

	return error;
}

//Flower state magnetization
BError SuperMesh::SetMagFlower(std::string meshName, int direction, DBL3 centre, double radius, double thickness)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMagFlower(direction, centre, radius, thickness);

	return error;
}

//Onion state magnetization
BError SuperMesh::SetMagOnion(std::string meshName, int direction, DBL3 centre, double radius1, double radius2, double thickness)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMagOnion(direction, centre, radius1, radius2, thickness);

	return error;
}

//Crosstie state magnetization
BError SuperMesh::SetMagCrosstie(std::string meshName, int direction, DBL3 centre, double radius, double thickness)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMagCrosstie(direction, centre, radius, thickness);

	return error;
}

//Invert magnetization direction in given mesh (must be magnetic)
BError SuperMesh::SetInvertedMag(std::string meshName, bool x, bool y, bool z)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetInvertedMag(x, y, z);

	return error;
}

//Mirror magnetization in given axis (literal x, y, or z) in given mesh (must be magnetic)
BError SuperMesh::SetMirroredMag(std::string meshName, std::string axis)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMirroredMag(axis);

	return error;
}

//Set random magentisation distribution in given mesh (must be magnetic)
BError SuperMesh::SetRandomMag(std::string meshName, int seed)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetRandomMag(seed);

	return error;
}

BError SuperMesh::SetRandomXYMag(std::string meshName, int seed)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetRandomXYMag(seed);

	return error;
}

BError SuperMesh::SetMagDomainWall(std::string meshName, std::string longitudinal, std::string transverse, double width, double position)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	auto set_component_from_literals = [](std::string literal) -> int {

		if (!literal.length()) return 0;

		std::vector<std::string> literals = { "-z", "-y", "-x", "", "x", "y", "z" };

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
BError SuperMesh::SetSkyrmion(std::string meshName, int orientation, int chirality, double diameter, DBL2 position)
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
BError SuperMesh::SetSkyrmionBloch(std::string meshName, int orientation, int chirality, double diameter, DBL2 position)
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

BError SuperMesh::SetField(std::string meshName, DBL3 field_cartesian)
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

//set prng_seed value in given mesh (all meshes if name empty)
BError SuperMesh::Set_PRNG_Seed(std::string meshName, unsigned seed)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->prng_seed = seed;
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->prng_seed = seed;
		}
	}

	error = UpdateConfiguration(UPDATECONFIG_PRNG);

	return error;
}

//Set uniform applied mechanical stress
BError SuperMesh::SetUniformStress(std::string meshName, DBL3 stress_cartesian)
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
BError SuperMesh::PrepareMovingMesh(std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set dipoles left and right, size 4 times the mesh length
	std::string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](std::string& name) {

		std::string adjusted_name = name;

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

	//6. couple ends to dipoles also
	Set_Coupled_To_Dipoles(true);

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - bloch wall
BError SuperMesh::PrepareMovingMesh_Bloch(std::string meshName)
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
	std::string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](std::string& name) {

		std::string adjusted_name = name;

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

	//6. couple ends to dipoles also
	Set_Coupled_To_Dipoles(true);

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - neel wall
BError SuperMesh::PrepareMovingMesh_Neel(std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set dipoles left and right, size 4 times the mesh length
	std::string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](std::string& name) {

		std::string adjusted_name = name;

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

	//6. couple ends to dipoles also
	Set_Coupled_To_Dipoles(true);

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//fully prepare moving mesh algorithm on named mesh (must be ferromagnetic) - skyrmion
BError SuperMesh::PrepareMovingMesh_Skyrmion(std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	//1. Set dipoles left and right, size 4 times the mesh length
	std::string left_dipole = "dip_left", right_dipole = "dip_right";

	auto adjust_name = [&](std::string& name) {

		std::string adjusted_name = name;

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

	//6. couple ends to dipoles also
	Set_Coupled_To_Dipoles(true);

	if (!error) error = UpdateConfiguration(UPDATECONFIG_GENERIC);

	return error;
}

//Clear moving mesh settings made by a prepare method
void SuperMesh::ClearMovingMesh(void)
{
	//1. Delete all dipole meshes which contain dip_left or dip_right in the name
	std::string left_dipole = "dip_left", right_dipole = "dip_right";

	std::vector<std::string> meshes_to_delete;

	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx]->GetMeshType() == MESH_DIPOLE) {

			std::string meshName = pMesh.get_key_from_index(idx);

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

	//6. uncouple ends from dipoles also
	Set_Coupled_To_Dipoles(false);

	UpdateConfiguration(UPDATECONFIG_GENERIC);
}

//set ferromagnetic mesh roughness refinement if Roughness module enabled in given mesh
BError SuperMesh::SetMeshRoughnessRefinement(std::string meshName, INT3 refine)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	if (!pMesh[meshName]->IsModuleSet(MOD_ROUGHNESS)) return error(BERROR_INCORRECTMODCONFIG);

	pMesh[meshName]->CallModuleMethod(&Roughness::set_refine, refine);

	return error;
}

//Set periodic boundary conditions for magnetization
//possible flags: x, y, z
BError SuperMesh::Set_PBC(std::string meshName, std::string flag, int images)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	auto set_mesh_pbc = [&](std::string meshName, BError& error) -> BError {

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
		else {

			//PBC in individual mesh without Demag module added : just apply PBC to differential operators, e.g. for exchange computation
			INT3 pbc_images = pMesh[meshName]->Get_Magnetic_PBC();

			if (flag == "x") pbc_images.x = images;
			else if (flag == "y") pbc_images.y = images;
			else if (flag == "z") pbc_images.z = images;
			else return error(BERROR_INCORRECTVALUE);

			pMesh[meshName]->Set_Magnetic_PBC(pbc_images);
		}

		return error;
	};

	if (meshName != superMeshHandle) {

		error = set_mesh_pbc(meshName, error);
	}
	else {

		//pbc setting for supermesh demag module if set
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			INT3 pbc_images = dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->Get_PBC();

			if (flag == "x") pbc_images.x = images;
			else if (flag == "y") pbc_images.y = images;
			else if (flag == "z") pbc_images.z = images;
			else return error(BERROR_INCORRECTVALUE);

			error = dynamic_cast<SDemag*>(pSMod(MODS_SDEMAG))->Set_PBC(pbc_images);
		}
		else {

			//if supermesh demag module not set then set pbc for all applicable meshes
			for (int idx = 0; idx < pMesh.size(); idx++) {

				error = set_mesh_pbc(pMesh.get_key_from_index(idx), error);
			}
		}
	}

	return error;
}

//set exchange coupling to neighboring meshes - an exchange-type module (i.e. inherit from ExchangeBase) must be enabled in the named mesh
BError SuperMesh::Set_ExchangeCoupledMeshes(bool status, std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetMeshExchangeCoupling(status);

	return error;
}

//Set/Get multilayered demag exclusion : will need to call UpdateConfiguration when the flag is changed, so the correct SDemag_Demag modules and related settings are set from the SDemag module.
BError SuperMesh::Set_Demag_Exclusion(bool exclude_from_multiconvdemag, std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) || !pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->Set_Demag_Exclusion(exclude_from_multiconvdemag);

	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	return error;
}

//set link_stochastic flag in named mesh, or all meshes if supermesh handle given
BError SuperMesh::SetLinkStochastic(bool link_stochastic, std::string meshName)
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

//set electric field VEC from a constant Jc value in named mesh
BError SuperMesh::SetEFromJcValue(DBL3 Jcvalue, std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (!pMesh[meshName]->EComputation_Enabled()) return error(BERROR_NOTRANSPORT);

	pMesh[meshName]->SetEFromJcValue(Jcvalue);

	return error;
}

//Set TMR type in named mesh (must be an insulator mesh, or leave blank to apply to all meshes)
BError SuperMesh::SetTMRType(std::string meshName, TMR_ TMR_type)
{
	BError error(__FUNCTION__);

	if (meshName == "") {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (pMesh[idx]->GetMeshType() == MESH_INSULATOR) {

				pMesh[idx]->SetTMRType(TMR_type);
			}
		}

		return error;
	}

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (pMesh[meshName]->GetMeshType() != MESH_INSULATOR) return error(BERROR_NOTINSULATOR);

	pMesh[meshName]->SetTMRType(TMR_type);

	return error;
}

//set tyext equation for RAp and RAap in insulator mesh with tmr module added
BError SuperMesh::SetTMR_BiasEquationParallel(std::string meshName, std::string equation_string)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (pMesh[meshName]->GetMeshType() != MESH_INSULATOR) return error(BERROR_NOTINSULATOR);
	if (!pMesh[meshName]->IsModuleSet(MOD_TMR)) return error(BERROR_NOTMR);

	pMesh[meshName]->CallModuleMethod(&TMR::SetBiasEquation_Parallel, equation_string);

	return error;
}

BError SuperMesh::SetTMR_BiasEquationAntiParallel(std::string meshName, std::string equation_string)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (pMesh[meshName]->GetMeshType() != MESH_INSULATOR) return error(BERROR_NOTINSULATOR);
	if (!pMesh[meshName]->IsModuleSet(MOD_TMR)) return error(BERROR_NOTMR);

	pMesh[meshName]->CallModuleMethod(&TMR::SetBiasEquation_AntiParallel, equation_string);

	return error;
}

//set text equation for TAMR conductivity in transport module
BError SuperMesh::Set_TAMR_Conductivity_Equation(std::string meshName, std::string equation_string)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (!pMesh[meshName]->IsModuleSet(MOD_TRANSPORT)) return error(BERROR_NOTRANSPORT);

	pMesh[meshName]->CallModuleMethod(&Transport::Set_TAMR_Conductivity_Equation, equation_string);

	return error;
}

//load field in supermesh (globalField) or in named mesh
BError SuperMesh::LoadOVF2Field(std::string meshName, std::string fileName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle) {

		//field in individual mesh
		if (pMesh[meshName]->Magnetism_Enabled() && pMesh[meshName]->IsModuleSet(MOD_ZEEMAN)) {

			pMesh[meshName]->CallModuleMethod(&ZeemanBase::SetFieldVEC_FromOVF2, fileName);
		}
		else return error(BERROR_INCORRECTMODCONFIG);
	}
	else {

		//global supermesh field
		OVF2 ovf2;
		error = ovf2.Read_OVF2_VEC(fileName, globalField);
		if (error) return error;

#if COMPILECUDA == 1
		if (pSMeshCUDA) {

			if (!(pSMeshCUDA->GetGlobalField())()->set_from_cpuvec(globalField)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
#endif

		//now that global field is set, must update configuration so each Zeeman module can uninitialize (on reinitialization correct settings will be made in Zeeman modules)
		error = UpdateConfiguration(UPDATECONFIG_SMESH_GLOBALFIELD);
	}

	return error;
}

BError SuperMesh::ClearGlobalField(void)
{
	BError error(__FUNCTION__);

	globalField.clear();

#if COMPILECUDA == 1
	if (pSMeshCUDA) {

		(pSMeshCUDA->GetGlobalField())()->clear();
	}
#endif

	error = UpdateConfiguration(UPDATECONFIG_SMESH_GLOBALFIELD);

	return error;
}

//shift globalField rectangle
void SuperMesh::ShiftGlobalField(DBL3 shift)
{
	globalField.shift_rect_start(shift);

#if COMPILECUDA == 1
	if (pSMeshCUDA) {

		(pSMeshCUDA->GetGlobalField())()->shift_rect_start(shift);
	}
#endif

	UpdateConfiguration(UPDATECONFIG_SMESH_GLOBALFIELD);
}

//--------------------------------------------------------- GET/SET PROPERTIES / VALUES at SuperMesh level

//search save data list (saveDataList) for given dataID set for given mesh. Return true if found and its rectangle is not Null (or equal to the mesh rect); else return false.
bool SuperMesh::IsOutputDataSet_withRect(int datumId, MeshBase* pmesh)
{
	//first get mesh name from mesh pointer	
	std::string meshName;
	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx] == pmesh) {

			meshName = pMesh.get_key_from_index(idx);
			break;
		}
	}
	
	if (!meshName.length()) return false;

	//now search output data for match
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		if (saveDataList[idx].datumId == (int)datumId && 
			saveDataList[idx].meshName == meshName &&
			!saveDataList[idx].rectangle.IsNull() &&
			saveDataList[idx].rectangle != pmesh->meshRect) return true;
	}

	return false;
}

//return true if data is set (with any rectangle)
bool SuperMesh::IsOutputDataSet(int datumId, MeshBase* pmesh)
{
	//first get mesh name from mesh pointer	
	std::string meshName;
	for (int idx = 0; idx < pMesh.size(); idx++) {

		if (pMesh[idx] == pmesh) {

			meshName = pMesh.get_key_from_index(idx);
			break;
		}
	}

	if (!meshName.length()) return false;

	//now search output data for match
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		if (saveDataList[idx].datumId == (int)datumId &&
			saveDataList[idx].meshName == meshName) return true;
	}

	return false;
}

//check if given stage is set
bool SuperMesh::IsStageSet(int stageType)
{
	for (int idx = 0; idx < simStages.size(); idx++) {
		
		if (simStages[idx].stage_type() == stageType) return true;
	}

	return false;
}