#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- MESH HANDLING - COMPONENTS

BError SuperMesh::AddMesh(string meshName, MESH_ meshType, Rect meshRect)
{
	BError error(__FUNCTION__);

	if (contains(meshName) || meshName == superMeshHandle) return error(BERROR_INCORRECTNAME);

	//when creating a new mesh the user supplies the mesh type, rectangle and name.
	//we also need a cellsize. The user can change the cellsize later but it's easier not to ask for it when the mesh is created - use a default value instead.
	//To choose default cellsize : use 5nm cubic cell up to a given number of maximum cells in each dimension. After this increase starting cellsize to keep to this maximum number of cells.
	//The user may be trying to create a large mesh for which there might not be enough memory with a 5nm cellsize, but intending to set a larger cellsize after mesh creation.

	DBL3 cellsize(DEFAULTCELLSIZE);

	INT3 cells = round(meshRect / cellsize);

	if (cells.x > MAXSTARTINGCELLS_X) cells.x = MAXSTARTINGCELLS_X;
	if (cells.y > MAXSTARTINGCELLS_Y) cells.y = MAXSTARTINGCELLS_Y;
	if (cells.z > MAXSTARTINGCELLS_Z) cells.z = MAXSTARTINGCELLS_Z;

	//adjusted cellsize which will result in an integer number of cells with upper limits set
	cellsize = meshRect / cells;

	switch (meshType) {

	case MESH_FERROMAGNETIC:
		pMesh.push_back(new FMesh(meshRect, cellsize, this), meshName);
		break;

	case MESH_ANTIFERROMAGNETIC:
		pMesh.push_back(new AFMesh(meshRect, cellsize, this), meshName);
		break;

	case MESH_DIAMAGNETIC:
		pMesh.push_back(new DiaMesh(meshRect, cellsize, this), meshName);
		break;

	case MESH_DIPOLE:
		pMesh.push_back(new DipoleMesh(meshRect, cellsize, this), meshName);
		break;

	case MESH_METAL:
		pMesh.push_back(new MetalMesh(meshRect, cellsize, this), meshName);
		break;

	case MESH_INSULATOR:
		pMesh.push_back(new InsulatorMesh(meshRect, cellsize, this), meshName);
		break;
	}

	//check the mesh was created correctly - if not, delete it
	error = pMesh.back()->Error_On_Create();

	if (!error) error = UpdateConfiguration(UPDATECONFIG_MESHADDED);
	else pMesh.pop_back();

	return error;
}

BError SuperMesh::DelMesh(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	//cannot delete last mesh
	if (pMesh.size() == 1) return error(BERROR_INCORRECTACTION);

	//delete allocated memory then erase it from list of meshes
	if(pMesh[meshName]) delete pMesh[meshName];
	pMesh.erase(meshName);

	if (activeMeshName == meshName) activeMeshName = pMesh.get_key_from_index(0);

	error = UpdateConfiguration(UPDATECONFIG_MESHDELETED);

	return error;
}

BError SuperMesh::RenameMesh(string oldName, string newName)
{
	BError error(__FUNCTION__);

	//must have old name and new name cannot be that of an existing mesh
	if (!contains(oldName) || contains(newName) || newName == superMeshHandle) return error(BERROR_INCORRECTNAME);

	//make changes
	pMesh.change_key(oldName, newName);

	if (activeMeshName == oldName) activeMeshName = newName;

	return error;
}

BError SuperMesh::SetMeshFocus(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle)
		activeMeshName = meshName;

	return error;
}

//set mesh rect for named mesh (any if applicable any other dependent meshes) and update dependent save data rects by calling the provided function.
BError SuperMesh::SetMeshRect(string meshName, Rect meshRect, std::function<void(string, Rect)> save_data_updater)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);
	
	//cannot set a plane, must be a proper 3D rect
	if (meshRect.IsPlane() || meshRect.e <= meshRect.s) return error(BERROR_INCORRECTVALUE);

	//snap mesh rectangle coordinates to nearest angstrom
	meshRect.snap(1e-10);

	if (!scale_rects) {

		Rect meshRect_old = pMesh[meshName]->GetMeshRect();

		//scale just the rectangle in named mesh
		error = pMesh[meshName]->SetMeshRect(meshRect);
		if (!error) save_data_updater(meshName, meshRect_old);
	}
	else {

		Rect meshRect_old = pMesh[meshName]->GetMeshRect();

		DBL3 old_sizes = meshRect_old.size();
		DBL3 new_sizes = meshRect.size();

		//scale all rectangles using this
		DBL3 ratios = new_sizes / old_sizes;

		//shift rectangles after scaling using this
		DBL3 origin = meshRect_old.s;
		DBL3 shift = meshRect.s - origin;

		//scale all mesh rectangles
		for (int idx = 0; idx < pMesh.size(); idx++) {

			Rect meshRect_original = pMesh[idx]->GetMeshRect();
			
			Rect meshRect_scaled_shifted = meshRect_original;
			meshRect_scaled_shifted.s = ((meshRect_scaled_shifted.s - origin) & ratios) + origin + shift;
			meshRect_scaled_shifted.e = ((meshRect_scaled_shifted.e - origin) & ratios) + origin + shift;

			//snap mesh rectangle coordinates to nearest angstrom
			meshRect_scaled_shifted.snap(1e-10);

			error = pMesh[idx]->SetMeshRect(meshRect_scaled_shifted);
			if (!error) save_data_updater(pMesh.get_key_from_index(idx), meshRect_original);
			else return error;
		}
	}

	return error;
}

//--------------------------------------------------------- MESH HANDLING - SHAPES

//copy all primary mesh data (magnetisation, elC, Temp, etc.) but do not change dimensions or discretisation
BError SuperMesh::copy_mesh_data(string meshName_from, string meshName_to)
{
	BError error(__FUNCTION__);

	if (!contains(meshName_from) || !contains(meshName_to)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName_to]->copy_mesh_data(*pMesh[meshName_from]);

	return error;
}

BError SuperMesh::delrect(string meshName, Rect rectangle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName]->delrect(rectangle);

	return error;
}

BError SuperMesh::setrect(string meshName, Rect rectangle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName]->setrect(rectangle);

	return error;
}

BError SuperMesh::resetrect(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName]->setrect(Rect(pMesh[meshName]->GetMeshDimensions()));

	return error;
}

//roughen mesh sides perpendicular to a named axis (axis = "x", "y", "z") to given depth (same units as h) with prng instantiated with given seed.
BError SuperMesh::RoughenMeshSides(string meshName, string axis, double depth, int seed)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	if ((axis != "x" && axis != "y" && axis != "z") || depth <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	error = pMesh[meshName]->RoughenMeshSides(axis, depth, seed);

	return error;
}

//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
//Rough both top and bottom if sides is empty, else it should be either -z or z.
BError SuperMesh::RoughenMeshSurfaces_Jagged(string meshName, double depth, double spacing, int seed, string sides)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (depth <= 0 || spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	error = pMesh[meshName]->RoughenMeshSurfaces_Jagged(depth, spacing, seed, sides);

	return error;
}

//clear roughness: set fine shape to coarse shape
BError SuperMesh::ClearMeshRoughness(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->CallModuleMethod(&Roughness::clear_roughness);

	return error;
}

//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
BError SuperMesh::GenerateGrains2D(string meshName, double spacing, int seed)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	error = pMesh[meshName]->GenerateGrains2D(spacing, seed);

	return error;
}

//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
BError SuperMesh::GenerateGrains3D(string meshName, double spacing, int seed)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	error = pMesh[meshName]->GenerateGrains3D(spacing, seed);

	return error;
}

//--------------------------------------------------------- MESH HANDLING - SETTINGS

BError SuperMesh::SetMagnetisationAngle(string meshName, double polar, double azim)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->SetMagnetisationAngle(polar, azim);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->SetMagnetisationAngle(polar, azim);
		}
	}

	return error;
}

BError SuperMesh::SetMagnetisationAngle_Rect(string meshName, double polar, double azim, Rect rectangle)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	Rect mesh_relrect = Rect(pMesh[meshName]->GetMeshDimensions());
	if (!mesh_relrect.contains(rectangle)) return error(BERROR_PARAMOUTOFBOUNDS);

	pMesh[meshName]->SetMagnetisationAngle(polar, azim, rectangle);

	return error;
}

//Invert magnetisation direction in given mesh (must be ferromagnetic)
BError SuperMesh::SetInvertedMagnetisation(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetInvertedMagnetisation();

	return error;
}

//Set random magentisation distribution in given mesh (must be ferromagnetic)
BError SuperMesh::SetRandomMagnetisation(string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	pMesh[meshName]->SetRandomMagnetisation();

	return error;
}

BError SuperMesh::SetMagnetisationDomainWall(string meshName, string longitudinal, string transverse, double width, double position)
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

	pMesh[meshName]->SetMagnetisationDomainWall(longi, trans, width, position);

	return error;
}

//Set Neel skyrmion in given mesh with chirality, diameter and relative x-y plane position
BError SuperMesh::SetSkyrmion(string meshName, int orientation, int chirality, double diameter, DBL2 position)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	DBL3 mesh_dim = pMesh[meshName]->GetMeshDimensions();

	Rect skyrmion_relrect = Rect(DBL3(position.x - diameter/2, position.y - diameter/2, 0.0), DBL3(position.x + diameter/2, position.y + diameter/2, mesh_dim.z));

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

		pMesh[meshName]->CallModuleMethod(&Zeeman::SetField, field_cartesian);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[idx]->CallModuleMethod(&Zeeman::SetField, field_cartesian);
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
	if (!error) error = SetMagnetisationAngle(left_dipole, 90.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagnetisationAngle(right_dipole, 90.0, 180.0);

	//2. Set moving mesh trigger, antisymmetric type and threshold
	odeSolver.SetMoveMeshTrigger(true, pMesh[meshName]->get_id());
	odeSolver.SetMoveMeshAntisymmetric(true);
	odeSolver.SetMoveMeshThreshold(MOVEMESH_ANTISYMMETRIC_THRESHOLD);

	//3. Set transverse domain wall structure along the x direction. Set the ends along x direction.
	SetMagnetisationDomainWall(meshName, "x", "y", pMesh[meshName]->GetMeshDimensions().x, 0.0);
	
	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	pMesh[meshName]->SetMagnetisationAngle(90, 0, Rect(DBL3(0), DBL3(MOVEMESH_ENDRATIO * dim.x, dim.y, dim.z)));
	pMesh[meshName]->SetMagnetisationAngle(90, 180, Rect(DBL3((1 - MOVEMESH_ENDRATIO) * dim.x, 0, 0), dim));

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
	odeSolver.SetMoveMeshTrigger(true, pMesh[meshName]->get_id());
	odeSolver.SetMoveMeshAntisymmetric(true);
	odeSolver.SetMoveMeshThreshold(MOVEMESH_ANTISYMMETRIC_THRESHOLD);

	//2. Set Bloch domain wall structure along the x direction. Set the ends along x direction.
	SetMagnetisationDomainWall(meshName, "z", "y", pMesh[meshName]->GetMeshDimensions().x, 0.0);

	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	pMesh[meshName]->SetMagnetisationAngle(0, 0, Rect(DBL3(0), DBL3(MOVEMESH_ENDRATIO * dim.x, dim.y, dim.z)));
	pMesh[meshName]->SetMagnetisationAngle(180, 0, Rect(DBL3((1 - MOVEMESH_ENDRATIO) * dim.x, 0, 0), dim));

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
	if (!error) error = SetMagnetisationAngle(left_dipole, 0.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagnetisationAngle(right_dipole, 180.0, 0.0);

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
	if (!error) error = SetMagnetisationAngle(left_dipole, 0.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagnetisationAngle(right_dipole, 180.0, 0.0);

	//2. Set moving mesh trigger, antisymmetric type and threshold
	odeSolver.SetMoveMeshTrigger(true, pMesh[meshName]->get_id());
	odeSolver.SetMoveMeshAntisymmetric(true);
	odeSolver.SetMoveMeshThreshold(MOVEMESH_ANTISYMMETRIC_THRESHOLD);

	//3. Set Bloch domain wall structure along the x direction. Set the ends along x direction.
	SetMagnetisationDomainWall(meshName, "z", "x", pMesh[meshName]->GetMeshDimensions().x, 0.0);

	DBL3 dim = pMesh[meshName]->GetMeshDimensions();
	pMesh[meshName]->SetMagnetisationAngle(0, 0, Rect(DBL3(0), DBL3(MOVEMESH_ENDRATIO * dim.x, dim.y, dim.z)));
	pMesh[meshName]->SetMagnetisationAngle(180, 0, Rect(DBL3((1 - MOVEMESH_ENDRATIO) * dim.x, 0, 0), dim));

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
	if (!error) error = SetMagnetisationAngle(left_dipole, 0.0, 0.0);
	if (!error) error = AddMesh(right_dipole, MESH_DIPOLE, meshRect_right);
	if (!error) error = SetMagnetisationAngle(right_dipole, 0.0, 0.0);

	//2. Set moving mesh trigger, symmetric type and threshold
	odeSolver.SetMoveMeshTrigger(true, pMesh[meshName]->get_id());
	odeSolver.SetMoveMeshAntisymmetric(false);
	odeSolver.SetMoveMeshThreshold(MOVEMESH_SYMMETRIC_THRESHOLD);

	//3. Set skyrmion in the centre amd up direction everywhere else
	pMesh[meshName]->SetMagnetisationAngle(0, 0);

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
	odeSolver.SetMoveMeshTrigger(false);

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

		if (flag == "x") error = pMesh[meshName]->Set_PBC_X(images);
		else if (flag == "y") error = pMesh[meshName]->Set_PBC_Y(images);
		else if (flag == "z") error = pMesh[meshName]->Set_PBC_Z(images);
		else return error(BERROR_INCORRECTVALUE);
	}
	else {

		//pbc setting for supermesh demag module
		if (IsSuperMeshModuleSet(MODS_SDEMAG)) {

			INT3 pbc_images = reinterpret_cast<SDemag*>(pSMod(MODS_SDEMAG))->Get_PBC();
			
			if (flag == "x") pbc_images.x = images;
			else if (flag == "y") pbc_images.y = images;
			else if (flag == "z") pbc_images.z = images;
			else return error(BERROR_INCORRECTVALUE);

			error = reinterpret_cast<SDemag*>(pSMod(MODS_SDEMAG))->Set_PBC(pbc_images);
		}
		else return error(BERROR_GPUERROR_CRIT);
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

		if (!pMesh[meshName]->MComputation_Enabled()) return error(BERROR_INCORRECTNAME);

		error = pMesh[meshName]->SetLinkStochastic(link_stochastic);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			if (pMesh[idx]->MComputation_Enabled()) error = pMesh[idx]->SetLinkStochastic(link_stochastic);
		}
	}

	return error;
}