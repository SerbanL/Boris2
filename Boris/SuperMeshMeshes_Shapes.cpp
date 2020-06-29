#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- MESH HANDLING - SHAPES

//set mesh rect for named mesh (any if applicable any other dependent meshes) and update dependent save data rects by calling the provided function.
BError SuperMesh::SetMeshRect(string meshName, Rect meshRect, std::function<void(string, Rect)> save_data_updater)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	//cannot set a plane, must be a proper 3D rect
	if (meshRect.IsPlane() || meshRect.e <= meshRect.s) return error(BERROR_INCORRECTVALUE);

	//snap mesh rectangle coordinates to this
	meshRect.snap(MINMESHSPACE);

	if (!scale_rects) {

		Rect meshRect_old = pMesh[meshName]->GetMeshRect();

		//scale just the rectangle in named mesh
		error = pMesh[meshName]->SetMeshRect(meshRect, true);
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

			//snap mesh rectangle coordinates to this
			meshRect_scaled_shifted.snap(MINMESHSPACE);

			error = pMesh[idx]->SetMeshRect(meshRect_scaled_shifted, true);
			if (!error) save_data_updater(pMesh.get_key_from_index(idx), meshRect_original);
			else return error;
		}
	}

	return error;
}

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