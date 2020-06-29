#include "stdafx.h"
#include "Atom_Mesh.h"

#include "SuperMesh.h"

//----------------------------------- OTHER MESH SHAPE CONTROL

BError Atom_Mesh::copy_mesh_data(MeshBase& copy_this)
{
	BError error(__FUNCTION__);

	Atom_Mesh* pcopy_this = dynamic_cast<Atom_Mesh*>(&copy_this);

	if (pcopy_this == nullptr || meshType != pcopy_this->GetMeshType()) return error(BERROR_INCORRECTVALUE);

#if COMPILECUDA == 1
	//if CUDA on then first transfer data to cpu as there may be a mismatch
	if (paMeshCUDA) {

		error = paMeshCUDA->copy_shapes_to_cpu();
	}
#endif

	//1a. shape atomic moments
	if (M1.linear_size() && pcopy_this->M1.linear_size()) {

		M1.copy_values(pcopy_this->M1);
	}

	//2. shape electrical conductivity
	if (elC.linear_size() && pcopy_this->elC.linear_size()) elC.copy_values(pcopy_this->elC);

	//3. shape temperature
	if (Temp.linear_size() && pcopy_this->Temp.linear_size()) Temp.copy_values(pcopy_this->Temp);

#if COMPILECUDA == 1
	//if CUDA on then load back to gpu
	if (paMeshCUDA) {

		error = paMeshCUDA->copy_shapes_from_cpu();
	}
#endif

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}

	return error;
}

//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
BError Atom_Mesh::applymask(double zDepth_m, string fileName, function<vector<unsigned char>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, double zDepth_m, string fileName, function<vector<unsigned char>(string, INT2)>& bitmap_loader) -> BError {

		BError error;

		INT3 cells = VEC_quantity.n;

		if (!VEC_quantity.apply_bitmap_mask(bitmap_loader(fileName, INT2(cells.x, cells.y)), zDepth_m))
			return error(BERROR_COULDNOTLOADFILE);

		return error;
	};

	error = change_mesh_shape(run_this, zDepth_m, fileName, bitmap_loader);

	return error;
}

//set cells to empty in given box (delete by setting entries to zero)
BError Atom_Mesh::delrect(Rect rectangle)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, Rect& rectangle) -> BError {

		VEC_quantity.delrect(rectangle);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	return error;
}

//set cells to non-empty in given box
BError Atom_Mesh::setrect(Rect rectangle)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, Rect& rectangle) -> BError {

		VEC_quantity.setrect(rectangle, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	return error;
}

//roughen mesh sides perpendicular to a named axis (axis = "x", "y", "z") to given depth (same units as h) with prng instantiated with given seed.
BError Atom_Mesh::RoughenMeshSides(string axis, double depth, unsigned seed)
{
	BError error(__FUNCTION__);

	if ((axis != "x" && axis != "y" && axis != "z") || depth <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, string axis, double depth, unsigned seed) -> BError {

		bool success = true;

		string side_m = "-z", side_p = "z";

		if (axis == "x") {

			side_m = "-x";
			side_p = "x";
		}

		if (axis == "y") {

			side_m = "-y";
			side_p = "y";
		}

		success &= VEC_VC_quantity.generate_roughside(side_m, depth, seed);
		success &= VEC_VC_quantity.generate_roughside(side_p, depth, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, axis, depth, seed);

	return error;
}

//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
//Rough both top and bottom if sides is empty, else it should be either -z or z.
BError Atom_Mesh::RoughenMeshSurfaces_Jagged(double depth, double spacing, unsigned seed, string sides)
{
	BError error(__FUNCTION__);

	if (depth <= 0 || spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double depth, double spacing, unsigned seed, string sides) -> BError {

		bool success = true;

		success &= VEC_VC_quantity.generate_jagged_surfaces(depth, spacing, seed, sides);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, depth, spacing, seed, sides);

	return error;
}

//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
BError Atom_Mesh::GenerateGrains2D(double spacing, unsigned seed)
{
	BError error(__FUNCTION__);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double spacing, unsigned seed) -> BError {

		bool success = true;

		success = VEC_VC_quantity.generate_Voronoi2D(spacing, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, spacing, seed);

	return error;
}

//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
BError Atom_Mesh::GenerateGrains3D(double spacing, unsigned seed)
{
	BError error(__FUNCTION__);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double spacing, unsigned seed) -> BError {

		bool success = true;

		success = VEC_VC_quantity.generate_Voronoi3D(spacing, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, spacing, seed);

	return error;
}