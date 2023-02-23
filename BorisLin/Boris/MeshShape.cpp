#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

//----------------------------------- OTHER MESH SHAPE CONTROL

BError Mesh::copy_mesh_data(MeshBase& copy_this)
{
	BError error(__FUNCTION__);

	Mesh* pcopy_this = dynamic_cast<Mesh*>(&copy_this);

	if (pcopy_this == nullptr || meshType != pcopy_this->GetMeshType()) return error(BERROR_INCORRECTVALUE);

#if COMPILECUDA == 1
	//if CUDA on then first transfer data to cpu as there may be a mismatch
	if (pMeshCUDA) {

		error = pcopy_this->pMeshCUDA->copy_shapes_to_cpu();
	}
#endif

	//1a. shape magnetization
	if (M.linear_size() && pcopy_this->M.linear_size()) {

		//if Roughness module is enabled then apply shape via the Roughness module instead
		if (IsModuleSet(MOD_ROUGHNESS) && pcopy_this->IsModuleSet(MOD_ROUGHNESS)) {

			error = dynamic_cast<Roughness*>(pMod(MOD_ROUGHNESS))->copy_roughness(dynamic_cast<Roughness*>(pcopy_this->pMod(MOD_ROUGHNESS)));

			if (error) return error;
		}
		else {

			M.copy_values(pcopy_this->M);

			//1b. shape magnetization in AF meshes
			if (M2.linear_size() && pcopy_this->M2.linear_size()) {

				M2.copy_values(pcopy_this->M2);
			}
		}
	}

	//2. shape electrical conductivity
	if (elC.linear_size() && pcopy_this->elC.linear_size()) elC.copy_values(pcopy_this->elC);

	//3. shape temperature
	if (Temp.linear_size() && pcopy_this->Temp.linear_size()) Temp.copy_values(pcopy_this->Temp);

#if COMPILECUDA == 1
	//if CUDA on then load back to gpu
	if (pMeshCUDA) {

		error = pMeshCUDA->copy_shapes_from_cpu();
	}
#endif

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}

	return error;
}

//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
BError Mesh::applymask(double zDepth_m, std::string fileName, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, double zDepth_m, std::string fileName, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader) -> BError {

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
BError Mesh::delrect(Rect rectangle)
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
BError Mesh::setrect(Rect rectangle)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, Rect& rectangle) -> BError {

		VEC_quantity.setrect(rectangle, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	return error;
}

//roughen mesh sides (side = "x", "y", "z", "-x", "-y", or "-z") to given depth (same units as h) with prng instantiated with given seed.
BError Mesh::RoughenMeshSides(std::string side, double depth, int seed)
{
	BError error(__FUNCTION__);

	if ((side != "x" && side != "y" && side != "z" && side != "-x" && side != "-y" && side != "-z") || depth <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, std::string side, double depth, int seed) -> BError {

		bool success = true;

		success &= VEC_VC_quantity.generate_roughside(side, depth, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, side, depth, seed);

	return error;
}

//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
//Rough both top and bottom if sides is empty, else it should be either -z or z.
BError Mesh::RoughenMeshSurfaces_Jagged(double depth, double spacing, int seed, std::string sides)
{
	BError error(__FUNCTION__);

	if (depth <= 0 || spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double depth, double spacing, int seed, std::string sides) -> BError {

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
BError Mesh::GenerateGrains2D(double spacing, int seed)
{
	BError error(__FUNCTION__);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double spacing, int seed) -> BError {

		bool success = true;

		success = VEC_VC_quantity.generate_Voronoi2D_Grains(spacing, seed);

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
BError Mesh::GenerateGrains3D(double spacing, int seed)
{
	BError error(__FUNCTION__);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double spacing, int seed) -> BError {

		bool success = true;

		success = VEC_VC_quantity.generate_Voronoi3D_Grains(spacing, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, spacing, seed);

	return error;
}

//Disk with dimensions (x, y diameters, thickness), centre position relative to mesh, rotation angles, number of repetitions along x, y, z (1, 1, 1 for no repetitions), displacement values if repetitions used
BError Mesh::shape_disk(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_disk(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//rectangle
BError Mesh::shape_rect(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_rect(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//triangle
BError Mesh::shape_triangle(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_triangle(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//prolate ellipsoid
BError Mesh::shape_ellipsoid(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_ellipsoid(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//pyramid
BError Mesh::shape_pyramid(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_pyramid(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//tetrahedron
BError Mesh::shape_tetrahedron(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_tetrahedron(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//cone
BError Mesh::shape_cone(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_cone(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//torus
BError Mesh::shape_torus(MeshShape shape)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, MeshShape shape) -> BError {

		VEC_quantity.shape_torus(shape, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shape);

	return error;
}

//general shape setting function, can set composite shape using combination of the above elementary shapes
BError Mesh::shape_set(std::vector<MeshShape> shapes)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, std::vector<MeshShape> shapes) -> BError {

		VEC_quantity.shape_set(shapes, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, shapes);

	return error;
}