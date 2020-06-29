#include "stdafx.h"
#include "Atom_Mesh.h"

#include "SuperMesh.h"

//----------------------------------- MESH INFO AND SIZE GET/SET METHODS

BError Atom_Mesh::SetMeshRect(Rect meshRect_, bool adjust_num_cells)
{
	BError error(__FUNCTION__);

	//cannot set a plane, must be a proper 3D rect
	if (meshRect_.IsPlane() || meshRect_.e <= meshRect_.s) return error(BERROR_INCORRECTVALUE);

	meshRect = meshRect_;

	auto adjust_default_cellsize = [](Rect meshRect, DBL3& cellsize, SZ3& cells, bool adjust_num_cells, bool min2 = false) {

		if (!cellsize.dim()) cellsize = DBL3(DEFAULTCELLSIZE);
		cells = round(meshRect / cellsize);

		if (adjust_num_cells) {

			if (cells.x > MAXSTARTINGCELLS_X) cells.x = MAXSTARTINGCELLS_X;
			if (cells.y > MAXSTARTINGCELLS_Y) cells.y = MAXSTARTINGCELLS_Y;
			if (cells.z > MAXSTARTINGCELLS_Z) cells.z = MAXSTARTINGCELLS_Z;
		}

		if (cells.x == 0) cells.x = 1;
		if (cells.y == 0) cells.y = 1;
		if (cells.z == 0) cells.z = 1;

		//minimum number of 2 cells in each dimension
		if (min2) {

			if (cells.x < 2) cells.x = 2;
			if (cells.y < 2) cells.y = 2;
			if (cells.z < 2) cells.z = 2;
		}

		//adjusted cellsize which will result in an integer number of cells with upper limits set
		cellsize = meshRect / cells;
	};

	//adjust cellsizes from their current values so they result in an integer number of cells with upper limits on number of cells in each dimension. The cellsizes can be manually adjusted later if needed.
	adjust_default_cellsize(meshRect, h, n, adjust_num_cells);
	adjust_default_cellsize(meshRect, h_dm, n_dm, adjust_num_cells);
	adjust_default_cellsize(meshRect, h_e, n_e, adjust_num_cells, true);
	adjust_default_cellsize(meshRect, h_t, n_t, adjust_num_cells, true);
	adjust_default_cellsize(meshRect, h_m, n_m, adjust_num_cells, true);

	error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHCHANGE);

	return error;
}

//magnetic properties
BError Atom_Mesh::SetMeshCellsize(DBL3 h_)
{
	BError error(__FUNCTION__);

	h = h_;

	//this sets the correct number of cells
	error = SetMeshRect(meshRect, false);

	return error;
}

//electrical conduction properties
BError Atom_Mesh::SetMeshECellsize(DBL3 h_e_)
{
	BError error(__FUNCTION__);

	h_e = h_e_;

	//this sets the correct number of cells
	error = SetMeshRect(meshRect, false);

	return error;
}

//thermal conduction properties
BError Atom_Mesh::SetMeshTCellsize(DBL3 h_t_)
{
	BError error(__FUNCTION__);

	h_t = h_t_;

	//this sets the correct number of cells
	error = SetMeshRect(meshRect, false);

	return error;
}

//Mechanical properties
BError Atom_Mesh::SetMeshMCellsize(DBL3 h_m_)
{
	BError error(__FUNCTION__);

	h_m = h_m_;

	//this sets the correct number of cells
	error = SetMeshRect(meshRect, false);

	return error;
}

//Set demagnetizing field evaluation macrocell size for atomistic meshes; may need to be adjusted if it doesn't result in an integer number of cells
BError Atom_Mesh::Set_Demag_Cellsize(DBL3 h_dm_)
{
	BError error(__FUNCTION__);

	h_dm = h_dm_;

	//this sets the correct number of cells
	error = SetMeshRect(meshRect, false);

	return error;
}