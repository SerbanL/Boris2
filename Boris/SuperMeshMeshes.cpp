#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- MESH HANDLING - COMPONENTS

BError SuperMesh::AddMesh(std::string meshName, MESH_ meshType, Rect meshRect)
{
	BError error(__FUNCTION__);

	if (contains(meshName) || meshName == superMeshHandle) return error(BERROR_INCORRECTNAME);

	//when creating a new mesh the user supplies the mesh type, rectangle and name.
	//we also need a cellsize. The user can change the cellsize later but it's easier not to ask for it when the mesh is created - use a default value instead.
	//To choose default cellsize : use 5nm cubic cell up to a given number of maximum cells in each dimension. After this increase starting cellsize to keep to this maximum number of cells.
	//The user may be trying to create a large mesh for which there might not be enough memory with a 5nm cellsize, but intending to set a larger cellsize after mesh creation.

	auto make_cellsize = [&](double default_cellsize, Rect meshRect) -> DBL3 {

		DBL3 cellsize(default_cellsize);

		INT3 cells = round(meshRect / cellsize);

		if (cells.x > MAXSTARTINGCELLS_X) cells.x = MAXSTARTINGCELLS_X;
		if (cells.y > MAXSTARTINGCELLS_Y) cells.y = MAXSTARTINGCELLS_Y;
		if (cells.z > MAXSTARTINGCELLS_Z) cells.z = MAXSTARTINGCELLS_Z;

		//adjusted cellsize which will result in an integer number of cells with upper limits set
		cellsize = meshRect / cells;

		return cellsize;
	};

	switch (meshType) {

	case MESH_FERROMAGNETIC:
		pMesh.push_back(new FMesh(meshRect, make_cellsize(DEFAULTCELLSIZE, meshRect), this), meshName);
		break;

	case MESH_ANTIFERROMAGNETIC:
		pMesh.push_back(new AFMesh(meshRect, make_cellsize(DEFAULTCELLSIZE, meshRect), this), meshName);
		break;

	case MESH_DIAMAGNETIC:
		pMesh.push_back(new DiaMesh(meshRect, make_cellsize(DEFAULTCELLSIZE, meshRect), this), meshName);
		break;

	case MESH_DIPOLE:
		pMesh.push_back(new DipoleMesh(meshRect, make_cellsize(DEFAULTCELLSIZE, meshRect), this), meshName);
		break;

	case MESH_METAL:
		pMesh.push_back(new MetalMesh(meshRect, make_cellsize(DEFAULTCELLSIZE, meshRect), this), meshName);
		break;

	case MESH_INSULATOR:
		pMesh.push_back(new InsulatorMesh(meshRect, make_cellsize(DEFAULTCELLSIZE, meshRect), this), meshName);
		break;

	case MESH_ATOM_CUBIC:
		pMesh.push_back(new Atom_Mesh_Cubic(meshRect, make_cellsize(DEFAULTATOMCELLSIZE, meshRect), this), meshName);
		break;
	}

	//check the mesh was created correctly - if not, delete it
	error = pMesh.back()->Error_On_Create();

	if (!error) error = UpdateConfiguration(UPDATECONFIG_MESHADDED);
	else pMesh.pop_back();

	return error;
}

BError SuperMesh::DelMesh(std::string meshName)
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

BError SuperMesh::RenameMesh(std::string oldName, std::string newName)
{
	BError error(__FUNCTION__);

	//must have old name and new name cannot be that of an existing mesh
	if (!contains(oldName) || contains(newName) || newName == superMeshHandle) return error(BERROR_INCORRECTNAME);

	//make changes
	pMesh.change_key(oldName, newName);

	if (activeMeshName == oldName) activeMeshName = newName;

	return error;
}

BError SuperMesh::SetMeshFocus(std::string meshName)
{
	BError error(__FUNCTION__);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName != superMeshHandle)
		activeMeshName = meshName;

	return error;
}