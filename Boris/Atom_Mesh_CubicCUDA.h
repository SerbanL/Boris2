#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "Atom_MeshCUDA.h"

#include "BorisCUDALib.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

class Atom_Mesh_Cubic;

//Store Mesh quantities as cu_obj managed cuda VECs
class Atom_Mesh_CubicCUDA :
	public Atom_MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	Atom_Mesh_Cubic *paMeshCubic;

public:

	//make this object by copying data from the Mesh holding this object
	Atom_Mesh_CubicCUDA(Atom_Mesh_Cubic* paMesh);

	~Atom_Mesh_CubicCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//get exchange_couple_to_meshes status flag from the cpu version
	bool GetMeshExchangeCoupling(void);

	//----------------------------------- VALUE GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	cuBReal GetTopologicalCharge(cuRect rectangle);

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void);
};

#else

class Atom_Mesh_CubicCUDA :
	public Atom_MeshCUDA
{

private:

public:

	//make this object by copying data from the Mesh holding this object
	Atom_Mesh_CubicCUDA(Atom_Mesh* paMesh) :
		Atom_MeshCUDA(paMesh)
	{}

	~Atom_Mesh_CubicCUDA() {}

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold) { return 0.0; }

	//----------------------------------- VALUE GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	cuBReal GetTopologicalCharge(cuRect rectangle) { return 0.0; }

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void) {}
};

#endif
#endif
