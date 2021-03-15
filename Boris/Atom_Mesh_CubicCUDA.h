#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "Atom_MeshCUDA.h"

#include "BorisCUDALib.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

class Atom_Mesh_Cubic;
class ManagedAtom_DiffEqCubicCUDA;

//Store Mesh quantities as cu_obj managed cuda VECs
class Atom_Mesh_CubicCUDA :
	public Atom_MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	Atom_Mesh_Cubic *paMeshCubic;

	// MONTE-CARLO DATA

	// Constrained MONTE-CARLO DATA

	//total moment along constrained direction
	cu_obj<cuBReal> cmc_M;

	//constraining direction
	cu_obj<cuReal3> cmc_n;

	//mc indices and shuffling auxiliary array : same as for the cpu version, but generate unsigned random numbers, not doubles, for most efficient sort-based shuffle
	cu_arr<unsigned> mc_indices_red, mc_indices_black;
	cu_arr<unsigned> mc_shuf_red, mc_shuf_black;

protected:

public:

	//make this object by copying data from the Mesh holding this object
	Atom_Mesh_CubicCUDA(Atom_Mesh_Cubic* paMesh);

	~Atom_Mesh_CubicCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//Take a Monte Carlo Metropolis step in this atomistic mesh
	cuBReal Iterate_MonteCarloCUDA_Classic(cuBReal mc_cone_angledeg, double target_acceptance_rate);

	//Take a constrained Monte Carlo Metropolis step in this atomistic mesh
	cuBReal Iterate_MonteCarloCUDA_Constrained(cuBReal mc_cone_angledeg, double target_acceptance_rate);

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold);

	//----------------------------------- OTHER CONTROL METHODS : implement pure virtual Atom_Mesh methods

	void Set_MonteCarlo_Constrained(DBL3 cmc_n_);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//get exchange_couple_to_meshes status flag from the cpu version
	bool GetMeshExchangeCoupling(void);

	//----------------------------------- VALUE GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	cuBReal GetTopologicalCharge(cuRect rectangle);

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void);

	//----------------------------------- ODE METHODS IN MAGNETIC MESH : Atom_Mesh_CubicCUDA.cu

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_dmdt(cuBox avBox);

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_mxdmdt(cuBox avBox);

	//-----------------------------------OBJECT GETTERS

	cu_obj<ManagedAtom_DiffEqCubicCUDA>& Get_ManagedAtom_DiffEqCUDA(void);
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

	//Take a Monte Carlo Metropolis step in this atomistic mesh
	cuBReal Iterate_MonteCarloCUDA_Classic(cuBReal mc_cone_angledeg, double target_acceptance_rate) { return 0.0; }

	//Take a constrained Monte Carlo Metropolis step in this atomistic mesh
	cuBReal Iterate_MonteCarloCUDA_Constrained(cuBReal mc_cone_angledeg, double target_acceptance_rate) { return 0.0; }

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold) { return 0.0; }

	//----------------------------------- VALUE GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	cuBReal GetTopologicalCharge(cuRect rectangle) { return 0.0; }

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void) {}

	//----------------------------------- ODE METHODS IN MAGNETIC MESH : Atom_Mesh_CubicCUDA.cu

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_dmdt(cuBox avBox) { return DBL3(); }

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_mxdmdt(cuBox avBox) { return DBL3(); }
};

#endif
#endif
