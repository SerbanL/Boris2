#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "MeshCUDA.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "BorisCUDALib.h"

class AFMesh;
class ManagedDiffEqAFMCUDA;

//Store Mesh quantities as cu_obj managed cuda VECs
class AFMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	AFMesh *pAFMesh;

public:

	//make this object by copying data from the Mesh holding this object
	AFMeshCUDA(AFMesh* pMesh);

	~AFMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//Take a Monte Carlo Metropolis step in this mesh
	cuBReal Iterate_MonteCarloCUDA_Classic(cuBReal mc_cone_angledeg, double target_acceptance_rate);

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//get exchange_couple_to_meshes status flag from the cpu version
	bool GetMeshExchangeCoupling(void);

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_AntiFerromagneticCUDA.cu

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_dmdt(cuBox avBox);
	cuReal3 Average_dmdt2(cuBox avBox);

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_mxdmdt(cuBox avBox);
	//for sub-lattice B
	cuReal3 Average_mxdmdt2(cuBox avBox);
	//mixed sub-lattices A and B
	cuReal3 Average_mxdm2dt(cuBox avBox);
	cuReal3 Average_m2xdmdt(cuBox avBox);

	//Save current magnetization in sM VECs (e.g. useful to reset dM / dt calculation)
	void SaveMagnetization(void);

	//-----------------------------------OBJECT GETTERS

	cu_obj<ManagedDiffEqAFMCUDA>& Get_ManagedDiffEqCUDA(void);
};

#else

class AFMeshCUDA :
	public MeshCUDA
{

private:

public:

	//make this object by copying data from the Mesh holding this object
	AFMeshCUDA(Mesh* pMesh) :
		MeshCUDA(pMesh)
	{}

	~AFMeshCUDA() {}

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_AntiFerromagneticCUDA.cu

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_dmdt(cuBox avBox) { return cuReal3(); }
	cuReal3 Average_dmdt2(cuBox avBox) { return cuReal3(); }

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_mxdmdt(cuBox avBox) { return cuReal3(); }
	cuReal3 Average_mxdmdt2(cuBox avBox) { return cuReal3(); }
	//mixed sub-lattices A and B
	cuReal3 Average_mxdm2dt(cuBox avBox) { return cuReal3(); }
	cuReal3 Average_m2xdmdt(cuBox avBox) { return cuReal3(); }
};

#endif
#endif
