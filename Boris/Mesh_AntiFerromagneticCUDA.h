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

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//get exchange_couple_to_meshes status flag from the cpu version
	bool GetMeshExchangeCoupling(void);

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_AntiFerromagneticCUDA.cu

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_dmdt(cuBox avBox);
	DBL3 Average_dmdt2(cuBox avBox);

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_mxdmdt(cuBox avBox);
	DBL3 Average_mxdmdt2(cuBox avBox);

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

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold) { return 0.0; }

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_AntiFerromagneticCUDA.cu

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_dmdt(cuBox avBox) { return DBL3(); }
	DBL3 Average_dmdt2(cuBox avBox) { return DBL3(); }

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	DBL3 Average_mxdmdt(cuBox avBox) { return DBL3(); }
	DBL3 Average_mxdmdt2(cuBox avBox) { return DBL3(); }
};

#endif
#endif
