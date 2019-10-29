#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshCUDA.h"

class FMesh;

//Store Mesh quantities as cu_obj managed cuda VECs
class FMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	FMesh *pFMesh;
	
public:

	//make this object by copying data from the Mesh holding this object
	FMeshCUDA(FMesh* pMesh);

	~FMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS
	
	//call when the mesh dimensions have changed - sets every quantity to the right dimensions
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	
	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//get exchange_couple_to_meshes status flag from the cpu version
	bool GetMeshExchangeCoupling(void);
};

#endif