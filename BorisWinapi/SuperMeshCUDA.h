#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshDisplayCUDA.h"

#include "SDemagCUDA.h"
#include "StrayFieldCUDA.h"
#include "STransportCUDA.h"
#include "SHeatCUDA.h"
#include "OerstedCUDA.h"

class SuperMesh;
class PhysQ;

class SuperMeshCUDA :
	public MeshDisplayCUDA
{

private:

	//pointer to cpu version of supermesh 
	SuperMesh* pSMesh;

	//the entire super-mesh rectangle (rectangle containing all meshes)
	Rect& sMeshRect;

	//ferromagnetic super-mesh dimensions (this is the super-mesh that encompasses all ferromagnetic meshes)
	SZ3& n_fm;

	//ferromagnetic super-mesh discretization cellsize
	DBL3& h_fm;

	//ferromagnetic super-mesh rectangle
	Rect& sMeshRect_fm;

	//electric super-mesh dimensions (this is the super-mesh that encompasses all meshes with electrical computations enabled)
	SZ3& n_e;

	//electric super-mesh discretization cellsize
	DBL3& h_e;

	//electric super-mesh rectangle
	Rect& sMeshRect_e;

public:

	//make this object by copying data from the Mesh holding this object
	SuperMeshCUDA(SuperMesh* pSMesh_);

	~SuperMeshCUDA();

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	vector<PhysQ> FetchOnScreenPhysicalQuantity(double detail_level);

	//----------------------------------- Getters

	//check if ODE solver needs spin accumulation solved
	bool SolveSpinCurrent(void);
};

#endif
