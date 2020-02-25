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
class Any;

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

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType);

	//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
	void PrepareDisplayedMeshValue(void);

	//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
	//Return an Any as the displayed quantity could be either a scalar or a vector.
	Any GetDisplayedMeshValue(DBL3 abs_pos);

	//return average value for currently displayed mesh quantity in the given relative rectangle
	Any GetAverageDisplayedMeshValue(Rect rel_rect);

	//----------------------------------- Getters

	//check if ODE solver needs spin accumulation solved
	bool SolveSpinCurrent(void);
};

#endif
