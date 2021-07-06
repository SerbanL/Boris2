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
struct MeshShape;

class SuperMeshCUDA :
	public MeshDisplayCUDA
{

private:

	//Auxiliary

	//storage for extracted mesh profiles
	cu_arr<cuBReal> profile_storage_sca;
	cu_arr<cuReal3> profile_storage_vec;

	//Others

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

	std::vector<PhysQ> FetchOnScreenPhysicalQuantity(double detail_level);

	//save the quantity currently displayed on screen for named mesh in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(std::string meshName, std::string fileName, std::string ovf2_dataType);

	//extract profile from named mesh, from currently display mesh quantity, but reading directly from the quantity
	//Displayed mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
	//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
	void GetPhysicalQuantityProfile(DBL3 start, DBL3 end, double step, DBL3 stencil, std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, std::string meshName, bool do_average, bool read_average);

	//return average value for currently displayed mesh quantity for named mesh in the given relative rectangle
	Any GetAverageDisplayedMeshValue(std::string meshName, Rect rel_rect, std::vector<MeshShape> shapes);

	//----------------------------------- Getters

	//check if ODE solver needs spin accumulation solved
	bool SolveSpinCurrent(void);

	//check disabled_transport_solver flag
	bool DisabledTransportSolver(void);
};

#endif
