#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include <vector>
#include <functional>
#include <string>

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"
#include "MeshDisplayCUDA.h"
#include "ManagedMeshCUDA.h"

#include "ZeemanCUDA.h"
#include "ExchangeCUDA.h"
#include "DMExchangeCUDA.h"
#include "iDMExchangeCUDA.h"
#include "AnisotropyCUDA.h"
#include "AnisotropyCubiCUDA.h"
#include "Demag_NCUDA.h"
#include "DemagCUDA.h"
#include "SDemagCUDA_Demag.h"
#include "TransportCUDA.h"
#include "HeatCUDA.h"
#include "SurfExchangeCUDA.h"
#include "SOTFieldCUDA.h"
#include "RoughnessCUDA.h"

using namespace std;

class Mesh;
class PhysQ;

//Store Mesh quantities as cu_obj managed cuda VECs
class MeshCUDA :
	public MeshParamsCUDA,
	public MeshDisplayCUDA
{

private:

	//pointer to cpu version of this mesh (base) - need it in the destructor. 
	//Still keep references to some Mesh data members here as we cannot use pMesh in .cu files (cannot have BorisLib.h in those compilation units - real headache, will need to fix this at some point somehow: problem is the nvcc compiler throws errors due to C++14 code in BorisLib)
	Mesh *pMesh;

	bool holder_mesh_destroyed = false;

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

public:

	//Managed Mesh
	cu_obj<ManagedMeshCUDA> cuMesh;

	//-----Mesh properties

	//the mesh rectangle in meters : this defines the mesh position in the world (in the super-mesh). cellsizes must divide the mesh rectangle, giving the number of cells in each dimension
	Rect& meshRect;

	//-----Ferromagnetic properties

	//number of cells (n.x, n.y, n.z)
	SZ3& n;

	//cellsizes (h.x, h.y, h.z)
	DBL3& h;

	//Magnetization
	cu_obj<cuVEC_VC<cuReal3>> M;

	//Additional magnetization used for antiferromagnetic meshes with 2 sub-lattice local approximation; exactly same dimensions as M
	cu_obj<cuVEC_VC<cuReal3>> M2;

	//effective field - sum total field of all the added modules
	cu_obj<cuVEC<cuReal3>> Heff;

	//Additional effective field used for antiferromagnetic meshes with 2 sub-lattice local approximation; exactly same dimensions as Heff
	cu_obj<cuVEC<cuReal3>> Heff2;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//number of cells for electrical properties
	SZ3& n_e;

	//cellsize for electrical properties
	DBL3& h_e;

	//electrical potential - on n_e, h_e mesh
	cu_obj<cuVEC_VC<cuBReal>> V;

	//electrical conductivity - on n_e, h_e mesh
	cu_obj<cuVEC_VC<cuBReal>> elC;

	//electric field - on n_e, h_e mesh
	cu_obj<cuVEC_VC<cuReal3>> E;

	//spin accumulation - on n_e, h_e mesh
	cu_obj<cuVEC_VC<cuReal3>> S;

	//-----Thermal conduction properties

	//number of cells for thermal properties
	SZ3& n_t;

	//cellsize for thermal properties
	DBL3& h_t;

	//temperature calculated by Heat module
	cu_obj<cuVEC_VC<cuBReal>> Temp;

	//-----Elastic properties

public:

	//make this object by copying data from the Mesh holding this object
	MeshCUDA(Mesh* pMesh);

	virtual ~MeshCUDA();

	//-------------------------- Error report / Management

	//obtain error_on_create from this mesh, as well as any set modules - return first error found
	BError Error_On_Create(void) { return error_on_create; }

	void Holder_Mesh_Destroyed(void) { holder_mesh_destroyed = true; }
	bool Holder_Mesh_Available(void) { return !holder_mesh_destroyed; }

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//call when the mesh dimensions have changed - sets every quantity to the right dimensions
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) = 0;
	
	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	PhysQ FetchOnScreenPhysicalQuantity(double detail_level);

	//----------------------------------- MESH INFO GET/SET METHODS

	int GetMeshType(void);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//magnetization dynamics computation enabled
	bool MComputation_Enabled(void);

	bool Magnetisation_Enabled(void);

	//electrical conduction computation enabled
	bool EComputation_Enabled(void);

	//thermal conduction computation enabled
	bool TComputation_Enabled(void);

	//check if interface conductance is enabled (for spin transport solver)
	bool GInterface_Enabled(void);

	//check if the ODECommon::available flag is true (ode step solved)
	bool CurrentTimeStepSolved(void);
	
	//check evaluation speedup flag in ODECommon
	int EvaluationSpeedup(void);

	//check in ODECommon the type of field update we need to do depending on the ODE evaluation step
	int Check_Step_Update(void);

	virtual bool GetMeshExchangeCoupling(void) { return false; }

	//----------------------------------- VALUE GETTERS

	//get average magnetisation in given rectangle (entire mesh if none specified)
	cuReal3 GetAverageMagnetisation(cuRect rectangle);
	cuReal3 GetAverageMagnetisation2(cuRect rectangle);

	cuBReal GetAverageElectricalPotential(cuRect rectangle);
	cuReal3 GetAverageSpinAccumulation(cuRect rectangle);
	cuBReal GetAverageElectricalConductivity(cuRect rectangle);

	cuBReal GetAverageTemperature(cuRect rectangle);

	//----------------------------------- MESH SHAPE CONTROL

	//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
	BError copy_shapes_from_cpu(void);

	//copy all meshes controlled using change_mesh_shape from gpu to cpu versions
	BError copy_shapes_to_cpu(void);
};

#endif