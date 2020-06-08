#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include <vector>
#include <functional>
#include <string>

#include "BorisCUDALib.h"

#include "CompileFlags.h"
#include "ErrorHandler.h"

#include "Atom_MeshParamsCUDA.h"
#include "MeshDisplayCUDA.h"
#include "ManagedAtom_MeshCUDA.h"

using namespace std;

class Atom_Mesh;
class PhysQ;
class Any;

//Store Atom_Mesh quantities as cu_obj managed cuda VECs
class Atom_MeshCUDA :
	public Atom_MeshParamsCUDA,
	public MeshDisplayCUDA
{
private:

	//pointer to cpu version of this mesh (base) - need it in the destructor. 
	//Still keep references to some Mesh data members here as we cannot use pMesh in .cu files (cannot have BorisLib.h in those compilation units - real headache, will need to fix this at some point somehow: problem is the nvcc compiler throws errors due to C++14 code in BorisLib)
	Atom_Mesh *paMesh;

	bool holder_mesh_destroyed = false;

protected:

	//auxiliary for computations
	cu_obj<cuBReal> aux_real;

	//auxiliary cuVEC for computations
	cu_obj<cuVEC<cuBReal>> aux_vec_sca;

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

public:

	//Managed Mesh
	cu_obj<ManagedAtom_MeshCUDA> cuaMesh;

	//-----Mesh properties

	//the mesh rectangle in meters : this defines the mesh position in the world (in the super-mesh). cellsizes must divide the mesh rectangle, giving the number of cells in each dimension
	Rect& meshRect;

	//-----Magnetic properties

	//number of cells (n.x, n.y, n.z)
	SZ3& n;

	//cellsizes (h.x, h.y, h.z)
	DBL3& h;

	//Atomic moments in units of Bohr magneton using double floating point precision (first sub-lattice : used in cubic, bcc, fcc, hcp)
	cu_obj<cuVEC_VC<cuReal3>> M1;

	//effective field (units of A/m : easier to integrate with micromagnetic meshes in a multiscale simulation this way) - sum total field of all the added modules (first sub-lattice : used in cubic, bcc, fcc, hcp)
	cu_obj<cuVEC<cuReal3>> Heff1;

	//-----Demagnetizing field macrocell : compute demagnetizing field or macroscopic dipole-dipole interaction using this cellsize

	//number of cells (n.x, n.y, n.z)
	SZ3& n_dm;

	//cellsize; doesn't have to be an integer number of atomistic cells
	DBL3& h_dm;

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

	//temperature calculated by Heat module (primary temperature, always used for 1-temperature model; for multi-temperature models in metals this is the itinerant electron temperature)
	cu_obj<cuVEC_VC<cuBReal>> Temp;

	//lattice temperature used in many-T models
	cu_obj<cuVEC_VC<cuBReal>> Temp_l;

	//-----Elastic properties

	//number of cells for mechanical properties
	SZ3& n_m;

	//cellsize for mechanical properties
	DBL3& h_m;

	//mechanical displacement vectors - on n_m, h_m mesh
	cu_obj<cuVEC_VC<cuReal3>> u_disp;

	//strain tensor (symmetric):
	//diagonal and off-diagonal components - on n_m, h_m mesh
	//xx, yy, zz
	cu_obj<cuVEC_VC<cuReal3>> strain_diag;
	//yz, xz, xy
	cu_obj<cuVEC_VC<cuReal3>> strain_odiag;

public:

	//------------------------CTOR/DTOR

	//make this object by copying data from the Mesh holding this object
	Atom_MeshCUDA(Atom_Mesh* paMesh);

	virtual ~Atom_MeshCUDA();

	//-------------------------- Error report / Management

	//obtain error_on_create from this mesh, as well as any set modules - return first error found
	BError Error_On_Create(void) { return error_on_create; }

	void Holder_Mesh_Destroyed(void) { holder_mesh_destroyed = true; }
	bool Holder_Mesh_Available(void) { return !holder_mesh_destroyed; }

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) = 0;

	//This is a "softer" version of UpdateConfiguration, which can be used any time and doesn't require the object to be Uninitialized; 
	//this will typically involve changing a value across multiple objects, thus better to call this method rather than try to remember which objects need the value changed.
	virtual void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) = 0;

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	PhysQ FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground);

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType);

	//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
	void PrepareDisplayedMeshValue(void);

	//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
	//Return an Any as the displayed quantity could be either a scalar or a vector.
	Any GetDisplayedMeshValue(DBL3 abs_pos);

	//return average value for currently displayed mesh quantity in the given relative rectangle
	Any GetAverageDisplayedMeshValue(Rect rel_rect);

	//copy aux_vec_sca in GPU memory to displayVEC in CPU memory
	void copy_aux_vec_sca(VEC<double>& displayVEC);

	//----------------------------------- MESH INFO GET/SET METHODS

	//NOTE : these type of methods are required here so we can access these properties in .cu files, where we cannot use a Mesh pointer
	//Using a Mesh pointer means having to include Mesh.h, thus also BorisLib.h -> this won't compile since nvcc doesn't recognise C++14 code in BorisLib
	//A future version of nvcc will make this practice redundant

	int GetMeshType(void);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//magnetization dynamics computation enabled
	bool MComputation_Enabled(void);

	bool Magnetism_Enabled(void);

	//electrical conduction computation enabled
	bool EComputation_Enabled(void);

	//thermal conduction computation enabled
	bool TComputation_Enabled(void);

	//mechanical computation enabled
	bool MechComputation_Enabled(void);

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

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	virtual cuBReal GetTopologicalCharge(cuRect rectangle) = 0;

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	virtual void Compute_TopoChargeDensity(void) = 0;

	cuBReal GetStageTime(void);
	int GetStageStep(void);

	//----------------------------------- MESH SHAPE CONTROL

	//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
	BError copy_shapes_from_cpu(void);

	//copy all meshes controlled using change_mesh_shape from gpu to cpu versions
	BError copy_shapes_to_cpu(void);

	//-----------------------------------OBJECT GETTERS

	cu_obj<ManagedAtom_DiffEq_CommonCUDA>& Get_ManagedAtom_DiffEq_CommonCUDA(void);
};

#endif
