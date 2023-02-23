#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include <string>

#include "BorisCUDALib.h"

#include "CompileFlags.h"
#include "ErrorHandler.h"

#include "MeshDisplayCUDA.h"

#include "MeshDefs.h"

class MeshBase;
class PhysQ;
class Any;

class MeshBaseCUDA
{
private:

	//pointer to cpu version of this mesh (base) - need it in the destructor. 
	MeshBase *pMeshBase;

protected:

	bool holder_mesh_destroyed = false;

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//--------Auxiliary

	//auxiliary for computations
	cu_obj<cuBReal> aux_real;
	cu_obj<cuReal3> aux_real3;
	cu_obj<size_t> aux_int;

	//auxiliary cuVEC for computations
	cu_obj<cuVEC<cuBReal>> aux_vec_sca;

	//storage for extracted mesh profiles
	cu_arr<cuBReal> profile_storage_sca;
	cu_arr<cuReal3> profile_storage_vec;

	//auxiliary VECs for computations
	cu_obj<cuVEC<cuReal3>> auxVEC_cuReal3;
	cu_obj<cuVEC<cuBReal>> auxVEC_cuBReal;

	// MONTE-CARLO DATA

	//random number generator - used by Monte Carlo methods
	cu_obj<cuBorisRand<>> prng;

	//last Monte-Carlo step acceptance probability
	cu_obj<cuBReal> mc_acceptance_rate;

	//save last acceptance rate
	double mc_acceptance_rate_last = 0.0;

	//don't need to compute the acceptance rate every single iteration, unless the acceptance rate is not within bounds : when this counter is zero perform reduction.
	int mc_acceptance_reduction_counter = 0;

	// Constrained MONTE-CARLO DATA

	//constraining direction
	cu_obj<cuReal3> cmc_n;

	//total moment along constrained direction
	cu_obj<cuBReal> cmc_M;

public:

	//-----Mesh properties

	//the mesh rectangle in meters : this defines the mesh position in the world (in the super-mesh). cellsizes must divide the mesh rectangle, giving the number of cells in each dimension
	Rect& meshRect;

	//-----Magnetic properties

	//number of cells (n.x, n.y, n.z)
	SZ3& n;

	//cellsizes (h.x, h.y, h.z)
	DBL3& h;

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

protected:

	//zero all single aux avalues
	void Zero_aux_values(void);

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS AUXILIARY

	//average into profile_storage_sca / profile_storage_vec
	void average_mesh_profile(cu_obj<cuVEC<cuBReal>>& cuvec_sca, int& num_profile_averages);
	void average_mesh_profile(cu_obj<cuVEC_VC<cuBReal>>& cuvecvc_sca, int& num_profile_averages);
	void average_mesh_profile(cu_obj<cuVEC<cuReal3>>& cuvec_vec, int& num_profile_averages);
	void average_mesh_profile(cu_obj<cuVEC_VC<cuReal3>>& cuvecvc_vec, int& num_profile_averages);

public:

	//------------------------CTOR/DTOR

	//make this object by copying data from the Mesh holding this object
	MeshBaseCUDA(MeshBase* pMeshBase_);

	virtual ~MeshBaseCUDA();

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

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	virtual cuBReal CheckMoveMesh(bool antisymmetric, double threshold) = 0;

	void Set_MonteCarlo_Constrained(DBL3 cmc_n_) { cmc_n.from_cpu(cmc_n_); }

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	virtual PhysQ FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground) = 0;

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	virtual BError SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType, MESHDISPLAY_ quantity) = 0;

	//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
	//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
	//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
	virtual void GetPhysicalQuantityProfile(
		DBL3 start, DBL3 end, double step, DBL3 stencil, 
		std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, 
		bool do_average, bool read_average, MESHDISPLAY_ quantity) = 0;

	//return average value for currently displayed mesh quantity in the given relative rectangle
	virtual Any GetAverageDisplayedMeshValue(Rect rel_rect, MESHDISPLAY_ quantity) = 0;

	//copy aux_vec_sca in GPU memory to displayVEC in CPU memory
	virtual void copy_aux_vec_sca(VEC<double>& displayVEC) = 0;

	//Get settings for module display data 
	//Return module displaying its effective field (MOD_ALL means total Heff)
	int Get_Module_Heff_Display(void);
	int Get_ActualModule_Heff_Display(void);
	//Return module displaying its energy density spatial variation (MOD_ERROR means none set)
	int Get_Module_Energy_Display(void);
	int Get_ActualModule_Energy_Display(void);

	//----------------------------------- MESH INFO GET/SET METHODS

	//NOTE : these type of methods are required here so we can access these properties in .cu files, where we cannot use a Mesh pointer
	//Using a Mesh pointer means having to include Mesh.h, thus also BorisLib.h -> this won't compile since nvcc doesn't recognise C++14 code in BorisLib
	//A future version of nvcc will make this practice redundant

	int GetMeshType(void);

	//search save data list (saveDataList) for given dataID set for this mesh. Return true if found and its rectangle is not Null; else return false.
	bool IsOutputDataSet_withRect(int datumId);
	//return true if data is set (with any rectangle)
	bool IsOutputDataSet(int datumId);

	//check if given stage is set
	bool IsStageSet(int stageType);

	//set computefields_if_MC flag on SuperMesh
	void Set_Force_MonteCarlo_ComputeFields(bool status);

	//others

	bool is_atomistic(void);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//magnetization dynamics computation enabled
	virtual bool MComputation_Enabled(void) = 0;

	virtual bool Magnetism_Enabled(void) = 0;

	//electrical conduction computation enabled
	virtual bool EComputation_Enabled(void) = 0;

	//thermal conduction computation enabled
	virtual bool TComputation_Enabled(void) = 0;

	//mechanical computation enabled
	virtual bool MechComputation_Enabled(void) = 0;

	//check if interface conductance is enabled (for spin transport solver)
	virtual bool GInterface_Enabled(void) = 0;

	virtual bool GetMeshExchangeCoupling(void) = 0;

	//----------------------------------- VALUE GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	virtual cuBReal GetTopologicalCharge(cuRect rectangle) = 0;

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	virtual void Compute_TopoChargeDensity(void) = 0;

	//check if the ODECommon::available flag is true (ode step solved)
	bool CurrentTimeStepSolved(void);

	//check evaluation speedup flag in ODECommon
	int GetEvaluationSpeedup(void);

	//check in ODECommon the type of field update we need to do depending on the ODE evaluation step
	int Check_Step_Update(void);

	//get total time with evaluation step resolution level
	cuBReal Get_EvalStep_Time(void);

	cuBReal GetStageTime(void);
	int GetStageStep(void);

	cuBReal GetTimeStep(void);

	//----------------------------------- MESH SHAPE CONTROL

	//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
	virtual BError copy_shapes_from_cpu(void) = 0;

	//copy all meshes controlled using change_mesh_shape from gpu to cpu versions
	virtual BError copy_shapes_to_cpu(void) = 0;
};

#endif
