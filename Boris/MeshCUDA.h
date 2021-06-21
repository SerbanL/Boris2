#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "MeshBaseCUDA.h"

#include "MeshParamsCUDA.h"
#include "ManagedMeshCUDA.h"

class Mesh;

//Store Mesh quantities as cu_obj managed cuda VECs
class MeshCUDA :
	public MeshBaseCUDA,
	public MeshParamsCUDA,
	public MeshDisplayCUDA
{
private:

	//pointer to cpu version of this mesh (base) - need it in the destructor. 
	//Still keep references to some Mesh data members here as we cannot use pMesh in .cu files (cannot have BorisLib.h in those compilation units - real headache, will need to fix this at some point somehow: problem is the nvcc compiler throws errors due to C++14 code in BorisLib)
	Mesh *pMesh;

public:

	//This points to data in ZeemanCUDA if set (else nullptr) - managed by ZeemanCUDA
	//pHa is nullptr if ZeemanCUDA not active
	cu_obj<cuReal3>* pHa = nullptr;

	//Managed Mesh
	cu_obj<ManagedMeshCUDA> cuMesh;

	//stores MOD_ identifiers of modules active in this mesh; used for Monte-Carlo energy methods so we can identify which contributions to include.
	cu_arr<int> cuModules;
	//number of elements stored in cuModules
	cu_obj<int> cuNumModules;

	//-----Magnetic properties

	//Magnetization
	cu_obj<cuVEC_VC<cuReal3>> M;

	//Additional magnetization used for antiferromagnetic meshes with 2 sub-lattice local approximation; exactly same dimensions as M
	cu_obj<cuVEC_VC<cuReal3>> M2;

	//effective field - sum total field of all the added modules
	cu_obj<cuVEC<cuReal3>> Heff;

	//Additional effective field used for antiferromagnetic meshes with 2 sub-lattice local approximation; exactly same dimensions as Heff
	cu_obj<cuVEC<cuReal3>> Heff2;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//In MeshBaseCUDA

	//-----Thermal conduction properties

	//In MeshBaseCUDA

	//-----Stochastic cellsize (VECs held in DiffEq)

	//number of cells for stochastic VECs
	SZ3& n_s;

	//cellsize for stochastic VECs
	DBL3& h_s;

	//link stochastic cellsize to magnetic cellsize by default (set this to false if you want to control h_s independently)
	bool& link_stochastic;

	//-----Elastic properties

	//In MeshBaseCUDA

public:

	//------------------------CTOR/DTOR

	//make this object by copying data from the Mesh holding this object
	MeshCUDA(Mesh* pMesh);

	virtual ~MeshCUDA();

	//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

	PhysQ FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground);

	//save the quantity currently displayed on screen in an ovf2 file using the specified format
	BError SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType);

	//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
	//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
	//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
	void GetPhysicalQuantityProfile(DBL3 start, DBL3 end, double step, DBL3 stencil, std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, bool do_average, bool read_average);

	//return average value for currently displayed mesh quantity in the given relative rectangle
	Any GetAverageDisplayedMeshValue(Rect rel_rect);

	//copy aux_vec_sca in GPU memory to displayVEC in CPU memory
	void copy_aux_vec_sca(VEC<double>& displayVEC);

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

	virtual bool GetMeshExchangeCoupling(void) { return false; }

	//----------------------------------- VALUE GETTERS

	//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
	cuBReal GetTopologicalCharge(cuRect rectangle);

	//compute topological charge density spatial dependence and have it available in aux_vec_sca
	//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
	void Compute_TopoChargeDensity(void);

	//----------------------------------- MESH SHAPE CONTROL

	//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
	BError copy_shapes_from_cpu(void);

	//copy all meshes controlled using change_mesh_shape from gpu to cpu versions
	BError copy_shapes_to_cpu(void);

	//-----------------------------------OBJECT GETTERS

	cu_obj<ManagedDiffEq_CommonCUDA>& Get_ManagedDiffEq_CommonCUDA(void);

	std::vector<DBL4>& get_tensorial_anisotropy(void);
	std::vector<DBL4>& get_tensorial_anisotropy2(void);
};

#endif