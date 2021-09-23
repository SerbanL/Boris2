#pragma once

#include "BorisLib.h"
#include "Modules.h"



class SuperMesh;

#ifdef MODULE_COMPILATION_SDEMAG

#include "Convolution.h"
#include "DemagKernel.h"
#include "DemagKernelCollection.h"

#include "SDemag_Demag.h"

#if COMPILECUDA == 1
#include "SDemagCUDA.h"
#include "SDemagCUDA_Demag.h"
#endif


////////////////////////////////////////////////////////////////////////////////////////////////
//
// Full demag on super-mesh, involving interpolating to/from multiple magnetic meshes

class SDemag :
	public Modules,
	public Convolution<SDemag, DemagKernel>,
	public ProgramState<SDemag, std::tuple<bool, SZ3, bool, int, INT3>, std::tuple<>>
{
#if COMPILECUDA == 1
	friend SDemagCUDA;
	friend SDemagCUDA_Demag;
#endif

	friend SDemag_Demag;

private:

	SuperMesh * pSMesh;

	//--------

	//super-mesh magnetization values used for computing demag field on the super-mesh
	VEC<DBL3> sm_Vals;

	//with supermesh convolution we need to know if any antiferromagnetic meshes are present so we can adjust the type of transfers used
	bool antiferromagnetic_meshes_present = false;

	//number of non-empty cells in sm_Vals
	int non_empty_cells;

	//total non-empty volume from all meshes participating in convolution
	double total_nonempty_volume = 0.0;

	//--------

	//the default method is to use the multilayered convolution, but you can revert to using convolution on the super-mesh directly
	//for convolution on the super-mesh to be exact, the super-mesh cells must either be empty or be covered fully by a single input mesh
	bool use_multilayered_convolution = true;

	//collection of all SDemag_Demag modules in individual (anti)ferromagnetic meshes - these modules are only created by this SDemag module
	std::vector<SDemag_Demag*> pSDemag_Demag;

	//collect FFT input spaces : after Forward FFT the ffts of M from the individual meshes will be found here
	//These are used as inputs to kernel multiplications. Same order as pSDemag_Demag.
	std::vector<VEC<ReIm3>*> FFT_Spaces_Input;

	//collection of rectangles of meshes, same ordering as for pSDemag_Demag and FFT_Spaces, used in multi-layered convolution
	//these are not necessarily the rectangles of the input M meshes, but are the rectangles of the transfer meshes (M -> transfer -> convolution)
	//for instance in 3D mode, all rectangles in multi-layered convolution must have same size
	//in 2D mode the rectangles can differ in thickness but must have the same xy size
	//thus in 3D mode find largest one and extend all the other rectangles to match (if possible try to have them overlapping in xy-plane projections so we can use kernel symmetries)
	//in 2D mode find largest xy dimension and extend all xy dimensions -> again try to overlap their xy-plane projections
	std::vector<Rect> Rect_collection;

	//demag kernels used for multilayered convolution, one collection per mesh/SDemag_Demag module. Don't recalculate redundant kernels in the collection.
	std::vector<DemagKernelCollection*> kernel_collection;

	//common discretisation number of cells
	//int multi-layered convolution, all layers must have exactly the same discretisation number of cells. This can be changed by the user but the default setting should suffice for the vast majority of cases.
	SZ3 n_common = SZ3(1);

	//keep calculating a default value for n when configuration changes, as opposed to the user having to set n_common
	bool use_default_n = true;

	//0 : don't use 2D
	//1 : always use 2D multi-layered convolution irrespective of how the individual meshes are discretised, i.e. each mesh will be taken as 2D
	//2 : each mesh will be layared along z direction as 2D layers. Thus this is a fully 2D layered convolution but more accurate than option 1.
	int force_2d_convolution = 0;

	//number of pbc images in each dimension (set to zero to disable).
	//There is also a copy of this in ConvolutionData inherited from Convolution - we need another copy here to detect changes
	//these pbc images are applicable for supermesh or multilayered convolution only - in the case of multilayered convolution all meshes have the same pbc conditions applied.
	INT3 demag_pbc_images = INT3();

	//Evaluation speedup mode data

	//times at which evaluations were done, used for extrapolation
	double time_demag1 = 0.0, time_demag2 = 0.0, time_demag3 = 0.0, time_demag4 = 0.0, time_demag5 = 0.0, time_demag6 = 0.0;

	int num_Hdemag_saved = 0;

private:

	//------------------ Helpers for multi-layered convolution control

	//when SDemag created, it needs to add one SDemag_Demag module to each (anti)ferromagnetic mesh (or multiple if the 2D layering option is enabled).
	BError Create_SDemag_Demag_Modules(void);

	//delete all SDemag_Demag modules - these are only created by SDemag
	void Destroy_SDemag_Demag_Modules(void);

	//make sure the pSDemag_Demag list is up to date : if any mismatches found return false
	bool Update_SDemag_Demag_List(void);

	//set default value for n_common : largest value from all (anti)ferromagnetic meshes
	void set_default_n_common(void);

	//get and adjust Rect_collection so the rectangles are matching for multi-layered convolution. n_common should be calculated already.
	//do this only before initialization
	void set_Rect_collection(void);

	//get maximum cellsize for multi-layered convolution (use it to normalize dimensions)
	double get_maximum_cellsize(void);

	//get convolution rectangle for the given SDemag_Demag module (remember this might not be the rectangle of M in that mesh, but an adjusted rectangle to make the convolution work)
	Rect get_convolution_rect(SDemag_Demag* demag_demag);

	//initialize transfer object for supermesh convolution
	BError Initialize_Mesh_Transfer(void);

	//Set PBC settings for M in all meshes
	BError Set_Magnetic_PBC(void);

public:

	SDemag(SuperMesh *pSMesh_);
	~SDemag();

	//-------------------Methods associated with saving/loading simulations

	void RepairObjectState(void) 
	{ 
		Uninitialize();

		//make sure convolution size and memory is allocated correctly
		UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
	}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	//if in multi-layered convolution mode this method uninitializes all SDemag_Demag modules
	//IMPORTANT : Only call this after UpdateConfiguration has been called (but do not put UpdateConfiguration in this method as it can also be called from within UpdateConfiguration)
	void UninitializeAll(void);

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Setters

	//change between demag calculation types : super-mesh (status = false) or multilayered (status = true)
	BError Set_Multilayered_Convolution(bool status);

	//enable multi-layered convolution and force it to 2D for all layers (1), or as 2D layers in each mesh (2). 0 to turn off any forced 2D layering.
	BError Set_2D_Multilayered_Convolution(int status);

	//set n_common for multi-layered convolution
	BError Set_n_common(SZ3 n);

	//set status for use_default_n
	BError Set_Default_n_status(bool status);

	//Set PBC images for supermesh demag
	BError Set_PBC(INT3 demag_pbc_images_);

	//-------------------Getters

	VEC<DBL3>& GetDemagField(void) 
	{ 
		return sm_Vals; 
	}

	//Get PBC images for supermesh demag
	INT3 Get_PBC(void) { return demag_pbc_images; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetDemagFieldCUDA(void) { return dynamic_cast<SDemagCUDA*>(pModuleCUDA)->GetDemagField(); }
#endif

	//getters for multi-layered convolution
	bool Get_Multilayered_Convolution_Status(void) { return use_multilayered_convolution; }
	int Get_2D_Multilayered_Convolution_Status(void) { return force_2d_convolution; }

	bool Use_Default_n_Status(void) { return use_default_n; }

	SZ3 Get_n_common(void) { return n_common; }
};

#else

class SDemag :
	public Modules
{

private:

	//super-mesh magnetization values used for computing demag field on the super-mesh
	VEC<DBL3> sm_Vals;

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>> sm_Vals_CUDA;
#endif

private:

public:

	SDemag(SuperMesh *pSMesh_) {}
	~SDemag() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Setters

	//change between demag calculation types : super-mesh (status = false) or multilayered (status = true)
	BError Set_Multilayered_Convolution(bool status) { return BError(); }

	//enable multi-layered convolution and force it to 2D for all layers
	BError Set_2D_Multilayered_Convolution(int status) { return BError(); }

	//set n_common for multi-layered convolution
	BError Set_n_common(SZ3 n) { return BError(); }

	//set status for use_default_n
	BError Set_Default_n_status(bool status) { return BError(); }

	BError Set_PBC(INT3 demag_pbc_images_) { return BError(); }

	//-------------------Getters

	VEC<DBL3>& GetDemagField(void) { return sm_Vals; }

	//Get PBC images for supermesh demag
	INT3 Get_PBC(void) { return INT3(); }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetDemagFieldCUDA(void) { return sm_Vals_CUDA; }
#endif

	//getters for multi-layered convolution
	bool Get_Multilayered_Convolution_Status(void) { return true; }
	int Get_2D_Multilayered_Convolution_Status(void) { return 0; }

	bool Use_Default_n_Status(void) { return true; }

	SZ3 Get_n_common(void) { return SZ3(); }
};

#endif