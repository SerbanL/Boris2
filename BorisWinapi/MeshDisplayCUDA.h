#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

//cannot include BorisLib.h (or even just VEC_VC.h and VEC.h) since it contains C++14 (ProgramState and Introspection in particular) code which nvcc doesn't compile.
//Since this header file inevitably gets included in .cu files (compiled with nvcc) then must forward declare VEC and VEC_VC here instead.
//To avoid incomplete type error then we must only declare pointers to these types in the class definition - otherwise how can the compiler know how much memory to allocate on object creation?
template <typename VType> class VEC;
template <typename VType> class VEC_VC;

//this class handles display-related conversions from cuVECs to VECs - take mesh data from gpu memory and prepare it in cpu memory so we can display it on screen
class MeshDisplayCUDA {

private:

	//-----Display cuVECs

	//extract quantity displayed in mesh to one of these at the current detail level, then copy over to cpu versions for display calculations

	//vectorial cuVEC_VC
	cu_obj<cuVEC_VC<cuReal3>> display_vec_vc_vec_CUDA;

	//vectorial cuVEC
	cu_obj<cuVEC<cuReal3>> display_vec_vec_CUDA;

	//scalar cuVEC_VC
	cu_obj<cuVEC_VC<cuBReal>> display_vec_vc_sca_CUDA;

	//scalar cuVEC
	cu_obj<cuVEC<cuBReal>> display_vec_sca_CUDA;

protected:

	//VECs for display - cuVEC data is extracted and placed in these
#if SINGLEPRECISION == 1
	//vectorial cuVEC_VC
	VEC_VC<FLT3>* pdisplay_vec_vc_vec;

	//vectorial cuVEC
	VEC<FLT3>* pdisplay_vec_vec;

	//scalar cuVEC_VC
	VEC_VC<float>* pdisplay_vec_vc_sca;

	//scalar cuVEC
	VEC<float>* pdisplay_vec_sca;
#else
	//vectorial cuVEC_VC
	VEC_VC<DBL3>* pdisplay_vec_vc_vec;

	//vectorial cuVEC
	VEC<DBL3>* pdisplay_vec_vec;

	//scalar cuVEC_VC
	VEC_VC<double>* pdisplay_vec_vc_sca;

	//scalar cuVEC
	VEC<double>* pdisplay_vec_sca;
#endif

protected:

	MeshDisplayCUDA(void);
	virtual ~MeshDisplayCUDA();

	//given type of mesh display quantity prepare appropriate display vec : 
	//extract from mesh display (in gpu memory) a coarser mesh display (still in gpu memory), then copy the coarse display to cpu memory ready to use further.
	//n_quantity and meshRect are the same as those in cu_obj_quantity : pass in values in cpu memory to speed up function execution
	bool prepare_display(SZ3 n_quantity, Rect meshRect, double detail_level, cu_obj<cuVEC_VC<cuReal3>>& cu_obj_quantity);
	bool prepare_display(SZ3 n_quantity, Rect meshRect, double detail_level, cu_obj<cuVEC_VC<cuBReal>>& cu_obj_quantity);
	bool prepare_display(SZ3 n_quantity, Rect meshRect, double detail_level, cu_obj<cuVEC<cuReal3>>& cu_obj_quantity);
	bool prepare_display(SZ3 n_quantity, Rect meshRect, double detail_level, cu_obj<cuVEC<cuBReal>>& cu_obj_quantity);
};

#endif