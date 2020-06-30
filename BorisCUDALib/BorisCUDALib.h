#pragma once

#include <cufft.h>
#include <cuda_runtime.h>

//Can include BorisCUDALib.h in both C++ code and CUDA C code files

#include "launchers.h"
#include "cuFuncs_Aux.h"
#include "CUDAError.h"
#include "alloc_cpy.h"

#include "cuObject.h"
#include "cu_prng.h"
#include "cuTypes.h"

#include "cuArray.h"
#include "cuArray_mng_gpu.h"
#include "cuArray_mng_cpu.h"
#include "cuArray_aux.h"
#include "cuArray_sizing.h"
#include "cuArray_transfer.h"


#include "TEquationCUDA.h"
#include "cuObj_Math_Special.h"

#include "cuFuncs_Math.h"

#include "cuVEC.h"
#include "cuVEC_mng.h"
#include "cuVEC_aux.h"
#include "cuVEC_oper.h"
#include "cuVEC_avg.h"
#include "cuVEC_MeshTransfer.h"

#include "cuVEC_VC.h"
#include "cuVEC_VC_mng.h"
#include "cuVEC_VC_aux.h"
#include "cuVEC_VC_flags.h"
#include "cuVEC_VC_shape.h"
#include "cuVEC_VC_cmbnd.h"
#include "cuVEC_VC_oper.h"
#include "cuVEC_VC_del.h"
#include "cuVEC_VC_diff2.h"
#include "cuVEC_VC_grad.h"
#include "cuVEC_VC_div.h"
#include "cuVEC_VC_curl.h"
#include "cuVEC_VC_ngbrsum.h"
