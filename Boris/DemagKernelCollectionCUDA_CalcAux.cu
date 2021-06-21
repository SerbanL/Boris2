#include "DemagKernelCollectionCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include <cuda_runtime.h>

//Auxiliary for kernel computations on the GPU

//--------------------------

//copy Re or Im parts of cuOut to cuIn

__global__ void cuOut_to_cuIn_Re_collection_kernel(size_t size, cufftDoubleReal* cuIn, cufftDoubleComplex* cuOut)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < size) {

		cuIn[idx] = cuOut[idx].x;
	}
}

void DemagKernelCollectionCUDA::cuOut_to_cuIn_Re(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut)
{
	cuOut_to_cuIn_Re_collection_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, cuIn, cuOut);
}

__global__ void cuOut_to_cuIn_Im_collection_kernel(size_t size, cufftDoubleReal* cuIn, cufftDoubleComplex* cuOut, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < size) {

		cuIn[idx] = cuOut[idx].y;
	}
}

void DemagKernelCollectionCUDA::cuOut_to_cuIn_Im(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut)
{
	cuOut_to_cuIn_Im_collection_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, cuIn, cuOut, cuN);
}

//--------------------------

#endif
#endif