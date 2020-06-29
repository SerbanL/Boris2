#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.cuh"

__global__ void Zero_Errors_kernel(cuBReal& max_error, cuBReal& max_value)
{
	if (threadIdx.x == 0) max_error = 0.0;
	else if (threadIdx.x == 1) max_value = 0.0;
}

void STransportCUDA::Zero_Errors(void)
{
	Zero_Errors_kernel <<< 1, CUDATHREADS >>> (max_error, max_value);
}

#endif

#endif