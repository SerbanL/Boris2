#include "MeshBaseCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

__global__ void Zero_aux_values_kernels(cuBReal& aux_real, cuReal3& aux_real3, size_t& aux_int)
{
	if (threadIdx.x == 0) aux_real = 0.0;
	if (threadIdx.x == 1) aux_real3 = cuReal3();
	if (threadIdx.x == 2) aux_int = 0.0;
}

//zero all single aux avalues
void MeshBaseCUDA::Zero_aux_values(void)
{
	Zero_aux_values_kernels <<< 1, CUDATHREADS >>> (aux_real, aux_real3, aux_int);
}

#endif

