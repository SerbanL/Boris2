#include "DemagKernelCUDA.h"

#if COMPILECUDA == 1

#if defined MODULE_COMPILATION_DEMAG || defined MODULE_COMPILATION_SDEMAG

#include <cuda_runtime.h>


//Auxiliary for kernel computations on the GPU

//--------------------------

//copy Re or Im parts of cuOut to cuIn

__global__ void cuOut_to_cuIn_Re_kernel(size_t size, cufftDoubleReal* cuIn, cufftDoubleComplex* cuOut)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < size) {

		cuIn[idx] = cuOut[idx].x;
	}
}

void DemagKernelCUDA::cuOut_to_cuIn_Re(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut)
{
	cuOut_to_cuIn_Re_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, cuIn, cuOut);
}

__global__ void cuOut_to_cuIn_Im_kernel(size_t size, cufftDoubleReal* cuIn, cufftDoubleComplex* cuOut, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < size) {

		cuIn[idx] = cuOut[idx].y;
	}
}

void DemagKernelCUDA::cuOut_to_cuIn_Im(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut)
{
	cuOut_to_cuIn_Im_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, cuIn, cuOut, cuN);
}

//--------------------------

//Copy Re parts of cuOut to Kdiag component (1: Kx, 2: Ky, 3: Kz)

__global__ void cuOut_to_Kdiagcomponent_kernel(cuVEC<cuReal3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		if (component == 1) Kdiag[idx].x = cuOut[idx].x;
		else if (component == 2) Kdiag[idx].y = cuOut[idx].x;
		else if (component == 3) Kdiag[idx].z = cuOut[idx].x;
	}
}

__global__ void cuOut_to_Kdiagcomponent_transpose_kernel(cuVEC<cuReal3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % (N.y / 2 + 1);
		int k = idx / ((N.x / 2 + 1) * (N.y / 2 + 1));

		if (component == 1) Kdiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].x = cuOut[idx].x;
		else if (component == 2) Kdiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].y = cuOut[idx].x;
		else if (component == 3) Kdiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].z = cuOut[idx].x;
	}
}

void DemagKernelCUDA::cuOut_to_Kdiagcomponent(cu_arr<cufftDoubleComplex>& cuOut, int component)
{
	if (!transpose_xy) {

		cuOut_to_Kdiagcomponent_kernel <<< ((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag, cuOut, cuN, component);
	}
	else {

		cuOut_to_Kdiagcomponent_transpose_kernel <<< ((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag, cuOut, cuN, component);
	}
}

//--------------------------

//Copy -Im parts of cuOut to K2D_odiag

__global__ void cuOut_to_K2D_odiag_kernel(cuVEC<cuBReal>& K2D_odiag, cufftDoubleComplex* cuOut, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1)) {

		K2D_odiag[idx] = -cuOut[idx].y;
	}
}

__global__ void cuOut_to_K2D_odiag_transpose_kernel(cuVEC<cuBReal>& K2D_odiag, cufftDoubleComplex* cuOut, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1)) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % (N.y / 2 + 1);

		K2D_odiag[j + i * (N.y / 2 + 1)] = -cuOut[idx].y;
	}
}

void DemagKernelCUDA::cuOut_to_K2D_odiag(cu_arr<cufftDoubleComplex>& cuOut)
{
	if (!transpose_xy) {

		cuOut_to_K2D_odiag_kernel <<< ((N.x / 2 + 1) * (N.y / 2 + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (K2D_odiag, cuOut, cuN);
	}
	else {

		cuOut_to_K2D_odiag_transpose_kernel <<< ((N.x / 2 + 1) * (N.y / 2 + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (K2D_odiag, cuOut, cuN);
	}
}

//--------------------------

__global__ void cuOut_to_Kodiagcomponent_kernel(cuVEC<cuReal3>& Kodiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		if (component == 1) Kodiag[idx].x = -cuOut[idx].x;
		else if (component == 2) Kodiag[idx].y = -cuOut[idx].y;
		else if (component == 3) Kodiag[idx].z = -cuOut[idx].y;
	}
}

__global__ void cuOut_to_Kodiagcomponent_transpose_kernel(cuVEC<cuReal3>& Kodiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % (N.y / 2 + 1);
		int k = idx / ((N.x / 2 + 1) * (N.y / 2 + 1));

		if (component == 1) Kodiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].x = -cuOut[idx].x;
		else if (component == 2) Kodiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].y = -cuOut[idx].y;
		else if (component == 3) Kodiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].z = -cuOut[idx].y;
	}
}

//Copy -(Re, Im, Im) parts of cuOut to Kodiag component (1: Kxy, 2: Kxz, 3: Kyz). Takes into account transpose_xy flag.
void DemagKernelCUDA::cuOut_to_Kodiagcomponent(cu_arr<cufftDoubleComplex>& cuOut, int component)
{
	if (!transpose_xy) {

		cuOut_to_Kodiagcomponent_kernel <<< ((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag, cuOut, cuN, component);
	}
	else {

		cuOut_to_Kodiagcomponent_transpose_kernel <<< ((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag, cuOut, cuN, component);
	}
}

#endif
#endif