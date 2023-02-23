#include "DemagKernelCollectionCUDA_KerType.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include <cuda_runtime.h>

//--------------------------

//Copy Re parts of cuOut to Kdiag component (1: Kx, 2: Ky, 3: Kz)

__global__ void Set_Kdiag_realcomponent_kernel(cuVEC<cuReal3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		if (component == 1) Kdiag[idx].x = cuOut[idx].x;
		else if (component == 2) Kdiag[idx].y = cuOut[idx].x;
		else if (component == 3) Kdiag[idx].z = cuOut[idx].x;
	}
}

__global__ void Set_Kdiag_realcomponent_transpose_kernel(cuVEC<cuReal3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
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

void cuKerType::Set_Kdiag_real(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_Kdiag_realcomponent_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (Kdiag_real, cuOut, N, component);
	}
	else {

		Set_Kdiag_realcomponent_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag_real, cuOut, N, component);
	}
}

//--------------------------

__global__ void Set_K2D_odiag_kernel(cuVEC<cuBReal>& K2D_odiag, cufftDoubleComplex* cuOut, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1)) {

		K2D_odiag[idx] = -cuOut[idx].y;
	}
}

__global__ void Set_K2D_odiag_transpose_kernel(cuVEC<cuBReal>& K2D_odiag, cufftDoubleComplex* cuOut, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1)) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % (N.y / 2 + 1);

		K2D_odiag[j + i * (N.y / 2 + 1)] = -cuOut[idx].y;
	}
}

//Copy -(Re, Im, Im) parts of cuOut to Kodiag component (1: Kxy, 2: Kxz, 3: Kyz). Takes into account transpose_xy flag.
void cuKerType::Set_K2D_odiag(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int transpose_xy)
{
	if (!transpose_xy) {

		Set_K2D_odiag_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (K2D_odiag, cuOut, N);
	}
	else {

		Set_K2D_odiag_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (K2D_odiag, cuOut, N);
	}
}

//--------------------------

//Copy -Im, Re, Im parts of cuOut to Kodiag component (1: Kxy, 2: Kxz, 3: Kyz)

__global__ void Set_Kodiag_real_2D_component_kernel(cuVEC<cuReal3>& Kodiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1)) {

		if (component == 1) Kodiag[idx].x = -cuOut[idx].y;
		else if (component == 2) Kodiag[idx].y = cuOut[idx].x;
		else if (component == 3) Kodiag[idx].z = cuOut[idx].y;
	}
}

__global__ void Set_Kodiag_real_2D_component_transpose_kernel(cuVEC<cuReal3>& Kodiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1)) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % (N.y / 2 + 1);

		if (component == 1) Kodiag[j + i * (N.y / 2 + 1)].x = -cuOut[idx].y;
		else if (component == 2) Kodiag[j + i * (N.y / 2 + 1)].y = cuOut[idx].x;
		else if (component == 3) Kodiag[j + i * (N.y / 2 + 1)].z = cuOut[idx].y;
	}
}

void cuKerType::Set_Kodiag_real_2D(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_Kodiag_real_2D_component_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_real, cuOut, N, component);
	}
	else {

		Set_Kodiag_real_2D_component_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_real, cuOut, N, component);
	}
}

//Copy -Re, -Im, -Im parts of cuOut to Kodiag component (1: Kxy, 2: Kxz, 3: Kyz)

__global__ void Set_Kodiag_real_3D_component_kernel(cuVEC<cuReal3>& Kodiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		if (component == 1) Kodiag[idx].x = -cuOut[idx].x;
		else if (component == 2) Kodiag[idx].y = -cuOut[idx].y;
		else if (component == 3) Kodiag[idx].z = -cuOut[idx].y;
	}
}

__global__ void Set_Kodiag_real_3D_component_transpose_kernel(cuVEC<cuReal3>& Kodiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
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

void cuKerType::Set_Kodiag_real_3D(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_Kodiag_real_3D_component_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_real, cuOut, N, component);
	}
	else {

		Set_Kodiag_real_3D_component_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_real, cuOut, N, component);
	}
}

//--------------------------

//Copy complex to complex for given kernel component

__global__ void Set_K_cmpl_kernel(cuVEC<cuReIm3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y * N.z) {

		if (component == 1) Kdiag[idx].x = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 2) Kdiag[idx].y = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 3) Kdiag[idx].z = cuReIm(cuOut[idx].x, cuOut[idx].y);
	}
}

__global__ void Set_K_cmpl_transpose_kernel(cuVEC<cuReIm3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y * N.z) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % N.y;
		int k = idx / ((N.x / 2 + 1) * N.y);

		if (component == 1) Kdiag[j + i * N.y + k * (N.x / 2 + 1) * N.y].x = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 2) Kdiag[j + i * N.y + k * (N.x / 2 + 1) * N.y].y = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 3) Kdiag[j + i * N.y + k * (N.x / 2 + 1) * N.y].z = cuReIm(cuOut[idx].x, cuOut[idx].y);
	}
}

void cuKerType::Set_Kdiag_cmpl(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_K_cmpl_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag_cmpl, cuOut, N, component);
	}
	else {

		Set_K_cmpl_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag_cmpl, cuOut, N, component);
	}
}

void cuKerType::Set_Kodiag_cmpl(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_K_cmpl_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_cmpl, cuOut, N, component);
	}
	else {

		Set_K_cmpl_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_cmpl, cuOut, N, component);
	}
}

//--------------------------

//Copy complex to complex for given kernel component

__global__ void Set_K_cmpl_reduced_kernel(cuVEC<cuReIm3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		if (component == 1) Kdiag[idx].x = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 2) Kdiag[idx].y = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 3) Kdiag[idx].z = cuReIm(cuOut[idx].x, cuOut[idx].y);
	}
}

__global__ void Set_K_cmpl_reduced_transpose_kernel(cuVEC<cuReIm3>& Kdiag, cufftDoubleComplex* cuOut, cuSZ3& N, int component)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1)) {

		int i = idx % (N.x / 2 + 1);
		int j = (idx / (N.x / 2 + 1)) % (N.y / 2 + 1);
		int k = idx / ((N.x / 2 + 1) * (N.y / 2 + 1));

		if (component == 1) Kdiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].x = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 2) Kdiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].y = cuReIm(cuOut[idx].x, cuOut[idx].y);
		else if (component == 3) Kdiag[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)].z = cuReIm(cuOut[idx].x, cuOut[idx].y);
	}
}

void cuKerType::Set_Kdiag_cmpl_reduced(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_K_cmpl_reduced_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag_cmpl, cuOut, N, component);
	}
	else {

		Set_K_cmpl_reduced_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kdiag_cmpl, cuOut, N, component);
	}
}

void cuKerType::Set_Kodiag_cmpl_reduced(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy)
{
	if (!transpose_xy) {

		Set_K_cmpl_reduced_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_cmpl, cuOut, N, component);
	}
	else {

		Set_K_cmpl_reduced_transpose_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Kodiag_cmpl, cuOut, N, component);
	}
}

#endif
#endif