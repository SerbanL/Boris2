#include "SDemagCUDA_KernelCollection.h"

//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include <cuda_runtime.h>

//---------------------------------------------------------------------

//2D

__global__ void cu_KernelMultiplication_Launcher_2D_Self(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < size && idx < (N.x / 2 + 1) * N.y) {

		kernel.cu_KernelMultiplication_2D_Self_Set(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
	}
}

__global__ void cu_KernelMultiplication_Launcher_2D_Self_transpose_xy(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < size && idx < (N.x / 2 + 1) * N.y) {

		kernel.cu_KernelMultiplication_2D_Self_transpose_xy_Set(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
	}
}

__global__ void cu_KernelMultiplication_Launcher_2D_zShifted(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	bool* inverse_zshifted,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < size && idx < (N.x / 2 + 1) * N.y) {

		if (inverse_zshifted[n]) {

			kernel.cu_KernelMultiplication_2D_inversezShifted_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
		}
		else {

			kernel.cu_KernelMultiplication_2D_zShifted_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
		}
	}
}

__global__ void cu_KernelMultiplication_Launcher_2D_zShifted_transpose_xy(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	bool* inverse_zshifted,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < size && idx < (N.x / 2 + 1) * N.y) {

		if (inverse_zshifted[n]) {

			kernel.cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
		}
		else {

			kernel.cu_KernelMultiplication_2D_zShifted_transpose_xy_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
		}
	}
}

__global__ void cu_KernelMultiplication_Launcher_2D_Regular(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < size && idx < (N.x / 2 + 1) * N.y) {

		kernel.cu_KernelMultiplication_2D_Regular_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
	}
}

//3D

__global__ void cu_KernelMultiplication_Launcher_3D_Self_transpose_xy(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y * N.z);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y * N.z);

	if (n < size && idx < (N.x / 2 + 1) * N.y * N.z) {

		kernel.cu_KernelMultiplication_3D_Self_transpose_xy_Set(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
	}
}

__global__ void cu_KernelMultiplication_Launcher_3D_zShifted_transpose_xy(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y * N.z);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y * N.z);

	if (n < size && idx < (N.x / 2 + 1) * N.y * N.z) {

		kernel.cu_KernelMultiplication_3D_zShifted_transpose_xy_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
	}
}

__global__ void cu_KernelMultiplication_Launcher_3D_Regular(
	cuKerType& kernel,
	cuComplex** Inx, cuComplex** Iny, cuComplex** Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	int& size)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	cuSZ3& N = kernel.N;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y * N.z);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y * N.z);

	if (n < size && idx < (N.x / 2 + 1) * N.y * N.z) {

		kernel.cu_KernelMultiplication_3D_Regular_Add(idx, Inx[n], Iny[n], Inz[n], Outx[n], Outy[n], Outz[n]);
	}
}

//---------------------------------------------------------------------

void KerTypeCollectionCUDA::Kernel_Multiplication_2D(bool transpose_xy)
{
	if (internal_demag) {

		if (transpose_xy) {

			cu_KernelMultiplication_Launcher_2D_Self_transpose_xy <<< ((N.x / 2 + 1) * N.y * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernel,
				InCol_x, InCol_y, InCol_z,
				OutCol_x, OutCol_y, OutCol_z,
				size);
		}
		else {

			cu_KernelMultiplication_Launcher_2D_Self <<< ((N.x / 2 + 1) * N.y * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernel,
				InCol_x, InCol_y, InCol_z,
				OutCol_x, OutCol_y, OutCol_z,
				size);
		}
	}
	else {

		if (zshifted) {

			if (transpose_xy) {

				cu_KernelMultiplication_Launcher_2D_zShifted_transpose_xy <<< ((N.x / 2 + 1) * N.y * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
					*kernel,
					InCol_x, InCol_y, InCol_z,
					OutCol_x, OutCol_y, OutCol_z,
					inverse_shifted,
					size);
			}
			else {

				cu_KernelMultiplication_Launcher_2D_zShifted <<< ((N.x / 2 + 1) * N.y * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
					*kernel,
					InCol_x, InCol_y, InCol_z,
					OutCol_x, OutCol_y, OutCol_z,
					inverse_shifted,
					size);
			}
		}
		else {

			cu_KernelMultiplication_Launcher_2D_Regular <<< ((N.x / 2 + 1) * N.y * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernel,
				InCol_x, InCol_y, InCol_z,
				OutCol_x, OutCol_y, OutCol_z,
				size);
		}
	}
}

void KerTypeCollectionCUDA::Kernel_Multiplication_3D(void)
{
	if (internal_demag) {

		cu_KernelMultiplication_Launcher_3D_Self_transpose_xy <<< ((N.x / 2 + 1) * N.y * N.z * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
			*kernel,
			InCol_x, InCol_y, InCol_z,
			OutCol_x, OutCol_y, OutCol_z,
			size);
	}
	else {

		if (zshifted) {

			cu_KernelMultiplication_Launcher_3D_zShifted_transpose_xy <<< ((N.x / 2 + 1) * N.y * N.z * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernel,
				InCol_x, InCol_y, InCol_z,
				OutCol_x, OutCol_y, OutCol_z,
				size);
		}
		else {

			cu_KernelMultiplication_Launcher_3D_Regular <<< ((N.x / 2 + 1) * N.y * N.z * size_cpu + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernel,
				InCol_x, InCol_y, InCol_z,
				OutCol_x, OutCol_y, OutCol_z,
				size);
		}
	}
}

#endif

#endif