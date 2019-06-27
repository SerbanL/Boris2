#include "DemagKernelCollectionCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include <cuda_runtime.h>

#include "DemagKernelCollectionCUDA_KerType.h"

//THESE ARE IN CURRENT USE

//-------------------------- STAND-ALONE VERSIONS FOR MULTIPLE INPUTS

//-------------------------- CONVOLUTION PRODUCT CUDA KERNELS

///////////////////////////////////////// 2D

//Self demag kernel multiplication - real and use symmetries
//N = (N.x/2 + 1, N.y, 1)
__global__ void cu_KernelMultiplication_2D_Self(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz, 
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		int j = (idx / (N.x / 2 + 1)) % N.y;

		if (j <= N.y / 2) {

			cuS2x[idx] = (kernel.Kdiag_real[idx].x  * FMx) + (kernel.K2D_odiag[idx] * FMy);
			cuS2y[idx] = (kernel.K2D_odiag[idx] * FMx) + (kernel.Kdiag_real[idx].y  * FMy);
			cuS2z[idx] = (kernel.Kdiag_real[idx].z  * FMz);
		}
		else {

			int i = idx % (N.x / 2 + 1);

			int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

			cuS2x[idx] = (kernel.Kdiag_real[ker_idx].x  * FMx) + (-kernel.K2D_odiag[ker_idx] * FMy);
			cuS2y[idx] = (-kernel.K2D_odiag[ker_idx] * FMx) + (kernel.Kdiag_real[ker_idx].y  * FMy);
			cuS2z[idx] = (kernel.Kdiag_real[ker_idx].z  * FMz);
		}
	}
}

//Self demag kernel multiplication - real and use symmetries
__global__ void cu_KernelMultiplication_2D_Self_transpose_xy(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1);

			cuS2x[idx] = (kernel.Kdiag_real[ker_idx].x  * FMx) + (kernel.K2D_odiag[ker_idx] * FMy);
			cuS2y[idx] = (kernel.K2D_odiag[ker_idx] * FMx) + (kernel.Kdiag_real[ker_idx].y  * FMy);
			cuS2z[idx] = (kernel.Kdiag_real[ker_idx].z  * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

			cuS2x[idx] = (kernel.Kdiag_real[ker_idx].x  * FMx) + (-kernel.K2D_odiag[ker_idx] * FMy);
			cuS2y[idx] = (-kernel.K2D_odiag[ker_idx] * FMx) + (kernel.Kdiag_real[ker_idx].y  * FMy);
			cuS2z[idx] = (kernel.Kdiag_real[ker_idx].z  * FMz);
		}
	}
}

//Self demag kernel multiplication - real and use symmetries
//N = (N.x/2 + 1, N.y, 1)
__global__ void cu_KernelMultiplication_2D_zShifted(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		cuVEC<cuReal3>& Kdiag = kernel.Kdiag_real;
		cuVEC<cuReal3>& Kodiag = kernel.Kodiag_real;

		int j = (idx / (N.x / 2 + 1)) % N.y;

		if (j <= N.y / 2) {

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[idx].x * FMx) + (Kodiag[idx].x * FMy) + !(Kodiag[idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[idx].x * FMx) + (Kdiag[idx].y * FMy) + !(Kodiag[idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] + !(Kodiag[idx].y * FMx) + !(Kodiag[idx].z * FMy) + (Kdiag[idx].z * FMz);
		}
		else {

			int i = idx % (N.x / 2 + 1);

			int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) - (Kodiag[ker_idx].x * FMy) + !(Kodiag[ker_idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] - (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) - !(Kodiag[ker_idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] + !(Kodiag[ker_idx].y * FMx) - !(Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
		}
	}
}

__global__ void cu_KernelMultiplication_2D_zShifted_transpose_xy(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		cuVEC<cuReal3>& Kdiag = kernel.Kdiag_real;
		cuVEC<cuReal3>& Kodiag = kernel.Kodiag_real;

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1);

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) + (Kodiag[ker_idx].x * FMy) + !(Kodiag[ker_idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + !(Kodiag[ker_idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] + !(Kodiag[ker_idx].y * FMx) + !(Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) - (Kodiag[ker_idx].x * FMy) + !(Kodiag[ker_idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] - (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) - !(Kodiag[ker_idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] + !(Kodiag[ker_idx].y * FMx) - !(Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
		}
	}
}

//Self demag kernel multiplication - real and use symmetries, adjusting for inverse z shift
//N = (N.x/2 + 1, N.y, 1)
__global__ void cu_KernelMultiplication_2D_inversezShifted(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		cuVEC<cuReal3>& Kdiag = kernel.Kdiag_real;
		cuVEC<cuReal3>& Kodiag = kernel.Kodiag_real;

		int j = (idx / (N.x / 2 + 1)) % N.y;

		if (j <= N.y / 2) {

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[idx].x * FMx) + (Kodiag[idx].x * FMy) - !(Kodiag[idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[idx].x * FMx) + (Kdiag[idx].y * FMy) - !(Kodiag[idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] - !(Kodiag[idx].y * FMx) - !(Kodiag[idx].z * FMy) + (Kdiag[idx].z * FMz);
		}
		else {

			int i = idx % (N.x / 2 + 1);

			int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) - (Kodiag[ker_idx].x * FMy) - !(Kodiag[ker_idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] - (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + !(Kodiag[ker_idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] - !(Kodiag[ker_idx].y * FMx) + !(Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
		}
	}
}

__global__ void cu_KernelMultiplication_2D_inversezShifted_transpose_xy(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		cuVEC<cuReal3>& Kdiag = kernel.Kdiag_real;
		cuVEC<cuReal3>& Kodiag = kernel.Kodiag_real;

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1);

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) + (Kodiag[ker_idx].x * FMy) - !(Kodiag[ker_idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) - !(Kodiag[ker_idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] - !(Kodiag[ker_idx].y * FMx) - !(Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

			cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) - (Kodiag[ker_idx].x * FMy) - !(Kodiag[ker_idx].y * FMz);
			cuS2y[idx] = (cuReIm)cuS2y[idx] - (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + !(Kodiag[ker_idx].z * FMz);
			cuS2z[idx] = (cuReIm)cuS2z[idx] - !(Kodiag[ker_idx].y * FMx) + !(Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
		}
	}
}

//Complex kernel multiplication with no symmetries used
//N = (N.x/2 + 1, N.y, 1)
__global__ void cu_KernelMultiplication_2D_Regular(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y) {

		cuVEC<cuReIm3>& Kdiag = kernel.Kdiag_cmpl;
		cuVEC<cuReIm3>& Kodiag = kernel.Kodiag_cmpl;

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[idx].x * FMx) + (Kodiag[idx].x * FMy) + (Kodiag[idx].y * FMz);
		cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[idx].x * FMx) + (Kdiag[idx].y * FMy) + (Kodiag[idx].z * FMz);
		cuS2z[idx] = (cuReIm)cuS2z[idx] + (Kodiag[idx].y * FMx) + (Kodiag[idx].z * FMy) + (Kdiag[idx].z * FMz);
	}
}

///////////////////////////////////////// 3D

//N = (N.x/2 + 1, N.y, N.z)
__global__ void cu_KernelMultiplication_3D_Self_transpose_xy(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz, 
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.z/2 and N.y/2 points
	//Kxy is even about N.z/2 and odd about N.y/2
	//Kxz is odd about N.z/2 and even about N.y/2
	//Kyz is odd about N.z/2 and odd about N.y/2

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y * N.z) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);
		int k = idx / ((N.x / 2 + 1) * N.y);

		cuVEC<cuReal3>& Kdiag = kernel.Kdiag_real;
		cuVEC<cuReal3>& Kodiag = kernel.Kodiag_real;

		if (k <= N.z / 2) {

			if (i <= N.y / 2) {

				int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuS2x[idx] = (Kdiag[ker_idx].x * FMx) + (Kodiag[ker_idx].x * FMy) + (Kodiag[ker_idx].y * FMz);
				cuS2y[idx] = (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + (Kodiag[ker_idx].z * FMz);
				cuS2z[idx] = (Kodiag[ker_idx].y * FMx) + (Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
			}
			else {

				int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuS2x[idx] = (Kdiag[ker_idx].x * FMx) + (-Kodiag[ker_idx].x * FMy) + (Kodiag[ker_idx].y * FMz);
				cuS2y[idx] = (-Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + (-Kodiag[ker_idx].z * FMz);
				cuS2z[idx] = (Kodiag[ker_idx].y * FMx) + (-Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
			}
		}
		else {

			if (i <= N.y / 2) {

				int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuS2x[idx] = (Kdiag[ker_idx].x * FMx) + (Kodiag[ker_idx].x * FMy) + (-Kodiag[ker_idx].y * FMz);
				cuS2y[idx] = (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + (-Kodiag[ker_idx].z * FMz);
				cuS2z[idx] = (-Kodiag[ker_idx].y * FMx) + (-Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
			}
			else {

				int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuS2x[idx] = (Kdiag[ker_idx].x * FMx) + (-Kodiag[ker_idx].x * FMy) + (-Kodiag[ker_idx].y * FMz);
				cuS2y[idx] = (-Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + (Kodiag[ker_idx].z * FMz);
				cuS2z[idx] = (-Kodiag[ker_idx].y * FMx) + (Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
			}
		}
	}
}

//z shifted for 3D : complex kernels, but use kernel symmetries
//N = (N.x/2 + 1, N.y, N.z)
__global__ void cu_KernelMultiplication_3D_zShifted_transpose_xy(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	//z shifted for 3D : can use kernels of reduced dimensions but must be complex
	//
	//Kxx : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
	//Kyy : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
	//Kzz : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
	//
	//Kxy : y - asymmetrical (-), z - Re part symmetrical  (+), Im part asymmetrical (-)
	//Kxz : y - symmetrical  (+), z - Re part asymmetrical (-), Im part symmetrical  (+)
	//Kyz : y - asymmetrical (-), z - Re part asymmetrical (-), Im part symmetrical  (+)

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y * N.z) {

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);
		int k = idx / ((N.x / 2 + 1) * N.y);

		cuVEC<cuReIm3>& Kdiag = kernel.Kdiag_cmpl;
		cuVEC<cuReIm3>& Kodiag = kernel.Kodiag_cmpl;

		if (k <= N.z / 2) {

			if (i <= N.y / 2) {

				int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);
				
				//lower z, lower y
				cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) + (Kodiag[ker_idx].x * FMy) + (Kodiag[ker_idx].y * FMz);
				cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) + (Kodiag[ker_idx].z * FMz);
				cuS2z[idx] = (cuReIm)cuS2z[idx] + (Kodiag[ker_idx].y * FMx) + (Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
			}
			else {

				int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				//lower z, upper y
				cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[ker_idx].x * FMx) - (Kodiag[ker_idx].x * FMy) + (Kodiag[ker_idx].y * FMz);
				cuS2y[idx] = (cuReIm)cuS2y[idx] - (Kodiag[ker_idx].x * FMx) + (Kdiag[ker_idx].y * FMy) - (Kodiag[ker_idx].z * FMz);
				cuS2z[idx] = (cuReIm)cuS2z[idx] + (Kodiag[ker_idx].y * FMx) - (Kodiag[ker_idx].z * FMy) + (Kdiag[ker_idx].z * FMz);
			}
		}
		else {

			if (i <= N.y / 2) {

				int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				//upper z, lower y
				cuS2x[idx] = (cuReIm)cuS2x[idx] + ((~Kdiag[ker_idx].x) * FMx) + ((~Kodiag[ker_idx].x) * FMy) - ((~Kodiag[ker_idx].y) * FMz);
				cuS2y[idx] = (cuReIm)cuS2y[idx] + ((~Kodiag[ker_idx].x) * FMx) + ((~Kdiag[ker_idx].y) * FMy) - ((~Kodiag[ker_idx].z) * FMz);
				cuS2z[idx] = (cuReIm)cuS2z[idx] - ((~Kodiag[ker_idx].y) * FMx) - ((~Kodiag[ker_idx].z) * FMy) + ((~Kdiag[ker_idx].z) * FMz);
			}
			else {

				int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				//upper z, upper y
				cuS2x[idx] = (cuReIm)cuS2x[idx] + ((~Kdiag[ker_idx].x) * FMx) - ((~Kodiag[ker_idx].x) * FMy) - ((~Kodiag[ker_idx].y) * FMz);
				cuS2y[idx] = (cuReIm)cuS2y[idx] - ((~Kodiag[ker_idx].x) * FMx) + ((~Kdiag[ker_idx].y) * FMy) + ((~Kodiag[ker_idx].z) * FMz);
				cuS2z[idx] = (cuReIm)cuS2z[idx] - ((~Kodiag[ker_idx].y) * FMx) + ((~Kodiag[ker_idx].z) * FMy) + ((~Kdiag[ker_idx].z) * FMz);
			}
		}
	}
}

//Complex kernel multiplication with no symmetries used
//N = (N.x/2 + 1, N.y, N.z)
__global__ void cu_KernelMultiplication_3D_Regular(
	cuKerType& kernel,
	cuComplex* cuSx, cuComplex* cuSy, cuComplex* cuSz,
	cuComplex* cuS2x, cuComplex* cuS2y, cuComplex* cuS2z,
	cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y * N.z) {

		cuVEC<cuReIm3>& Kdiag = kernel.Kdiag_cmpl;
		cuVEC<cuReIm3>& Kodiag = kernel.Kodiag_cmpl;

		cuReIm FMx = cuSx[idx];
		cuReIm FMy = cuSy[idx];
		cuReIm FMz = cuSz[idx];

		cuS2x[idx] = (cuReIm)cuS2x[idx] + (Kdiag[idx].x * FMx) + (Kodiag[idx].x * FMy) + (Kodiag[idx].y * FMz);
		cuS2y[idx] = (cuReIm)cuS2y[idx] + (Kodiag[idx].x * FMx) + (Kdiag[idx].y * FMy) + (Kodiag[idx].z * FMz);
		cuS2z[idx] = (cuReIm)cuS2z[idx] + (Kodiag[idx].y * FMx) + (Kodiag[idx].z * FMy) + (Kdiag[idx].z * FMz);
	}
}

//-------------------------- RUN-TIME KERNEL MULTIPLICATION - MULTIPLE INPUTS TO SINGLE OUTPUT (TESTING ONLY)

void DemagKernelCollectionCUDA::KernelMultiplication_2D(
	std::vector<cu_arr<cuComplex>*>& Incol_x, std::vector<cu_arr<cuComplex>*>& Incol_y, std::vector<cu_arr<cuComplex>*>& Incol_z,
	cu_arr<cuComplex>& Out_x, cu_arr<cuComplex>& Out_y, cu_arr<cuComplex>& Out_z)
{
	//first compute the self contribution -> this sets Out
	if (transpose_xy) {

		cu_KernelMultiplication_2D_Self_transpose_xy <<< ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
			*kernels[self_contribution_index],
			*Incol_x[self_contribution_index], *Incol_y[self_contribution_index], *Incol_z[self_contribution_index],
			Out_x, Out_y, Out_z, cuN);
	}
	else {

		cu_KernelMultiplication_2D_Self << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
			*kernels[self_contribution_index],
			*Incol_x[self_contribution_index], *Incol_y[self_contribution_index], *Incol_z[self_contribution_index],
			Out_x, Out_y, Out_z, cuN);
	}

	//the rest add to Out
	for (int mesh_index = 0; mesh_index < Incol_x.size(); mesh_index++) {

		if (mesh_index == self_contribution_index) continue;

		//z-shifted : use symmetries
		else if (zshifted[mesh_index]) {

			//inverse : adjust signs
			if (inverse_shifted[mesh_index]) {

				if (transpose_xy) {

					cu_KernelMultiplication_2D_inversezShifted_transpose_xy << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
						*kernels[mesh_index],
						*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
						Out_x, Out_y, Out_z, cuN);
				}
				else {
					
					cu_KernelMultiplication_2D_inversezShifted << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
						*kernels[mesh_index],
						*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
						Out_x, Out_y, Out_z, cuN);
				}
			}
			//z-shifted regular
			else {

				if (transpose_xy) {

					cu_KernelMultiplication_2D_zShifted_transpose_xy << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
						*kernels[mesh_index],
						*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
						Out_x, Out_y, Out_z, cuN);
				}
				else {

					cu_KernelMultiplication_2D_zShifted << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
						*kernels[mesh_index],
						*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
						Out_x, Out_y, Out_z, cuN);
				}
			}
		}

		//now compute the other contributions by adding to Out : general kernel multiplication without any symmetries used
		else {
			
			cu_KernelMultiplication_2D_Regular <<< ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernels[mesh_index],
				*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
				Out_x, Out_y, Out_z, cuN);
		}
	}
}

void DemagKernelCollectionCUDA::KernelMultiplication_3D(
	std::vector<cu_arr<cuComplex>*>& Incol_x, std::vector<cu_arr<cuComplex>*>& Incol_y, std::vector<cu_arr<cuComplex>*>& Incol_z,
	cu_arr<cuComplex>& Out_x, cu_arr<cuComplex>& Out_y, cu_arr<cuComplex>& Out_z)
{
	//transpose_xy always true in 3D
	
	//first compute the self contribution -> this sets Out
	cu_KernelMultiplication_3D_Self_transpose_xy <<< ((N.x / 2 + 1)*N.y*N.z + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
		*kernels[self_contribution_index],
		*Incol_x[self_contribution_index], *Incol_y[self_contribution_index], *Incol_z[self_contribution_index],
		Out_x, Out_y, Out_z, cuN);
	
	//now compute the other contribution by adding to Out
	for (int mesh_index = 0; mesh_index < Incol_x.size(); mesh_index++) {

		if (mesh_index == self_contribution_index) continue;

		//z-shifted : use symmetries
		else if (zshifted[mesh_index]) {
			
			cu_KernelMultiplication_3D_zShifted_transpose_xy <<< ((N.x / 2 + 1)*N.y*N.z + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernels[mesh_index],
				*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
				Out_x, Out_y, Out_z, cuN);
		}

		//now compute the other contributions by adding to Out : general kernel multiplication without any symmetries used
		else {
			
			cu_KernelMultiplication_3D_Regular <<< ((N.x / 2 + 1)*N.y*N.z + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				*kernels[mesh_index],
				*Incol_x[mesh_index], *Incol_y[mesh_index], *Incol_z[mesh_index],
				Out_x, Out_y, Out_z, cuN);
		}
	}
}

//-------------------------- RUN-TIME KERNEL MULTIPLICATION - MULTIPLE OUTPUTS FROM SINGLE INPUT

//-------------------------- Kernels

//2D

__global__ void cu_KernelMultiplication_Launcher_2D_Set(
	cuKerType* kernels,
	cuComplex* Inx, cuComplex* Iny, cuComplex* Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	bool* inverse_zshifted,
	cuSZ3& N, int& num_outputs, bool& transpose_xy)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < num_outputs && idx < (N.x / 2 + 1) * N.y) {

		//now determine the device multiplication function to call

		if (kernels[n].internal_demag) {

			if (transpose_xy) {

				kernels[n].cu_KernelMultiplication_2D_Self_transpose_xy_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
			else {

				kernels[n].cu_KernelMultiplication_2D_Self_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
		}
		else {

			if (kernels[n].zshifted) {

				if (inverse_zshifted[n]) {

					if (transpose_xy) {

						kernels[n].cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
					else {

						kernels[n].cu_KernelMultiplication_2D_inversezShifted_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
				}
				else {

					if (transpose_xy) {

						kernels[n].cu_KernelMultiplication_2D_zShifted_transpose_xy_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
					else {

						kernels[n].cu_KernelMultiplication_2D_zShifted_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
				}
			}
			else {

				kernels[n].cu_KernelMultiplication_2D_Regular_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
		}
	}
}

__global__ void cu_KernelMultiplication_Launcher_2D_Add(
	cuKerType* kernels,
	cuComplex* Inx, cuComplex* Iny, cuComplex* Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	bool* inverse_zshifted,
	cuSZ3& N, int& num_outputs, bool& transpose_xy)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y);

	if (n < num_outputs && idx < (N.x / 2 + 1) * N.y) {

		//now determine the device multiplication function to call

		if (kernels[n].internal_demag) {

			if (transpose_xy) {

				kernels[n].cu_KernelMultiplication_2D_Self_transpose_xy_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
			else {

				kernels[n].cu_KernelMultiplication_2D_Self_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
		}
		else {

			if (kernels[n].zshifted) {

				if (inverse_zshifted[n]) {

					if (transpose_xy) {

						kernels[n].cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
					else {

						kernels[n].cu_KernelMultiplication_2D_inversezShifted_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
				}
				else {

					if (transpose_xy) {

						kernels[n].cu_KernelMultiplication_2D_zShifted_transpose_xy_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
					else {

						kernels[n].cu_KernelMultiplication_2D_zShifted_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
					}
				}
			}
			else {

				kernels[n].cu_KernelMultiplication_2D_Regular_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
		}
	}
}

//-------------------------- Kernels

//3D

__global__ void cu_KernelMultiplication_Launcher_3D_Set(
	cuKerType* kernels,
	cuComplex* Inx, cuComplex* Iny, cuComplex* Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	cuSZ3& N, int& num_outputs)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y * N.z);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y * N.z);

	if (n < num_outputs && idx < (N.x / 2 + 1) * N.y * N.z) {

		//now determine the device multiplication function to call

		if (kernels[n].internal_demag) {

			kernels[n].cu_KernelMultiplication_3D_Self_transpose_xy_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
		}
		else {

			if (kernels[n].zshifted) {

				kernels[n].cu_KernelMultiplication_3D_zShifted_transpose_xy_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
			else {

				kernels[n].cu_KernelMultiplication_3D_Regular_Set(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
		}
	}
}

__global__ void cu_KernelMultiplication_Launcher_3D_Add(
	cuKerType* kernels,
	cuComplex* Inx, cuComplex* Iny, cuComplex* Inz,
	cuComplex** Outx, cuComplex** Outy, cuComplex** Outz,
	cuSZ3& N, int& num_outputs)
{
	int idx_big = blockDim.x * blockIdx.x + threadIdx.x;

	//index of output arrays
	int n = idx_big / ((N.x / 2 + 1) * N.y * N.z);

	//index in indexed output arrays (and also index in input)
	int idx = idx_big % ((N.x / 2 + 1) * N.y * N.z);

	if (n < num_outputs && idx < (N.x / 2 + 1) * N.y * N.z) {

		//now determine the device multiplication function to call

		if (kernels[n].internal_demag) {

			kernels[n].cu_KernelMultiplication_3D_Self_transpose_xy_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
		}
		else {

			if (kernels[n].zshifted) {

				kernels[n].cu_KernelMultiplication_3D_zShifted_transpose_xy_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
			else {

				kernels[n].cu_KernelMultiplication_3D_Regular_Add(idx, Inx, Iny, Inz, Outx[n], Outy[n], Outz[n]);
			}
		}
	}
}

//-------------------------- Kernel Launchers

//2D

void DemagKernelCollectionCUDA::KernelMultiplication_2D_Set(
	cu_arr<cuComplex>& In_x, cu_arr<cuComplex>& In_y, cu_arr<cuComplex>& In_z,
	cu_arr<cuComplex*>& Outcol_x, cu_arr<cuComplex*>& Outcol_y, cu_arr<cuComplex*>& Outcol_z)
{
	cu_KernelMultiplication_Launcher_2D_Set <<< ((N.x / 2 + 1) * N.y * kernels.size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
		kernels_gpu,
		In_x, In_y, In_z,
		Outcol_x, Outcol_y, Outcol_z,
		inverse_shifted_gpu,
		cuN, num_kernels, transpose_xy_gpu);
}

void DemagKernelCollectionCUDA::KernelMultiplication_2D_Add(
	cu_arr<cuComplex>& In_x, cu_arr<cuComplex>& In_y, cu_arr<cuComplex>& In_z,
	cu_arr<cuComplex*>& Outcol_x, cu_arr<cuComplex*>& Outcol_y, cu_arr<cuComplex*>& Outcol_z)
{
	cu_KernelMultiplication_Launcher_2D_Add <<< ((N.x / 2 + 1) * N.y * kernels.size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
		kernels_gpu,
		In_x, In_y, In_z,
		Outcol_x, Outcol_y, Outcol_z,
		inverse_shifted_gpu,
		cuN, num_kernels, transpose_xy_gpu);
}

//-------------------------- Kernel Launchers

//3D

void DemagKernelCollectionCUDA::KernelMultiplication_3D_Set(
	cu_arr<cuComplex>& In_x, cu_arr<cuComplex>& In_y, cu_arr<cuComplex>& In_z,
	cu_arr<cuComplex*>& Outcol_x, cu_arr<cuComplex*>& Outcol_y, cu_arr<cuComplex*>& Outcol_z)
{
	cu_KernelMultiplication_Launcher_3D_Set <<< ((N.x / 2 + 1) * N.y * N.z * kernels.size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
		kernels_gpu,
		In_x, In_y, In_z,
		Outcol_x, Outcol_y, Outcol_z,
		cuN, num_kernels);
}

void DemagKernelCollectionCUDA::KernelMultiplication_3D_Add(
	cu_arr<cuComplex>& In_x, cu_arr<cuComplex>& In_y, cu_arr<cuComplex>& In_z,
	cu_arr<cuComplex*>& Outcol_x, cu_arr<cuComplex*>& Outcol_y, cu_arr<cuComplex*>& Outcol_z)
{
	cu_KernelMultiplication_Launcher_3D_Add <<< ((N.x / 2 + 1) * N.y * N.z * kernels.size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
		kernels_gpu,
		In_x, In_y, In_z,
		Outcol_x, Outcol_y, Outcol_z,
		cuN, num_kernels);
}

#endif

#endif