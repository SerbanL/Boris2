#include "DemagKernelCollectionCUDA_KerType.h"

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include <cuda_runtime.h>

//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE

///////////////////////////////////////// 2D

__device__ void cuKerType::cu_KernelMultiplication_2D_Self_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	int j = (idx / (N.x / 2 + 1)) % N.y;

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (j <= N.y / 2) {

		Outx[idx] = (Kdiag_real[idx].x  * FMx) + (K2D_odiag[idx] * FMy);
		Outy[idx] = (K2D_odiag[idx] * FMx) + (Kdiag_real[idx].y  * FMy);
		Outz[idx] = (Kdiag_real[idx].z  * FMz);
	}
	else {

		int i = idx % (N.x / 2 + 1);

		int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x  * FMx) + (-K2D_odiag[ker_idx] * FMy);
		Outy[idx] = (-K2D_odiag[ker_idx] * FMx) + (Kdiag_real[ker_idx].y  * FMy);
		Outz[idx] = (Kdiag_real[ker_idx].z  * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_Self_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	int j = (idx / (N.x / 2 + 1)) % N.y;

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (j <= N.y / 2) {

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[idx].x  * FMx) + (K2D_odiag[idx] * FMy);
		Outy[idx] = (cuReIm)Outy[idx] + (K2D_odiag[idx] * FMx) + (Kdiag_real[idx].y  * FMy);
		Outz[idx] = (cuReIm)Outz[idx] + (Kdiag_real[idx].z  * FMz);
	}
	else {

		int i = idx % (N.x / 2 + 1);

		int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x  * FMx) + (-K2D_odiag[ker_idx] * FMy);
		Outy[idx] = (cuReIm)Outy[idx] + (-K2D_odiag[ker_idx] * FMx) + (Kdiag_real[ker_idx].y  * FMy);
		Outz[idx] = (cuReIm)Outz[idx] + (Kdiag_real[ker_idx].z  * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_Self_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (i <= N.y / 2) {

		int ker_idx = i + j * (N.y / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x  * FMx) + (K2D_odiag[ker_idx] * FMy);
		Outy[idx] = (K2D_odiag[ker_idx] * FMx) + (Kdiag_real[ker_idx].y  * FMy);
		Outz[idx] = (Kdiag_real[ker_idx].z  * FMz);
	}
	else {

		int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x  * FMx) + (-K2D_odiag[ker_idx] * FMy);
		Outy[idx] = (-K2D_odiag[ker_idx] * FMx) + (Kdiag_real[ker_idx].y  * FMy);
		Outz[idx] = (Kdiag_real[ker_idx].z  * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_Self_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (i <= N.y / 2) {

		int ker_idx = i + j * (N.y / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x  * FMx) + (K2D_odiag[ker_idx] * FMy);
		Outy[idx] = (cuReIm)Outy[idx] + (K2D_odiag[ker_idx] * FMx) + (Kdiag_real[ker_idx].y  * FMy);
		Outz[idx] = (cuReIm)Outz[idx] + (Kdiag_real[ker_idx].z  * FMz);
	}
	else {

		int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x  * FMx) + (-K2D_odiag[ker_idx] * FMy);
		Outy[idx] = (cuReIm)Outy[idx] + (-K2D_odiag[ker_idx] * FMx) + (Kdiag_real[ker_idx].y  * FMy);
		Outz[idx] = (cuReIm)Outz[idx] + (Kdiag_real[ker_idx].z  * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_zShifted_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int j = (idx / (N.x / 2 + 1)) % N.y;

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (j <= N.y / 2) {

		Outx[idx] = (Kdiag_real[idx].x * FMx) + (Kodiag_real[idx].x * FMy) + !(Kodiag_real[idx].y * FMz);
		Outy[idx] = (Kodiag_real[idx].x * FMx) + (Kdiag_real[idx].y * FMy) + !(Kodiag_real[idx].z * FMz);
		Outz[idx] = !(Kodiag_real[idx].y * FMx) + !(Kodiag_real[idx].z * FMy) + (Kdiag_real[idx].z * FMz);
	}
	else {

		int i = idx % (N.x / 2 + 1);

		int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) + !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = -1.0 * (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) - !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = !(Kodiag_real[ker_idx].y * FMx) - !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_zShifted_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int j = (idx / (N.x / 2 + 1)) % N.y;

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (j <= N.y / 2) {

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[idx].x * FMx) + (Kodiag_real[idx].x * FMy) + !(Kodiag_real[idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_real[idx].x * FMx) + (Kdiag_real[idx].y * FMy) + !(Kodiag_real[idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] + !(Kodiag_real[idx].y * FMx) + !(Kodiag_real[idx].z * FMy) + (Kdiag_real[idx].z * FMz);
	}
	else {

		int i = idx % (N.x / 2 + 1);

		int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) + !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] - (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) - !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] + !(Kodiag_real[ker_idx].y * FMx) - !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_zShifted_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (i <= N.y / 2) {

		int ker_idx = i + j * (N.y / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) + !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = !(Kodiag_real[ker_idx].y * FMx) + !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
	else {

		int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) + !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = -1.0 * (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) - !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = !(Kodiag_real[ker_idx].y * FMx) - !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_zShifted_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (i <= N.y / 2) {

		int ker_idx = i + j * (N.y / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) + !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] + !(Kodiag_real[ker_idx].y * FMx) + !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
	else {

		int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) + !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] - (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) - !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] + !(Kodiag_real[ker_idx].y * FMx) - !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_inversezShifted_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int j = (idx / (N.x / 2 + 1)) % N.y;

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (j <= N.y / 2) {

		Outx[idx] = (Kdiag_real[idx].x * FMx) + (Kodiag_real[idx].x * FMy) - !(Kodiag_real[idx].y * FMz);
		Outy[idx] = (Kodiag_real[idx].x * FMx) + (Kdiag_real[idx].y * FMy) - !(Kodiag_real[idx].z * FMz);
		Outz[idx] = -1.0 * !(Kodiag_real[idx].y * FMx) - !(Kodiag_real[idx].z * FMy) + (Kdiag_real[idx].z * FMz);
	}
	else {

		int i = idx % (N.x / 2 + 1);

		int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) - !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = -1.0 * (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = -1.0 * !(Kodiag_real[ker_idx].y * FMx) + !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_inversezShifted_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int j = (idx / (N.x / 2 + 1)) % N.y;

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (j <= N.y / 2) {

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[idx].x * FMx) + (Kodiag_real[idx].x * FMy) - !(Kodiag_real[idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_real[idx].x * FMx) + (Kdiag_real[idx].y * FMy) - !(Kodiag_real[idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] - !(Kodiag_real[idx].y * FMx) - !(Kodiag_real[idx].z * FMy) + (Kdiag_real[idx].z * FMz);
	}
	else {

		int i = idx % (N.x / 2 + 1);

		int ker_idx = i + (N.y - j) * (N.x / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) - !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] - (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] - !(Kodiag_real[ker_idx].y * FMx) + !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (i <= N.y / 2) {

		int ker_idx = i + j * (N.y / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) - !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) - !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = -1.0 * !(Kodiag_real[ker_idx].y * FMx) - !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
	else {

		int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

		Outx[idx] = (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) - !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = -1.0 * (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = -1.0 * !(Kodiag_real[ker_idx].y * FMx) + !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (i <= N.y / 2) {

		int ker_idx = i + j * (N.y / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) - !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) - !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] - !(Kodiag_real[ker_idx].y * FMx) - !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
	else {

		int ker_idx = (N.y - i) + j * (N.y / 2 + 1);

		Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) - (Kodiag_real[ker_idx].x * FMy) - !(Kodiag_real[ker_idx].y * FMz);
		Outy[idx] = (cuReIm)Outy[idx] - (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + !(Kodiag_real[ker_idx].z * FMz);
		Outz[idx] = (cuReIm)Outz[idx] - !(Kodiag_real[ker_idx].y * FMx) + !(Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
	}
}

__device__ void cuKerType::cu_KernelMultiplication_2D_Regular_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	Outx[idx] = (Kdiag_cmpl[idx].x * FMx) + (Kodiag_cmpl[idx].x * FMy) + (Kodiag_cmpl[idx].y * FMz);
	Outy[idx] = (Kodiag_cmpl[idx].x * FMx) + (Kdiag_cmpl[idx].y * FMy) + (Kodiag_cmpl[idx].z * FMz);
	Outz[idx] = (Kodiag_cmpl[idx].y * FMx) + (Kodiag_cmpl[idx].z * FMy) + (Kdiag_cmpl[idx].z * FMz);
}

__device__ void cuKerType::cu_KernelMultiplication_2D_Regular_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_cmpl[idx].x * FMx) + (Kodiag_cmpl[idx].x * FMy) + (Kodiag_cmpl[idx].y * FMz);
	Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_cmpl[idx].x * FMx) + (Kdiag_cmpl[idx].y * FMy) + (Kodiag_cmpl[idx].z * FMz);
	Outz[idx] = (cuReIm)Outz[idx] + (Kodiag_cmpl[idx].y * FMx) + (Kodiag_cmpl[idx].z * FMy) + (Kdiag_cmpl[idx].z * FMz);
}

///////////////////////////////////////// 3D

__device__ void cuKerType::cu_KernelMultiplication_3D_Self_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.z/2 and N.y/2 points
	//Kxy is even about N.z/2 and odd about N.y/2
	//Kxz is odd about N.z/2 and even about N.y/2
	//Kyz is odd about N.z/2 and odd about N.y/2

	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);
	int k = idx / ((N.x / 2 + 1) * N.y);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (k <= N.z / 2) {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) + (Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (Kodiag_real[ker_idx].y * FMx) + (Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (Kdiag_real[ker_idx].x * FMx) + (-Kodiag_real[ker_idx].x * FMy) + (Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (-Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (-Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (Kodiag_real[ker_idx].y * FMx) + (-Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
	}
	else {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) + (-Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (-Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (-Kodiag_real[ker_idx].y * FMx) + (-Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (Kdiag_real[ker_idx].x * FMx) + (-Kodiag_real[ker_idx].x * FMy) + (-Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (-Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (-Kodiag_real[ker_idx].y * FMx) + (Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
	}
}

__device__ void cuKerType::cu_KernelMultiplication_3D_Self_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);
	int k = idx / ((N.x / 2 + 1) * N.y);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (k <= N.z / 2) {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) + (Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (cuReIm)Outz[idx] + (Kodiag_real[ker_idx].y * FMx) + (Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) + (-Kodiag_real[ker_idx].x * FMy) + (Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (cuReIm)Outy[idx] + (-Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (-Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (cuReIm)Outz[idx] + (Kodiag_real[ker_idx].y * FMx) + (-Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
	}
	else {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) + (Kodiag_real[ker_idx].x * FMy) + (-Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (-Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (cuReIm)Outz[idx] + (-Kodiag_real[ker_idx].y * FMx) + (-Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_real[ker_idx].x * FMx) + (-Kodiag_real[ker_idx].x * FMy) + (-Kodiag_real[ker_idx].y * FMz);
			Outy[idx] = (cuReIm)Outy[idx] + (-Kodiag_real[ker_idx].x * FMx) + (Kdiag_real[ker_idx].y * FMy) + (Kodiag_real[ker_idx].z * FMz);
			Outz[idx] = (cuReIm)Outz[idx] + (-Kodiag_real[ker_idx].y * FMx) + (Kodiag_real[ker_idx].z * FMy) + (Kdiag_real[ker_idx].z * FMz);
		}
	}
}

__device__ void cuKerType::cu_KernelMultiplication_3D_zShifted_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
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

	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);
	int k = idx / ((N.x / 2 + 1) * N.y);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (k <= N.z / 2) {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			//lower z, lower y
			Outx[idx] = (Kdiag_cmpl[ker_idx].x * FMx) + (Kodiag_cmpl[ker_idx].x * FMy) + (Kodiag_cmpl[ker_idx].y * FMz);
			Outy[idx] = (Kodiag_cmpl[ker_idx].x * FMx) + (Kdiag_cmpl[ker_idx].y * FMy) + (Kodiag_cmpl[ker_idx].z * FMz);
			Outz[idx] = (Kodiag_cmpl[ker_idx].y * FMx) + (Kodiag_cmpl[ker_idx].z * FMy) + (Kdiag_cmpl[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			//lower z, upper y
			Outx[idx] = (Kdiag_cmpl[ker_idx].x * FMx) - (Kodiag_cmpl[ker_idx].x * FMy) + (Kodiag_cmpl[ker_idx].y * FMz);
			Outy[idx] = -1.0 * (Kodiag_cmpl[ker_idx].x * FMx) + (Kdiag_cmpl[ker_idx].y * FMy) - (Kodiag_cmpl[ker_idx].z * FMz);
			Outz[idx] = (Kodiag_cmpl[ker_idx].y * FMx) - (Kodiag_cmpl[ker_idx].z * FMy) + (Kdiag_cmpl[ker_idx].z * FMz);
		}
	}
	else {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			//upper z, lower y
			Outx[idx] = ((~Kdiag_cmpl[ker_idx].x) * FMx) + ((~Kodiag_cmpl[ker_idx].x) * FMy) - ((~Kodiag_cmpl[ker_idx].y) * FMz);
			Outy[idx] = ((~Kodiag_cmpl[ker_idx].x) * FMx) + ((~Kdiag_cmpl[ker_idx].y) * FMy) - ((~Kodiag_cmpl[ker_idx].z) * FMz);
			Outz[idx] = -1.0 * ((~Kodiag_cmpl[ker_idx].y) * FMx) - ((~Kodiag_cmpl[ker_idx].z) * FMy) + ((~Kdiag_cmpl[ker_idx].z) * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			//upper z, upper y
			Outx[idx] = ((~Kdiag_cmpl[ker_idx].x) * FMx) - ((~Kodiag_cmpl[ker_idx].x) * FMy) - ((~Kodiag_cmpl[ker_idx].y) * FMz);
			Outy[idx] = -1.0 * ((~Kodiag_cmpl[ker_idx].x) * FMx) + ((~Kdiag_cmpl[ker_idx].y) * FMy) + ((~Kodiag_cmpl[ker_idx].z) * FMz);
			Outz[idx] = -1.0 * ((~Kodiag_cmpl[ker_idx].y) * FMx) + ((~Kodiag_cmpl[ker_idx].z) * FMy) + ((~Kdiag_cmpl[ker_idx].z) * FMz);
		}
	}
}

__device__ void cuKerType::cu_KernelMultiplication_3D_zShifted_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
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

	int i = idx % N.y;
	int j = (idx / N.y) % (N.x / 2 + 1);
	int k = idx / ((N.x / 2 + 1) * N.y);

	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	if (k <= N.z / 2) {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			//lower z, lower y
			Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_cmpl[ker_idx].x * FMx) + (Kodiag_cmpl[ker_idx].x * FMy) + (Kodiag_cmpl[ker_idx].y * FMz);
			Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_cmpl[ker_idx].x * FMx) + (Kdiag_cmpl[ker_idx].y * FMy) + (Kodiag_cmpl[ker_idx].z * FMz);
			Outz[idx] = (cuReIm)Outz[idx] + (Kodiag_cmpl[ker_idx].y * FMx) + (Kodiag_cmpl[ker_idx].z * FMy) + (Kdiag_cmpl[ker_idx].z * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			//lower z, upper y
			Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_cmpl[ker_idx].x * FMx) - (Kodiag_cmpl[ker_idx].x * FMy) + (Kodiag_cmpl[ker_idx].y * FMz);
			Outy[idx] = (cuReIm)Outy[idx] - (Kodiag_cmpl[ker_idx].x * FMx) + (Kdiag_cmpl[ker_idx].y * FMy) - (Kodiag_cmpl[ker_idx].z * FMz);
			Outz[idx] = (cuReIm)Outz[idx] + (Kodiag_cmpl[ker_idx].y * FMx) - (Kodiag_cmpl[ker_idx].z * FMy) + (Kdiag_cmpl[ker_idx].z * FMz);
		}
	}
	else {

		if (i <= N.y / 2) {

			int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			//upper z, lower y
			Outx[idx] = (cuReIm)Outx[idx] + ((~Kdiag_cmpl[ker_idx].x) * FMx) + ((~Kodiag_cmpl[ker_idx].x) * FMy) - ((~Kodiag_cmpl[ker_idx].y) * FMz);
			Outy[idx] = (cuReIm)Outy[idx] + ((~Kodiag_cmpl[ker_idx].x) * FMx) + ((~Kdiag_cmpl[ker_idx].y) * FMy) - ((~Kodiag_cmpl[ker_idx].z) * FMz);
			Outz[idx] = (cuReIm)Outz[idx] - ((~Kodiag_cmpl[ker_idx].y) * FMx) - ((~Kodiag_cmpl[ker_idx].z) * FMy) + ((~Kdiag_cmpl[ker_idx].z) * FMz);
		}
		else {

			int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			//upper z, upper y
			Outx[idx] = (cuReIm)Outx[idx] + ((~Kdiag_cmpl[ker_idx].x) * FMx) - ((~Kodiag_cmpl[ker_idx].x) * FMy) - ((~Kodiag_cmpl[ker_idx].y) * FMz);
			Outy[idx] = (cuReIm)Outy[idx] - ((~Kodiag_cmpl[ker_idx].x) * FMx) + ((~Kdiag_cmpl[ker_idx].y) * FMy) + ((~Kodiag_cmpl[ker_idx].z) * FMz);
			Outz[idx] = (cuReIm)Outz[idx] - ((~Kodiag_cmpl[ker_idx].y) * FMx) + ((~Kodiag_cmpl[ker_idx].z) * FMy) + ((~Kdiag_cmpl[ker_idx].z) * FMz);
		}
	}
}

__device__ void cuKerType::cu_KernelMultiplication_3D_Regular_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	Outx[idx] = (Kdiag_cmpl[idx].x * FMx) + (Kodiag_cmpl[idx].x * FMy) + (Kodiag_cmpl[idx].y * FMz);
	Outy[idx] = (Kodiag_cmpl[idx].x * FMx) + (Kdiag_cmpl[idx].y * FMy) + (Kodiag_cmpl[idx].z * FMz);
	Outz[idx] = (Kodiag_cmpl[idx].y * FMx) + (Kodiag_cmpl[idx].z * FMy) + (Kdiag_cmpl[idx].z * FMz);
}

__device__ void cuKerType::cu_KernelMultiplication_3D_Regular_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz)
{
	cuReIm FMx = Inx[idx];
	cuReIm FMy = Iny[idx];
	cuReIm FMz = Inz[idx];

	Outx[idx] = (cuReIm)Outx[idx] + (Kdiag_cmpl[idx].x * FMx) + (Kodiag_cmpl[idx].x * FMy) + (Kodiag_cmpl[idx].y * FMz);
	Outy[idx] = (cuReIm)Outy[idx] + (Kodiag_cmpl[idx].x * FMx) + (Kdiag_cmpl[idx].y * FMy) + (Kodiag_cmpl[idx].z * FMz);
	Outz[idx] = (cuReIm)Outz[idx] + (Kodiag_cmpl[idx].y * FMx) + (Kodiag_cmpl[idx].z * FMy) + (Kdiag_cmpl[idx].z * FMz);
}

#endif

#endif
