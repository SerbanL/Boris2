#include "OerstedKernelCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_OERSTED

//-------------------------- CONVOLUTION PRODUCT CUDA KERNELS

__global__ void cu_Oersted_ConvProd_3D_transpose_xy(cuVEC<cuReal3>& KOe, cuBComplex* cuSx, cuBComplex* cuSy, cuBComplex* cuSz, cuSZ3& N)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//Kxy is odd about N.z/2 and even about N.y/2
	//Kxz is even about N.z/2 and odd about N.y/2
	//Kyz is even about N.z/2 and even about N.y/2

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (N.x / 2 + 1) * N.y * N.z) {

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);
		int k = idx / ((N.x / 2 + 1) * N.y);

		if (k <= N.z / 2) {

			if (i <= N.y / 2) {

				cuReIm FMx = cuSx[idx];
				cuReIm FMy = cuSy[idx];
				cuReIm FMz = cuSz[idx];

				int ker_idx = i + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuSx[idx] = !((KOe[ker_idx].x * FMy) + (KOe[ker_idx].y * FMz));
				cuSy[idx] = !((-KOe[ker_idx].x * FMx) + (KOe[ker_idx].z * FMz));
				cuSz[idx] = !((-KOe[ker_idx].y * FMx) + (-KOe[ker_idx].z * FMy));
			}
			else {

				cuReIm FMx = cuSx[idx];
				cuReIm FMy = cuSy[idx];
				cuReIm FMz = cuSz[idx];

				int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuSx[idx] = !((KOe[ker_idx].x * FMy) + (-KOe[ker_idx].y * FMz));
				cuSy[idx] = !((-KOe[ker_idx].x * FMx) + (KOe[ker_idx].z * FMz));
				cuSz[idx] = !((KOe[ker_idx].y * FMx) + (-KOe[ker_idx].z * FMy));
			}
		}
		else {

			if (i <= N.y / 2) {

				cuReIm FMx = cuSx[idx];
				cuReIm FMy = cuSy[idx];
				cuReIm FMz = cuSz[idx];

				int ker_idx = i + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuSx[idx] = !((-KOe[ker_idx].x * FMy) + (KOe[ker_idx].y * FMz));
				cuSy[idx] = !((KOe[ker_idx].x * FMx) + (KOe[ker_idx].z * FMz));
				cuSz[idx] = !((-KOe[ker_idx].y * FMx) + (-KOe[ker_idx].z * FMy));
			}
			else {

				cuReIm FMx = cuSx[idx];
				cuReIm FMy = cuSy[idx];
				cuReIm FMz = cuSz[idx];

				int ker_idx = (N.y - i) + j * (N.y / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				cuSx[idx] = !((-KOe[ker_idx].x * FMy) + (-KOe[ker_idx].y * FMz));
				cuSy[idx] = !((KOe[ker_idx].x * FMx) + (KOe[ker_idx].z * FMz));
				cuSz[idx] = !((KOe[ker_idx].y * FMx) + (-KOe[ker_idx].z * FMy));
			}
		}
	}
}

//N = (N.x/2 + 1, N.y, 4)
//xy is transposed
__global__ void cu_Oersted_ConvProd_q2D_4_transpose_xy(cuVEC<cuReal3>& KOe, cuBComplex* cuSx, cuBComplex* cuSy, cuBComplex* cuSz, cuSZ3& N)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//Kxy is odd about N.z/2 and even about N.y/2
	//Kxz is even about N.z/2 and odd about N.y/2
	//Kyz is even about N.z/2 and even about N.y/2
	
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	//N.z = 4, and this kernel was called with (N.x/2 + 1) * N.y points: handle all z points in one go
	int planecount = (N.x / 2 + 1) * N.y;

	//kernels packed into planes of (N.y / 2 + 1) * (N.x / 2 + 1) size
	int kerplanecount = (N.x / 2 + 1) * (N.y / 2 + 1);

	if (idx < planecount) {

		//the z-axis points (the others are zero)
		cuReIm3 a = cuReIm3(cuSx[idx], cuSy[idx], cuSz[idx]);
		cuReIm3 b = cuReIm3(cuSx[idx + planecount], cuSy[idx + planecount], cuSz[idx + planecount]);

		//forward z-axis fft
		//NOTE: cuda fft uses -i for the forward fft and +i for the inverse fft.
		//The kernels are purely imaginary so must use same convention here: you would get a wrong result if you used the convention of +i forward and -i inverse.
		cuReIm3 X0 = a + b;
		cuReIm3 X1 = a - !b;
		cuReIm3 X2 = a - b;
		cuReIm3 X3 = a + !b;

		//kernel multiplication
		cuReIm3 F0, F1, F2, F3;

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		if (i <= N.y / 2) {

			int ker_baseidx = i + j * (N.y / 2 + 1);

			F0.x = !((KOe[ker_baseidx].x * X0.y) + (KOe[ker_baseidx].y * X0.z));
			F0.y = !((-KOe[ker_baseidx].x * X0.x) + (KOe[ker_baseidx].z * X0.z));
			F0.z = !((-KOe[ker_baseidx].y * X0.x) + (-KOe[ker_baseidx].z * X0.y));

			F1.x = !((KOe[ker_baseidx + kerplanecount].x * X1.y) + (KOe[ker_baseidx + kerplanecount].y * X1.z));
			F1.y = !((-KOe[ker_baseidx + kerplanecount].x * X1.x) + (KOe[ker_baseidx + kerplanecount].z * X1.z));
			F1.z = !((-KOe[ker_baseidx + kerplanecount].y * X1.x) + (-KOe[ker_baseidx + kerplanecount].z * X1.y));

			F2.x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X2.y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X2.z));
			F2.y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X2.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X2.z));
			F2.z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X2.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X2.y));

			F3.x = !((-KOe[ker_baseidx + kerplanecount].x * X3.y) + (KOe[ker_baseidx + kerplanecount].y * X3.z));
			F3.y = !((KOe[ker_baseidx + kerplanecount].x * X3.x) + (KOe[ker_baseidx + kerplanecount].z * X3.z));
			F3.z = !((-KOe[ker_baseidx + kerplanecount].y * X3.x) + (-KOe[ker_baseidx + kerplanecount].z * X3.y));
		}
		else {

			int ker_baseidx = (N.y - i) + j * (N.y / 2 + 1);

			F0.x = !((KOe[ker_baseidx].x * X0.y) + (-KOe[ker_baseidx].y * X0.z));
			F0.y = !((-KOe[ker_baseidx].x * X0.x) + (KOe[ker_baseidx].z * X0.z));
			F0.z = !((KOe[ker_baseidx].y * X0.x) + (-KOe[ker_baseidx].z * X0.y));

			F1.x = !((KOe[ker_baseidx + kerplanecount].x * X1.y) + (-KOe[ker_baseidx + kerplanecount].y * X1.z));
			F1.y = !((-KOe[ker_baseidx + kerplanecount].x * X1.x) + (KOe[ker_baseidx + kerplanecount].z * X1.z));
			F1.z = !((KOe[ker_baseidx + kerplanecount].y * X1.x) + (-KOe[ker_baseidx + kerplanecount].z * X1.y));

			F2.x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X2.y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X2.z));
			F2.y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X2.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X2.z));
			F2.z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X2.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X2.y));

			F3.x = !((-KOe[ker_baseidx + kerplanecount].x * X3.y) + (-KOe[ker_baseidx + kerplanecount].y * X3.z));
			F3.y = !((KOe[ker_baseidx + kerplanecount].x * X3.x) + (KOe[ker_baseidx + kerplanecount].z * X3.z));
			F3.z = !((KOe[ker_baseidx + kerplanecount].y * X3.x) + (-KOe[ker_baseidx + kerplanecount].z * X3.y));
		}

		//inverse z-axis fft (but without division by 4). Also only keep first 2 points

		cuSx[idx] = F0.x + F1.x + F2.x + F3.x;
		cuSy[idx] = F0.y + F1.y + F2.y + F3.y;
		cuSz[idx] = F0.z + F1.z + F2.z + F3.z;

		cuReIm3 F1c = !F1;
		cuReIm3 F3c = !F3;

		cuSx[idx + planecount] = F0.x + F1c.x + F2.x - F3c.x;
		cuSy[idx + planecount] = F0.y + F1c.y + F2.y - F3c.y;
		cuSz[idx + planecount] = F0.z + F1c.z + F2.z - F3c.z;
	}
}

//N = (N.x/2 + 1, N.y, 8)
//xy is transposed
__global__ void cu_Oersted_ConvProd_q2D_8_transpose_xy(cuVEC<cuReal3>& KOe, cuBComplex* cuSx, cuBComplex* cuSy, cuBComplex* cuSz, cuSZ3& N)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//Kxy is odd about N.z/2 and even about N.y/2
	//Kxz is even about N.z/2 and odd about N.y/2
	//Kyz is even about N.z/2 and even about N.y/2

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	//N.z = 8, and this kernel was called with (N.x/2 + 1) * N.y points: handle all z points in one go
	int planecount = (N.x / 2 + 1) * N.y;

	//kernels packed into planes of (N.y / 2 + 1) * (N.x / 2 + 1) size
	int kerplanecount = (N.x / 2 + 1) * (N.y / 2 + 1);

	if (idx < planecount) {

#define a (cuBReal)0.7071067811865

		//the z-axis points (the others are zero)
		cuReIm3 x0 = cuReIm3(cuSx[idx], cuSy[idx], cuSz[idx]);
		cuReIm3 x1 = cuReIm3(cuSx[idx + planecount], cuSy[idx + planecount], cuSz[idx + planecount]);
		cuReIm3 x2 = cuReIm3(cuSx[idx + 2 * planecount], cuSy[idx + 2 * planecount], cuSz[idx + 2 * planecount]);
		cuReIm3 x3 = cuReIm3(cuSx[idx + 3 * planecount], cuSy[idx + 3 * planecount], cuSz[idx + 3 * planecount]);

		//Radix-4 step
		cuReIm3 X0 = x0 + x2;
		cuReIm3 X2 = x0 - x2;
		cuReIm3 X4 = x0 - !x2;
		cuReIm3 X6 = x0 + !x2;

		cuReIm3 X1 = x1 + x3;
		cuReIm3 X3 = !(x3 - x1);
		cuReIm3 X5 = (x1 - !x3) * cuReIm(a, -a);
		cuReIm3 X7 = (x1 + !x3) * cuReIm(-a, -a);

		//Radix-2 step
		cuReIm3 temp = X0 - X1;
		X0 = X0 + X1;
		X1 = temp;

		temp = X2 - X3;
		X2 = X2 + X3;
		X3 = temp;

		temp = X4 - X5;
		X4 = X4 + X5;
		X5 = temp;

		temp = X6 - X7;
		X6 = X6 + X7;
		X7 = temp;

		//data set in shuffled order:
		//X0, X4, X2, X6, X1, X5, X3, X7

		cuReIm3 F0, F1, F2, F3, F4, F5, F6, F7;

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		if (i <= N.y / 2) {

			int ker_baseidx = i + j * (N.y / 2 + 1);

			F0.x = !((KOe[ker_baseidx].x * X0.y) + (KOe[ker_baseidx].y * X0.z));
			F0.y = !((-KOe[ker_baseidx].x * X0.x) + (KOe[ker_baseidx].z * X0.z));
			F0.z = !((-KOe[ker_baseidx].y * X0.x) + (-KOe[ker_baseidx].z * X0.y));

			F4.x = !((KOe[ker_baseidx + kerplanecount].x * X4.y) + (KOe[ker_baseidx + kerplanecount].y * X4.z));
			F4.y = !((-KOe[ker_baseidx + kerplanecount].x * X4.x) + (KOe[ker_baseidx + kerplanecount].z * X4.z));
			F4.z = !((-KOe[ker_baseidx + kerplanecount].y * X4.x) + (-KOe[ker_baseidx + kerplanecount].z * X4.y));

			F2.x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X2.y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X2.z));
			F2.y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X2.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X2.z));
			F2.z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X2.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X2.y));

			F6.x = !((KOe[ker_baseidx + 3 * kerplanecount].x * X6.y) + (KOe[ker_baseidx + 3 * kerplanecount].y * X6.z));
			F6.y = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X6.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X6.z));
			F6.z = !((-KOe[ker_baseidx + 3 * kerplanecount].y * X6.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X6.y));

			F1.x = !((KOe[ker_baseidx + 4 * kerplanecount].x * X1.y) + (KOe[ker_baseidx + 4 * kerplanecount].y * X1.z));
			F1.y = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X1.x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X1.z));
			F1.z = !((-KOe[ker_baseidx + 4 * kerplanecount].y * X1.x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X1.y));

			F5.x = !((-KOe[ker_baseidx + 3* kerplanecount].x * X5.y) + (KOe[ker_baseidx + 3 * kerplanecount].y * X5.z));
			F5.y = !((KOe[ker_baseidx + 3 * kerplanecount].x * X5.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X5.z));
			F5.z = !((-KOe[ker_baseidx + 3 * kerplanecount].y * X5.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X5.y));

			F3.x = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X3.y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X3.z));
			F3.y = !((KOe[ker_baseidx + 2 * kerplanecount].x * X3.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X3.z));
			F3.z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X3.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X3.y));

			F7.x = !((-KOe[ker_baseidx + kerplanecount].x * X7.y) + (KOe[ker_baseidx + kerplanecount].y * X7.z));
			F7.y = !((KOe[ker_baseidx + kerplanecount].x * X7.x) + (KOe[ker_baseidx + kerplanecount].z * X7.z));
			F7.z = !((-KOe[ker_baseidx + kerplanecount].y * X7.x) + (-KOe[ker_baseidx + kerplanecount].z * X7.y));
		}
		else {

			int ker_baseidx = (N.y - i) + j * (N.y / 2 + 1);

			F0.x = !((KOe[ker_baseidx].x * X0.y) + (-KOe[ker_baseidx].y * X0.z));
			F0.y = !((-KOe[ker_baseidx].x * X0.x) + (KOe[ker_baseidx].z * X0.z));
			F0.z = !((KOe[ker_baseidx].y * X0.x) + (-KOe[ker_baseidx].z * X0.y));

			F4.x = !((KOe[ker_baseidx + kerplanecount].x * X4.y) + (-KOe[ker_baseidx + kerplanecount].y * X4.z));
			F4.y = !((-KOe[ker_baseidx + kerplanecount].x * X4.x) + (KOe[ker_baseidx + kerplanecount].z * X4.z));
			F4.z = !((KOe[ker_baseidx + kerplanecount].y * X4.x) + (-KOe[ker_baseidx + kerplanecount].z * X4.y));

			F2.x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X2.y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X2.z));
			F2.y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X2.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X2.z));
			F2.z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X2.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X2.y));

			F6.x = !((KOe[ker_baseidx + 3 * kerplanecount].x * X6.y) + (-KOe[ker_baseidx + 3 * kerplanecount].y * X6.z));
			F6.y = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X6.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X6.z));
			F6.z = !((KOe[ker_baseidx + 3 * kerplanecount].y * X6.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X6.y));

			F1.x = !((KOe[ker_baseidx + 4 * kerplanecount].x * X1.y) + (-KOe[ker_baseidx + 4 * kerplanecount].y * X1.z));
			F1.y = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X1.x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X1.z));
			F1.z = !((KOe[ker_baseidx + 4 * kerplanecount].y * X1.x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X1.y));

			F5.x = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X5.y) + (-KOe[ker_baseidx + 3 * kerplanecount].y * X5.z));
			F5.y = !((KOe[ker_baseidx + 3 * kerplanecount].x * X5.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X5.z));
			F5.z = !((KOe[ker_baseidx + 3 * kerplanecount].y * X5.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X5.y));

			F3.x = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X3.y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X3.z));
			F3.y = !((KOe[ker_baseidx + 2 * kerplanecount].x * X3.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X3.z));
			F3.z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X3.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X3.y));

			F7.x = !((-KOe[ker_baseidx + kerplanecount].x * X7.y) + (-KOe[ker_baseidx + kerplanecount].y * X7.z));
			F7.y = !((KOe[ker_baseidx + kerplanecount].x * X7.x) + (KOe[ker_baseidx + kerplanecount].z * X7.z));
			F7.z = !((KOe[ker_baseidx + kerplanecount].y * X7.x) + (-KOe[ker_baseidx + kerplanecount].z * X7.y));
		}

		//inverse z-axis fft (but without division by 8). Also only keep first 4 points.

		//Radix-2 step
		X0 = F0 + F1;
		X1 = F0 - F1;

		X2 = F2 + F3;
		X3 = F2 - F3;

		X4 = F4 + F5;
		X5 = F4 - F5;

		X6 = F6 + F7;
		X7 = F6 - F7;

		//Radix-4 step
		cuReIm3 t0 = X0 + X2;
		cuReIm3 t1 = X0 - X2;
		cuReIm3 t2 = X4 + X6;
		cuReIm3 t3 = !(X6 - X4);

		X0 = (t0 + t2);
		X2 = (t1 - t3);

		t0 = X1 + !X3;
		t1 = X1 - !X3;
		t2 = X5 * cuReIm(a, a) + X7 * cuReIm(-a, a);
		t3 = X7 * cuReIm(-a, -a) - X5 * cuReIm(-a, a);

		X1 = (t0 + t2);
		X3 = (t1 - t3);

		cuSx[idx] = X0.x;
		cuSy[idx] = X0.y;
		cuSz[idx] = X0.z;

		cuSx[idx + planecount] = X1.x;
		cuSy[idx + planecount] = X1.y;
		cuSz[idx + planecount] = X1.z;

		cuSx[idx + 2 * planecount] = X2.x;
		cuSy[idx + 2 * planecount] = X2.y;
		cuSz[idx + 2 * planecount] = X2.z;

		cuSx[idx + 3 * planecount] = X3.x;
		cuSy[idx + 3 * planecount] = X3.y;
		cuSz[idx + 3 * planecount] = X3.z;

#undef a
	}
}

//N = (N.x/2 + 1, N.y, 16)
//xy is transposed
__global__ void cu_Oersted_ConvProd_q2D_16_transpose_xy(cuVEC<cuReal3>& KOe, cuBComplex* cuSx, cuBComplex* cuSy, cuBComplex* cuSz, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	//N.z = 16, and this kernel was called with (N.x/2 + 1) * N.y points: handle all z points in one go
	int planecount = (N.x / 2 + 1) * N.y;

	//kernels packed into planes of (N.y / 2 + 1) * (N.x / 2 + 1) size
	int kerplanecount = (N.x / 2 + 1) * (N.y / 2 + 1);

	if (idx < planecount) {

		//the z-axis points (the others are zero)
		cuReIm3 x0 = cuReIm3(cuSx[idx], cuSy[idx], cuSz[idx]);
		cuReIm3 x1 = cuReIm3(cuSx[idx + planecount], cuSy[idx + planecount], cuSz[idx + planecount]);
		cuReIm3 x2 = cuReIm3(cuSx[idx + 2 * planecount], cuSy[idx + 2 * planecount], cuSz[idx + 2 * planecount]);
		cuReIm3 x3 = cuReIm3(cuSx[idx + 3 * planecount], cuSy[idx + 3 * planecount], cuSz[idx + 3 * planecount]);
		cuReIm3 x4 = cuReIm3(cuSx[idx + 4 * planecount], cuSy[idx + 4 * planecount], cuSz[idx + 4 * planecount]);
		cuReIm3 x5 = cuReIm3(cuSx[idx + 5 * planecount], cuSy[idx + 5 * planecount], cuSz[idx + 5 * planecount]);
		cuReIm3 x6 = cuReIm3(cuSx[idx + 6 * planecount], cuSy[idx + 6 * planecount], cuSz[idx + 6 * planecount]);
		cuReIm3 x7 = cuReIm3(cuSx[idx + 7 * planecount], cuSy[idx + 7 * planecount], cuSz[idx + 7 * planecount]);

#define a	(cuBReal)9.238795325113E-01
#define b	(cuBReal)3.826834323651E-01
#define c	(cuBReal)7.071067811865E-01

		//First stage
		cuReIm3 X0 = x0 + x4;
		cuReIm3 X4 = x0 - x4;
		cuReIm3 X8 = x0 - !x4;
		cuReIm3 X12 = x0 + !x4;

		cuReIm3 X1 = x1 + x5;
		cuReIm3 X5 = (x1 - x5) * cuReIm(c, -c);
		cuReIm3 X9 = (x1 - !x5) * cuReIm(a, -b);
		cuReIm3 X13 = (x1 + !x5) * cuReIm(b, -a);

		cuReIm3 X2 = x2 + x6;
		cuReIm3 X6 = !(x6 - x2);
		cuReIm3 X10 = (x2 - !x6) * cuReIm(c, -c);
		cuReIm3 X14 = (x2 + !x6) * cuReIm(-c, -c);

		cuReIm3 X3 = x3 + x7;
		cuReIm3 X7 = (x3 - x7) * cuReIm(-c, -c);
		cuReIm3 X11 = (x3 - !x7) * cuReIm(b, -a);
		cuReIm3 X15 = (x3 + !x7) * cuReIm(-a, b);

		//Second stage
		cuReIm3 t0 = X0 + X2;
		cuReIm3 t1 = X0 - X2;
		cuReIm3 t2 = X1 + X3;
		cuReIm3 t3 = !(X3 - X1);

		X0 = t0 + t2;
		X1 = t0 - t2;
		X2 = t1 + t3;
		X3 = t1 - t3;

		t0 = X4 + X6;
		t1 = X4 - X6;
		t2 = X5 + X7;
		t3 = !(X7 - X5);

		X4 = t0 + t2;
		X5 = t0 - t2;
		X6 = t1 + t3;
		X7 = t1 - t3;

		t0 = X8 + X10;
		t1 = X8 - X10;
		t2 = X9 + X11;
		t3 = !(X11 - X9);

		X8 = t0 + t2;
		X9 = t0 - t2;
		X10 = t1 + t3;
		X11 = t1 - t3;

		t0 = X12 + X14;
		t1 = X12 - X14;
		t2 = X13 + X15;
		t3 = !(X15 - X13);

		X12 = t0 + t2;
		X13 = t0 - t2;
		X14 = t1 + t3;
		X15 = t1 - t3;

		//output is shuffled now, i.e. it is ordered as:
		//X0, X8, X4, X12, X2, X10, X6, X14, X1, X9, X5, X13, X3, X11, X7, X15

		cuReIm3 F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15;

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		if (i <= N.y / 2) {

			int ker_baseidx = i + j * (N.y / 2 + 1);

			F0.x = !((KOe[ker_baseidx].x * X0.y) + (KOe[ker_baseidx].y * X0.z));
			F0.y = !((-KOe[ker_baseidx].x * X0.x) + (KOe[ker_baseidx].z * X0.z));
			F0.z = !((-KOe[ker_baseidx].y * X0.x) + (-KOe[ker_baseidx].z * X0.y));

			F8.x = !((KOe[ker_baseidx + kerplanecount].x * X8.y) + (KOe[ker_baseidx + kerplanecount].y * X8.z));
			F8.y = !((-KOe[ker_baseidx + kerplanecount].x * X8.x) + (KOe[ker_baseidx + kerplanecount].z * X8.z));
			F8.z = !((-KOe[ker_baseidx + kerplanecount].y * X8.x) + (-KOe[ker_baseidx + kerplanecount].z * X8.y));

			F4.x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X4.y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X4.z));
			F4.y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X4.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X4.z));
			F4.z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X4.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X4.y));

			F12.x = !((KOe[ker_baseidx + 3 * kerplanecount].x * X12.y) + (KOe[ker_baseidx + 3 * kerplanecount].y * X12.z));
			F12.y = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X12.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X12.z));
			F12.z = !((-KOe[ker_baseidx + 3 * kerplanecount].y * X12.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X12.y));

			F2.x = !((KOe[ker_baseidx + 4 * kerplanecount].x * X2.y) + (KOe[ker_baseidx + 4 * kerplanecount].y * X2.z));
			F2.y = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X2.x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X2.z));
			F2.z = !((-KOe[ker_baseidx + 4 * kerplanecount].y * X2.x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X2.y));
			
			F10.x = !((KOe[ker_baseidx + 5 * kerplanecount].x * X10.y) + (KOe[ker_baseidx + 5 * kerplanecount].y * X10.z));
			F10.y = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X10.x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X10.z));
			F10.z = !((-KOe[ker_baseidx + 5 * kerplanecount].y * X10.x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X10.y));

			F6.x = !((KOe[ker_baseidx + 6 * kerplanecount].x * X6.y) + (KOe[ker_baseidx + 6 * kerplanecount].y * X6.z));
			F6.y = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X6.x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X6.z));
			F6.z = !((-KOe[ker_baseidx + 6 * kerplanecount].y * X6.x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X6.y));

			F14.x = !((KOe[ker_baseidx + 7 * kerplanecount].x * X14.y) + (KOe[ker_baseidx + 7 * kerplanecount].y * X14.z));
			F14.y = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X14.x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X14.z));
			F14.z = !((-KOe[ker_baseidx + 7 * kerplanecount].y * X14.x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X14.y));

			F1.x = !((KOe[ker_baseidx + 8 * kerplanecount].x * X1.y) + (KOe[ker_baseidx + 8 * kerplanecount].y * X1.z));
			F1.y = !((-KOe[ker_baseidx + 8 * kerplanecount].x * X1.x) + (KOe[ker_baseidx + 8 * kerplanecount].z * X1.z));
			F1.z = !((-KOe[ker_baseidx + 8 * kerplanecount].y * X1.x) + (-KOe[ker_baseidx + 8 * kerplanecount].z * X1.y));

			F9.x = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X9.y) + (KOe[ker_baseidx + 7 * kerplanecount].y * X9.z));
			F9.y = !((KOe[ker_baseidx + 7 * kerplanecount].x * X9.x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X9.z));
			F9.z = !((-KOe[ker_baseidx + 7 * kerplanecount].y * X9.x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X9.y));

			F5.x = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X5.y) + (KOe[ker_baseidx + 6 * kerplanecount].y * X5.z));
			F5.y = !((KOe[ker_baseidx + 6 * kerplanecount].x * X5.x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X5.z));
			F5.z = !((-KOe[ker_baseidx + 6 * kerplanecount].y * X5.x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X5.y));

			F13.x = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X13.y) + (KOe[ker_baseidx + 5 * kerplanecount].y * X13.z));
			F13.y = !((KOe[ker_baseidx + 5 * kerplanecount].x * X13.x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X13.z));
			F13.z = !((-KOe[ker_baseidx + 5 * kerplanecount].y * X13.x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X13.y));

			F3.x = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X3.y) + (KOe[ker_baseidx + 4 * kerplanecount].y * X3.z));
			F3.y = !((KOe[ker_baseidx + 4 * kerplanecount].x * X3.x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X3.z));
			F3.z = !((-KOe[ker_baseidx + 4 * kerplanecount].y * X3.x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X3.y));

			F11.x = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X11.y) + (KOe[ker_baseidx + 3 * kerplanecount].y * X11.z));
			F11.y = !((KOe[ker_baseidx + 3 * kerplanecount].x * X11.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X11.z));
			F11.z = !((-KOe[ker_baseidx + 3 * kerplanecount].y * X11.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X11.y));

			F7.x = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X7.y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X7.z));
			F7.y = !((KOe[ker_baseidx + 2 * kerplanecount].x * X7.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X7.z));
			F7.z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X7.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X7.y));

			F15.x = !((-KOe[ker_baseidx + kerplanecount].x * X15.y) + (KOe[ker_baseidx + kerplanecount].y * X15.z));
			F15.y = !((KOe[ker_baseidx + kerplanecount].x * X15.x) + (KOe[ker_baseidx + kerplanecount].z * X15.z));
			F15.z = !((-KOe[ker_baseidx + kerplanecount].y * X15.x) + (-KOe[ker_baseidx + kerplanecount].z * X15.y));
		}
		else {

			int ker_baseidx = (N.y - i) + j * (N.y / 2 + 1);

			F0.x = !((KOe[ker_baseidx].x * X0.y) + (-KOe[ker_baseidx].y * X0.z));
			F0.y = !((-KOe[ker_baseidx].x * X0.x) + (KOe[ker_baseidx].z * X0.z));
			F0.z = !((KOe[ker_baseidx].y * X0.x) + (-KOe[ker_baseidx].z * X0.y));

			F8.x = !((KOe[ker_baseidx + kerplanecount].x * X8.y) + (-KOe[ker_baseidx + kerplanecount].y * X8.z));
			F8.y = !((-KOe[ker_baseidx + kerplanecount].x * X8.x) + (KOe[ker_baseidx + kerplanecount].z * X8.z));
			F8.z = !((KOe[ker_baseidx + kerplanecount].y * X8.x) + (-KOe[ker_baseidx + kerplanecount].z * X8.y));

			F4.x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X4.y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X4.z));
			F4.y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X4.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X4.z));
			F4.z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X4.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X4.y));

			F12.x = !((KOe[ker_baseidx + 3 * kerplanecount].x * X12.y) + (-KOe[ker_baseidx + 3 * kerplanecount].y * X12.z));
			F12.y = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X12.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X12.z));
			F12.z = !((KOe[ker_baseidx + 3 * kerplanecount].y * X12.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X12.y));

			F2.x = !((KOe[ker_baseidx + 4 * kerplanecount].x * X2.y) + (-KOe[ker_baseidx + 4 * kerplanecount].y * X2.z));
			F2.y = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X2.x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X2.z));
			F2.z = !((KOe[ker_baseidx + 4 * kerplanecount].y * X2.x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X2.y));

			F10.x = !((KOe[ker_baseidx + 5 * kerplanecount].x * X10.y) + (-KOe[ker_baseidx + 5 * kerplanecount].y * X10.z));
			F10.y = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X10.x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X10.z));
			F10.z = !((KOe[ker_baseidx + 5 * kerplanecount].y * X10.x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X10.y));

			F6.x = !((KOe[ker_baseidx + 6 * kerplanecount].x * X6.y) + (-KOe[ker_baseidx + 6 * kerplanecount].y * X6.z));
			F6.y = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X6.x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X6.z));
			F6.z = !((KOe[ker_baseidx + 6 * kerplanecount].y * X6.x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X6.y));

			F14.x = !((KOe[ker_baseidx + 7 * kerplanecount].x * X14.y) + (-KOe[ker_baseidx + 7 * kerplanecount].y * X14.z));
			F14.y = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X14.x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X14.z));
			F14.z = !((KOe[ker_baseidx + 7 * kerplanecount].y * X14.x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X14.y));

			F1.x = !((KOe[ker_baseidx + 8 * kerplanecount].x * X1.y) + (-KOe[ker_baseidx + 8 * kerplanecount].y * X1.z));
			F1.y = !((-KOe[ker_baseidx + 8 * kerplanecount].x * X1.x) + (KOe[ker_baseidx + 8 * kerplanecount].z * X1.z));
			F1.z = !((KOe[ker_baseidx + 8 * kerplanecount].y * X1.x) + (-KOe[ker_baseidx + 8 * kerplanecount].z * X1.y));

			F9.x = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X9.y) + (-KOe[ker_baseidx + 7 * kerplanecount].y * X9.z));
			F9.y = !((KOe[ker_baseidx + 7 * kerplanecount].x * X9.x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X9.z));
			F9.z = !((KOe[ker_baseidx + 7 * kerplanecount].y * X9.x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X9.y));

			F5.x = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X5.y) + (-KOe[ker_baseidx + 6 * kerplanecount].y * X5.z));
			F5.y = !((KOe[ker_baseidx + 6 * kerplanecount].x * X5.x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X5.z));
			F5.z = !((KOe[ker_baseidx + 6 * kerplanecount].y * X5.x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X5.y));

			F13.x = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X13.y) + (-KOe[ker_baseidx + 5 * kerplanecount].y * X13.z));
			F13.y = !((KOe[ker_baseidx + 5 * kerplanecount].x * X13.x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X13.z));
			F13.z = !((KOe[ker_baseidx + 5 * kerplanecount].y * X13.x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X13.y));

			F3.x = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X3.y) + (-KOe[ker_baseidx + 4 * kerplanecount].y * X3.z));
			F3.y = !((KOe[ker_baseidx + 4 * kerplanecount].x * X3.x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X3.z));
			F3.z = !((KOe[ker_baseidx + 4 * kerplanecount].y * X3.x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X3.y));

			F11.x = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X11.y) + (-KOe[ker_baseidx + 3 * kerplanecount].y * X11.z));
			F11.y = !((KOe[ker_baseidx + 3 * kerplanecount].x * X11.x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X11.z));
			F11.z = !((KOe[ker_baseidx + 3 * kerplanecount].y * X11.x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X11.y));

			F7.x = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X7.y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X7.z));
			F7.y = !((KOe[ker_baseidx + 2 * kerplanecount].x * X7.x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X7.z));
			F7.z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X7.x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X7.y));

			F15.x = !((-KOe[ker_baseidx + kerplanecount].x * X15.y) + (-KOe[ker_baseidx + kerplanecount].y * X15.z));
			F15.y = !((KOe[ker_baseidx + kerplanecount].x * X15.x) + (KOe[ker_baseidx + kerplanecount].z * X15.z));
			F15.z = !((KOe[ker_baseidx + kerplanecount].y * X15.x) + (-KOe[ker_baseidx + kerplanecount].z * X15.y));
		}

		//inverse z-axis fft (but without division by 16). Also only keep first 8 points.

		//First stage
		t0 = F0 + F1;
		t1 = F0 - F1;
		t2 = F2 + F3;
		t3 = !(F3 - F2);

		X0 = t0 + t2;
		X1 = t1 - t3;
		X2 = t0 - t2;
		X3 = t1 + t3;

		t0 = F4 + F5;
		t1 = F4 - F5;
		t2 = F6 + F7;
		t3 = !(F7 - F6);

		X4 = t0 + t2;
		X5 = t1 - t3;
		X6 = t0 - t2;
		X7 = t1 + t3;

		t0 = F8 + F9;
		t1 = F8 - F9;
		t2 = F10 + F11;
		t3 = !(F11 - F10);

		X8 = t0 + t2;
		X9 = t1 - t3;
		X10 = t0 - t2;
		X11 = t1 + t3;

		t0 = F12 + F13;
		t1 = F12 - F13;
		t2 = F14 + F15;
		t3 = !(F15 - F14);

		X12 = t0 + t2;
		X13 = t1 - t3;
		X14 = t0 - t2;
		X15 = t1 + t3;

		//Second stage

		t0 = X0 + X4;
		t1 = X0 - X4;
		t2 = X8 + X12;
		t3 = !(X12 - X8);

		X0 = t0 + t2;
		X4 = t1 - t3;

		t0 = X1 + X5 * cuReIm(c, c);
		t1 = X1 - X5 * cuReIm(c, c);
		t2 = X9 * cuReIm(a, b) + X13 * cuReIm(b, a);
		t3 = (X13 * cuReIm(-a, b) - X9 * cuReIm(-b, a));

		X1 = t0 + t2;
		X5 = t1 - t3;

		t0 = X2 + !X6;
		t1 = X2 - !X6;
		t2 = X10 * cuReIm(c, c) + X14 * cuReIm(-c, c);
		t3 = (X14 * cuReIm(-c, -c) - X10 * cuReIm(-c, c));

		X2 = t0 + t2;
		X6 = t1 - t3;

		t0 = X3 + X7 * cuReIm(-c, c);
		t1 = X3 - X7 * cuReIm(-c, c);
		t2 = X11 * cuReIm(b, a) + X15 * cuReIm(-a, -b);
		t3 = (X15 * cuReIm(b, -a) - X11 * cuReIm(-a, b));

		X3 = t0 + t2;
		X7 = t1 - t3;

		cuSx[idx] = X0.x;
		cuSy[idx] = X0.y;
		cuSz[idx] = X0.z;
		cuSx[idx + 4 * planecount] = X4.x;
		cuSy[idx + 4 * planecount] = X4.y;
		cuSz[idx + 4 * planecount] = X4.z;

		cuSx[idx + 1 * planecount] = X1.x;
		cuSy[idx + 1 * planecount] = X1.y;
		cuSz[idx + 1 * planecount] = X1.z;
		cuSx[idx + 5 * planecount] = X5.x;
		cuSy[idx + 5 * planecount] = X5.y;
		cuSz[idx + 5 * planecount] = X5.z;

		cuSx[idx + 2 * planecount] = X2.x;
		cuSy[idx + 2 * planecount] = X2.y;
		cuSz[idx + 2 * planecount] = X2.z;
		cuSx[idx + 6 * planecount] = X6.x;
		cuSy[idx + 6 * planecount] = X6.y;
		cuSz[idx + 6 * planecount] = X6.z;

		cuSx[idx + 3 * planecount] = X3.x;
		cuSy[idx + 3 * planecount] = X3.y;
		cuSz[idx + 3 * planecount] = X3.z;
		cuSx[idx + 7 * planecount] = X7.x;
		cuSy[idx + 7 * planecount] = X7.y;
		cuSz[idx + 7 * planecount] = X7.z;

#undef a
#undef b
#undef c
	}
}

//N = (N.x/2 + 1, N.y, 32)
//xy is transposed
__global__ void cu_Oersted_ConvProd_q2D_32_transpose_xy(cuVEC<cuReal3>& KOe, cuBComplex* cuSx, cuBComplex* cuSy, cuBComplex* cuSz, cuSZ3& N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	//N.z = 32, and this kernel was called with (N.x/2 + 1) * N.y points: handle all z points in one go
	int planecount = (N.x / 2 + 1) * N.y;

	//kernels packed into planes of (N.y / 2 + 1) * (N.x / 2 + 1) size
	int kerplanecount = (N.x / 2 + 1) * (N.y / 2 + 1);

	if (idx < planecount) {

		//input data
#define x(n)	(cuReIm3(cuSx[idx + (n) * planecount], cuSy[idx + (n) * planecount], cuSz[idx + (n) * planecount]))

		//no performance gain to be had from setting these as X0, X1, ... etc.
		//unrolling loops does make a slight difference though - probably last case for which you want to unroll loops
		cuReIm3 X[32];

		cuReIm3 t0, t1, t2, t3;

		//input stage

#define a	(cuBReal)0.980785280403230
#define b	(cuBReal)0.195090322016128
#define c	(cuBReal)0.923879532511287
#define d	(cuBReal)0.382683432365090
#define e	(cuBReal)0.831469612302545
#define f	(cuBReal)0.555570233019602
#define g	(cuBReal)0.707106781186548

		//j = 0
		X[0] = (x(0) + x(8));
		X[8] = (x(0) - x(8));
		X[16] = (x(0) - !x(8));
		X[24] = (x(0) + !x(8));

		//j = 1
		X[1] = (x(1) + x(9));
		X[9] = (x(1) - x(9)) * cuReIm(c, -d);
		X[17] = (x(1) - !x(9)) * cuReIm(a, -b);
		X[25] = (x(1) + !x(9)) * cuReIm(e, -f);

		//j = 2
		X[2] = (x(2) + x(10));
		X[10] = (x(2) - x(10)) * cuReIm(g, -g);
		X[18] = (x(2) - !x(10)) * cuReIm(c, -d);
		X[26] = (x(2) + !x(10)) * cuReIm(d, -c);

		//j = 3
		X[3] = (x(3) + x(11));
		X[11] = (x(3) - x(11)) * cuReIm(d, -c);
		X[19] = (x(3) - !x(11)) * cuReIm(e, -f);
		X[27] = (x(3) + !x(11)) * cuReIm(-b, -a);

		//j = 4
		X[4] = (x(4) + x(12));
		X[12] = !(x(12) - x(4));
		X[20] = (x(4) - !x(12)) * cuReIm(g, -g);
		X[28] = (x(4) + !x(12)) * cuReIm(-g, -g);

		//j = 5
		X[5] = (x(5) + x(13));
		X[13] = (x(5) - x(13)) * cuReIm(-d, -c);
		X[21] = (x(5) - !x(13)) * cuReIm(f, -e);
		X[29] = (x(5) + !x(13)) * cuReIm(-a, -b);

		//j = 6
		X[6] = (x(6) + x(14));
		X[14] = (x(6) - x(14)) * cuReIm(-g, -g);
		X[22] = (x(6) - !x(14)) * cuReIm(d, -c);
		X[30] = (x(6) + !x(14)) * cuReIm(-c, d);

		//j = 7
		X[7] = (x(7) + x(15));
		X[15] = (x(7) - x(15)) * cuReIm(-c, -d);
		X[23] = (x(7) - !x(15)) * cuReIm(b, -a);
		X[31] = (x(7) + !x(15)) * cuReIm(-f, e);

#undef x

		//final radix4 stage

		//j = 0
		t0 = (X[0] + X[4]);
		t1 = (X[0] - X[4]);
		t2 = (X[2] + X[6]);
		t3 = !(X[6] - X[2]);

		X[0] = (t0 + t2);
		X[2] = (t0 - t2);
		X[4] = (t1 + t3);
		X[6] = (t1 - t3);

		t0 = (X[8] + X[12]);
		t1 = (X[8] - X[12]);
		t2 = (X[10] + X[14]);
		t3 = !(X[14] - X[10]);

		X[8] = (t0 + t2);
		X[10] = (t0 - t2);
		X[12] = (t1 + t3);
		X[14] = (t1 - t3);

		t0 = (X[16] + X[20]);
		t1 = (X[16] - X[20]);
		t2 = (X[18] + X[22]);
		t3 = !(X[22] - X[18]);

		X[16] = (t0 + t2);
		X[18] = (t0 - t2);
		X[20] = (t1 + t3);
		X[22] = (t1 - t3);

		t0 = (X[24] + X[28]);
		t1 = (X[24] - X[28]);
		t2 = (X[26] + X[30]);
		t3 = !(X[30] - X[26]);

		X[24] = (t0 + t2);
		X[26] = (t0 - t2);
		X[28] = (t1 + t3);
		X[30] = (t1 - t3);

		//j = 1
		t0 = (X[1] + X[5]);
		t1 = (X[1] - X[5]);
		t2 = (X[3] + X[7]);
		t3 = !(X[7] - X[3]);

		X[1] = (t0 + t2);
		X[3] = !(t2 - t0);
		X[5] = (t1 + t3) * cuReIm(g, -g);
		X[7] = (t1 - t3) * cuReIm(-g, -g);

		t0 = (X[9] + X[13]);
		t1 = (X[9] - X[13]);
		t2 = (X[11] + X[15]);
		t3 = !(X[15] - X[11]);

		X[9] = (t0 + t2);
		X[11] = !(t2 - t0);
		X[13] = (t1 + t3) * cuReIm(g, -g);
		X[15] = (t1 - t3) * cuReIm(-g, -g);

		t0 = (X[17] + X[21]);
		t1 = (X[17] - X[21]);
		t2 = (X[19] + X[23]);
		t3 = !(X[23] - X[19]);

		X[17] = (t0 + t2);
		X[19] = !(t2 - t0);
		X[21] = (t1 + t3) * cuReIm(g, -g);
		X[23] = (t1 - t3) * cuReIm(-g, -g);

		t0 = (X[25] + X[29]);
		t1 = (X[25] - X[29]);
		t2 = (X[27] + X[31]);
		t3 = !(X[31] - X[27]);

		X[25] = (t0 + t2);
		X[27] = !(t2 - t0);
		X[29] = (t1 + t3) * cuReIm(g, -g);
		X[31] = (t1 - t3) * cuReIm(-g, -g);

		//radix-2 step to finish
		t0 = X[0] - X[1];
		X[0] = X[0] + X[1];
		X[1] = t0;

		t0 = X[2] - X[3];
		X[2] = X[2] + X[3];
		X[3] = t0;

		t0 = X[4] - X[5];
		X[4] = X[4] + X[5];
		X[5] = t0;

		t0 = X[6] - X[7];
		X[6] = X[6] + X[7];
		X[7] = t0;

		t0 = X[8] - X[9];
		X[8] = X[8] + X[9];
		X[9] = t0;

		t0 = X[10] - X[11];
		X[10] = X[10] + X[11];
		X[11] = t0;

		t0 = X[12] - X[13];
		X[12] = X[12] + X[13];
		X[13] = t0;

		t0 = X[14] - X[15];
		X[14] = X[14] + X[15];
		X[15] = t0;

		t0 = X[16] - X[17];
		X[16] = X[16] + X[17];
		X[17] = t0;

		t0 = X[18] - X[19];
		X[18] = X[18] + X[19];
		X[19] = t0;

		t0 = X[20] - X[21];
		X[20] = X[20] + X[21];
		X[21] = t0;

		t0 = X[22] - X[23];
		X[22] = X[22] + X[23];
		X[23] = t0;

		t0 = X[24] - X[25];
		X[24] = X[24] + X[25];
		X[25] = t0;

		t0 = X[26] - X[27];
		X[26] = X[26] + X[27];
		X[27] = t0;

		t0 = X[28] - X[29];
		X[28] = X[28] + X[29];
		X[29] = t0;

		t0 = X[30] - X[31];
		X[30] = X[30] + X[31];
		X[31] = t0;

		//output is shuffled now, i.e. it is ordered as:
		//0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31

		int i = idx % N.y;
		int j = (idx / N.y) % (N.x / 2 + 1);

		cuReIm3 F[32];

		if (i <= N.y / 2) {

			int ker_baseidx = i + j * (N.y / 2 + 1);

			F[0].x = !((KOe[ker_baseidx].x * X[0].y) + (KOe[ker_baseidx].y * X[0].z));
			F[0].y = !((-KOe[ker_baseidx].x * X[0].x) + (KOe[ker_baseidx].z * X[0].z));
			F[0].z = !((-KOe[ker_baseidx].y * X[0].x) + (-KOe[ker_baseidx].z * X[0].y));

			F[16].x = !((KOe[ker_baseidx + kerplanecount].x * X[16].y) + (KOe[ker_baseidx + kerplanecount].y * X[16].z));
			F[16].y = !((-KOe[ker_baseidx + kerplanecount].x * X[16].x) + (KOe[ker_baseidx + kerplanecount].z * X[16].z));
			F[16].z = !((-KOe[ker_baseidx + kerplanecount].y * X[16].x) + (-KOe[ker_baseidx + kerplanecount].z * X[16].y));

			F[8].x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X[8].y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X[8].z));
			F[8].y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X[8].x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X[8].z));
			F[8].z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X[8].x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X[8].y));

			F[24].x = !((KOe[ker_baseidx + 3 * kerplanecount].x * X[24].y) + (KOe[ker_baseidx + 3 * kerplanecount].y * X[24].z));
			F[24].y = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X[24].x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X[24].z));
			F[24].z = !((-KOe[ker_baseidx + 3 * kerplanecount].y * X[24].x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X[24].y));

			F[4].x = !((KOe[ker_baseidx + 4 * kerplanecount].x * X[4].y) + (KOe[ker_baseidx + 4 * kerplanecount].y * X[4].z));
			F[4].y = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X[4].x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X[4].z));
			F[4].z = !((-KOe[ker_baseidx + 4 * kerplanecount].y * X[4].x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X[4].y));

			F[20].x = !((KOe[ker_baseidx + 5 * kerplanecount].x * X[20].y) + (KOe[ker_baseidx + 5 * kerplanecount].y * X[20].z));
			F[20].y = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X[20].x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X[20].z));
			F[20].z = !((-KOe[ker_baseidx + 5 * kerplanecount].y * X[20].x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X[20].y));

			F[12].x = !((KOe[ker_baseidx + 6 * kerplanecount].x * X[12].y) + (KOe[ker_baseidx + 6 * kerplanecount].y * X[12].z));
			F[12].y = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X[12].x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X[12].z));
			F[12].z = !((-KOe[ker_baseidx + 6 * kerplanecount].y * X[12].x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X[12].y));

			F[28].x = !((KOe[ker_baseidx + 7 * kerplanecount].x * X[28].y) + (KOe[ker_baseidx + 7 * kerplanecount].y * X[28].z));
			F[28].y = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X[28].x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X[28].z));
			F[28].z = !((-KOe[ker_baseidx + 7 * kerplanecount].y * X[28].x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X[28].y));

			F[2].x = !((KOe[ker_baseidx + 8 * kerplanecount].x * X[2].y) + (KOe[ker_baseidx + 8 * kerplanecount].y * X[2].z));
			F[2].y = !((-KOe[ker_baseidx + 8 * kerplanecount].x * X[2].x) + (KOe[ker_baseidx + 8 * kerplanecount].z * X[2].z));
			F[2].z = !((-KOe[ker_baseidx + 8 * kerplanecount].y * X[2].x) + (-KOe[ker_baseidx + 8 * kerplanecount].z * X[2].y));

			F[18].x = !((KOe[ker_baseidx + 9 * kerplanecount].x * X[18].y) + (KOe[ker_baseidx + 9 * kerplanecount].y * X[18].z));
			F[18].y = !((-KOe[ker_baseidx + 9 * kerplanecount].x * X[18].x) + (KOe[ker_baseidx + 9 * kerplanecount].z * X[18].z));
			F[18].z = !((-KOe[ker_baseidx + 9 * kerplanecount].y * X[18].x) + (-KOe[ker_baseidx + 9 * kerplanecount].z * X[18].y));

			F[10].x = !((KOe[ker_baseidx + 10 * kerplanecount].x * X[10].y) + (KOe[ker_baseidx + 10 * kerplanecount].y * X[10].z));
			F[10].y = !((-KOe[ker_baseidx + 10 * kerplanecount].x * X[10].x) + (KOe[ker_baseidx + 10 * kerplanecount].z * X[10].z));
			F[10].z = !((-KOe[ker_baseidx + 10 * kerplanecount].y * X[10].x) + (-KOe[ker_baseidx + 10 * kerplanecount].z * X[10].y));

			F[26].x = !((KOe[ker_baseidx + 11 * kerplanecount].x * X[26].y) + (KOe[ker_baseidx + 11 * kerplanecount].y * X[26].z));
			F[26].y = !((-KOe[ker_baseidx + 11 * kerplanecount].x * X[26].x) + (KOe[ker_baseidx + 11 * kerplanecount].z * X[26].z));
			F[26].z = !((-KOe[ker_baseidx + 11 * kerplanecount].y * X[26].x) + (-KOe[ker_baseidx + 11 * kerplanecount].z * X[26].y));

			F[6].x = !((KOe[ker_baseidx + 12 * kerplanecount].x * X[6].y) + (KOe[ker_baseidx + 12 * kerplanecount].y * X[6].z));
			F[6].y = !((-KOe[ker_baseidx + 12 * kerplanecount].x * X[6].x) + (KOe[ker_baseidx + 12 * kerplanecount].z * X[6].z));
			F[6].z = !((-KOe[ker_baseidx + 12 * kerplanecount].y * X[6].x) + (-KOe[ker_baseidx + 12 * kerplanecount].z * X[6].y));

			F[22].x = !((KOe[ker_baseidx + 13 * kerplanecount].x * X[22].y) + (KOe[ker_baseidx + 13 * kerplanecount].y * X[22].z));
			F[22].y = !((-KOe[ker_baseidx + 13 * kerplanecount].x * X[22].x) + (KOe[ker_baseidx + 13 * kerplanecount].z * X[22].z));
			F[22].z = !((-KOe[ker_baseidx + 13 * kerplanecount].y * X[22].x) + (-KOe[ker_baseidx + 13 * kerplanecount].z * X[22].y));

			F[14].x = !((KOe[ker_baseidx + 14 * kerplanecount].x * X[14].y) + (KOe[ker_baseidx + 14 * kerplanecount].y * X[14].z));
			F[14].y = !((-KOe[ker_baseidx + 14 * kerplanecount].x * X[14].x) + (KOe[ker_baseidx + 14 * kerplanecount].z * X[14].z));
			F[14].z = !((-KOe[ker_baseidx + 14 * kerplanecount].y * X[14].x) + (-KOe[ker_baseidx + 14 * kerplanecount].z * X[14].y));

			F[30].x = !((KOe[ker_baseidx + 15 * kerplanecount].x * X[30].y) + (KOe[ker_baseidx + 15 * kerplanecount].y * X[30].z));
			F[30].y = !((-KOe[ker_baseidx + 15 * kerplanecount].x * X[30].x) + (KOe[ker_baseidx + 15 * kerplanecount].z * X[30].z));
			F[30].z = !((-KOe[ker_baseidx + 15 * kerplanecount].y * X[30].x) + (-KOe[ker_baseidx + 15 * kerplanecount].z * X[30].y));

			F[1].x = !((KOe[ker_baseidx + 16 * kerplanecount].x * X[1].y) + (KOe[ker_baseidx + 16 * kerplanecount].y * X[1].z));
			F[1].y = !((-KOe[ker_baseidx + 16 * kerplanecount].x * X[1].x) + (KOe[ker_baseidx + 16 * kerplanecount].z * X[1].z));
			F[1].z = !((-KOe[ker_baseidx + 16 * kerplanecount].y * X[1].x) + (-KOe[ker_baseidx + 16 * kerplanecount].z * X[1].y));

			F[17].x = !((-KOe[ker_baseidx + 15 * kerplanecount].x * X[17].y) + (KOe[ker_baseidx + 15 * kerplanecount].y * X[17].z));
			F[17].y = !((KOe[ker_baseidx + 15 * kerplanecount].x * X[17].x) + (KOe[ker_baseidx + 15 * kerplanecount].z * X[17].z));
			F[17].z = !((-KOe[ker_baseidx + 15 * kerplanecount].y * X[17].x) + (-KOe[ker_baseidx + 15 * kerplanecount].z * X[17].y));

			F[9].x = !((-KOe[ker_baseidx + 14 * kerplanecount].x * X[9].y) + (KOe[ker_baseidx + 14 * kerplanecount].y * X[9].z));
			F[9].y = !((KOe[ker_baseidx + 14 * kerplanecount].x * X[9].x) + (KOe[ker_baseidx + 14 * kerplanecount].z * X[9].z));
			F[9].z = !((-KOe[ker_baseidx + 14 * kerplanecount].y * X[9].x) + (-KOe[ker_baseidx + 14 * kerplanecount].z * X[9].y));

			F[25].x = !((-KOe[ker_baseidx + 13 * kerplanecount].x * X[25].y) + (KOe[ker_baseidx + 13 * kerplanecount].y * X[25].z));
			F[25].y = !((KOe[ker_baseidx + 13 * kerplanecount].x * X[25].x) + (KOe[ker_baseidx + 13 * kerplanecount].z * X[25].z));
			F[25].z = !((-KOe[ker_baseidx + 13 * kerplanecount].y * X[25].x) + (-KOe[ker_baseidx + 13 * kerplanecount].z * X[25].y));

			F[5].x = !((-KOe[ker_baseidx + 12 * kerplanecount].x * X[5].y) + (KOe[ker_baseidx + 12 * kerplanecount].y * X[5].z));
			F[5].y = !((KOe[ker_baseidx + 12 * kerplanecount].x * X[5].x) + (KOe[ker_baseidx + 12 * kerplanecount].z * X[5].z));
			F[5].z = !((-KOe[ker_baseidx + 12 * kerplanecount].y * X[5].x) + (-KOe[ker_baseidx + 12* kerplanecount].z * X[5].y));

			F[21].x = !((-KOe[ker_baseidx + 11 * kerplanecount].x * X[21].y) + (KOe[ker_baseidx + 11 * kerplanecount].y * X[21].z));
			F[21].y = !((KOe[ker_baseidx + 11 * kerplanecount].x * X[21].x) + (KOe[ker_baseidx + 11 * kerplanecount].z * X[21].z));
			F[21].z = !((-KOe[ker_baseidx + 11 * kerplanecount].y * X[21].x) + (-KOe[ker_baseidx + 11 * kerplanecount].z * X[21].y));

			F[13].x = !((-KOe[ker_baseidx + 10 * kerplanecount].x * X[13].y) + (KOe[ker_baseidx + 10 * kerplanecount].y * X[13].z));
			F[13].y = !((KOe[ker_baseidx + 10 * kerplanecount].x * X[13].x) + (KOe[ker_baseidx + 10 * kerplanecount].z * X[13].z));
			F[13].z = !((-KOe[ker_baseidx + 10 * kerplanecount].y * X[13].x) + (-KOe[ker_baseidx + 10 * kerplanecount].z * X[13].y));

			F[29].x = !((-KOe[ker_baseidx + 9 * kerplanecount].x * X[29].y) + (KOe[ker_baseidx + 9 * kerplanecount].y * X[29].z));
			F[29].y = !((KOe[ker_baseidx + 9 * kerplanecount].x * X[29].x) + (KOe[ker_baseidx + 9 * kerplanecount].z * X[29].z));
			F[29].z = !((-KOe[ker_baseidx + 9 * kerplanecount].y * X[29].x) + (-KOe[ker_baseidx + 9 * kerplanecount].z * X[29].y));

			F[3].x = !((-KOe[ker_baseidx + 8 * kerplanecount].x * X[3].y) + (KOe[ker_baseidx + 8 * kerplanecount].y * X[3].z));
			F[3].y = !((KOe[ker_baseidx + 8 * kerplanecount].x * X[3].x) + (KOe[ker_baseidx + 8 * kerplanecount].z * X[3].z));
			F[3].z = !((-KOe[ker_baseidx + 8 * kerplanecount].y * X[3].x) + (-KOe[ker_baseidx + 8 * kerplanecount].z * X[3].y));

			F[19].x = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X[19].y) + (KOe[ker_baseidx + 7 * kerplanecount].y * X[19].z));
			F[19].y = !((KOe[ker_baseidx + 7 * kerplanecount].x * X[19].x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X[19].z));
			F[19].z = !((-KOe[ker_baseidx + 7 * kerplanecount].y * X[19].x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X[19].y));

			F[11].x = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X[11].y) + (KOe[ker_baseidx + 6 * kerplanecount].y * X[11].z));
			F[11].y = !((KOe[ker_baseidx + 6 * kerplanecount].x * X[11].x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X[11].z));
			F[11].z = !((-KOe[ker_baseidx + 6 * kerplanecount].y * X[11].x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X[11].y));

			F[27].x = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X[27].y) + (KOe[ker_baseidx + 5 * kerplanecount].y * X[27].z));
			F[27].y = !((KOe[ker_baseidx + 5 * kerplanecount].x * X[27].x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X[27].z));
			F[27].z = !((-KOe[ker_baseidx + 5 * kerplanecount].y * X[27].x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X[27].y));

			F[7].x = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X[7].y) + (KOe[ker_baseidx + 4 * kerplanecount].y * X[7].z));
			F[7].y = !((KOe[ker_baseidx + 4 * kerplanecount].x * X[7].x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X[7].z));
			F[7].z = !((-KOe[ker_baseidx + 4 * kerplanecount].y * X[7].x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X[7].y));

			F[23].x = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X[23].y) + (KOe[ker_baseidx + 3 * kerplanecount].y * X[23].z));
			F[23].y = !((KOe[ker_baseidx + 3 * kerplanecount].x * X[23].x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X[23].z));
			F[23].z = !((-KOe[ker_baseidx + 3 * kerplanecount].y * X[23].x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X[23].y));

			F[15].x = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X[15].y) + (KOe[ker_baseidx + 2 * kerplanecount].y * X[15].z));
			F[15].y = !((KOe[ker_baseidx + 2 * kerplanecount].x * X[15].x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X[15].z));
			F[15].z = !((-KOe[ker_baseidx + 2 * kerplanecount].y * X[15].x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X[15].y));

			F[31].x = !((-KOe[ker_baseidx + kerplanecount].x * X[31].y) + (KOe[ker_baseidx + kerplanecount].y * X[31].z));
			F[31].y = !((KOe[ker_baseidx + kerplanecount].x * X[31].x) + (KOe[ker_baseidx + kerplanecount].z * X[31].z));
			F[31].z = !((-KOe[ker_baseidx + kerplanecount].y * X[31].x) + (-KOe[ker_baseidx + kerplanecount].z * X[31].y));
		}
		else {

			int ker_baseidx = (N.y - i) + j * (N.y / 2 + 1);

			F[0].x = !((KOe[ker_baseidx].x * X[0].y) + (-KOe[ker_baseidx].y * X[0].z));
			F[0].y = !((-KOe[ker_baseidx].x * X[0].x) + (KOe[ker_baseidx].z * X[0].z));
			F[0].z = !((KOe[ker_baseidx].y * X[0].x) + (-KOe[ker_baseidx].z * X[0].y));

			F[16].x = !((KOe[ker_baseidx + kerplanecount].x * X[16].y) + (-KOe[ker_baseidx + kerplanecount].y * X[16].z));
			F[16].y = !((-KOe[ker_baseidx + kerplanecount].x * X[16].x) + (KOe[ker_baseidx + kerplanecount].z * X[16].z));
			F[16].z = !((KOe[ker_baseidx + kerplanecount].y * X[16].x) + (-KOe[ker_baseidx + kerplanecount].z * X[16].y));

			F[8].x = !((KOe[ker_baseidx + 2 * kerplanecount].x * X[8].y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X[8].z));
			F[8].y = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X[8].x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X[8].z));
			F[8].z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X[8].x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X[8].y));

			F[24].x = !((KOe[ker_baseidx + 3 * kerplanecount].x * X[24].y) + (-KOe[ker_baseidx + 3 * kerplanecount].y * X[24].z));
			F[24].y = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X[24].x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X[24].z));
			F[24].z = !((KOe[ker_baseidx + 3 * kerplanecount].y * X[24].x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X[24].y));

			F[4].x = !((KOe[ker_baseidx + 4 * kerplanecount].x * X[4].y) + (-KOe[ker_baseidx + 4 * kerplanecount].y * X[4].z));
			F[4].y = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X[4].x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X[4].z));
			F[4].z = !((KOe[ker_baseidx + 4 * kerplanecount].y * X[4].x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X[4].y));

			F[20].x = !((KOe[ker_baseidx + 5 * kerplanecount].x * X[20].y) + (-KOe[ker_baseidx + 5 * kerplanecount].y * X[20].z));
			F[20].y = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X[20].x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X[20].z));
			F[20].z = !((KOe[ker_baseidx + 5 * kerplanecount].y * X[20].x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X[20].y));

			F[12].x = !((KOe[ker_baseidx + 6 * kerplanecount].x * X[12].y) + (-KOe[ker_baseidx + 6 * kerplanecount].y * X[12].z));
			F[12].y = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X[12].x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X[12].z));
			F[12].z = !((KOe[ker_baseidx + 6 * kerplanecount].y * X[12].x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X[12].y));

			F[28].x = !((KOe[ker_baseidx + 7 * kerplanecount].x * X[28].y) + (-KOe[ker_baseidx + 7 * kerplanecount].y * X[28].z));
			F[28].y = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X[28].x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X[28].z));
			F[28].z = !((KOe[ker_baseidx + 7 * kerplanecount].y * X[28].x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X[28].y));

			F[2].x = !((KOe[ker_baseidx + 8 * kerplanecount].x * X[2].y) + (-KOe[ker_baseidx + 8 * kerplanecount].y * X[2].z));
			F[2].y = !((-KOe[ker_baseidx + 8 * kerplanecount].x * X[2].x) + (KOe[ker_baseidx + 8 * kerplanecount].z * X[2].z));
			F[2].z = !((KOe[ker_baseidx + 8 * kerplanecount].y * X[2].x) + (-KOe[ker_baseidx + 8 * kerplanecount].z * X[2].y));

			F[18].x = !((KOe[ker_baseidx + 9 * kerplanecount].x * X[18].y) + (-KOe[ker_baseidx + 9 * kerplanecount].y * X[18].z));
			F[18].y = !((-KOe[ker_baseidx + 9 * kerplanecount].x * X[18].x) + (KOe[ker_baseidx + 9 * kerplanecount].z * X[18].z));
			F[18].z = !((KOe[ker_baseidx + 9 * kerplanecount].y * X[18].x) + (-KOe[ker_baseidx + 9 * kerplanecount].z * X[18].y));

			F[10].x = !((KOe[ker_baseidx + 10 * kerplanecount].x * X[10].y) + (-KOe[ker_baseidx + 10 * kerplanecount].y * X[10].z));
			F[10].y = !((-KOe[ker_baseidx + 10 * kerplanecount].x * X[10].x) + (KOe[ker_baseidx + 10 * kerplanecount].z * X[10].z));
			F[10].z = !((KOe[ker_baseidx + 10 * kerplanecount].y * X[10].x) + (-KOe[ker_baseidx + 10 * kerplanecount].z * X[10].y));

			F[26].x = !((KOe[ker_baseidx + 11 * kerplanecount].x * X[26].y) + (-KOe[ker_baseidx + 11 * kerplanecount].y * X[26].z));
			F[26].y = !((-KOe[ker_baseidx + 11 * kerplanecount].x * X[26].x) + (KOe[ker_baseidx + 11 * kerplanecount].z * X[26].z));
			F[26].z = !((KOe[ker_baseidx + 11 * kerplanecount].y * X[26].x) + (-KOe[ker_baseidx + 11 * kerplanecount].z * X[26].y));

			F[6].x = !((KOe[ker_baseidx + 12 * kerplanecount].x * X[6].y) + (-KOe[ker_baseidx + 12 * kerplanecount].y * X[6].z));
			F[6].y = !((-KOe[ker_baseidx + 12 * kerplanecount].x * X[6].x) + (KOe[ker_baseidx + 12 * kerplanecount].z * X[6].z));
			F[6].z = !((-KOe[ker_baseidx + 12 * kerplanecount].y * X[6].x) + (-KOe[ker_baseidx + 12 * kerplanecount].z * X[6].y));

			F[22].x = !((KOe[ker_baseidx + 13 * kerplanecount].x * X[22].y) + (KOe[ker_baseidx + 13 * kerplanecount].y * X[22].z));
			F[22].y = !((-KOe[ker_baseidx + 13 * kerplanecount].x * X[22].x) + (KOe[ker_baseidx + 13 * kerplanecount].z * X[22].z));
			F[22].z = !((KOe[ker_baseidx + 13 * kerplanecount].y * X[22].x) + (-KOe[ker_baseidx + 13 * kerplanecount].z * X[22].y));

			F[14].x = !((KOe[ker_baseidx + 14 * kerplanecount].x * X[14].y) + (-KOe[ker_baseidx + 14 * kerplanecount].y * X[14].z));
			F[14].y = !((-KOe[ker_baseidx + 14 * kerplanecount].x * X[14].x) + (KOe[ker_baseidx + 14 * kerplanecount].z * X[14].z));
			F[14].z = !((KOe[ker_baseidx + 14 * kerplanecount].y * X[14].x) + (-KOe[ker_baseidx + 14 * kerplanecount].z * X[14].y));

			F[30].x = !((KOe[ker_baseidx + 15 * kerplanecount].x * X[30].y) + (-KOe[ker_baseidx + 15 * kerplanecount].y * X[30].z));
			F[30].y = !((-KOe[ker_baseidx + 15 * kerplanecount].x * X[30].x) + (KOe[ker_baseidx + 15 * kerplanecount].z * X[30].z));
			F[30].z = !((KOe[ker_baseidx + 15 * kerplanecount].y * X[30].x) + (-KOe[ker_baseidx + 15 * kerplanecount].z * X[30].y));

			F[1].x = !((KOe[ker_baseidx + 16 * kerplanecount].x * X[1].y) + (-KOe[ker_baseidx + 16 * kerplanecount].y * X[1].z));
			F[1].y = !((-KOe[ker_baseidx + 16 * kerplanecount].x * X[1].x) + (KOe[ker_baseidx + 16 * kerplanecount].z * X[1].z));
			F[1].z = !((KOe[ker_baseidx + 16 * kerplanecount].y * X[1].x) + (-KOe[ker_baseidx + 16 * kerplanecount].z * X[1].y));

			F[17].x = !((-KOe[ker_baseidx + 15 * kerplanecount].x * X[17].y) + (-KOe[ker_baseidx + 15 * kerplanecount].y * X[17].z));
			F[17].y = !((KOe[ker_baseidx + 15 * kerplanecount].x * X[17].x) + (KOe[ker_baseidx + 15 * kerplanecount].z * X[17].z));
			F[17].z = !((KOe[ker_baseidx + 15 * kerplanecount].y * X[17].x) + (-KOe[ker_baseidx + 15 * kerplanecount].z * X[17].y));

			F[9].x = !((-KOe[ker_baseidx + 14 * kerplanecount].x * X[9].y) + (-KOe[ker_baseidx + 14 * kerplanecount].y * X[9].z));
			F[9].y = !((KOe[ker_baseidx + 14 * kerplanecount].x * X[9].x) + (KOe[ker_baseidx + 14 * kerplanecount].z * X[9].z));
			F[9].z = !((KOe[ker_baseidx + 14 * kerplanecount].y * X[9].x) + (-KOe[ker_baseidx + 14 * kerplanecount].z * X[9].y));

			F[25].x = !((-KOe[ker_baseidx + 13 * kerplanecount].x * X[25].y) + (-KOe[ker_baseidx + 13 * kerplanecount].y * X[25].z));
			F[25].y = !((KOe[ker_baseidx + 13 * kerplanecount].x * X[25].x) + (KOe[ker_baseidx + 13 * kerplanecount].z * X[25].z));
			F[25].z = !((KOe[ker_baseidx + 13 * kerplanecount].y * X[25].x) + (-KOe[ker_baseidx + 13 * kerplanecount].z * X[25].y));

			F[5].x = !((-KOe[ker_baseidx + 12 * kerplanecount].x * X[5].y) + (-KOe[ker_baseidx + 12 * kerplanecount].y * X[5].z));
			F[5].y = !((KOe[ker_baseidx + 12 * kerplanecount].x * X[5].x) + (KOe[ker_baseidx + 12 * kerplanecount].z * X[5].z));
			F[5].z = !((KOe[ker_baseidx + 12 * kerplanecount].y * X[5].x) + (-KOe[ker_baseidx + 12 * kerplanecount].z * X[5].y));

			F[21].x = !((-KOe[ker_baseidx + 11 * kerplanecount].x * X[21].y) + (-KOe[ker_baseidx + 11 * kerplanecount].y * X[21].z));
			F[21].y = !((KOe[ker_baseidx + 11 * kerplanecount].x * X[21].x) + (KOe[ker_baseidx + 11 * kerplanecount].z * X[21].z));
			F[21].z = !((KOe[ker_baseidx + 11 * kerplanecount].y * X[21].x) + (-KOe[ker_baseidx + 11 * kerplanecount].z * X[21].y));

			F[13].x = !((-KOe[ker_baseidx + 10 * kerplanecount].x * X[13].y) + (-KOe[ker_baseidx + 10 * kerplanecount].y * X[13].z));
			F[13].y = !((KOe[ker_baseidx + 10 * kerplanecount].x * X[13].x) + (KOe[ker_baseidx + 10 * kerplanecount].z * X[13].z));
			F[13].z = !((-KOe[ker_baseidx + 10 * kerplanecount].y * X[13].x) + (-KOe[ker_baseidx + 10 * kerplanecount].z * X[13].y));

			F[29].x = !((-KOe[ker_baseidx + 9 * kerplanecount].x * X[29].y) + (KOe[ker_baseidx + 9 * kerplanecount].y * X[29].z));
			F[29].y = !((KOe[ker_baseidx + 9 * kerplanecount].x * X[29].x) + (KOe[ker_baseidx + 9 * kerplanecount].z * X[29].z));
			F[29].z = !((KOe[ker_baseidx + 9 * kerplanecount].y * X[29].x) + (-KOe[ker_baseidx + 9 * kerplanecount].z * X[29].y));

			F[3].x = !((-KOe[ker_baseidx + 8 * kerplanecount].x * X[3].y) + (-KOe[ker_baseidx + 8 * kerplanecount].y * X[3].z));
			F[3].y = !((KOe[ker_baseidx + 8 * kerplanecount].x * X[3].x) + (KOe[ker_baseidx + 8 * kerplanecount].z * X[3].z));
			F[3].z = !((KOe[ker_baseidx + 8 * kerplanecount].y * X[3].x) + (-KOe[ker_baseidx + 8 * kerplanecount].z * X[3].y));

			F[19].x = !((-KOe[ker_baseidx + 7 * kerplanecount].x * X[19].y) + (-KOe[ker_baseidx + 7 * kerplanecount].y * X[19].z));
			F[19].y = !((KOe[ker_baseidx + 7 * kerplanecount].x * X[19].x) + (KOe[ker_baseidx + 7 * kerplanecount].z * X[19].z));
			F[19].z = !((KOe[ker_baseidx + 7 * kerplanecount].y * X[19].x) + (-KOe[ker_baseidx + 7 * kerplanecount].z * X[19].y));

			F[11].x = !((-KOe[ker_baseidx + 6 * kerplanecount].x * X[11].y) + (-KOe[ker_baseidx + 6 * kerplanecount].y * X[11].z));
			F[11].y = !((KOe[ker_baseidx + 6 * kerplanecount].x * X[11].x) + (KOe[ker_baseidx + 6 * kerplanecount].z * X[11].z));
			F[11].z = !((KOe[ker_baseidx + 6 * kerplanecount].y * X[11].x) + (-KOe[ker_baseidx + 6 * kerplanecount].z * X[11].y));

			F[27].x = !((-KOe[ker_baseidx + 5 * kerplanecount].x * X[27].y) + (-KOe[ker_baseidx + 5 * kerplanecount].y * X[27].z));
			F[27].y = !((KOe[ker_baseidx + 5 * kerplanecount].x * X[27].x) + (KOe[ker_baseidx + 5 * kerplanecount].z * X[27].z));
			F[27].z = !((KOe[ker_baseidx + 5 * kerplanecount].y * X[27].x) + (-KOe[ker_baseidx + 5 * kerplanecount].z * X[27].y));

			F[7].x = !((-KOe[ker_baseidx + 4 * kerplanecount].x * X[7].y) + (-KOe[ker_baseidx + 4 * kerplanecount].y * X[7].z));
			F[7].y = !((KOe[ker_baseidx + 4 * kerplanecount].x * X[7].x) + (KOe[ker_baseidx + 4 * kerplanecount].z * X[7].z));
			F[7].z = !((KOe[ker_baseidx + 4 * kerplanecount].y * X[7].x) + (-KOe[ker_baseidx + 4 * kerplanecount].z * X[7].y));

			F[23].x = !((-KOe[ker_baseidx + 3 * kerplanecount].x * X[23].y) + (-KOe[ker_baseidx + 3 * kerplanecount].y * X[23].z));
			F[23].y = !((KOe[ker_baseidx + 3 * kerplanecount].x * X[23].x) + (KOe[ker_baseidx + 3 * kerplanecount].z * X[23].z));
			F[23].z = !((KOe[ker_baseidx + 3 * kerplanecount].y * X[23].x) + (-KOe[ker_baseidx + 3 * kerplanecount].z * X[23].y));

			F[15].x = !((-KOe[ker_baseidx + 2 * kerplanecount].x * X[15].y) + (-KOe[ker_baseidx + 2 * kerplanecount].y * X[15].z));
			F[15].y = !((KOe[ker_baseidx + 2 * kerplanecount].x * X[15].x) + (KOe[ker_baseidx + 2 * kerplanecount].z * X[15].z));
			F[15].z = !((KOe[ker_baseidx + 2 * kerplanecount].y * X[15].x) + (-KOe[ker_baseidx + 2 * kerplanecount].z * X[15].y));

			F[31].x = !((-KOe[ker_baseidx + kerplanecount].x * X[31].y) + (-KOe[ker_baseidx + kerplanecount].y * X[31].z));
			F[31].y = !((KOe[ker_baseidx + kerplanecount].x * X[31].x) + (KOe[ker_baseidx + kerplanecount].z * X[31].z));
			F[31].z = !((KOe[ker_baseidx + kerplanecount].y * X[31].x) + (-KOe[ker_baseidx + kerplanecount].z * X[31].y));
		}

		//inverse z-axis fft (but without division by 32). Also only keep first 16 points.

		//radix-2 stage to start
		X[0] = F[0] + F[1];
		X[1] = F[0] - F[1];

		X[2] = F[2] + F[3];
		X[3] = F[2] - F[3];

		X[4] = F[4] + F[5];
		X[5] = F[4] - F[5];

		X[6] = F[6] + F[7];
		X[7] = F[6] - F[7];

		X[8] = F[8] + F[9];
		X[9] = F[8] - F[9];

		X[10] = F[10] + F[11];
		X[11] = F[10] - F[11];

		X[12] = F[12] + F[13];
		X[13] = F[12] - F[13];

		X[14] = F[14] + F[15];
		X[15] = F[14] - F[15];

		X[16] = F[16] + F[17];
		X[17] = F[16] - F[17];

		X[18] = F[18] + F[19];
		X[19] = F[18] - F[19];

		X[20] = F[20] + F[21];
		X[21] = F[20] - F[21];

		X[22] = F[22] + F[23];
		X[23] = F[22] - F[23];

		X[24] = F[24] + F[25];
		X[25] = F[24] - F[25];

		X[26] = F[26] + F[27];
		X[27] = F[26] - F[27];

		X[28] = F[28] + F[29];
		X[29] = F[28] - F[29];

		X[30] = F[30] + F[31];
		X[31] = F[30] - F[31];

		//First radix-4 stage

		//j = 0 (no multiplications)
		t0 = (X[0] + X[2]);
		t1 = (X[0] - X[2]);
		t2 = (X[4] + X[6]);
		t3 = !(X[6] - X[4]);

		X[0] = t0 + t2;
		X[2] = t1 - t3;
		X[4] = t0 - t2;
		X[6] = t1 + t3;

		t0 = (X[8] + X[10]);
		t1 = (X[8] - X[10]);
		t2 = (X[12] + X[14]);
		t3 = !(X[14] - X[12]);

		X[8] = t0 + t2;
		X[10] = t1 - t3;
		X[12] = t0 - t2;
		X[14] = t1 + t3;

		t0 = (X[16] + X[18]);
		t1 = (X[16] - X[18]);
		t2 = (X[20] + X[22]);
		t3 = !(X[22] - X[20]);

		X[16] = t0 + t2;
		X[18] = t1 - t3;
		X[20] = t0 - t2;
		X[22] = t1 + t3;

		t0 = (X[24] + X[26]);
		t1 = (X[24] - X[26]);
		t2 = (X[28] + X[30]);
		t3 = !(X[30] - X[28]);

		X[24] = t0 + t2;
		X[26] = t1 - t3;
		X[28] = t0 - t2;
		X[30] = t1 + t3;

		//j = 1
		t0 = (X[1] + !X[3]);
		t1 = (X[1] - !X[3]);
		t2 = (X[5] * cuReIm(g, g) + X[7] * cuReIm(-g, g));
		t3 = (X[7] * cuReIm(-g, -g) - X[5] * cuReIm(-g, g));

		X[1] = t0 + t2;
		X[3] = t1 - t3;
		X[5] = t0 - t2;
		X[7] = t1 + t3;

		t0 = (X[9] + !X[11]);
		t1 = (X[9] - !X[11]);
		t2 = (X[13] * cuReIm(g, g) + X[15] * cuReIm(-g, g));
		t3 = (X[15] * cuReIm(-g, -g) - X[13] * cuReIm(-g, g));

		X[9] = t0 + t2;
		X[11] = t1 - t3;
		X[13] = t0 - t2;
		X[15] = t1 + t3;

		t0 = (X[17] + !X[19]);
		t1 = (X[17] - !X[19]);
		t2 = (X[21] * cuReIm(g, g) + X[23] * cuReIm(-g, g));
		t3 = (X[23] * cuReIm(-g, -g) - X[21] * cuReIm(-g, g));

		X[17] = t0 + t2;
		X[19] = t1 - t3;
		X[21] = t0 - t2;
		X[23] = t1 + t3;

		t0 = (X[25] + !X[27]);
		t1 = (X[25] - !X[27]);
		t2 = (X[29] * cuReIm(g, g) + X[31] * cuReIm(-g, g));
		t3 = (X[31] * cuReIm(-g, -g) - X[29] * cuReIm(-g, g));

		X[25] = t0 + t2;
		X[27] = t1 - t3;
		X[29] = t0 - t2;
		X[31] = t1 + t3;

		//Output radix-4 stage (truncated output)
		//j = 0
		t0 = (X[0] + X[8]);
		t1 = (X[0] - X[8]);
		t2 = (X[16] + X[24]);
		t3 = !(X[24] - X[16]);

		cuReIm3 l = t0 + t2;
		cuReIm3 h = t1 - t3;

		cuSx[idx] = l.x;
		cuSy[idx] = l.y;
		cuSz[idx] = l.z;
		cuSx[idx + 8 * planecount] = h.x;
		cuSy[idx + 8 * planecount] = h.y;
		cuSz[idx + 8 * planecount] = h.z;

		//j = 1
		t0 = (X[1] + X[9] * cuReIm(c, d));
		t1 = (X[1] - X[9] * cuReIm(c, d));
		t2 = (X[17] * cuReIm(a, b) + X[25] * cuReIm(e, f));
		t3 = (X[25] * cuReIm(-f, e) - X[17] * cuReIm(-b, a));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + planecount] = l.x;
		cuSy[idx + planecount] = l.y;
		cuSz[idx + planecount] = l.z;
		cuSx[idx + 9 * planecount] = h.x;
		cuSy[idx + 9 * planecount] = h.y;
		cuSz[idx + 9 * planecount] = h.z;

		//j = 2
		t0 = (X[2] + X[10] * cuReIm(g, g));
		t1 = (X[2] - X[10] * cuReIm(g, g));
		t2 = (X[18] * cuReIm(c, d) + X[26] * cuReIm(d, c));
		t3 = (X[26] * cuReIm(-c, d) - X[18] * cuReIm(-d, c));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + 2 * planecount] = l.x;
		cuSy[idx + 2 * planecount] = l.y;
		cuSz[idx + 2 * planecount] = l.z;
		cuSx[idx + 10 * planecount] = h.x;
		cuSy[idx + 10 * planecount] = h.y;
		cuSz[idx + 10 * planecount] = h.z;

		//j = 3
		t0 = (X[3] + X[11] * cuReIm(d, c));
		t1 = (X[3] - X[11] * cuReIm(d, c));
		t2 = (X[19] * cuReIm(e, f) + X[27] * cuReIm(-b, a));
		t3 = (X[27] * cuReIm(-a, -b) - X[19] * cuReIm(-f, e));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + 3 * planecount] = l.x;
		cuSy[idx + 3 * planecount] = l.y;
		cuSz[idx + 3 * planecount] = l.z;
		cuSx[idx + 11 * planecount] = h.x;
		cuSy[idx + 11 * planecount] = h.y;
		cuSz[idx + 11 * planecount] = h.z;

		//j = 4
		t0 = (X[4] + !X[12]);
		t1 = (X[4] - !X[12]);
		t2 = (X[20] * cuReIm(g, g) + X[28] * cuReIm(-g, g));
		t3 = (X[28] * cuReIm(-g, -g) - X[20] * cuReIm(-g, g));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + 4 * planecount] = l.x;
		cuSy[idx + 4 * planecount] = l.y;
		cuSz[idx + 4 * planecount] = l.z;
		cuSx[idx + 12 * planecount] = h.x;
		cuSy[idx + 12 * planecount] = h.y;
		cuSz[idx + 12 * planecount] = h.z;

		//j = 5
		t0 = (X[5] + X[13] * cuReIm(-d, c));
		t1 = (X[5] - X[13] * cuReIm(-d, c));
		t2 = (X[21] * cuReIm(f, e) + X[29] * cuReIm(-a, b));
		t3 = (X[29] * cuReIm(-b, -a) - X[21] * cuReIm(-e, f));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + 5 * planecount] = l.x;
		cuSy[idx + 5 * planecount] = l.y;
		cuSz[idx + 5 * planecount] = l.z;
		cuSx[idx + 13 * planecount] = h.x;
		cuSy[idx + 13 * planecount] = h.y;
		cuSz[idx + 13 * planecount] = h.z;

		//j = 6
		t0 = (X[6] + X[14] * cuReIm(-g, g));
		t1 = (X[6] - X[14] * cuReIm(-g, g));
		t2 = (X[22] * cuReIm(d, c) + X[30] * cuReIm(-c, -d));
		t3 = (X[30] * cuReIm(d, -c) - X[22] * cuReIm(-c, d));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + 6 * planecount] = l.x;
		cuSy[idx + 6 * planecount] = l.y;
		cuSz[idx + 6 * planecount] = l.z;
		cuSx[idx + 14 * planecount] = h.x;
		cuSy[idx + 14 * planecount] = h.y;
		cuSz[idx + 14 * planecount] = h.z;

		//j = 7
		t0 = (X[7] + X[15] * cuReIm(-c, d));
		t1 = (X[7] - X[15] * cuReIm(-c, d));
		t2 = (X[23] * cuReIm(b, a) + X[31] * cuReIm(-f, -e));
		t3 = (X[31] * cuReIm(e, -f) - X[23] * cuReIm(-a, b));

		l = t0 + t2;
		h = t1 - t3;

		cuSx[idx + 7 * planecount] = l.x;
		cuSy[idx + 7 * planecount] = l.y;
		cuSz[idx + 7 * planecount] = l.z;
		cuSx[idx + 15 * planecount] = h.x;
		cuSy[idx + 15 * planecount] = h.y;
		cuSz[idx + 15 * planecount] = h.z;

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
#undef g
	}
}

//-------------------------- RUN-TIME KERNEL MULTIPLICATION

void OerstedKernelCUDA::KernelMultiplication_2D(void)
{
	
}

void OerstedKernelCUDA::KernelMultiplication_3D(void)
{
	cu_Oersted_ConvProd_3D_transpose_xy <<< ((N.x / 2 + 1)*N.y*N.z + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (KOe, cuS_x, cuS_y, cuS_z, cuN);
}

//Kernel multiplication in quasi-2D mode : z-axis fft / kernel multiplication / z-axis ifft rolled into one (but do not divide by N for the ifft)
void OerstedKernelCUDA::KernelMultiplication_q2D(int q2D_level)
{
	//N.z = 4, n.z = 2

	switch (q2D_level)
	{
		//N.z = 4, n.z = 2
	case 4:
		cu_Oersted_ConvProd_q2D_4_transpose_xy << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (KOe, cuS_x, cuS_y, cuS_z, cuN);
		break;

		//N.z = 8, n.z = 3, 4
	case 8:
		cu_Oersted_ConvProd_q2D_8_transpose_xy << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (KOe, cuS_x, cuS_y, cuS_z, cuN);
		break;

		//N.z = 16, n.z = 5, 6, 7, 8
	case 16:
		cu_Oersted_ConvProd_q2D_16_transpose_xy << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (KOe, cuS_x, cuS_y, cuS_z, cuN);
		break;

		//N.z = 32, n.z = 9, 10, ..., 16
	case 32:
		cu_Oersted_ConvProd_q2D_32_transpose_xy << < ((N.x / 2 + 1)*N.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (KOe, cuS_x, cuS_y, cuS_z, cuN);
		break;

		//higher values not handled in q2D mode as they are slower than full 3D mode
	}
}

#endif

#endif