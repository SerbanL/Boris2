#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>
#include "DemagTFunc_LCUDA.cuh"

//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

//REGULAR VERSIONS FOR INTERNAL FIELD

//2D

__global__ void CalcTens2D_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);

	if (idx < i0lim * j0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;

		cuDBL3 val;

		//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
		//There is significant redundancy remaining, could be optimized further.
		int i = abs(image.i * N.x + i0);
		int j = abs(image.j * N.y + j0);
		int k = abs(image.k);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {
						
			//D12, D13, D23
			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {
						
			val = cuDBL3(
				Ldia(i, j, k, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
				Ldia(j, i, k, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
				Ldia(k, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
		}
		
		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D11[i0 + j0 * N.x] += val.x;
		D22[i0 + j0 * N.x] += val.y;
		D33[i0 + j0 * N.x] += val.z;
		
		if (!pbc_images.x) {

			if (i0) {
				D11[((N.x - i0) % N.x) + j0 * N.x] += val.x;
				D22[((N.x - i0) % N.x) + j0 * N.x] += val.y;
				D33[((N.x - i0) % N.x) + j0 * N.x] += val.z;
			}

			if (!pbc_images.y) {

				if (i0 && j0) {
					D11[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += val.x;
					D22[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += val.y;
					D33[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += val.z;
				}
			}
		}

		if (!pbc_images.y) {

			if (j0) {
				D11[i0 + ((N.y - j0) % N.y) * N.x] += val.x;
				D22[i0 + ((N.y - j0) % N.y) * N.x] += val.y;
				D33[i0 + ((N.y - j0) % N.y) * N.x] += val.z;
			}
		}
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Ldia_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	int i0lim = (x_images ? N.x : N.x / 2);
	int j0lim = (y_images ? N.y : N.y / 2);

	for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
		for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
			for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

				CalcTens2D_Ldia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D11, D22, D33,
					f_vals_xx, f_vals_yy, f_vals_zz,
					N, hRatios, sign, asymptotic_distance,
					demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
					cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));

			}
		}
	}
}

__global__ void CalcTens2D_Lodia_PBC_kernel(
	double* D12,
	cuVEC<double>& g_vals_xy,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);

	if (idx < i0lim * j0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;

		double val = 0.0;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = image.k;

		//use modulus of indexes and adjust for tensor element signs based on symmetries
		int sign_i = cu_get_sign(i);
		int sign_j = cu_get_sign(j);
		i = abs(i);
		j = abs(j);
		k = abs(k);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

			//D12, D13, D23
			val = demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign * sign_i * sign_j;
		}
		else {

			val = Lodia(i, j, k, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign * sign_i * sign_j;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[i0 + j0 * N.x] += val;

		if (!pbc_images.x) {

			if (i0) D12[((N.x - i0) % N.x) + j0 * N.x] += -val;

			if (!pbc_images.y) {

				if (i0 && j0) D12[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += val;
			}
		}

		if (!pbc_images.y) {

			if (j0) D12[i0 + ((N.y - j0) % N.y) * N.x] += -val;
		}
	}
}

void DemagTFuncCUDA::CalcTens2D_Lodia_PBC(
	cu_arr<double>& D12,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	int i0lim = (x_images ? N.x : N.x / 2);
	int j0lim = (y_images ? N.y : N.y / 2);

	for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
		for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
			for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

				CalcTens2D_Lodia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D12,
					g_vals_xy,
					N, hRatios, sign, asymptotic_distance,
					demagAsymptoticOffDiag_xy,
					cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
			}
		}
	}
}

//3D

__global__ void CalcTens3D_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);
	int k0lim = (pbc_images.z ? N.z : N.z / 2);

	if (idx < i0lim * j0lim * k0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;
		int k0 = idx / (i0lim * j0lim);

		cuDBL3 val;

		//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
		//There is significant redundancy remaining, could be optimized further.
		int i = abs(image.i * N.x + i0);
		int j = abs(image.j * N.y + j0);
		int k = abs(image.k * N.z + k0);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

			//D12, D13, D23
			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {
						
			val = cuDBL3(
				Ldia(i, j, k, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
				Ldia(j, i, k, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
				Ldia(k, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
		}
		

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D11[i0 + j0 * N.x + k0 * N.x*N.y] += val.x;
		D22[i0 + j0 * N.x + k0 * N.x*N.y] += val.y;
		D33[i0 + j0 * N.x + k0 * N.x*N.y] += val.z;

		if (!pbc_images.x) {

			if (i0) {
				D11[((N.x - i0) % N.x) + j0 * N.x + k0 * N.x*N.y] += val.x;
				D22[((N.x - i0) % N.x) + j0 * N.x + k0 * N.x*N.y] += val.y;
				D33[((N.x - i0) % N.x) + j0 * N.x + k0 * N.x*N.y] += val.z;
			}

			if (!pbc_images.y) {

				if (i0 && j0) {
					D11[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.x;
					D22[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.y;
					D33[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.z;
				}

				if (!pbc_images.z) {

					if (i0 && j0 && k0) {
						D11[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.x;
						D22[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.y;
						D33[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.z;
					}
				}
			}

			if (!pbc_images.z) {

				if (i0 && k0) {
					D11[((N.x - i0) % N.x) + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.x;
					D22[((N.x - i0) % N.x) + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.y;
					D33[((N.x - i0) % N.x) + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.z;
				}
			}
		}

		if (!pbc_images.y) {

			if (j0) {
				D11[i0 + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.x;
				D22[i0 + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.y;
				D33[i0 + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.z;
			}

			if (!pbc_images.z) {

				if (j0 && k0) {
					D11[i0 + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.x;
					D22[i0 + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.y;
					D33[i0 + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.z;
				}
			}
		}
		
		if (!pbc_images.z) {

			if (k0) {
				D11[i0 + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.x;
				D22[i0 + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.y;
				D33[i0 + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.z;
			}
		}
	}
}

void DemagTFuncCUDA::CalcTens3D_Ldia_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	int i0lim = (x_images ? N.x : N.x / 2);
	int j0lim = (y_images ? N.y : N.y / 2);
	int k0lim = (z_images ? N.z : N.z / 2);
	
	for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
		for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
			for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

				CalcTens3D_Ldia_PBC_kernel <<< (i0lim*j0lim*k0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D11, D22, D33,
					f_vals_xx, f_vals_yy, f_vals_zz,
					N, hRatios, sign, asymptotic_distance,
					demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
					cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
			}
		}
	}
}

__global__ void CalcTens3D_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);
	int k0lim = (pbc_images.z ? N.z : N.z / 2);

	if (idx < i0lim * j0lim * k0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;
		int k0 = idx / (i0lim * j0lim);

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = image.k * N.z + k0;

		//use modulus of indexes and adjust for tensor element signs based on symmetries
		int sign_i = cu_get_sign(i);
		int sign_j = cu_get_sign(j);
		int sign_k = cu_get_sign(k);
		i = abs(i);
		j = abs(j);
		k = abs(k);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

			//D12, D13, D23
			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z, j * hRatios.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z, i * hRatios.x) * sign)
				& cuDBL3(sign_i * sign_j, sign_i * sign_k, sign_j * sign_k);
		}
		else {
					
			val = cuDBL3(
				Lodia(i, j, k, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
				Lodia(i, k, j, hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
				Lodia(j, k, i, hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign)
				& cuDBL3(sign_i * sign_j, sign_i * sign_k, sign_j * sign_k);
		}
		
		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[i0 + j0 * N.x + k0 * N.x*N.y] += val.x;
		D13[i0 + j0 * N.x + k0 * N.x*N.y] += val.y;
		D23[i0 + j0 * N.x + k0 * N.x*N.y] += val.z;

		if (!pbc_images.x) {

			if (i0) {
				D12[((N.x - i0) % N.x) + j0 * N.x + k0 * N.x*N.y] += -val.x;
				D13[((N.x - i0) % N.x) + j0 * N.x + k0 * N.x*N.y] += -val.y;
				D23[((N.x - i0) % N.x) + j0 * N.x + k0 * N.x*N.y] += val.z;
			}

			if (!pbc_images.y) {

				if (i0 && j0) {
					D12[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.x;
					D13[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += -val.y;
					D23[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += -val.z;
				}

				if (!pbc_images.z) {

					if (i0 && j0 && k0) {
						D12[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.x;
						D13[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.y;
						D23[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.z;
					}
				}
			}

			if (!pbc_images.z) {

				if (i0 && k0) {
					D12[((N.x - i0) % N.x) + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += -val.x;
					D13[((N.x - i0) % N.x) + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.y;
					D23[((N.x - i0) % N.x) + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += -val.z;
				}
			}
		}

		if (!pbc_images.y) {

			if (j0) {
				D12[i0 + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += -val.x;
				D13[i0 + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += val.y;
				D23[i0 + ((N.y - j0) % N.y) * N.x + k0 * N.x*N.y] += -val.z;
			}

			if (!pbc_images.z) {

				if (j0 && k0) {
					D12[i0 + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += -val.x;
					D13[i0 + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += -val.y;
					D23[i0 + ((N.y - j0) % N.y) * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.z;
				}
			}
		}

		if (!pbc_images.z) {

			if (k0) {
				D12[i0 + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += val.x;
				D13[i0 + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += -val.y;
				D23[i0 + j0 * N.x + ((N.z - k0) % N.z) * N.x*N.y] += -val.z;
			}
		}
	}
}

void DemagTFuncCUDA::CalcTens3D_Lodia_PBC(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	int i0lim = (x_images ? N.x : N.x / 2);
	int j0lim = (y_images ? N.y : N.y / 2);
	int k0lim = (z_images ? N.z : N.z / 2);

	for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
		for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
			for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

				CalcTens3D_Lodia_PBC_kernel <<< (i0lim*j0lim*k0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D12, D13, D23,
					g_vals_xy, g_vals_xz, g_vals_yz,
					N, hRatios, sign, asymptotic_distance,
					demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz,
					cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
			}
		}
	}
}

#endif
