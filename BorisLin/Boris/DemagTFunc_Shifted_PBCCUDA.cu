#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>
#include "DemagTFunc_LCUDA.cuh"

//SHIFTED VERSIONS FOR STRAY FIELD

//2D

__global__ void CalcTens2D_zShifted_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
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

		double ks = cu_mod(shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && 
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_zshifted_xx_yy(i, j, 0, (pbc_images.x ? N.x : N.x / 2), (pbc_images.y ? N.y : N.y / 2), 1, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
				Ldia_zshifted_xx_yy(j, i, 0, (pbc_images.y ? N.y : N.y / 2), (pbc_images.x ? N.x : N.x / 2), 1, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
				Ldia_zshifted_zz(0, j, i, 1, (pbc_images.y ? N.y : N.y / 2), (pbc_images.x ? N.x : N.x / 2), hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
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

__global__ void CalcTens2D_Shifted_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : 2 * (N.x / 2) - 1);
	int j0lim = (pbc_images.y ? N.y : 2 * (N.y / 2) - 1);

	if (idx < i0lim * j0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;

		if (!pbc_images.x) i0 += -N.x / 2 + 1;
		if (!pbc_images.y) j0 += -N.y / 2 + 1;

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = image.k;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(k + (shift.z / hRatios.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
				Ldia_single(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z, hRatios.y, hRatios.x, hRatios.z),
				Ldia_single(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x, hRatios.z, hRatios.y, hRatios.x)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D11[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.x;
		D22[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.y;
		D33[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Ldia_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y) && !z_images) {

		int i0lim = (x_images ? N.x : N.x / 2);
		int j0lim = (y_images ? N.y : N.y / 2);

		//z shift only, so use full xy plane symmetries. no pbc along z.

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

				CalcTens2D_zShifted_Ldia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D11, D22, D33,
					f_vals_xx, f_vals_yy, f_vals_zz,
					N, hRatios, shift, sign, asymptotic_distance,
					demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
					cuINT3(x_images, y_images, 0), cuINT3(i_img, j_img, 0));
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible

		int i0lim = (x_images ? N.x : 2 * (N.x / 2) - 1);
		int j0lim = (y_images ? N.y : 2 * (N.y / 2) - 1);

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
				for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

					CalcTens2D_Shifted_Ldia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(D11, D22, D33,
						N, hRatios, shift, sign, asymptotic_distance,
						demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
						cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
				}
			}
		}
	}
}

__global__ void CalcTens2D_zShifted_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);

	if (idx < i0lim * j0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;

		//use modulus of indexes and adjust for tensor element signs based on symmetries
		int sign_i = cu_get_sign(i);
		int sign_j = cu_get_sign(j);
		i = abs(i);
		j = abs(j);

		double ks = cu_mod(shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && 
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
			//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
			//D12, D13, D23
			val += cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, shift.z, j * hRatios.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, shift.z, i * hRatios.x) * sign)
				& cuDBL3(sign_i * sign_j, sign_i, sign_j);
		}
		else {

			val = cuDBL3(
				Lodia_xy_zshifted(i, j, 0, (pbc_images.x ? N.x : N.x / 2), (pbc_images.y ? N.y : N.y / 2), 1, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
				Lodia_xz_yz_zshifted(i, 0, j, (pbc_images.x ? N.x : N.x / 2), 1, (pbc_images.y ? N.y : N.y / 2), hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
				Lodia_xz_yz_zshifted(j, 0, i, (pbc_images.y ? N.y : N.y / 2), 1, (pbc_images.x ? N.x : N.x / 2), hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign)
				& cuDBL3(sign_i * sign_j, sign_i, sign_j);
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[i0 + j0 * N.x] += val.x;
		D13[i0 + j0 * N.x] += val.y;
		D23[i0 + j0 * N.x] += val.z;

		if (!pbc_images.x) {

			if (i0) {

				D12[((N.x - i0) % N.x) + j0 * N.x] += -val.x;
				D13[((N.x - i0) % N.x) + j0 * N.x] += -val.y;
				D23[((N.x - i0) % N.x) + j0 * N.x] += val.z;
			}

			if (!pbc_images.y) {

				if (i0 && j0) {

					D12[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += val.x;
					D13[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += -val.y;
					D23[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] += -val.z;
				}
			}
		}

		if (!pbc_images.y) {

			if (j0) {

				D12[i0 + ((N.y - j0) % N.y) * N.x] += -val.x;
				D13[i0 + ((N.y - j0) % N.y) * N.x] += val.y;
				D23[i0 + ((N.y - j0) % N.y) * N.x] += -val.z;
			}
		}
	}
}

__global__ void CalcTens2D_Shifted_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : 2 * (N.x / 2) - 1);
	int j0lim = (pbc_images.y ? N.y : 2 * (N.y / 2) - 1);

	if (idx < i0lim * j0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;

		if (!pbc_images.x) i0 += -N.x / 2 + 1;
		if (!pbc_images.y) j0 += -N.y / 2 + 1;

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = image.k;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(k + (shift.z / hRatios.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z),
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y),
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x)) * sign;
		}
		else {

			val = cuDBL3(
				Lodia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
				Lodia_single(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y, hRatios.x, hRatios.z, hRatios.y),
				Lodia_single(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x, hRatios.y, hRatios.z, hRatios.x)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.x;
		D13[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.y;
		D23[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Lodia_PBC(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y) && !z_images) {

		int i0lim = (x_images ? N.x : N.x / 2);
		int j0lim = (y_images ? N.y : N.y / 2);

		//z shift only, so use full xy plane symmetries. no pbc along z.

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

				CalcTens2D_zShifted_Lodia_PBC_kernel << < (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
					(D12, D13, D23,
					g_vals_xy, g_vals_xz, g_vals_yz,
					N, hRatios, shift, sign, asymptotic_distance,
					demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz,
					cuINT3(x_images, y_images, 0), cuINT3(i_img, j_img, 0));
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible

		int i0lim = (x_images ? N.x : 2 * (N.x / 2) - 1);
		int j0lim = (y_images ? N.y : 2 * (N.y / 2) - 1);

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
				for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

					CalcTens2D_Shifted_Lodia_PBC_kernel << < (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
						(D12, D13, D23,
						N, hRatios, shift, sign, asymptotic_distance,
						demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz,
						cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
				}
			}
		}
	}
}

//3D

__global__ void CalcTens3D_zShifted_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);

	if (idx < i0lim * j0lim * (2 * (N.z/2) - 1)) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;
		int k0 = (idx / (i0lim * j0lim)) - N.z/2 + 1;

		cuDBL3 val;

		//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
		//There is significant redundancy remaining, could be optimized further.
		int i = abs(image.i * N.x + i0);
		int j = abs(image.j * N.y + j0);
		int k = k0;

		double ks = cu_mod(k + shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_zshifted_xx_yy(i, j, k, (pbc_images.x ? N.x : N.x / 2), (pbc_images.y ? N.y : N.y / 2), N.z / 2, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
				Ldia_zshifted_xx_yy(j, i, k, (pbc_images.y ? N.y : N.y / 2), (pbc_images.x ? N.x : N.x / 2), N.z / 2, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
				Ldia_zshifted_zz(k, j, i, N.z / 2, (pbc_images.y ? N.y : N.y / 2), (pbc_images.x ? N.x : N.x / 2), hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D11[i0 + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
		D22[i0 + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
		D33[i0 + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;

		if (!pbc_images.x) {

			if (i0) {
				D11[((N.x - i0) % N.x) + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
				D22[((N.x - i0) % N.x) + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
				D33[((N.x - i0) % N.x) + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;
			}

			if (!pbc_images.y) {

				if (i0 && j0) {
					D11[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
					D22[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
					D33[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;
				}
			}
		}

		if (!pbc_images.y) {

			if (j0) {
				D11[i0 + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
				D22[i0 + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
				D33[i0 + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;
			}
		}
	}
}

__global__ void CalcTens3D_Shifted_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : 2 * (N.x / 2) - 1);
	int j0lim = (pbc_images.y ? N.y : 2 * (N.y / 2) - 1);
	int k0lim = (pbc_images.z ? N.z : 2 * (N.z / 2) - 1);

	if (idx < i0lim * j0lim * k0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;
		int k0 = idx / (i0lim * j0lim);

		if (!pbc_images.x) i0 += -N.x / 2 + 1;
		if (!pbc_images.y) j0 += -N.y / 2 + 1;
		if (!pbc_images.z) k0 += -N.z / 2 + 1;

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = image.k * N.z + k0;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(k + (shift.z / hRatios.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
				Ldia_single(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z, hRatios.y, hRatios.x, hRatios.z),
				Ldia_single(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x, hRatios.z, hRatios.y, hRatios.x)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D11[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
		D22[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
		D33[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens3D_Shifted_Ldia_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y) && !z_images) {

		int i0lim = (x_images ? N.x : N.x / 2);
		int j0lim = (y_images ? N.y : N.y / 2);

		//z shift only, so use full xy plane symmetries. no pbc along z.

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

				CalcTens3D_zShifted_Ldia_PBC_kernel <<< (i0lim * j0lim * (2*(N.z/2) - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D11, D22, D33,
					f_vals_xx, f_vals_yy, f_vals_zz,
					N, hRatios, shift, sign, asymptotic_distance,
					demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
					cuINT3(x_images, y_images, 0), cuINT3(i_img, j_img, 0));
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible

		int i0lim = (x_images ? N.x : 2 * (N.x / 2) - 1);
		int j0lim = (y_images ? N.y : 2 * (N.y / 2) - 1);
		int k0lim = (z_images ? N.z : 2 * (N.z / 2) - 1);

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
				for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

					CalcTens3D_Shifted_Ldia_PBC_kernel <<< (i0lim*j0lim*k0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(D11, D22, D33,
						N, hRatios, shift, sign, asymptotic_distance,
						demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
						cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
				}
			}
		}
	}
}

__global__ void CalcTens3D_zShifted_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : N.x / 2);
	int j0lim = (pbc_images.y ? N.y : N.y / 2);

	if (idx < i0lim * j0lim * (2 * (N.z/2) - 1)) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;
		int k0 = (idx / (i0lim * j0lim)) - N.z/2 + 1;

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = k0;

		//use modulus of indexes and adjust for tensor element signs based on symmetries
		int sign_i = cu_get_sign(i);
		int sign_j = cu_get_sign(j);
		i = abs(i);
		j = abs(j);

		double ks = cu_mod(k + shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
			//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
			//D12, D13, D23
			val += cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z + shift.z, j * hRatios.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z + shift.z, i * hRatios.x) * sign)
				& cuDBL3(sign_i * sign_j, sign_i, sign_j);
		}
		else {

			val = cuDBL3(
				Lodia_xy_zshifted(i, j, k, (pbc_images.x ? N.x : N.x / 2), (pbc_images.y ? N.y : N.y / 2), N.z / 2, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
				Lodia_xz_yz_zshifted(i, k, j, (pbc_images.x ? N.x : N.x / 2), N.z / 2, (pbc_images.y ? N.y : N.y / 2), hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
				Lodia_xz_yz_zshifted(j, k, i, (pbc_images.y ? N.y : N.y / 2), N.z / 2, (pbc_images.x ? N.x : N.x / 2), hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign)
				& cuDBL3(sign_i * sign_j, sign_i, sign_j);
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[i0 + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
		D13[i0 + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
		D23[i0 + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;

		if (!pbc_images.x) {

			if (i0) {

				D12[((N.x - i0) % N.x) + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += -val.x;
				D13[((N.x - i0) % N.x) + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += -val.y;
				D23[((N.x - i0) % N.x) + j0 * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;
			}

			if (!pbc_images.y) {

				if (i0 && j0) {

					D12[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
					D13[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += -val.y;
					D23[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += -val.z;
				}
			}
		}

		if (!pbc_images.y) {

			if (j0) {

				D12[i0 + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += -val.x;
				D13[i0 + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
				D23[i0 + ((N.y - j0) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += -val.z;
			}
		}
	}
}

__global__ void CalcTens3D_Shifted_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz,
	cuINT3 pbc_images, cuINT3 image)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int i0lim = (pbc_images.x ? N.x : 2 * (N.x / 2) - 1);
	int j0lim = (pbc_images.y ? N.y : 2 * (N.y / 2) - 1);
	int k0lim = (pbc_images.z ? N.z : 2 * (N.z / 2) - 1);

	if (idx < i0lim * j0lim * k0lim) {

		int i0 = idx % i0lim;
		int j0 = (idx / i0lim) % j0lim;
		int k0 = idx / (i0lim * j0lim);

		if (!pbc_images.x) i0 += -N.x / 2 + 1;
		if (!pbc_images.y) j0 += -N.y / 2 + 1;
		if (!pbc_images.z) k0 += -N.z / 2 + 1;

		cuDBL3 val;

		int i = image.i * N.x + i0;
		int j = image.j * N.y + j0;
		int k = image.k * N.z + k0;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(k + (shift.z / hRatios.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z),
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y),
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x)) * sign;
		}
		else {

			val = cuDBL3(
				Lodia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
				Lodia_single(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y, hRatios.x, hRatios.z, hRatios.y),
				Lodia_single(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x, hRatios.y, hRatios.z, hRatios.x)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.x;
		D13[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.y;
		D23[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x + ((k0 + N.z) % N.z) * N.x*N.y] += val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens3D_Shifted_Lodia_PBC(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y) && !z_images) {

		int i0lim = (x_images ? N.x : N.x / 2);
		int j0lim = (y_images ? N.y : N.y / 2);

		//z shift only, so use full xy plane symmetries. no pbc along z.

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
				
				CalcTens3D_zShifted_Lodia_PBC_kernel <<< (i0lim * j0lim * (2 * (N.z/2) - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D12, D13, D23,
					g_vals_xy, g_vals_xz, g_vals_yz,
					N, hRatios, shift, sign, asymptotic_distance,
					demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz,
					cuINT3(x_images, y_images, 0), cuINT3(i_img, j_img, 0));
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible

		int i0lim = (x_images ? N.x : 2 * (N.x / 2) - 1);
		int j0lim = (y_images ? N.y : 2 * (N.y / 2) - 1);
		int k0lim = (z_images ? N.z : 2 * (N.z / 2) - 1);

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
				for (int k_img = -z_images; k_img < z_images + 1; k_img++) {
					
					CalcTens3D_Shifted_Lodia_PBC_kernel <<< (i0lim*j0lim*k0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(D12, D13, D23,
						N, hRatios, shift, sign, asymptotic_distance,
						demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz,
						cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
				}
			}
		}
	}
}

#endif
