#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>
#include "DemagTFunc_LCUDA.cuh"

//IRREGULAR AND SHIFTED VERSIONS FOR STRAY FIELD

//2D

__global__ void CalcTens2D_zShifted_Irregular_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuVEC<double>& f_vals_xx_del, cuVEC<double>& f_vals_yy_del, cuVEC<double>& f_vals_zz_del,
	cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
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

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && 
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(cu_mod(shift.z / d.z))) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

			//D12, D13, D23
			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x, j * d.y, shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y, i * d.x, shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * d.y, i * d.x) * sign) * s.z;
		}
		else {

			val = cuDBL3(
				Ldia_zshifted_irregular_xx_yy(i, j, 0, (pbc_images.x ? N.x : N.x / 2), (pbc_images.y ? N.y : N.y / 2), 1, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
				Ldia_zshifted_irregular_xx_yy(j, i, 0, (pbc_images.y ? N.y : N.y / 2), (pbc_images.x ? N.x : N.x / 2), 1, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
				Ldia_zshifted_irregular_zz(0, j, i, 1, (pbc_images.y ? N.y : N.y / 2), (pbc_images.x ? N.x : N.x / 2), d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;
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

__global__ void CalcTens2D_Shifted_Irregular_Ldia_PBC_kernel(
	double* D11, double* D22, double* D33,
	cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
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

		double is = cu_mod(i + (shift.x / d.x));
		double js = cu_mod(j + (shift.y / d.y));
		double ks = cu_mod(k + (shift.z / d.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y + shift.y, i * d.x + shift.x, k * d.z + shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * d.z + shift.z, j * d.y + shift.y, i * d.x + shift.x) * sign) * s.z;
		}
		else {
			
			val = cuDBL3(
				Ldia_shifted_irregular_xx_yy_single(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z, d.x, d.y, s.z, d.z),
				Ldia_shifted_irregular_xx_yy_single(j * d.y + shift.y, i * d.x + shift.x, k * d.z + shift.z, d.y, d.x, s.z, d.z),
				Ldia_shifted_irregular_zz_single(k * d.z + shift.z, j * d.y + shift.y, i * d.x + shift.x, s.z, d.z, d.y, d.x)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D11[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.x;
		D22[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.y;
		D33[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Irregular_Ldia_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y) && !z_images) {

		int i0lim = (x_images ? N.x : N.x / 2);
		int j0lim = (y_images ? N.y : N.y / 2);

		//z shift only, so use full xy plane symmetries. no pbc along z.

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

				CalcTens2D_zShifted_Irregular_Ldia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D11, D22, D33,
					f_vals_xx, f_vals_yy, f_vals_zz,
					f_vals_xx_del, f_vals_yy_del, f_vals_zz_del,
					N, s, d, shift, sign, asymptotic_distance,
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

					CalcTens2D_Shifted_Irregular_Ldia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(D11, D22, D33,
						N, s, d, shift, sign, asymptotic_distance,
						demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz,
						cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
				}
			}
		}
	}
}

__global__ void CalcTens2D_zShifted_Irregular_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuVEC<double>& g_vals_xy_del, cuVEC<double>& g_vals_xz_del, cuVEC<double>& g_vals_yz_del,
	cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
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

		//apply asymptotic equations?
		//apply asymptotic equations?
		if (asymptotic_distance > 0 && 
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(cu_mod(shift.z / d.z))) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

			//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
			//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
			//D12, D13, D23
			val += (cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x, j * d.y, shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x, shift.z, j * d.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y, shift.z, i * d.x) * sign)
				& cuDBL3(sign_i * sign_j, sign_i, sign_j)) * s.z;
		}
		else {

			val = cuDBL3(
				Lodia_zshifted_irregular_xy(i, j, 0, (pbc_images.x ? N.x : N.x / 2), (pbc_images.y ? N.y : N.y / 2), 1, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del) * sign,
				Lodia_zshifted_irregular_xz_yz(i, 0, j, (pbc_images.x ? N.x : N.x / 2), 1, (pbc_images.y ? N.y : N.y / 2), d.x*d.y*d.z, g_vals_xz, g_vals_xz_del) * sign,
				Lodia_zshifted_irregular_xz_yz(j, 0, i, (pbc_images.y ? N.y : N.y / 2), 1, (pbc_images.x ? N.x : N.x / 2), d.x*d.y*d.z, g_vals_yz, g_vals_yz_del) * sign)
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

__global__ void CalcTens2D_Shifted_Irregular_Lodia_PBC_kernel(
	double* D12, double* D13, double* D23,
	cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
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

		double is = cu_mod(i + (shift.x / d.x));
		double js = cu_mod(j + (shift.y / d.y));
		double ks = cu_mod(k + (shift.z / d.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
			int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z),
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x + shift.x, k * d.z + shift.z, j * d.y + shift.y),
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y + shift.y, k * d.z + shift.z, i * d.x + shift.x)) * sign * s.z;
		}
		else {

			val = cuDBL3(
				Lodia_shifted_irregular_xy_single(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z, d.x, d.y, s.z, d.z),
				Lodia_shifted_irregular_xz_yz_single(i * d.x + shift.x, k * d.z + shift.z, j * d.y + shift.y, d.x, s.z, d.z, d.y),
				Lodia_shifted_irregular_xz_yz_single(j * d.y + shift.y, k * d.z + shift.z, i * d.x + shift.x, d.y, s.z, d.z, d.x)) * sign;
		}

		//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
		D12[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.x;
		D13[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.y;
		D23[((i0 + N.x) % N.x) + ((j0 + N.y) % N.y) * N.x] += val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Irregular_Lodia_PBC(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y) && !z_images) {

		int i0lim = (x_images ? N.x : N.x / 2);
		int j0lim = (y_images ? N.y : N.y / 2);

		//z shift only, so use full xy plane symmetries. no pbc along z.

		for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
			for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

				CalcTens2D_zShifted_Irregular_Lodia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(D12, D13, D23,
					g_vals_xy, g_vals_xz, g_vals_yz,
					g_vals_xy_del, g_vals_xz_del, g_vals_yz_del,
					N, s, d, shift, sign, asymptotic_distance,
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

					CalcTens2D_Shifted_Irregular_Lodia_PBC_kernel <<< (i0lim*j0lim + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(D12, D13, D23,
						N, s, d, shift, sign, asymptotic_distance,
						demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz,
						cuINT3(x_images, y_images, z_images), cuINT3(i_img, j_img, k_img));
				}
			}
		}
	}
}

#endif
