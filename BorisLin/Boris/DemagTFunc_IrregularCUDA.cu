#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>
#include "DemagTFunc_LCUDA.cuh"

//IRREGULAR AND SHIFTED VERSIONS FOR STRAY FIELD

//2D

//z shift only, so use full xy plane symmetries
__global__ void CalcTens2D_zShifted_Irregular_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuVEC<double>& f_vals_xx_del, cuVEC<double>& f_vals_yy_del, cuVEC<double>& f_vals_zz_del,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;

		cuDBL3 val;

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && 
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(cu_mod(shift.z / d.z))) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x, j * d.y, shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y, i * d.x, shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * d.y, i * d.x) * sign) * s.z;
		}
		else {

			val = cuDBL3(
				Ldia_zshifted_irregular_xx_yy(i, j, 0, n.x, n.y, n.z, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
				Ldia_zshifted_irregular_xx_yy(j, i, 0, n.y, n.x, n.z, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
				Ldia_zshifted_irregular_zz(0, j, i, n.z, n.y, n.x, d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;
		}

		D11[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.x;
		D22[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.y;
		D33[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.z;

		D11[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.x;
		D22[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.y;
		D33[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.z;

		D11[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.x;
		D22[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.y;
		D33[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.z;

		D11[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.x;
		D22[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.y;
		D33[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.z;

		D11[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.x;
		D22[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.y;
		D33[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.z;
	}
}

//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)
__global__ void CalcTens2D_Shifted_Irregular_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuVEC<double>& f_vals_xx_del, cuVEC<double>& f_vals_yy_del, cuVEC<double>& f_vals_zz_del,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (2 * n.y - 1) * (2 * n.x - 1)) {

		int i = (idx % (2 * n.x - 1)) - n.x + 1;
		int j = ((idx / (2 * n.x - 1)) % (2 * n.y - 1)) - n.y + 1;

		cuDBL3 val;

		double is = cu_mod(i + (shift.x / d.x));
		double js = cu_mod(j + (shift.y / d.y));
		double ks = cu_mod(shift.z / d.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
			//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x + shift.x, j * d.y + shift.y, shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y + shift.y, i * d.x + shift.x, shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * d.y + shift.y, i * d.x + shift.x) * sign) * s.z;
		}
		else {

			val = cuDBL3(
				Ldia_shifted_irregular_xx_yy(i, j, 0, n.x, n.y, n.z, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
				Ldia_shifted_irregular_xx_yy(j, i, 0, n.y, n.x, n.z, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
				Ldia_shifted_irregular_zz(0, j, i, n.z, n.y, n.x, d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;
		}

		D11[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.x;
		D22[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.y;
		D33[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Irregular_Ldia(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

		CalcTens2D_zShifted_Irregular_Ldia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D11, D22, D33,
			f_vals_xx, f_vals_yy, f_vals_zz,
			f_vals_xx_del, f_vals_yy_del, f_vals_zz_del,
			n, N, s, d, shift, sign, asymptotic_distance,
			demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

		CalcTens2D_Shifted_Irregular_Ldia_kernel <<< ((2 * n.y - 1) * (2 * n.x - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D11, D22, D33,
			f_vals_xx, f_vals_yy, f_vals_zz,
			f_vals_xx_del, f_vals_yy_del, f_vals_zz_del,
			n, N, s, d, shift, sign, asymptotic_distance,
			demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
	}
}

//z shift only, so use full xy plane symmetries
__global__ void CalcTens2D_zShifted_Irregular_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuVEC<double>& g_vals_xy_del, cuVEC<double>& g_vals_xz_del, cuVEC<double>& g_vals_yz_del,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;

		cuDBL3 val;

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && 
			(i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(cu_mod(shift.z / d.z))) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

			//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
			//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
			//D12, D13, D23
			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x, j * d.y, shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x, shift.z, j * d.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y, shift.z, i * d.x) * sign) * s.z;
		}
		else {

			val = cuDBL3(
				Lodia_zshifted_irregular_xy(i, j, 0, n.x, n.y, n.z, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del),
				Lodia_zshifted_irregular_xz_yz(i, 0, j, n.x, n.z, n.y, d.x*d.y*d.z, g_vals_xz, g_vals_xz_del),
				Lodia_zshifted_irregular_xz_yz(j, 0, i, n.y, n.z, n.x, d.x*d.y*d.z, g_vals_yz, g_vals_yz_del)) * sign;
		}

		D12[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.x;
		D13[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.y;
		D23[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.z;

		D12[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = -val.x;
		D13[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = -val.y;
		D23[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val.z;

		D12[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = -val.x;
		D13[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.y;
		D23[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = -val.z;

		D12[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val.x;
		D13[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = -val.y;
		D23[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = -val.z;
	}
}

//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)
__global__ void CalcTens2D_Shifted_Irregular_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuVEC<double>& g_vals_xy_del, cuVEC<double>& g_vals_xz_del, cuVEC<double>& g_vals_yz_del,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (2 * n.y - 1) * (2 * n.x - 1)) {

		int i = (idx % (2 * n.x - 1)) - n.x + 1;
		int j = ((idx / (2 * n.x - 1)) % (2 * n.y - 1)) - n.y + 1;

		cuDBL3 val;

		double is = cu_mod(i + (shift.x / d.x));
		double js = cu_mod(j + (shift.y / d.y));
		double ks = cu_mod(shift.z / d.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
			int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
			//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
			//D12, D13, D23
			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x + shift.x, j * d.y + shift.y, shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x + shift.x, shift.z, j * d.y + shift.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y + shift.y, shift.z, i * d.x + shift.x) * sign) * s.z;
		}
		else {

			val = cuDBL3(
				Lodia_shifted_irregular_xy(i, j, 0, n.x, n.y, n.z, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del),
				Lodia_shifted_irregular_xz_yz(i, 0, j, n.x, n.z, n.y, d.x*d.y*d.z, g_vals_xz, g_vals_xz_del),
				Lodia_shifted_irregular_xz_yz(j, 0, i, n.y, n.z, n.x, d.x*d.y*d.z, g_vals_yz, g_vals_yz_del)) * sign;
		}

		D12[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.x;
		D13[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.y;
		D23[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Irregular_Lodia(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

		CalcTens2D_zShifted_Irregular_Lodia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D12, D13, D23,
			g_vals_xy, g_vals_xz, g_vals_yz,
			g_vals_xy_del, g_vals_xz_del, g_vals_yz_del,
			n, N, s, d, shift, sign, asymptotic_distance,
			demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

		CalcTens2D_Shifted_Irregular_Lodia_kernel <<< ((2 * n.y - 1) * (2 * n.x - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D12, D13, D23,
			g_vals_xy, g_vals_xz, g_vals_yz,
			g_vals_xy_del, g_vals_xz_del, g_vals_yz_del,
			n, N, s, d, shift, sign, asymptotic_distance,
			demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
	}
}

#endif