#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>
#include "DemagTFunc_LCUDA.cuh"

//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

//SHIFTED VERSIONS FOR STRAY FIELD

//2D

//z shift only, so use full xy plane symmetries
__global__ void CalcTens2D_zShifted_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;

		cuDBL3 val;

		double ks = cu_mod(shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_zshifted_xx_yy(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
				Ldia_zshifted_xx_yy(j, i, 0, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
				Ldia_zshifted_zz(0, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
		}

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
__global__ void CalcTens2D_Shifted_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (2 * n.y - 1) * (2 * n.x - 1)) {

		int i = (idx % (2 * n.x - 1)) - n.x + 1;
		int j = ((idx / (2 * n.x - 1)) % (2 * n.y - 1)) - n.y + 1;

		cuDBL3 val;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y + shift.y, i * hRatios.x + shift.x, shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_shifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
				Ldia_shifted(j, i, 0, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
				Ldia_shifted(0, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
		}

		D11[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.x;
		D22[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.y;
		D33[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Ldia(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

		CalcTens2D_zShifted_Ldia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D11, D22, D33,
				f_vals_xx, f_vals_yy, f_vals_zz,
				n, N, hRatios, shift, sign, asymptotic_distance,
				demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

		CalcTens2D_Shifted_Ldia_kernel <<< ((2 * n.y - 1) * (2 * n.x - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D11, D22, D33,
			f_vals_xx, f_vals_yy, f_vals_zz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
	}
}

//z shift only, so use full xy plane symmetries
__global__ void CalcTens2D_zShifted_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;

		cuDBL3 val;

		double ks = cu_mod(shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, shift.z, j * hRatios.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, shift.z, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Lodia_xy_zshifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
				Lodia_xz_yz_zshifted(i, 0, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
				Lodia_xz_yz_zshifted(j, 0, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;
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
__global__ void CalcTens2D_Shifted_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (2 * n.y - 1) * (2 * n.x - 1)) {

		int i = (idx % (2 * n.x - 1)) - n.x + 1;
		int j = ((idx / (2 * n.x - 1)) % (2 * n.y - 1)) - n.y + 1;

		cuDBL3 val;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x + shift.x, shift.z, j * hRatios.y + shift.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y + shift.y, shift.z, i * hRatios.x + shift.x) * sign);
		}
		else {

			val = cuDBL3(
				Lodia_shifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
				Lodia_shifted(i, 0, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
				Lodia_shifted(j, 0, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;
		}

		D12[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.x;
		D13[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.y;
		D23[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x] = val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Shifted_Lodia(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

		CalcTens2D_zShifted_Lodia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D12, D13, D23,
			g_vals_xy, g_vals_xz, g_vals_yz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

		CalcTens2D_Shifted_Lodia_kernel <<< ((2 * n.y - 1) * (2 * n.x - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D12, D13, D23,
			g_vals_xy, g_vals_xz, g_vals_yz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
	}
}

//3D

//z shift only, so use full xy plane symmetries
__global__ void CalcTens3D_zShifted_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.x * n.y * (2 * n.z - 1)) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y) - n.z + 1;

		cuDBL3 val;

		double ks = cu_mod(k + shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z + shift.z) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia_zshifted_xx_yy(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
				Ldia_zshifted_xx_yy(j, i, k, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
				Ldia_zshifted_zz(k, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
		}

		D11[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;
	}
}

//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)
__global__ void CalcTens3D_Shifted_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (2 * n.y - 1) * (2 * n.x - 1) * (2 * n.z - 1)) {

		int i = (idx % (2 * n.x - 1)) - n.x + 1;
		int j = ((idx / (2 * n.x - 1)) % (2 * n.y - 1)) - n.y + 1;
		int k = idx / ((2 * n.x - 1) * (2 * n.y - 1)) - n.z + 1;

		cuDBL3 val;

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
				Ldia_shifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
				Ldia_shifted(j, i, k, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
				Ldia_shifted(k, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
		}

		D11[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens3D_Shifted_Ldia(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

		CalcTens3D_zShifted_Ldia_kernel <<< (n.x * n.y * (2 * n.z - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D11, D22, D33,
			f_vals_xx, f_vals_yy, f_vals_zz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

		CalcTens3D_Shifted_Ldia_kernel <<< ((2 * n.y - 1) * (2 * n.x - 1) * (2 * n.z - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D11, D22, D33,
			f_vals_xx, f_vals_yy, f_vals_zz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
	}
}

//z shift only, so use full xy plane symmetries
__global__ void CalcTens3D_zShifted_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.x * n.y * (2 * n.z - 1)) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y) - n.z + 1;

		cuDBL3 val;

		double ks = cu_mod(k + shift.z / hRatios.z);

		//apply asymptotic equations?
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance || int(cu_floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z + shift.z, j * hRatios.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z + shift.z, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Lodia_xy_zshifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
				Lodia_xz_yz_zshifted(i, k, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
				Lodia_xz_yz_zshifted(j, k, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;
		}

		D12[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D13[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D23[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D12[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.x;
		D13[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.y;
		D23[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D12[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.x;
		D13[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D23[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.z;

		D12[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D13[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.y;
		D23[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.z;
	}
}

//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)
__global__ void CalcTens3D_Shifted_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (2 * n.y - 1) * (2 * n.x - 1) * (2 * n.z - 1)) {

		int i = (idx % (2 * n.x - 1)) - n.x + 1;
		int j = ((idx / (2 * n.x - 1)) % (2 * n.y - 1)) - n.y + 1;
		int k = idx / ((2 * n.x - 1) * (2 * n.y - 1)) - n.z + 1;

		cuDBL3 val;

		double is = cu_mod(i + (shift.x / hRatios.x));
		double js = cu_mod(j + (shift.y / hRatios.y));
		double ks = cu_mod(k + (shift.z / hRatios.z));

		//apply asymptotic equations?
		if (asymptotic_distance > 0 &&
			(int(cu_floor_epsilon(is)) >= asymptotic_distance || int(cu_floor_epsilon(js)) >= asymptotic_distance || int(cu_floor_epsilon(ks)) >= asymptotic_distance ||
				int(cu_floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x) * sign);
		}
		else {

			val = cuDBL3(
				Lodia_shifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
				Lodia_shifted(i, k, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
				Lodia_shifted(j, k, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;
		}

		D12[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D13[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D23[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;
	}
}

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens3D_Shifted_Lodia(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance)
{
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

		CalcTens3D_zShifted_Lodia_kernel <<< (n.x * n.y * (2 * n.z - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D12, D13, D23,
			g_vals_xy, g_vals_xz, g_vals_yz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

		CalcTens3D_Shifted_Lodia_kernel <<< ((2 * n.y - 1) * (2 * n.x - 1) * (2 * n.z - 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(D12, D13, D23,
			g_vals_xy, g_vals_xz, g_vals_yz,
			n, N, hRatios, shift, sign, asymptotic_distance,
			demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
	}
}

#endif