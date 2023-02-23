#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>
#include "DemagTFunc_LCUDA.cuh"

//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

//REGULAR VERSIONS FOR INTERNAL FIELD

//2D

__global__ void CalcTens2D_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;

		cuDBL3 val;

		//apply asymptotic equations?
		//Must include the separate i, j checks otherwise the i*i + j*j expression can overflow the 4 byte integer range
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || i * i + j * j >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, 0) * sign,
				demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, 0) * sign,
				demagAsymptoticDiag_zz.AsymptoticLdia(0, j * hRatios.y, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Ldia(i, j, 0, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
				Ldia(j, i, 0, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
				Ldia(0, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
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

//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
void DemagTFuncCUDA::CalcTens2D_Ldia(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance)
{
	CalcTens2D_Ldia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(D11, D22, D33,
		f_vals_xx, f_vals_yy, f_vals_zz,
		n, N, hRatios, sign, asymptotic_distance,
		demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
}

__global__ void CalcTens2D_Lodia_kernel(
	double* D12,
	cuVEC<double>& g_vals_xy,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;

		double val = 0.0;

		//apply asymptotic equations?
		//Must include the separate i, j checks otherwise the i*i + j*j expression can overflow the 4 byte integer range
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || i * i + j * j >= asymptotic_distance * asymptotic_distance)) {

			val = demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, 0) * sign;
		}
		else {

			val = Lodia(i, j, 0, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign;
		}

		D12[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = val;
		D12[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x] = -1 * val;
		D12[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = -1 * val;
		D12[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x] = val;
	}
}

void DemagTFuncCUDA::CalcTens2D_Lodia(
	cu_arr<double>& D12,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance)
{
	CalcTens2D_Lodia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(D12,
		g_vals_xy,
		n, N, hRatios, sign, asymptotic_distance,
		demagAsymptoticOffDiag_xy);
}

//3D

__global__ void CalcTens3D_Ldia_kernel(
	double* D11, double* D22, double* D33,
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticDiagCUDA& demagAsymptoticDiag_xx, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_yy, DemagAsymptoticDiagCUDA& demagAsymptoticDiag_zz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		cuDBL3 val;

		//apply asymptotic equations?
		//Must include the separate i, j, k checks otherwise the i*i + j*j + k*k expression can overflow the 4 byte integer range
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

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

		D11[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.z;

		D11[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.x;
		D22[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.y;
		D33[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.z;
	}
}

__global__ void CalcTens3D_Lodia_kernel(
	double* D12, double* D13, double* D23,
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
	DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xy, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_xz, DemagAsymptoticOffDiagCUDA& demagAsymptoticOffDiag_yz)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < n.dim()) {

		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		cuDBL3 val;

		//apply asymptotic equations?
		//Must include the separate i, j, k checks otherwise the i*i + j*j + k*k expression can overflow the 4 byte integer range
		if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

			val = cuDBL3(
				demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
				demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z, j * hRatios.y) * sign,
				demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z, i * hRatios.x) * sign);
		}
		else {

			val = cuDBL3(
				Lodia(i, j, k, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
				Lodia(i, k, j, hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
				Lodia(j, k, i, hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign);
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
		
		D12[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.x;
		D13[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = -val.y;
		D23[((i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = -val.z;
		
		D12[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = val.x;
		D13[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.y;
		D23[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((k + N.z) % N.z)*N.x*N.y] = -val.z;
		
		D12[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = -val.x;
		D13[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.y;
		D23[((-i + N.x) % N.x) + ((j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = -val.z;
		
		D12[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = -val.x;
		D13[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = -val.y;
		D23[((i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.z;
		
		D12[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.x;
		D13[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.y;
		D23[((-i + N.x) % N.x) + ((-j + N.y) % N.y)*N.x + ((-k + N.z) % N.z)*N.x*N.y] = val.z;
	}
}

//3D
void DemagTFuncCUDA::CalcTens3D_Ldia(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance)
{
	CalcTens3D_Ldia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(D11, D22, D33,
		f_vals_xx, f_vals_yy, f_vals_zz,
		n, N, hRatios, sign, asymptotic_distance,
		demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz);
}

void DemagTFuncCUDA::CalcTens3D_Lodia(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance)
{
	CalcTens3D_Lodia_kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(D12, D13, D23,
		g_vals_xy, g_vals_xz, g_vals_yz,
		n, N, hRatios, sign, asymptotic_distance,
		demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz);
}

#endif