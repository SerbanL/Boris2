#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include "DemagTFunc_fgCUDA.cuh"

#include <cuda_runtime.h>

//REGULAR

//--------------------- fill_f_vals

//--------------------- CUBE CELL

__global__ void fill_f_vals_xx_cubecell_kernel(cuVEC<double>& f_vals_xx, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = f_vals_xx.n.x - 2;
	int ny_dist = f_vals_xx.n.y - 2;
	int nz_dist = f_vals_xx.n.z - 2;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - 1;

		f_vals_xx[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z);
	}
}

__global__ void fill_f_vals_yy_zz_cubecell_kernel(cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = f_vals_xx.n.x - 2;
	int ny_dist = f_vals_xx.n.y - 2;
	int nz_dist = f_vals_xx.n.z - 2;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - 1;

		if (j <= nx_dist && i <= ny_dist) {

			//we already have this value : get it
			f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f_vals_xx[cuINT3(j + 1, i + 1, k + 1)];
		}
		else f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z);

		if (k <= nx_dist && i <= nz_dist) {

			//we already have this value : get it
			f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f_vals_xx[cuINT3(k + 1, j + 1, i + 1)];
		}
		else f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f(k * hRatios.z, j * hRatios.y, i * hRatios.x);
	}
}

//---------------------	SQUARE CELL

__global__ void fill_f_vals_xx_zz_squarecell_kernel(cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_zz, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = f_vals_xx.n.x - 2;
	int ny_dist = f_vals_xx.n.y - 2;
	int nz_dist = f_vals_xx.n.z - 2;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - 1;

		f_vals_xx[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z);
		f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f(k * hRatios.z, j * hRatios.y, i * hRatios.x);
	}
}

__global__ void fill_f_vals_yy_squarecell_kernel(cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = f_vals_xx.n.x - 2;
	int ny_dist = f_vals_xx.n.y - 2;
	int nz_dist = f_vals_xx.n.z - 2;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - 1;

		if (j <= nx_dist && i <= ny_dist) {

			//we already have this value : get it
			f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f_vals_xx[cuINT3(j + 1, i + 1, k + 1)];
		}
		else f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z);
	}
}

//---------------------	GENERAL CELL

__global__ void fill_f_vals_xx_yy_zz_generalcell_kernel(cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = f_vals_xx.n.x - 2;
	int ny_dist = f_vals_xx.n.y - 2;
	int nz_dist = f_vals_xx.n.z - 2;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - 1;

		f_vals_xx[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z);
		f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z);
		f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f(k * hRatios.z, j * hRatios.y, i * hRatios.x);
	}
}

bool DemagTFuncCUDA::fill_f_vals(cuINT3 n, cuDBL3 hRatios, int asymptotic_distance)
{
	int nx_dist = n.x, ny_dist = n.y, nz_dist = n.z;

	if (asymptotic_distance > 0) {

		nx_dist = minimum(n.x, asymptotic_distance);
		ny_dist = minimum(n.y, asymptotic_distance);
		nz_dist = minimum(n.z, asymptotic_distance);
	}

	if (!f_vals_xx()->assign(cuSZ3(nx_dist + 2, ny_dist + 2, nz_dist + 2), 0.0)) return false;
	if (!f_vals_yy()->assign(cuSZ3(ny_dist + 2, nx_dist + 2, nz_dist + 2), 0.0)) return false;
	if (!f_vals_zz()->assign(cuSZ3(nz_dist + 2, ny_dist + 2, nx_dist + 2), 0.0)) return false;

	if (cuIsE(hRatios.x, hRatios.y) && cuIsE(hRatios.y, hRatios.z)) {

		//cubic cell : re-use calculated f values as much as possible

		//first calculate xx values
		fill_f_vals_xx_cubecell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, hRatios);

		//now calculate yy and zz values only if we don't have them already from previous step
		fill_f_vals_yy_zz_cubecell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, f_vals_yy, f_vals_zz, hRatios);
	}
	else if (cuIsE(hRatios.x, hRatios.y)) {

		//square cell : re-use calculated f values as much as possible

		//first calculate xx and zz values
		fill_f_vals_xx_zz_squarecell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, f_vals_zz, hRatios);

		//now calculate yy values only if we don't have them already from previous step
		fill_f_vals_yy_squarecell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, f_vals_yy, hRatios);
	}
	else {

		//general case : calculate xx, yy, and zz values separately
		fill_f_vals_xx_yy_zz_generalcell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, f_vals_yy, f_vals_zz, hRatios);
	}

	return true;
}

//--------------------- fill_g2D_vals

__global__ void fill_g2D_vals_xy_kernel(cuVEC<double>& g_vals_xy, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = g_vals_xy.n.x - 2;
	int ny_dist = g_vals_xy.n.y - 2;
	int nz_dist = g_vals_xy.n.z - 2;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - 1;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - 1;
		int k = (idx / (g_vals_xy.n.x * g_vals_xy.n.y)) - 1;

		g_vals_xy[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z);
	}
}

bool DemagTFuncCUDA::fill_g2D_vals(cuINT3 n, cuDBL3 hRatios, int asymptotic_distance)
{
	int nx_dist = n.x, ny_dist = n.y, nz_dist = n.z;

	if (asymptotic_distance > 0) {

		nx_dist = minimum(n.x, asymptotic_distance);
		ny_dist = minimum(n.y, asymptotic_distance);
		nz_dist = minimum(n.z, asymptotic_distance);
	}

	if (!g_vals_xy()->assign(cuSZ3(nx_dist + 2, ny_dist + 2, nz_dist + 2), 0.0)) return false;

	fill_g2D_vals_xy_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (g_vals_xy, hRatios);
	
	return true;
}

//--------------------- fill_g_vals

//--------------------- CUBE CELL

__global__ void fill_g_vals_xy_cubecell_kernel(cuVEC<double>& g_vals_xy, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = g_vals_xy.n.x - 2;
	int ny_dist = g_vals_xy.n.y - 2;
	int nz_dist = g_vals_xy.n.z - 2;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - 1;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - 1;
		int k = (idx / (g_vals_xy.n.x * g_vals_xy.n.y)) - 1;

		g_vals_xy[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z);
	}
}

__global__ void fill_g_vals_xz_yz_cubecell_kernel(cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = g_vals_xy.n.x - 2;
	int ny_dist = g_vals_xy.n.y - 2;
	int nz_dist = g_vals_xy.n.z - 2;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - 1;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - 1;
		int k = (idx / (g_vals_xy.n.x * g_vals_xy.n.y)) - 1;

		if (k <= ny_dist && j <= nz_dist) {

			//we already have this value : get it
			g_vals_xz[(i + 1) + (k + 1) * (nx_dist + 2) + (j + 1) * (nx_dist + 2) * (nz_dist + 2)] = g_vals_xy[cuINT3(i + 1, k + 1, j + 1)];
		}
		else g_vals_xz[(i + 1) + (k + 1) * (nx_dist + 2) + (j + 1) * (nx_dist + 2) * (nz_dist + 2)] = g(i * hRatios.x, k * hRatios.z, j * hRatios.y);


		if (j <= nx_dist && k <= ny_dist && i <= nz_dist) {

			//we already have this value : get it
			g_vals_yz[(j + 1) + (k + 1) * (ny_dist + 2) + (i + 1) * (ny_dist + 2) * (nz_dist + 2)] = g_vals_xy[cuINT3(j + 1, k + 1, i + 1)];
		}
		else g_vals_yz[(j + 1) + (k + 1) * (ny_dist + 2) + (i + 1) * (ny_dist + 2) * (nz_dist + 2)] = g(j * hRatios.y, k * hRatios.z, i * hRatios.x);
	}
}

//--------------------- GENERAL CELL

__global__ void fill_g_vals_xy_xz_yz_generalcell_kernel(cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz, cuDBL3 hRatios)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int nx_dist = g_vals_xy.n.x - 2;
	int ny_dist = g_vals_xy.n.y - 2;
	int nz_dist = g_vals_xy.n.z - 2;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - 1;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - 1;
		int k = (idx / (g_vals_xy.n.x * g_vals_xy.n.y)) - 1;

		g_vals_xy[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z);
		g_vals_xz[(i + 1) + (k + 1) * (nx_dist + 2) + (j + 1) * (nx_dist + 2) * (nz_dist + 2)] = g(i * hRatios.x, k * hRatios.z, j * hRatios.y);
		g_vals_yz[(j + 1) + (k + 1) * (ny_dist + 2) + (i + 1) * (ny_dist + 2) * (nz_dist + 2)] = g(j * hRatios.y, k * hRatios.z, i * hRatios.x);
	}
}

bool DemagTFuncCUDA::fill_g_vals(cuINT3 n, cuDBL3 hRatios, int asymptotic_distance)
{
	int nx_dist = n.x, ny_dist = n.y, nz_dist = n.z;

	if (asymptotic_distance > 0) {

		nx_dist = minimum(n.x, asymptotic_distance);
		ny_dist = minimum(n.y, asymptotic_distance);
		nz_dist = minimum(n.z, asymptotic_distance);
	}

	if (!g_vals_xy()->assign(cuSZ3(nx_dist + 2, ny_dist + 2, nz_dist + 2), 0.0)) return false;
	if (!g_vals_xz()->assign(cuSZ3(nx_dist + 2, nz_dist + 2, ny_dist + 2), 0.0)) return false;
	if (!g_vals_yz()->assign(cuSZ3(ny_dist + 2, nz_dist + 2, nx_dist + 2), 0.0)) return false;

	if (cuIsE(hRatios.x, hRatios.y) && cuIsE(hRatios.y, hRatios.z)) {

		//cubic cell : re-use calculated f values as much as possible
		
		//first calculate xy values
		fill_g_vals_xy_cubecell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (g_vals_xy, hRatios);

		//now calculate xy and xz values only if we don't have them from previous step
		fill_g_vals_xz_yz_cubecell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (g_vals_xy, g_vals_xz, g_vals_yz, hRatios);
	}
	else {

		//general case : calculate xy, xz, and yz values separately
		fill_g_vals_xy_xz_yz_generalcell_kernel <<< ((nx_dist + 2) * (ny_dist + 2) * (nz_dist + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (g_vals_xy, g_vals_xz, g_vals_yz, hRatios);
	}

	return true;
}

//SHIFTED

__global__ void fill_f_vals_shifted_kernel(cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz, cuINT3 n, cuDBL3 hRatios, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - n.x;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - n.y;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - n.z;

		f_vals_xx[(i + n.x) + (j + n.y) * (2 * n.x + 1) + (k + n.z) * (2 * n.x + 1) * (2 * n.y + 1)] = f(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z);
		f_vals_yy[(j + n.y) + (i + n.x) * (2 * n.y + 1) + (k + n.z) * (2 * n.y + 1) * (2 * n.x + 1)] = f(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z);
		f_vals_zz[(k + n.z) + (j + n.y) * (2 * n.z + 1) + (i + n.x) * (2 * n.z + 1) * (2 * n.y + 1)] = f(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x);
	}
}

//shifted versions for stray field
bool DemagTFuncCUDA::fill_f_vals_shifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift)
{
	//shifted versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!f_vals_xx()->assign(cuSZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy()->assign(cuSZ3(2 * n.y + 1, 2 * n.x + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz()->assign(cuSZ3(2 * n.z + 1, 2 * n.y + 1, 2 * n.x + 1), 0.0)) return false;

	fill_f_vals_shifted_kernel <<< ((2 * n.x + 1) * (2 * n.y + 1) * (2 * n.z + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, f_vals_yy, f_vals_zz, n, hRatios, shift);

	return true;
}

__global__ void fill_g_vals_shifted_kernel(cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz, cuINT3 n, cuDBL3 hRatios, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - n.x;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - n.y;
		int k = (idx / (g_vals_xy.n.x * g_vals_xy.n.y)) - n.z;

		g_vals_xy[(i + n.x) + (j + n.y) * (2 * n.x + 1) + (k + n.z) * (2 * n.x + 1) * (2 * n.y + 1)] = g(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z);
		g_vals_xz[(i + n.x) + (k + n.z) * (2 * n.x + 1) + (j + n.y) * (2 * n.x + 1) * (2 * n.z + 1)] = g(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y);
		g_vals_yz[(j + n.y) + (k + n.z) * (2 * n.y + 1) + (i + n.x) * (2 * n.y + 1) * (2 * n.z + 1)] = g(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x);
	}
}

bool DemagTFuncCUDA::fill_g_vals_shifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift)
{
	//shifted versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!g_vals_xy()->assign(cuSZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz()->assign(cuSZ3(2 * n.x + 1, 2 * n.z + 1, 2 * n.y + 1), 0.0)) return false;
	if (!g_vals_yz()->assign(cuSZ3(2 * n.y + 1, 2 * n.z + 1, 2 * n.x + 1), 0.0)) return false;

	fill_g_vals_shifted_kernel <<< ((2 * n.x + 1) * (2 * n.y + 1) * (2 * n.z + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (g_vals_xy, g_vals_xz, g_vals_yz, n, hRatios, shift);

	return true;
}

//Z-SHIFTED

__global__ void fill_f_vals_zshifted_kernel(cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz, cuINT3 n, cuDBL3 hRatios, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;
		int k = (idx / (f_vals_xx.n.x * f_vals_xx.n.y)) - n.z;

		f_vals_xx[cuINT3(i + 1, j + 1, k + n.z)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z);
		f_vals_yy[cuINT3(j + 1, i + 1, k + n.z)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z + shift.z);
		f_vals_zz[cuINT3(k + n.z, j + 1, i + 1)] = f(k * hRatios.z + shift.z, j * hRatios.y, i * hRatios.x);
	}
}

//shifted versions for stray field -> in particular for z shift only
bool DemagTFuncCUDA::fill_f_vals_zshifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift, int asymptotic_distance)
{
	//for z shift only can use xy plane symmetries

	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!f_vals_xx()->assign(cuSZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy()->assign(cuSZ3(n.y + 2, n.x + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz()->assign(cuSZ3(2 * n.z + 1, n.y + 2, n.x + 2), 0.0)) return false;

	fill_f_vals_zshifted_kernel <<< ((n.x + 2) * (n.y + 2) * (2 * n.z + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (f_vals_xx, f_vals_yy, f_vals_zz, n, hRatios, shift);

	return true;
}

__global__ void fill_g_vals_zshifted_kernel(cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz, cuINT3 n, cuDBL3 hRatios, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - 1;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - 1;
		int k = (idx / (g_vals_xy.n.x * g_vals_xy.n.y)) - n.z;

		g_vals_xy[cuINT3(i + 1, j + 1, k + n.z)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z);
		g_vals_xz[cuINT3(i + 1, k + n.z, j + 1)] = g(i * hRatios.x, k * hRatios.z + shift.z, j * hRatios.y);
		g_vals_yz[cuINT3(j + 1, k + n.z, i + 1)] = g(j * hRatios.y, k * hRatios.z + shift.z, i * hRatios.x);
	}
}

bool DemagTFuncCUDA::fill_g_vals_zshifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift, int asymptotic_distance)
{
	//for z shift only can use xy plane symmetries

	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!g_vals_xy()->assign(cuSZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz()->assign(cuSZ3(n.x + 2, 2 * n.z + 1, n.y + 2), 0.0)) return false;
	if (!g_vals_yz()->assign(cuSZ3(n.y + 2, 2 * n.z + 1, n.x + 2), 0.0)) return false;

	fill_g_vals_zshifted_kernel <<< ((n.x + 2) * (n.y + 2) * (2 * n.z + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (g_vals_xy, g_vals_xz, g_vals_yz, n, hRatios, shift);

	return true;
}

//SHIFTED IRREGULAR

__global__ void fill_f_vals_shifted_irregular_kernel(
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz, 
	cuVEC<double>& f_vals_xx_del, cuVEC<double>& f_vals_yy_del, cuVEC<double>& f_vals_zz_del,
	cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < f_vals_xx.linear_size()) {

		int i = (idx % f_vals_xx.n.x) - n.x;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - n.y;

		double del = h_src.z - h_dst.z;

		//xx, yy : +dz, -sz, -del
		//zz : +sx, -dx, +del

		//k = -1
		f_vals_xx[(i + n.x) + (j + n.y) * (2 * n.x + 1)] = f(i * h_dst.x + shift.x, j * h_dst.y + shift.y, -h_src.z + shift.z);
		f_vals_yy[(j + n.y) + (i + n.x) * (2 * n.y + 1)] = f(j * h_dst.y + shift.y, i * h_dst.x + shift.x, -h_src.z + shift.z);
		f_vals_zz[(j + n.y) * (2 * n.z + 1) + (i + n.x) * (2 * n.z + 1) * (2 * n.y + 1)] = f(-h_dst.z + shift.z, j * h_dst.y + shift.y, i * h_dst.x + shift.x);

		//k = 0
		f_vals_xx[(i + n.x) + (j + n.y) * (2 * n.x + 1) + (2 * n.x + 1) * (2 * n.y + 1)] = f(i * h_dst.x + shift.x, j * h_dst.y + shift.y, shift.z);
		f_vals_yy[(j + n.y) + (i + n.x) * (2 * n.y + 1) + (2 * n.y + 1) * (2 * n.x + 1)] = f(j * h_dst.y + shift.y, i * h_dst.x + shift.x, shift.z);
		f_vals_zz[1 + (j + n.y) * (2 * n.z + 1) + (i + n.x) * (2 * n.z + 1) * (2 * n.y + 1)] = f(shift.z, j * h_dst.y + shift.y, i * h_dst.x + shift.x);

		//k = 1
		f_vals_xx[(i + n.x) + (j + n.y) * (2 * n.x + 1) + 2 * (2 * n.x + 1) * (2 * n.y + 1)] = f(i * h_dst.x + shift.x, j * h_dst.y + shift.y, h_dst.z + shift.z);
		f_vals_yy[(j + n.y) + (i + n.x) * (2 * n.y + 1) + 2 * (2 * n.y + 1) * (2 * n.x + 1)] = f(j * h_dst.y + shift.y, i * h_dst.x + shift.x, h_dst.z + shift.z);
		f_vals_zz[2 + (j + n.y) * (2 * n.z + 1) + (i + n.x) * (2 * n.z + 1) * (2 * n.y + 1)] = f(h_src.z + shift.z, j * h_dst.y + shift.y, i * h_dst.x + shift.x);

		//del versions
		f_vals_xx_del[(i + n.x) + (j + n.y) * (2 * n.x + 1)] = f(i * h_dst.x + shift.x, j * h_dst.y + shift.y, -del + shift.z);
		f_vals_yy_del[(j + n.y) + (i + n.x) * (2 * n.y + 1)] = f(j * h_dst.y + shift.y, i * h_dst.x + shift.x, -del + shift.z);
		f_vals_zz_del[(j + n.y) + (i + n.x) * (2 * n.y + 1)] = f(+del + shift.z, j * h_dst.y + shift.y, i * h_dst.x + shift.x);
	}
}

//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz)
bool DemagTFuncCUDA::fill_f_vals_shifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift)
{
	//shifted (and irregular) versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!f_vals_xx()->assign(cuSZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy()->assign(cuSZ3(2 * n.y + 1, 2 * n.x + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz()->assign(cuSZ3(2 * n.z + 1, 2 * n.y + 1, 2 * n.x + 1), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!f_vals_xx_del()->assign(cuSZ3(2 * n.x + 1, 2 * n.y + 1, 1), 0.0)) return false;
	if (!f_vals_yy_del()->assign(cuSZ3(2 * n.y + 1, 2 * n.x + 1, 1), 0.0)) return false;
	if (!f_vals_zz_del()->assign(cuSZ3(1, 2 * n.y + 1, 2 * n.x + 1), 0.0)) return false;

	fill_f_vals_shifted_irregular_kernel <<< ((2 * n.x + 1) * (2 * n.y + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
		(f_vals_xx, f_vals_yy, f_vals_zz,
		f_vals_xx_del, f_vals_yy_del, f_vals_zz_del,
		n, h_src, h_dst, shift);

	return true;
}

__global__ void fill_g_vals_shifted_irregular_kernel(
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuVEC<double>& g_vals_xy_del, cuVEC<double>& g_vals_xz_del, cuVEC<double>& g_vals_yz_del,
	cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < g_vals_xy.linear_size()) {

		int i = (idx % g_vals_xy.n.x) - n.x;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - n.y;

		double del = h_src.z - h_dst.z;

		//xy : +sz, -dz, +del
		//xz, yz : +sy, -dy, +del

		//k = -1
		g_vals_xy[(i + n.x) + (j + n.y) * (2 * n.x + 1)] = g(i * h_dst.x + shift.x, j * h_dst.y + shift.y, -h_dst.z + shift.z);
		g_vals_xz[(i + n.x) + (j + n.y) * (2 * n.x + 1) * (2 * n.z + 1)] = g(i * h_dst.x + shift.x, -h_dst.z + shift.z, j * h_dst.y + shift.y);
		g_vals_yz[(j + n.y) + (i + n.x) * (2 * n.y + 1) * (2 * n.z + 1)] = g(j * h_dst.y + shift.y, -h_dst.z + shift.z, i * h_dst.x + shift.x);

		//k = 0
		g_vals_xy[(i + n.x) + (j + n.y) * (2 * n.x + 1) + (2 * n.x + 1) * (2 * n.y + 1)] = g(i * h_dst.x + shift.x, j * h_dst.y + shift.y, shift.z);
		g_vals_xz[(i + n.x) + (2 * n.x + 1) + (j + n.y) * (2 * n.x + 1) * (2 * n.z + 1)] = g(i * h_dst.x + shift.x, shift.z, j * h_dst.y + shift.y);
		g_vals_yz[(j + n.y) + (2 * n.y + 1) + (i + n.x) * (2 * n.y + 1) * (2 * n.z + 1)] = g(j * h_dst.y + shift.y, shift.z, i * h_dst.x + shift.x);

		//k = 1
		g_vals_xy[(i + n.x) + (j + n.y) * (2 * n.x + 1) + 2 * (2 * n.x + 1) * (2 * n.y + 1)] = g(i * h_dst.x + shift.x, j * h_dst.y + shift.y, h_src.z + shift.z);
		g_vals_xz[(i + n.x) + 2 * (2 * n.x + 1) + (j + n.y) * (2 * n.x + 1) * (2 * n.z + 1)] = g(i * h_dst.x + shift.x, h_src.z + shift.z, j * h_dst.y + shift.y);
		g_vals_yz[(j + n.y) + 2 * (2 * n.y + 1) + (i + n.x) * (2 * n.y + 1) * (2 * n.z + 1)] = g(j * h_dst.y + shift.y, h_src.z + shift.z, i * h_dst.x + shift.x);

		//del versions
		g_vals_xy_del[(i + n.x) + (j + n.y) * (2 * n.x + 1)] = g(i * h_dst.x + shift.x, j * h_dst.y + shift.y, del + shift.z);
		g_vals_xz_del[(i + n.x) + (j + n.y) * (2 * n.x + 1)] = g(i * h_dst.x + shift.x, del + shift.z, j * h_dst.y + shift.y);
		g_vals_yz_del[(j + n.y) + (i + n.x) * (2 * n.y + 1)] = g(j * h_dst.y + shift.y, del + shift.z, i * h_dst.x + shift.x);
	}
}

//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz)
bool DemagTFuncCUDA::fill_g_vals_shifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift)
{
	//shifted (and irregular) versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!g_vals_xy()->assign(cuSZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz()->assign(cuSZ3(2 * n.x + 1, 2 * n.z + 1, 2 * n.y + 1), 0.0)) return false;
	if (!g_vals_yz()->assign(cuSZ3(2 * n.y + 1, 2 * n.z + 1, 2 * n.x + 1), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!g_vals_xy_del()->assign(cuSZ3(2 * n.x + 1, 2 * n.y + 1, 1), 0.0)) return false;
	if (!g_vals_xz_del()->assign(cuSZ3(2 * n.x + 1, 1, 2 * n.y + 1), 0.0)) return false;
	if (!g_vals_yz_del()->assign(cuSZ3(2 * n.y + 1, 1, 2 * n.x + 1), 0.0)) return false;

	fill_g_vals_shifted_irregular_kernel <<< ((2 * n.x + 1) * (2 * n.y + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(g_vals_xy, g_vals_xz, g_vals_yz,
		g_vals_xy_del, g_vals_xz_del, g_vals_yz_del,
		n, h_src, h_dst, shift);

	return true;
}

//Z-SHIFTED IRREGULAR

__global__ void fill_f_vals_zshifted_irregular_kernel(
	cuVEC<double>& f_vals_xx, cuVEC<double>& f_vals_yy, cuVEC<double>& f_vals_zz,
	cuVEC<double>& f_vals_xx_del, cuVEC<double>& f_vals_yy_del, cuVEC<double>& f_vals_zz_del,
	cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (n.x + 2) * (n.y + 2)) {

		int i = (idx % f_vals_xx.n.x) - 1;
		int j = ((idx / f_vals_xx.n.x) % f_vals_xx.n.y) - 1;

		double del = h_src.z - h_dst.z;

		//xx, yy : +dz, -sz, -del
		//zz : +sx, -dx, +del

		//k = -1

		f_vals_xx[cuINT3(i + 1, j + 1, 0)] = f(i * h_dst.x, j * h_dst.y, -h_src.z + shift.z);
		f_vals_yy[cuINT3(j + 1, i + 1, 0)] = f(j * h_dst.y, i * h_dst.x, -h_src.z + shift.z);
		f_vals_zz[cuINT3(0, j + 1, i + 1)] = f(-h_dst.z + shift.z, j * h_dst.y, i * h_dst.x);

		//k = 0
		f_vals_xx[cuINT3(i + 1, j + 1, 1)] = f(i * h_dst.x, j * h_dst.y, shift.z);
		f_vals_yy[cuINT3(j + 1, i + 1, 1)] = f(j * h_dst.y, i * h_dst.x, shift.z);
		f_vals_zz[cuINT3(1, j + 1, i + 1)] = f(shift.z, j * h_dst.y, i * h_dst.x);

		//k = 1
		f_vals_xx[cuINT3(i + 1, j + 1, 2)] = f(i * h_dst.x, j * h_dst.y, h_dst.z + shift.z);
		f_vals_yy[cuINT3(j + 1, i + 1, 2)] = f(j * h_dst.y, i * h_dst.x, h_dst.z + shift.z);
		f_vals_zz[cuINT3(2, j + 1, i + 1)] = f(h_src.z + shift.z, j * h_dst.y, i * h_dst.x);

		//del versions
		f_vals_xx_del[cuINT3(i + 1, j + 1, 0)] = f(i * h_dst.x, j * h_dst.y, -del + shift.z);
		f_vals_yy_del[cuINT3(j + 1, i + 1, 0)] = f(j * h_dst.y, i * h_dst.x, -del + shift.z);
		f_vals_zz_del[cuINT3(0, j + 1, i + 1)] = f(+del + shift.z, j * h_dst.y, i * h_dst.x);
	}
}

//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz) -> in particular for z shift only
bool DemagTFuncCUDA::fill_f_vals_zshifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift, int asymptotic_distance)
{
	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!f_vals_xx()->assign(cuSZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy()->assign(cuSZ3(n.y + 2, n.x + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz()->assign(cuSZ3(2 * n.z + 1, n.y + 2, n.x + 2), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!f_vals_xx_del()->assign(cuSZ3(n.x + 2, n.y + 2, 1), 0.0)) return false;
	if (!f_vals_yy_del()->assign(cuSZ3(n.y + 2, n.x + 2, 1), 0.0)) return false;
	if (!f_vals_zz_del()->assign(cuSZ3(1, n.y + 2, n.x + 2), 0.0)) return false;

	fill_f_vals_zshifted_irregular_kernel <<< ((n.x + 2) * (n.y + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(f_vals_xx, f_vals_yy, f_vals_zz,
		f_vals_xx_del, f_vals_yy_del, f_vals_zz_del,
		n, h_src, h_dst, shift);

	return true;
}

__global__ void fill_g_vals_zshifted_irregular_kernel(
	cuVEC<double>& g_vals_xy, cuVEC<double>& g_vals_xz, cuVEC<double>& g_vals_yz,
	cuVEC<double>& g_vals_xy_del, cuVEC<double>& g_vals_xz_del, cuVEC<double>& g_vals_yz_del,
	cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < (n.x + 2) * (n.y + 2)) {

		int i = (idx % g_vals_xy.n.x) - 1;
		int j = ((idx / g_vals_xy.n.x) % g_vals_xy.n.y) - 1;

		double del = h_src.z - h_dst.z;

		//xy : +sz, -dz, +del
			//xz, yz : +sy, -dy, +del

			//k = -1
		g_vals_xy[cuINT3(i + 1, j + 1, 0)] = g(i * h_dst.x, j * h_dst.y, -h_dst.z + shift.z);
		g_vals_xz[cuINT3(i + 1, 0, j + 1)] = g(i * h_dst.x, -h_dst.z + shift.z, j * h_dst.y);
		g_vals_yz[cuINT3(j + 1, 0, i + 1)] = g(j * h_dst.y, -h_dst.z + shift.z, i * h_dst.x);

		//k = 0
		g_vals_xy[cuINT3(i + 1, j + 1, 1)] = g(i * h_dst.x, j * h_dst.y, shift.z);
		g_vals_xz[cuINT3(i + 1, 1, j + 1)] = g(i * h_dst.x, shift.z, j * h_dst.y);
		g_vals_yz[cuINT3(j + 1, 1, i + 1)] = g(j * h_dst.y, shift.z, i * h_dst.x);

		//k = 1
		g_vals_xy[cuINT3(i + 1, j + 1, 2)] = g(i * h_dst.x, j * h_dst.y, h_src.z + shift.z);
		g_vals_xz[cuINT3(i + 1, 2, j + 1)] = g(i * h_dst.x, h_src.z + shift.z, j * h_dst.y);
		g_vals_yz[cuINT3(j + 1, 2, i + 1)] = g(j * h_dst.y, h_src.z + shift.z, i * h_dst.x);

		//del versions
		g_vals_xy_del[cuINT3(i + 1, j + 1, 0)] = g(i * h_dst.x, j * h_dst.y, del + shift.z);
		g_vals_xz_del[cuINT3(i + 1, 0, j + 1)] = g(i * h_dst.x, del + shift.z, j * h_dst.y);
		g_vals_yz_del[cuINT3(j + 1, 0, i + 1)] = g(j * h_dst.y, del + shift.z, i * h_dst.x);
	}
}

bool DemagTFuncCUDA::fill_g_vals_zshifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift, int asymptotic_distance)
{
	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!g_vals_xy()->assign(cuSZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz()->assign(cuSZ3(n.x + 2, 2 * n.z + 1, n.y + 2), 0.0)) return false;
	if (!g_vals_yz()->assign(cuSZ3(n.y + 2, 2 * n.z + 1, n.x + 2), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!g_vals_xy_del()->assign(cuSZ3(n.x + 2, n.y + 2, 1), 0.0)) return false;
	if (!g_vals_xz_del()->assign(cuSZ3(n.x + 2, 1, n.y + 2), 0.0)) return false;
	if (!g_vals_yz_del()->assign(cuSZ3(n.y + 2, 1, n.x + 2), 0.0)) return false;

	fill_g_vals_zshifted_irregular_kernel <<< ((n.x + 2) * (n.y + 2) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(g_vals_xy, g_vals_xz, g_vals_yz,
		g_vals_xy_del, g_vals_xz_del, g_vals_yz_del,
		n, h_src, h_dst, shift);

	return true;
}

#endif
