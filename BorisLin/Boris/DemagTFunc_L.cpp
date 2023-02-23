#include "stdafx.h"
#include "DemagTFunc.h"

//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

//REGULAR VERSIONS FOR INTERNAL FIELD

double DemagTFunc::Ldia(int i, int j, int k, double hx, double hy, double hz, VEC<double>& f_vals)
{
	//From Newell's paper, this is what the combination of F, F1, and F2 functions reduces to in terms of the f integrals
	//Note in Newell's paper we F2(x, y, z) = f(x, y, z) - f(x, 0, z) - f(x, y, 0) - f(x, 0, 0).
	//We can drop the f(x, 0, z), f(x, y, 0), f(x, 0, 0) since they don't affect the overall sum below, as they are cancelled out

	//this function should only be called with 0 <= i < n.x, 0 <= j < n.y, 0 <= k < n.z

	//to read values from f_vals increment i, j, k by 1
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (n.x + 2), (n.y + 2), (n.z + 2)

	i++; j++; k++;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * f_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * f_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * f_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * f_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * f_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * f_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * f_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * f_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * f_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * f_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * f_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * f_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * f_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * f_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * f_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * f_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

double DemagTFunc::Lodia(int i, int j, int k, double hx, double hy, double hz, VEC<double>& g_vals)
{
	//From Newell's paper, this is what the combination of G, G1, and G2 functions reduces to in terms of the g integrals
	//Note in Newell's paper we G2(x, y, z) = g(x, y, z) - g(x, y, 0).
	//We can drop the g(x, y, 0) since it doesn't affect the overall sum below, as it's cancelled out (e.g. 8-4-4=0, -4+2+2=0, 2-1-1=0, etc)

	//this function should only be called with 0 <= i < n.x, 0 <= j < n.y, 0 <= k < n.z

	//to read values from g_vals increment i, j, k by 1
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (n.x + 2), (n.y + 2), (n.z + 2)

	i++; j++; k++;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * g_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * g_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * g_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * g_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * g_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * g_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * g_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * g_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * g_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * g_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * g_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * g_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * g_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * g_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * g_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * g_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * g_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

//SHIFTED VERSIONS FOR STRAY FIELD

double DemagTFunc::Ldia_shifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& f_vals)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * f_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * f_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * f_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * f_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * f_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * f_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * f_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * f_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * f_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * f_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * f_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * f_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * f_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * f_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * f_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * f_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

double DemagTFunc::Lodia_shifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& g_vals)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * g_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * g_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * g_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * g_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * g_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * g_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * g_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * g_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * g_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * g_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * g_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * g_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * g_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * g_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * g_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * g_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * g_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

//special versions for z shift only
double DemagTFunc::Ldia_zshifted_xx_yy(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& f_vals)
{
	i++; j++; k += nz;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * f_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * f_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * f_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * f_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * f_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * f_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * f_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * f_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * f_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * f_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * f_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * f_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * f_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * f_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * f_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * f_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

double DemagTFunc::Ldia_zshifted_zz(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& f_vals)
{
	i += nx; j++; k++;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * f_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * f_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * f_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * f_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * f_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * f_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * f_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * f_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * f_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * f_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * f_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * f_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * f_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * f_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * f_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * f_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

double DemagTFunc::Lodia_xy_zshifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& g_vals)
{
	i++; j++; k += nz;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * g_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * g_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * g_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * g_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * g_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * g_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * g_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * g_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * g_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * g_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * g_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * g_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * g_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * g_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * g_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * g_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * g_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

double DemagTFunc::Lodia_xz_yz_zshifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& g_vals)
{
	i++; j += ny; k++;

	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * g_vals[INT3(i, j, k)];

	main_sum[tn][1] = -4 * g_vals[INT3(i + 1, j, k)];
	main_sum[tn][2] = -4 * g_vals[INT3(i - 1, j, k)];
	main_sum[tn][3] = -4 * g_vals[INT3(i, j + 1, k)];
	main_sum[tn][4] = -4 * g_vals[INT3(i, j - 1, k)];
	main_sum[tn][5] = -4 * g_vals[INT3(i, j, k + 1)];
	main_sum[tn][6] = -4 * g_vals[INT3(i, j, k - 1)];

	main_sum[tn][7] = +2 * g_vals[INT3(i - 1, j - 1, k)];
	main_sum[tn][8] = +2 * g_vals[INT3(i - 1, j + 1, k)];
	main_sum[tn][9] = +2 * g_vals[INT3(i + 1, j - 1, k)];
	main_sum[tn][10] = +2 * g_vals[INT3(i + 1, j + 1, k)];

	main_sum[tn][11] = +2 * g_vals[INT3(i - 1, j, k - 1)];
	main_sum[tn][12] = +2 * g_vals[INT3(i - 1, j, k + 1)];
	main_sum[tn][13] = +2 * g_vals[INT3(i + 1, j, k - 1)];
	main_sum[tn][14] = +2 * g_vals[INT3(i + 1, j, k + 1)];

	main_sum[tn][15] = +2 * g_vals[INT3(i, j - 1, k - 1)];
	main_sum[tn][16] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	main_sum[tn][17] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	main_sum[tn][18] = +2 * g_vals[INT3(i, j + 1, k + 1)];

	main_sum[tn][19] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	main_sum[tn][20] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	main_sum[tn][21] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	main_sum[tn][22] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	main_sum[tn][23] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	main_sum[tn][24] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	main_sum[tn][25] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	main_sum[tn][26] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

//IRREGULAR AND SHIFTED VERSIONS FOR STRAY FIELD

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
double DemagTFunc::Ldia_shifted_irregular_xx_yy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * f_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * f_vals_del[INT3(i, j, 0)];

	irregular_sum[tn][2] = -2 * f_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -2 * f_vals_del[INT3(i + 1, j, 0)];
	irregular_sum[tn][4] = -2 * f_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][5] = -2 * f_vals_del[INT3(i - 1, j, 0)];
	irregular_sum[tn][6] = -2 * f_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][7] = -2 * f_vals_del[INT3(i, j + 1, 0)];
	irregular_sum[tn][8] = -2 * f_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][9] = -2 * f_vals_del[INT3(i, j - 1, 0)];
	irregular_sum[tn][10] = -4 * f_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][11] = -4 * f_vals[INT3(i, j, k - 1)];

	irregular_sum[tn][12] = +f_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +f_vals_del[INT3(i + 1, j + 1, 0)];
	irregular_sum[tn][14] = +f_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][15] = +f_vals_del[INT3(i + 1, j - 1, 0)];
	irregular_sum[tn][16] = +f_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][17] = +f_vals_del[INT3(i - 1, j + 1, 0)];
	irregular_sum[tn][18] = +f_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][19] = +f_vals_del[INT3(i - 1, j - 1, 0)];
	irregular_sum[tn][20] = +2 * f_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][21] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][22] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][23] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][24] = +2 * f_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][25] = +2 * f_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][26] = +2 * f_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][27] = +2 * f_vals[INT3(i, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
double DemagTFunc::Ldia_shifted_irregular_zz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * f_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * f_vals_del[INT3(0, j, k)];

	irregular_sum[tn][2] = -4 * f_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -4 * f_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][4] = -2 * f_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][5] = -2 * f_vals_del[INT3(0, j, k + 1)];
	irregular_sum[tn][6] = -2 * f_vals[INT3(i, j, k - 1)];
	irregular_sum[tn][7] = -2 * f_vals_del[INT3(0, j, k - 1)];
	irregular_sum[tn][8] = -2 * f_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][9] = -2 * f_vals_del[INT3(0, j + 1, k)];
	irregular_sum[tn][10] = -2 * f_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][11] = -2 * f_vals_del[INT3(0, j - 1, k)];

	irregular_sum[tn][12] = +2 * f_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +2 * f_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][14] = +2 * f_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][15] = +2 * f_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][16] = +2 * f_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][17] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][18] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][19] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][20] = +f_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][21] = +f_vals_del[INT3(0, j + 1, k + 1)];
	irregular_sum[tn][22] = +f_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][23] = +f_vals_del[INT3(0, j + 1, k - 1)];
	irregular_sum[tn][24] = +f_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][25] = +f_vals_del[INT3(0, j - 1, k + 1)];
	irregular_sum[tn][26] = +f_vals[INT3(i, j - 1, k - 1)];
	irregular_sum[tn][27] = +f_vals_del[INT3(0, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
double DemagTFunc::Lodia_shifted_irregular_xy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * g_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * g_vals_del[INT3(i, j, 0)];

	irregular_sum[tn][2] = -2 * g_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -2 * g_vals_del[INT3(i + 1, j, 0)];
	irregular_sum[tn][4] = -2 * g_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][5] = -2 * g_vals_del[INT3(i - 1, j, 0)];
	irregular_sum[tn][6] = -2 * g_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][7] = -2 * g_vals_del[INT3(i, j + 1, 0)];
	irregular_sum[tn][8] = -2 * g_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][9] = -2 * g_vals_del[INT3(i, j - 1, 0)];
	irregular_sum[tn][10] = -4 * g_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][11] = -4 * g_vals[INT3(i, j, k - 1)];

	irregular_sum[tn][12] = +g_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +g_vals_del[INT3(i + 1, j + 1, 0)];
	irregular_sum[tn][14] = +g_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][15] = +g_vals_del[INT3(i + 1, j - 1, 0)];
	irregular_sum[tn][16] = +g_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][17] = +g_vals_del[INT3(i - 1, j + 1, 0)];
	irregular_sum[tn][18] = +g_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][19] = +g_vals_del[INT3(i - 1, j - 1, 0)];
	irregular_sum[tn][20] = +2 * g_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][21] = +2 * g_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][22] = +2 * g_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][23] = +2 * g_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][24] = +2 * g_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][25] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][26] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][27] = +2 * g_vals[INT3(i, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
double DemagTFunc::Lodia_shifted_irregular_xz_yz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * g_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * g_vals_del[INT3(i, 0, k)];

	irregular_sum[tn][2] = -2 * g_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -2 * g_vals_del[INT3(i + 1, 0, k)];
	irregular_sum[tn][4] = -2 * g_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][5] = -2 * g_vals_del[INT3(i - 1, 0, k)];
	irregular_sum[tn][6] = -4 * g_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][7] = -4 * g_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][8] = -2 * g_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][9] = -2 * g_vals_del[INT3(i, 0, k + 1)];
	irregular_sum[tn][10] = -2 * g_vals[INT3(i, j, k - 1)];
	irregular_sum[tn][11] = -2 * g_vals_del[INT3(i, 0, k - 1)];

	irregular_sum[tn][12] = +2 * g_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +2 * g_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][14] = +2 * g_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][15] = +2 * g_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][16] = +g_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][17] = +g_vals_del[INT3(i + 1, 0, k + 1)];
	irregular_sum[tn][18] = +g_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][19] = +g_vals_del[INT3(i + 1, 0, k - 1)];
	irregular_sum[tn][20] = +g_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][21] = +g_vals_del[INT3(i - 1, 0, k + 1)];
	irregular_sum[tn][22] = +g_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][23] = +g_vals_del[INT3(i - 1, 0, k - 1)];
	irregular_sum[tn][24] = +2 * g_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][25] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][26] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][27] = +2 * g_vals[INT3(i, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
double DemagTFunc::Ldia_zshifted_irregular_xx_yy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i++; j++; k += nz;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * f_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * f_vals_del[INT3(i, j, 0)];

	irregular_sum[tn][2] = -2 * f_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -2 * f_vals_del[INT3(i + 1, j, 0)];
	irregular_sum[tn][4] = -2 * f_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][5] = -2 * f_vals_del[INT3(i - 1, j, 0)];
	irregular_sum[tn][6] = -2 * f_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][7] = -2 * f_vals_del[INT3(i, j + 1, 0)];
	irregular_sum[tn][8] = -2 * f_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][9] = -2 * f_vals_del[INT3(i, j - 1, 0)];
	irregular_sum[tn][10] = -4 * f_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][11] = -4 * f_vals[INT3(i, j, k - 1)];

	irregular_sum[tn][12] = +f_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +f_vals_del[INT3(i + 1, j + 1, 0)];
	irregular_sum[tn][14] = +f_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][15] = +f_vals_del[INT3(i + 1, j - 1, 0)];
	irregular_sum[tn][16] = +f_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][17] = +f_vals_del[INT3(i - 1, j + 1, 0)];
	irregular_sum[tn][18] = +f_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][19] = +f_vals_del[INT3(i - 1, j - 1, 0)];
	irregular_sum[tn][20] = +2 * f_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][21] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][22] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][23] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][24] = +2 * f_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][25] = +2 * f_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][26] = +2 * f_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][27] = +2 * f_vals[INT3(i, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
double DemagTFunc::Ldia_zshifted_irregular_zz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j++; k++;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * f_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * f_vals_del[INT3(0, j, k)];

	irregular_sum[tn][2] = -4 * f_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -4 * f_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][4] = -2 * f_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][5] = -2 * f_vals_del[INT3(0, j, k + 1)];
	irregular_sum[tn][6] = -2 * f_vals[INT3(i, j, k - 1)];
	irregular_sum[tn][7] = -2 * f_vals_del[INT3(0, j, k - 1)];
	irregular_sum[tn][8] = -2 * f_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][9] = -2 * f_vals_del[INT3(0, j + 1, k)];
	irregular_sum[tn][10] = -2 * f_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][11] = -2 * f_vals_del[INT3(0, j - 1, k)];

	irregular_sum[tn][12] = +2 * f_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +2 * f_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][14] = +2 * f_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][15] = +2 * f_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][16] = +2 * f_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][17] = +2 * f_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][18] = +2 * f_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][19] = +2 * f_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][20] = +f_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][21] = +f_vals_del[INT3(0, j + 1, k + 1)];
	irregular_sum[tn][22] = +f_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][23] = +f_vals_del[INT3(0, j + 1, k - 1)];
	irregular_sum[tn][24] = +f_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][25] = +f_vals_del[INT3(0, j - 1, k + 1)];
	irregular_sum[tn][26] = +f_vals[INT3(i, j - 1, k - 1)];
	irregular_sum[tn][27] = +f_vals_del[INT3(0, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * f_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * f_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * f_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * f_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * f_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * f_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * f_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * f_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
double DemagTFunc::Lodia_zshifted_irregular_xy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i++; j++; k += nz;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * g_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * g_vals_del[INT3(i, j, 0)];

	irregular_sum[tn][2] = -2 * g_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -2 * g_vals_del[INT3(i + 1, j, 0)];
	irregular_sum[tn][4] = -2 * g_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][5] = -2 * g_vals_del[INT3(i - 1, j, 0)];
	irregular_sum[tn][6] = -2 * g_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][7] = -2 * g_vals_del[INT3(i, j + 1, 0)];
	irregular_sum[tn][8] = -2 * g_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][9] = -2 * g_vals_del[INT3(i, j - 1, 0)];
	irregular_sum[tn][10] = -4 * g_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][11] = -4 * g_vals[INT3(i, j, k - 1)];

	irregular_sum[tn][12] = +g_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +g_vals_del[INT3(i + 1, j + 1, 0)];
	irregular_sum[tn][14] = +g_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][15] = +g_vals_del[INT3(i + 1, j - 1, 0)];
	irregular_sum[tn][16] = +g_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][17] = +g_vals_del[INT3(i - 1, j + 1, 0)];
	irregular_sum[tn][18] = +g_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][19] = +g_vals_del[INT3(i - 1, j - 1, 0)];
	irregular_sum[tn][20] = +2 * g_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][21] = +2 * g_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][22] = +2 * g_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][23] = +2 * g_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][24] = +2 * g_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][25] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][26] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][27] = +2 * g_vals[INT3(i, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
double DemagTFunc::Lodia_zshifted_irregular_xz_yz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i++; j += ny; k++;

	int tn = omp_get_thread_num();

	irregular_sum[tn][0] = +4 * g_vals[INT3(i, j, k)];
	irregular_sum[tn][1] = +4 * g_vals_del[INT3(i, 0, k)];

	irregular_sum[tn][2] = -2 * g_vals[INT3(i + 1, j, k)];
	irregular_sum[tn][3] = -2 * g_vals_del[INT3(i + 1, 0, k)];
	irregular_sum[tn][4] = -2 * g_vals[INT3(i - 1, j, k)];
	irregular_sum[tn][5] = -2 * g_vals_del[INT3(i - 1, 0, k)];
	irregular_sum[tn][6] = -4 * g_vals[INT3(i, j + 1, k)];
	irregular_sum[tn][7] = -4 * g_vals[INT3(i, j - 1, k)];
	irregular_sum[tn][8] = -2 * g_vals[INT3(i, j, k + 1)];
	irregular_sum[tn][9] = -2 * g_vals_del[INT3(i, 0, k + 1)];
	irregular_sum[tn][10] = -2 * g_vals[INT3(i, j, k - 1)];
	irregular_sum[tn][11] = -2 * g_vals_del[INT3(i, 0, k - 1)];

	irregular_sum[tn][12] = +2 * g_vals[INT3(i + 1, j + 1, k)];
	irregular_sum[tn][13] = +2 * g_vals[INT3(i + 1, j - 1, k)];
	irregular_sum[tn][14] = +2 * g_vals[INT3(i - 1, j + 1, k)];
	irregular_sum[tn][15] = +2 * g_vals[INT3(i - 1, j - 1, k)];
	irregular_sum[tn][16] = +g_vals[INT3(i + 1, j, k + 1)];
	irregular_sum[tn][17] = +g_vals_del[INT3(i + 1, 0, k + 1)];
	irregular_sum[tn][18] = +g_vals[INT3(i + 1, j, k - 1)];
	irregular_sum[tn][19] = +g_vals_del[INT3(i + 1, 0, k - 1)];
	irregular_sum[tn][20] = +g_vals[INT3(i - 1, j, k + 1)];
	irregular_sum[tn][21] = +g_vals_del[INT3(i - 1, 0, k + 1)];
	irregular_sum[tn][22] = +g_vals[INT3(i - 1, j, k - 1)];
	irregular_sum[tn][23] = +g_vals_del[INT3(i - 1, 0, k - 1)];
	irregular_sum[tn][24] = +2 * g_vals[INT3(i, j + 1, k + 1)];
	irregular_sum[tn][25] = +2 * g_vals[INT3(i, j + 1, k - 1)];
	irregular_sum[tn][26] = +2 * g_vals[INT3(i, j - 1, k + 1)];
	irregular_sum[tn][27] = +2 * g_vals[INT3(i, j - 1, k - 1)];

	irregular_sum[tn][28] = -1 * g_vals[INT3(i - 1, j - 1, k - 1)];
	irregular_sum[tn][29] = -1 * g_vals[INT3(i - 1, j - 1, k + 1)];
	irregular_sum[tn][30] = -1 * g_vals[INT3(i - 1, j + 1, k - 1)];
	irregular_sum[tn][31] = -1 * g_vals[INT3(i + 1, j - 1, k - 1)];
	irregular_sum[tn][32] = -1 * g_vals[INT3(i - 1, j + 1, k + 1)];
	irregular_sum[tn][33] = -1 * g_vals[INT3(i + 1, j - 1, k + 1)];
	irregular_sum[tn][34] = -1 * g_vals[INT3(i + 1, j + 1, k - 1)];
	irregular_sum[tn][35] = -1 * g_vals[INT3(i + 1, j + 1, k + 1)];

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * tau);
}