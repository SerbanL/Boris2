#include "stdafx.h"
#include "DemagTFunc.h"

//TEST METHODS

//get diagonal component at given position for given cellsizes
double DemagTFunc::Ldia_single(double x, double y, double z, double hx, double hy, double hz)
{
	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * f(x, y, z);

	main_sum[tn][1] = -4 * f(x + hx, y, z);
	main_sum[tn][2] = -4 * f(x - hx, y, z);
	main_sum[tn][3] = -4 * f(x, y + hy, z);
	main_sum[tn][4] = -4 * f(x, y - hy, z);
	main_sum[tn][5] = -4 * f(x, y, z + hz);
	main_sum[tn][6] = -4 * f(x, y, z - hz);

	main_sum[tn][7] = +2 * f(x - hx, y - hy, z);
	main_sum[tn][8] = +2 * f(x - hx, y + hy, z);
	main_sum[tn][9] = +2 * f(x + hx, y - hy, z);
	main_sum[tn][10] = +2 * f(x + hx, y + hy, z);

	main_sum[tn][11] = +2 * f(x - hx, y, z - hz);
	main_sum[tn][12] = +2 * f(x - hx, y, z + hz);
	main_sum[tn][13] = +2 * f(x + hx, y, z - hz);
	main_sum[tn][14] = +2 * f(x + hx, y, z + hz);

	main_sum[tn][15] = +2 * f(x, y - hy, z - hz);
	main_sum[tn][16] = +2 * f(x, y - hy, z + hz);
	main_sum[tn][17] = +2 * f(x, y + hy, z - hz);
	main_sum[tn][18] = +2 * f(x, y + hy, z + hz);

	main_sum[tn][19] = -1 * f(x - hx, y - hy, z - hz);
	main_sum[tn][20] = -1 * f(x - hx, y - hy, z + hz);
	main_sum[tn][21] = -1 * f(x - hx, y + hy, z - hz);
	main_sum[tn][22] = -1 * f(x + hx, y - hy, z - hz);
	main_sum[tn][23] = -1 * f(x - hx, y + hy, z + hz);
	main_sum[tn][24] = -1 * f(x + hx, y - hy, z + hz);
	main_sum[tn][25] = -1 * f(x + hx, y + hy, z - hz);
	main_sum[tn][26] = -1 * f(x + hx, y + hy, z + hz);

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

//get off-diagonal component at given position for given cellsizes
double DemagTFunc::Lodia_single(double x, double y, double z, double hx, double hy, double hz)
{
	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * g(x, y, z);

	main_sum[tn][1] = -4 * g(x + hx, y, z);
	main_sum[tn][2] = -4 * g(x - hx, y, z);
	main_sum[tn][3] = -4 * g(x, y + hy, z);
	main_sum[tn][4] = -4 * g(x, y - hy, z);
	main_sum[tn][5] = -4 * g(x, y, z + hz);
	main_sum[tn][6] = -4 * g(x, y, z - hz);

	main_sum[tn][7] = +2 * g(x - hx, y - hy, z);
	main_sum[tn][8] = +2 * g(x - hx, y + hy, z);
	main_sum[tn][9] = +2 * g(x + hx, y - hy, z);
	main_sum[tn][10] = +2 * g(x + hx, y + hy, z);

	main_sum[tn][11] = +2 * g(x - hx, y, z - hz);
	main_sum[tn][12] = +2 * g(x - hx, y, z + hz);
	main_sum[tn][13] = +2 * g(x + hx, y, z - hz);
	main_sum[tn][14] = +2 * g(x + hx, y, z + hz);

	main_sum[tn][15] = +2 * g(x, y - hy, z - hz);
	main_sum[tn][16] = +2 * g(x, y - hy, z + hz);
	main_sum[tn][17] = +2 * g(x, y + hy, z - hz);
	main_sum[tn][18] = +2 * g(x, y + hy, z + hz);

	main_sum[tn][19] = -1 * g(x - hx, y - hy, z - hz);
	main_sum[tn][20] = -1 * g(x - hx, y - hy, z + hz);
	main_sum[tn][21] = -1 * g(x - hx, y + hy, z - hz);
	main_sum[tn][22] = -1 * g(x + hx, y - hy, z - hz);
	main_sum[tn][23] = -1 * g(x - hx, y + hy, z + hz);
	main_sum[tn][24] = -1 * g(x + hx, y - hy, z + hz);
	main_sum[tn][25] = -1 * g(x + hx, y + hy, z - hz);
	main_sum[tn][26] = -1 * g(x + hx, y + hy, z + hz);

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

//all diagonal components - single value computation version
DBL3 DemagTFunc::Ldia_single(DBL3 dist, DBL3 h, bool minus)
{
	int sign = 1;
	if (minus) sign = -1;

	return DBL3(
		Ldia_single(dist.x, dist.y, dist.z, h.x, h.y, h.z),
		Ldia_single(dist.y, dist.x, dist.z, h.y, h.x, h.z),
		Ldia_single(dist.z, dist.y, dist.x, h.z, h.y, h.x)) * sign;
}

//all off-diagonal components - single value computation version
DBL3 DemagTFunc::Lodia_single(DBL3 dist, DBL3 h, bool minus)
{
	int sign = 1;
	if (minus) sign = -1;

	return DBL3(
		Lodia_single(dist.x, dist.y, dist.z, h.x, h.y, h.z),
		Lodia_single(dist.x, dist.z, dist.y, h.x, h.z, h.y),
		Lodia_single(dist.y, dist.z, dist.x, h.y, h.z, h.x)) * sign;
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
double DemagTFunc::Ldia_shifted_irregular_xx_yy_single(double x, double y, double z, double hx, double hy, double sz, double dz)
{
	int tn = omp_get_thread_num();

	double del = sz - dz;

	irregular_sum[tn][0] = +4 * f(x, y, z);
	irregular_sum[tn][1] = +4 * f(x, y, z - del);

	irregular_sum[tn][2] = -2 * f(x + hx, y, z);
	irregular_sum[tn][3] = -2 * f(x + hx, y, z - del);
	irregular_sum[tn][4] = -2 * f(x - hx, y, z);
	irregular_sum[tn][5] = -2 * f(x - hx, y, z - del);
	irregular_sum[tn][6] = -2 * f(x, y + hy, z);
	irregular_sum[tn][7] = -2 * f(x, y + hy, z - del);
	irregular_sum[tn][8] = -2 * f(x, y - hy, z);
	irregular_sum[tn][9] = -2 * f(x, y - hy, z - del);
	irregular_sum[tn][10] = -4 * f(x, y, z + dz);
	irregular_sum[tn][11] = -4 * f(x, y, z - sz);

	irregular_sum[tn][12] = +f(x + hx, y + hy, z);
	irregular_sum[tn][13] = +f(x + hx, y + hy, z - del);
	irregular_sum[tn][14] = +f(x + hx, y - hy, z);
	irregular_sum[tn][15] = +f(x + hx, y - hy, z - del);
	irregular_sum[tn][16] = +f(x - hx, y + hy, z);
	irregular_sum[tn][17] = +f(x - hx, y + hy, z - del);
	irregular_sum[tn][18] = +f(x - hx, y - hy, z);
	irregular_sum[tn][19] = +f(x - hx, y - hy, z - del);
	irregular_sum[tn][20] = +2 * f(x + hx, y, z + dz);
	irregular_sum[tn][21] = +2 * f(x + hx, y, z - sz);
	irregular_sum[tn][22] = +2 * f(x - hx, y, z + dz);
	irregular_sum[tn][23] = +2 * f(x - hx, y, z - sz);
	irregular_sum[tn][24] = +2 * f(x, y + hy, z + dz);
	irregular_sum[tn][25] = +2 * f(x, y + hy, z - sz);
	irregular_sum[tn][26] = +2 * f(x, y - hy, z + dz);
	irregular_sum[tn][27] = +2 * f(x, y - hy, z - sz);

	irregular_sum[tn][28] = -1 * f(x - hx, y - hy, z - sz);
	irregular_sum[tn][29] = -1 * f(x - hx, y - hy, z + dz);
	irregular_sum[tn][30] = -1 * f(x - hx, y + hy, z - sz);
	irregular_sum[tn][31] = -1 * f(x + hx, y - hy, z - sz);
	irregular_sum[tn][32] = -1 * f(x - hx, y + hy, z + dz);
	irregular_sum[tn][33] = -1 * f(x + hx, y - hy, z + dz);
	irregular_sum[tn][34] = -1 * f(x + hx, y + hy, z - sz);
	irregular_sum[tn][35] = -1 * f(x + hx, y + hy, z + dz);

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * hx * hy * dz);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
double DemagTFunc::Ldia_shifted_irregular_zz_single(double x, double y, double z, double sx, double dx, double hy, double hz)
{
	int tn = omp_get_thread_num();

	double del = sx - dx;

	irregular_sum[tn][0] = +4 * f(x, y, z);
	irregular_sum[tn][1] = +4 * f(x + del, y, z);

	irregular_sum[tn][2] = -4 * f(x + sx, y, z);
	irregular_sum[tn][3] = -4 * f(x - dx, y, z);
	irregular_sum[tn][4] = -2 * f(x, y, z + hz);
	irregular_sum[tn][5] = -2 * f(x + del, y, z + hz);
	irregular_sum[tn][6] = -2 * f(x, y, z - hz);
	irregular_sum[tn][7] = -2 * f(x + del, y, z - hz);
	irregular_sum[tn][8] = -2 * f(x, y + hy, z);
	irregular_sum[tn][9] = -2 * f(x + del, y + hy, z);
	irregular_sum[tn][10] = -2 * f(x, y - hy, z);
	irregular_sum[tn][11] = -2 * f(x + del, y - hy, z);

	irregular_sum[tn][12] = +2 * f(x + sx, y + hy, z);
	irregular_sum[tn][13] = +2 * f(x + sx, y - hy, z);
	irregular_sum[tn][14] = +2 * f(x - dx, y + hy, z);
	irregular_sum[tn][15] = +2 * f(x - dx, y - hy, z);
	irregular_sum[tn][16] = +2 * f(x + sx, y, z + hz);
	irregular_sum[tn][17] = +2 * f(x + sx, y, z - hz);
	irregular_sum[tn][18] = +2 * f(x - dx, y, z + hz);
	irregular_sum[tn][19] = +2 * f(x - dx, y, z - hz);
	irregular_sum[tn][20] = +f(x, y + hy, z + hz);
	irregular_sum[tn][21] = +f(x + del, y + hy, z + hz);
	irregular_sum[tn][22] = +f(x, y + hy, z - hz);
	irregular_sum[tn][23] = +f(x + del, y + hy, z - hz);
	irregular_sum[tn][24] = +f(x, y - hy, z + hz);
	irregular_sum[tn][25] = +f(x + del, y - hy, z + hz);
	irregular_sum[tn][26] = +f(x, y - hy, z - hz);
	irregular_sum[tn][27] = +f(x + del, y - hy, z - hz);
	
	irregular_sum[tn][28] = -1 * f(x - dx, y - hy, z - hz);
	irregular_sum[tn][29] = -1 * f(x - dx, y - hy, z + hz);
	irregular_sum[tn][30] = -1 * f(x - dx, y + hy, z - hz);
	irregular_sum[tn][31] = -1 * f(x + sx, y - hy, z - hz);
	irregular_sum[tn][32] = -1 * f(x - dx, y + hy, z + hz);
	irregular_sum[tn][33] = -1 * f(x + sx, y - hy, z + hz);
	irregular_sum[tn][34] = -1 * f(x + sx, y + hy, z - hz);
	irregular_sum[tn][35] = -1 * f(x + sx, y + hy, z + hz);

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * dx * hy * hz);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
double DemagTFunc::Lodia_shifted_irregular_xy_single(double x, double y, double z, double hx, double hy, double sz, double dz)
{
	int tn = omp_get_thread_num();

	double del = sz - dz;

	irregular_sum[tn][0] = +4 * g(x, y, z);
	irregular_sum[tn][1] = +4 * g(x, y, z + del);

	irregular_sum[tn][2] = -2 * g(x + hx, y, z);
	irregular_sum[tn][3] = -2 * g(x + hx, y, z + del);
	irregular_sum[tn][4] = -2 * g(x - hx, y, z);
	irregular_sum[tn][5] = -2 * g(x - hx, y, z + del);
	irregular_sum[tn][6] = -2 * g(x, y + hy, z);
	irregular_sum[tn][7] = -2 * g(x, y + hy, z + del);
	irregular_sum[tn][8] = -2 * g(x, y - hy, z);
	irregular_sum[tn][9] = -2 * g(x, y - hy, z + del);
	irregular_sum[tn][10] = -4 * g(x, y, z + sz);
	irregular_sum[tn][11] = -4 * g(x, y, z - dz);

	irregular_sum[tn][12] = +g(x + hx, y + hy, z);
	irregular_sum[tn][13] = +g(x + hx, y + hy, z + del);
	irregular_sum[tn][14] = +g(x + hx, y - hy, z);
	irregular_sum[tn][15] = +g(x + hx, y - hy, z + del);
	irregular_sum[tn][16] = +g(x - hx, y + hy, z);
	irregular_sum[tn][17] = +g(x - hx, y + hy, z + del);
	irregular_sum[tn][18] = +g(x - hx, y - hy, z);
	irregular_sum[tn][19] = +g(x - hx, y - hy, z + del);
	irregular_sum[tn][20] = +2 * g(x + hx, y, z + sz);
	irregular_sum[tn][21] = +2 * g(x + hx, y, z - dz);
	irregular_sum[tn][22] = +2 * g(x - hx, y, z + sz);
	irregular_sum[tn][23] = +2 * g(x - hx, y, z - dz);
	irregular_sum[tn][24] = +2 * g(x, y + hy, z + sz);
	irregular_sum[tn][25] = +2 * g(x, y + hy, z - dz);
	irregular_sum[tn][26] = +2 * g(x, y - hy, z + sz);
	irregular_sum[tn][27] = +2 * g(x, y - hy, z - dz);

	irregular_sum[tn][28] = -1 * g(x - hx, y - hy, z - dz);
	irregular_sum[tn][29] = -1 * g(x - hx, y - hy, z + sz);
	irregular_sum[tn][30] = -1 * g(x - hx, y + hy, z - dz);
	irregular_sum[tn][31] = -1 * g(x + hx, y - hy, z - dz);
	irregular_sum[tn][32] = -1 * g(x - hx, y + hy, z + sz);
	irregular_sum[tn][33] = -1 * g(x + hx, y - hy, z + sz);
	irregular_sum[tn][34] = -1 * g(x + hx, y + hy, z - dz);
	irregular_sum[tn][35] = -1 * g(x + hx, y + hy, z + sz);

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * hx * hy * dz);
}

//off-diagonal components for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
double DemagTFunc::Lodia_shifted_irregular_xz_yz_single(double x, double y, double z, double hx, double sy, double dy, double hz)
{
	int tn = omp_get_thread_num();

	double del = sy - dy;

	irregular_sum[tn][0] = +4 * g(x, y, z);
	irregular_sum[tn][1] = +4 * g(x, y + del, z);

	irregular_sum[tn][2] = -2 * g(x + hx, y, z);
	irregular_sum[tn][3] = -2 * g(x + hx, y + del, z);
	irregular_sum[tn][4] = -2 * g(x - hx, y, z);
	irregular_sum[tn][5] = -2 * g(x - hx, y + del, z);
	irregular_sum[tn][6] = -4 * g(x, y + sy, z);
	irregular_sum[tn][7] = -4 * g(x, y - dy, z);
	irregular_sum[tn][8] = -2 * g(x, y, z + hz);
	irregular_sum[tn][9] = -2 * g(x, y + del, z + hz);
	irregular_sum[tn][10] = -2 * g(x, y, z - hz);
	irregular_sum[tn][11] = -2 * g(x, y + del, z - hz);

	irregular_sum[tn][12] = +2 * g(x + hx, y + sy, z);
	irregular_sum[tn][13] = +2 * g(x + hx, y - dy, z);
	irregular_sum[tn][14] = +2 * g(x - hx, y + sy, z);
	irregular_sum[tn][15] = +2 * g(x - hx, y - dy, z);
	irregular_sum[tn][16] = +g(x + hx, y, z + hz);
	irregular_sum[tn][17] = +g(x + hx, y + del, z + hz);
	irregular_sum[tn][18] = +g(x + hx, y, z - hz);
	irregular_sum[tn][19] = +g(x + hx, y + del, z - hz);
	irregular_sum[tn][20] = +g(x - hx, y, z + hz);
	irregular_sum[tn][21] = +g(x - hx, y + del, z + hz);
	irregular_sum[tn][22] = +g(x - hx, y, z - hz);
	irregular_sum[tn][23] = +g(x - hx, y + del, z - hz);
	irregular_sum[tn][24] = +2 * g(x, y + sy, z + hz);
	irregular_sum[tn][25] = +2 * g(x, y + sy, z - hz);
	irregular_sum[tn][26] = +2 * g(x, y - dy, z + hz);
	irregular_sum[tn][27] = +2 * g(x, y - dy, z - hz);

	irregular_sum[tn][28] = -1 * g(x - hx, y - dy, z - hz);
	irregular_sum[tn][29] = -1 * g(x - hx, y - dy, z + hz);
	irregular_sum[tn][30] = -1 * g(x - hx, y + sy, z - hz);
	irregular_sum[tn][31] = -1 * g(x + hx, y - dy, z - hz);
	irregular_sum[tn][32] = -1 * g(x - hx, y + sy, z + hz);
	irregular_sum[tn][33] = -1 * g(x + hx, y - dy, z + hz);
	irregular_sum[tn][34] = -1 * g(x + hx, y + sy, z - hz);
	irregular_sum[tn][35] = -1 * g(x + hx, y + sy, z + hz);

	return sum_KahanNeumaier(irregular_sum[tn]) / (4 * PI * hx * dy * hz);
}