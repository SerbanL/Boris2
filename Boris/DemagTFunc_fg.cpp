#include "stdafx.h"
#include "DemagTFunc.h"

//---------------------f and g VECTORS COMPUTATION

//Basis f and g functions - these are correct
//The checks in them are absolutely essential, otherwise you'll run into NaNs sooner or later, especially in multi-layered convolution.
//Moreover don't check a double for zero, use the IsZ or IsNZ functions. If you do, guaranteed eventually the check will fail and result in NaNs.

double DemagTFunc::f(double x, double y, double z)
{
	int tn = omp_get_thread_num();

	double xsq = x * x;
	double ysq = y * y;
	double zsq = z * z;

	double Rxy_sq = xsq + ysq;
	double Rxz_sq = xsq + zsq;

	double R = sqrt(xsq + ysq + zsq);
	
	f_sum[tn][0] = (2.0 * xsq - ysq - zsq) * R / 6.0;

	f_sum[tn][1] = 0;
	if (IsNZ(Rxz_sq)) {

		double arg = 2.0 * y * (y + R) / Rxz_sq;
		if (arg > -1.0) f_sum[tn][1] = y * (zsq - xsq) * log1p(arg) / 4.0;
	}

	f_sum[tn][2] = 0;
	if (IsNZ(Rxy_sq)) {

		double arg = 2.0 * z * (z + R) / Rxy_sq;
		if (arg > -1.0) f_sum[tn][2] = z * (ysq - xsq) * log1p(arg) / 4.0;
	}

	if (IsNZ(x)) f_sum[tn][3] = -x*y*z*atan(y*z / (x*R));
	else f_sum[tn][3] = 0;

	return sum_KahanNeumaier(f_sum[tn]);
}

double DemagTFunc::g(double x, double y, double z)
{
	int tn = omp_get_thread_num();

	double xsq = x * x;
	double ysq = y * y;
	double zsq = z * z;

	double Rxy_sq = xsq + ysq;
	double Rxz_sq = xsq + zsq;
	double Ryz_sq = ysq + zsq;

	double R = sqrt(xsq + ysq + zsq);

	g_sum[tn][0] = -x * y * R / 3.0;
	
	g_sum[tn][1] = 0;
	if (IsNZ(Rxy_sq)) {

		double arg = 2.0 * z * (z + R) / Rxy_sq;
		if (arg > -1.0) g_sum[tn][1] = x*y*z*log1p(arg) / 2.0;
	}

	g_sum[tn][2] = 0;
	if (IsNZ(Ryz_sq)) {

		double arg = 2.0 * x * (x + R) / Ryz_sq;
		if (arg > -1.0) g_sum[tn][2] = y * (3.0 * zsq - ysq) * log1p(arg) / 12.0;
	}

	g_sum[tn][3] = 0;
	if (IsNZ(Rxz_sq)) {

		double arg = 2.0 * y * (y + R) / Rxz_sq;
		if (arg > -1.0) g_sum[tn][3] = x * (3.0 * zsq - xsq) * log1p(arg) / 12.0;
	}

	if (IsNZ(z)) g_sum[tn][4] = -zsq * z * atan(x*y / (z*R)) / 6.0;
	else g_sum[tn][4] = 0;

	if (IsNZ(y)) g_sum[tn][5] = -ysq * z * atan(x*z / (y*R)) / 2.0;
	else g_sum[tn][5] = 0;

	if (IsNZ(x)) g_sum[tn][6] = -xsq * z * atan(y*z / (x*R)) / 2.0;
	else g_sum[tn][6] = 0;

	return sum_KahanNeumaier(g_sum[tn]);
}

bool DemagTFunc::fill_f_vals(INT3 n, DBL3 hRatios, int asymptotic_distance)
{
	int nx_dist = n.x, ny_dist = n.y, nz_dist = n.z;

	if (asymptotic_distance > 0) {

		nx_dist = minimum(n.x, asymptotic_distance);
		ny_dist = minimum(n.y, asymptotic_distance);
		nz_dist = minimum(n.z, asymptotic_distance);
	}

	if (!f_vals_xx.assign(SZ3(nx_dist + 2, ny_dist + 2, nz_dist + 2), 0.0)) return false;
	if (!f_vals_yy.assign(SZ3(ny_dist + 2, nx_dist + 2, nz_dist + 2), 0.0)) return false;
	if (!f_vals_zz.assign(SZ3(nz_dist + 2, ny_dist + 2, nx_dist + 2), 0.0)) return false;

	if (IsE(hRatios.x, hRatios.y) && IsE(hRatios.y, hRatios.z)) {

		//cubic cell : re-use calculated f values as much as possible

		//first calculate xx values
#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					f_vals_xx[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z);
				}
			}
		}

		//now calculate yy and zz values only if we don't have them already from previous step

#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					if (j <= nx_dist && i <= ny_dist) {

						//we already have this value : get it
						f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f_vals_xx[INT3(j + 1, i + 1, k + 1)];
					}
					else f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z);

					if (k <= nx_dist && i <= nz_dist) {

						//we already have this value : get it
						f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f_vals_xx[INT3(k + 1, j + 1, i + 1)];
					}
					else f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f(k * hRatios.z, j * hRatios.y, i * hRatios.x);
				}
			}
		}
	}
	else if (IsE(hRatios.x, hRatios.y)) {

		//square cell : re-use calculated f values as much as possible

		//first calculate xx and zz values
#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					f_vals_xx[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z);
					f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f(k * hRatios.z, j * hRatios.y, i * hRatios.x);
				}
			}
		}

		//now calculate yy values only if we don't have them already from previous step

#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					if (j <= nx_dist && i <= ny_dist) {

						//we already have this value : get it
						f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f_vals_xx[INT3(j + 1, i + 1, k + 1)];
					}
					else f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z);
				}
			}
		}
	}
	else {

		//general case : calculate xx, yy, and zz values separately

#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					f_vals_xx[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z);
					f_vals_yy[(j + 1) + (i + 1) * (ny_dist + 2) + (k + 1) * (ny_dist + 2) * (nx_dist + 2)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z);
					f_vals_zz[(k + 1) + (j + 1) * (nz_dist + 2) + (i + 1) * (nz_dist + 2) * (ny_dist + 2)] = f(k * hRatios.z, j * hRatios.y, i * hRatios.x);
				}
			}
		}
	}

	return true;
}

bool DemagTFunc::fill_g_vals(INT3 n, DBL3 hRatios, int asymptotic_distance)
{
	int nx_dist = n.x, ny_dist = n.y, nz_dist = n.z;

	if (asymptotic_distance > 0) {

		nx_dist = minimum(n.x, asymptotic_distance);
		ny_dist = minimum(n.y, asymptotic_distance);
		nz_dist = minimum(n.z, asymptotic_distance);
	}

	if (!g_vals_xy.assign(SZ3(nx_dist + 2, ny_dist + 2, nz_dist + 2), 0.0)) return false;
	if (!g_vals_xz.assign(SZ3(nx_dist + 2, nz_dist + 2, ny_dist + 2), 0.0)) return false;
	if (!g_vals_yz.assign(SZ3(ny_dist + 2, nz_dist + 2, nx_dist + 2), 0.0)) return false;

	if (IsE(hRatios.x, hRatios.y) && IsE(hRatios.y, hRatios.z)) {

		//cubic cell : re-use calculated f values as much as possible

		//first calculate xy values
#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					g_vals_xy[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z);
				}
			}
		}

		//now calculate xy and xz values only if we don't have them from previous step

#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					if (k <= ny_dist && j <= nz_dist) {

						//we already have this value : get it
						g_vals_xz[(i + 1) + (k + 1) * (nx_dist + 2) + (j + 1) * (nx_dist + 2) * (nz_dist + 2)] = g_vals_xy[INT3(i + 1, k + 1, j + 1)];
					}
					else g_vals_xz[(i + 1) + (k + 1) * (nx_dist + 2) + (j + 1) * (nx_dist + 2) * (nz_dist + 2)] = g(i * hRatios.x, k * hRatios.z, j * hRatios.y);


					if (j <= nx_dist && k <= ny_dist && i <= nz_dist) {

						//we already have this value : get it
						g_vals_yz[(j + 1) + (k + 1) * (ny_dist + 2) + (i + 1) * (ny_dist + 2) * (nz_dist + 2)] = g_vals_xy[INT3(j + 1, k + 1, i + 1)];
					}
					else g_vals_yz[(j + 1) + (k + 1) * (ny_dist + 2) + (i + 1) * (ny_dist + 2) * (nz_dist + 2)] = g(j * hRatios.y, k * hRatios.z, i * hRatios.x);
				}
			}
		}
	}
	else {

		//general case : calculate xy, xz, and yz values separately

#pragma omp parallel for
		for (int j = -1; j <= ny_dist; j++) {
			for (int k = -1; k <= nz_dist; k++) {
				for (int i = -1; i <= nx_dist; i++) {

					g_vals_xy[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z);
					g_vals_xz[(i + 1) + (k + 1) * (nx_dist + 2) + (j + 1) * (nx_dist + 2) * (nz_dist + 2)] = g(i * hRatios.x, k * hRatios.z, j * hRatios.y);
					g_vals_yz[(j + 1) + (k + 1) * (ny_dist + 2) + (i + 1) * (ny_dist + 2) * (nz_dist + 2)] = g(j * hRatios.y, k * hRatios.z, i * hRatios.x);
				}
			}
		}
	}

	return true;
}

bool DemagTFunc::fill_g2D_vals(INT3 n, DBL3 hRatios, int asymptotic_distance)
{
	int nx_dist = n.x, ny_dist = n.y, nz_dist = n.z;

	if (asymptotic_distance > 0) {

		nx_dist = minimum(n.x, asymptotic_distance);
		ny_dist = minimum(n.y, asymptotic_distance);
		nz_dist = minimum(n.z, asymptotic_distance);
	}

	if (!g_vals_xy.assign(SZ3(nx_dist + 2, ny_dist + 2, nz_dist + 2), 0.0)) return false;

#pragma omp parallel for
	for (int j = -1; j <= ny_dist; j++) {
		for (int k = -1; k <= nz_dist; k++) {
			for (int i = -1; i <= nx_dist; i++) {

				g_vals_xy[(i + 1) + (j + 1) * (nx_dist + 2) + (k + 1) * (nx_dist + 2) * (ny_dist + 2)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z);
			}
		}
	}

	return true;
}

//shifted versions for stray field
bool DemagTFunc::fill_f_vals_shifted(INT3 n, DBL3 hRatios, DBL3 shift)
{
	//shifted versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!f_vals_xx.assign(SZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy.assign(SZ3(2 * n.y + 1, 2 * n.x + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz.assign(SZ3(2 * n.z + 1, 2 * n.y + 1, 2 * n.x + 1), 0.0)) return false;

#pragma omp parallel for
	for (int j = -n.y; j <= n.y; j++) {
		for (int k = -n.z; k <= n.z; k++) {
			for (int i = -n.x; i <= n.x; i++) {

				f_vals_xx[(i + n.x) + (j + n.y) * (2 * n.x + 1) + (k + n.z) * (2 * n.x + 1) * (2 * n.y + 1)] = f(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z);
				f_vals_yy[(j + n.y) + (i + n.x) * (2 * n.y + 1) + (k + n.z) * (2 * n.y + 1) * (2 * n.x + 1)] = f(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z);
				f_vals_zz[(k + n.z) + (j + n.y) * (2 * n.z + 1) + (i + n.x) * (2 * n.z + 1) * (2 * n.y + 1)] = f(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x);
			}
		}
	}

	return true;
}

bool DemagTFunc::fill_g_vals_shifted(INT3 n, DBL3 hRatios, DBL3 shift)
{
	//shifted versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!g_vals_xy.assign(SZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz.assign(SZ3(2 * n.x + 1, 2 * n.z + 1, 2 * n.y + 1), 0.0)) return false;
	if (!g_vals_yz.assign(SZ3(2 * n.y + 1, 2 * n.z + 1, 2 * n.x + 1), 0.0)) return false;

#pragma omp parallel for
	for (int j = -n.y; j <= n.y; j++) {
		for (int k = -n.z; k <= n.z; k++) {
			for (int i = -n.x; i <= n.x; i++) {

				g_vals_xy[(i + n.x) + (j + n.y) * (2 * n.x + 1) + (k + n.z) * (2 * n.x + 1) * (2 * n.y + 1)] = g(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z);
				g_vals_xz[(i + n.x) + (k + n.z) * (2 * n.x + 1) + (j + n.y) * (2 * n.x + 1) * (2 * n.z + 1)] = g(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y);
				g_vals_yz[(j + n.y) + (k + n.z) * (2 * n.y + 1) + (i + n.x) * (2 * n.y + 1) * (2 * n.z + 1)] = g(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x);
			}
		}
	}

	return true;
}

//shifted versions for stray field -> in particular for z shift only
bool DemagTFunc::fill_f_vals_zshifted(INT3 n, DBL3 hRatios, DBL3 shift, int asymptotic_distance)
{
	//for z shift only can use xy plane symmetries

	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!f_vals_xx.assign(SZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy.assign(SZ3(n.y + 2, n.x + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz.assign(SZ3(2 * n.z + 1, n.y + 2, n.x + 2), 0.0)) return false;

#pragma omp parallel for
	for (int j = -1; j <= n.y; j++) {
		for (int k = -n.z; k <= n.z; k++) {
			for (int i = -1; i <= n.x; i++) {

				f_vals_xx[INT3(i + 1, j + 1, k + n.z)] = f(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z);
				f_vals_yy[INT3(j + 1, i + 1, k + n.z)] = f(j * hRatios.y, i * hRatios.x, k * hRatios.z + shift.z);
				f_vals_zz[INT3(k + n.z, j + 1, i + 1)] = f(k * hRatios.z + shift.z, j * hRatios.y, i * hRatios.x);
			}
		}
	}

	return true;
}

bool DemagTFunc::fill_g_vals_zshifted(INT3 n, DBL3 hRatios, DBL3 shift, int asymptotic_distance)
{
	//for z shift only can use xy plane symmetries

	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!g_vals_xy.assign(SZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz.assign(SZ3(n.x + 2, 2 * n.z + 1, n.y + 2), 0.0)) return false;
	if (!g_vals_yz.assign(SZ3(n.y + 2, 2 * n.z + 1, n.x + 2), 0.0)) return false;

#pragma omp parallel for
	for (int j = -1; j <= n.y; j++) {
		for (int k = -n.z; k <= n.z; k++) {
			for (int i = -1; i <= n.x; i++) {
				
				g_vals_xy[INT3(i + 1, j + 1, k + n.z)] = g(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z);
				g_vals_xz[INT3(i + 1, k + n.z, j + 1)] = g(i * hRatios.x, k * hRatios.z + shift.z, j * hRatios.y);
				g_vals_yz[INT3(j + 1, k + n.z, i + 1)] = g(j * hRatios.y, k * hRatios.z + shift.z, i * hRatios.x);
			}
		}
	}

	return true;
}

//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz)
bool DemagTFunc::fill_f_vals_shifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift)
{
	//shifted (and irregular) versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!f_vals_xx.assign(SZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy.assign(SZ3(2 * n.y + 1, 2 * n.x + 1, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz.assign(SZ3(2 * n.z + 1, 2 * n.y + 1, 2 * n.x + 1), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!f_vals_xx_del.assign(SZ3(2 * n.x + 1, 2 * n.y + 1, 1), 0.0)) return false;
	if (!f_vals_yy_del.assign(SZ3(2 * n.y + 1, 2 * n.x + 1, 1), 0.0)) return false;
	if (!f_vals_zz_del.assign(SZ3(1, 2 * n.y + 1, 2 * n.x + 1), 0.0)) return false;

	double del = h_src.z - h_dst.z;

#pragma omp parallel for
	for (int j = -n.y; j <= n.y; j++) {
		for (int i = -n.x; i <= n.x; i++) {

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

	return true;
}

//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz)
bool DemagTFunc::fill_g_vals_shifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift)
{
	//shifted (and irregular) versions cannot use symmetries in general so fill from -n to n inclusive, hence 2n + 1 in each dimension -> tensor elements range in [-n + 1, n)

	if (!g_vals_xy.assign(SZ3(2 * n.x + 1, 2 * n.y + 1, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz.assign(SZ3(2 * n.x + 1, 2 * n.z + 1, 2 * n.y + 1), 0.0)) return false;
	if (!g_vals_yz.assign(SZ3(2 * n.y + 1, 2 * n.z + 1, 2 * n.x + 1), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!g_vals_xy_del.assign(SZ3(2 * n.x + 1, 2 * n.y + 1, 1), 0.0)) return false;
	if (!g_vals_xz_del.assign(SZ3(2 * n.x + 1, 1, 2 * n.y + 1), 0.0)) return false;
	if (!g_vals_yz_del.assign(SZ3(2 * n.y + 1, 1, 2 * n.x + 1), 0.0)) return false;

	double del = h_src.z - h_dst.z;

#pragma omp parallel for
	for (int j = -n.y; j <= n.y; j++) {
		for (int i = -n.x; i <= n.x; i++) {

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

	return true;
}

//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz) -> in particular for z shift only
bool DemagTFunc::fill_f_vals_zshifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift, int asymptotic_distance)
{
	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!f_vals_xx.assign(SZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_yy.assign(SZ3(n.y + 2, n.x + 2, 2 * n.z + 1), 0.0)) return false;
	if (!f_vals_zz.assign(SZ3(2 * n.z + 1, n.y + 2, n.x + 2), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!f_vals_xx_del.assign(SZ3(n.x + 2, n.y + 2, 1), 0.0)) return false;
	if (!f_vals_yy_del.assign(SZ3(n.y + 2, n.x + 2, 1), 0.0)) return false;
	if (!f_vals_zz_del.assign(SZ3(1, n.y + 2, n.x + 2), 0.0)) return false;

	double del = h_src.z - h_dst.z;

#pragma omp parallel for
	for (int j = -1; j <= n.y; j++) {
		for (int i = -1; i <= n.x; i++) {

			//xx, yy : +dz, -sz, -del
			//zz : +sx, -dx, +del

			//k = -1
			
			f_vals_xx[INT3(i + 1, j + 1, 0)] = f(i * h_dst.x, j * h_dst.y, -h_src.z + shift.z);
			f_vals_yy[INT3(j + 1, i + 1, 0)] = f(j * h_dst.y, i * h_dst.x, -h_src.z + shift.z);
			f_vals_zz[INT3(0, j + 1, i + 1)] = f(-h_dst.z + shift.z, j * h_dst.y, i * h_dst.x);

			//k = 0
			f_vals_xx[INT3(i + 1, j + 1, 1)] = f(i * h_dst.x, j * h_dst.y, shift.z);
			f_vals_yy[INT3(j + 1, i + 1, 1)] = f(j * h_dst.y, i * h_dst.x, shift.z);
			f_vals_zz[INT3(1, j + 1, i + 1)] = f(shift.z, j * h_dst.y, i * h_dst.x);

			//k = 1
			f_vals_xx[INT3(i + 1, j + 1, 2)] = f(i * h_dst.x, j * h_dst.y, h_dst.z + shift.z);
			f_vals_yy[INT3(j + 1, i + 1, 2)] = f(j * h_dst.y, i * h_dst.x, h_dst.z + shift.z);
			f_vals_zz[INT3(2, j + 1, i + 1)] = f(h_src.z + shift.z, j * h_dst.y, i * h_dst.x);

			//del versions
			f_vals_xx_del[INT3(i + 1, j + 1, 0)] = f(i * h_dst.x, j * h_dst.y, -del + shift.z);
			f_vals_yy_del[INT3(j + 1, i + 1, 0)] = f(j * h_dst.y, i * h_dst.x, -del + shift.z);
			f_vals_zz_del[INT3(0, j + 1, i + 1)] = f(+del + shift.z, j * h_dst.y, i * h_dst.x);
		}
	}

	return true;
}

bool DemagTFunc::fill_g_vals_zshifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift, int asymptotic_distance)
{
	if (asymptotic_distance > 0) {

		n.x = minimum(n.x, asymptotic_distance);
		n.y = minimum(n.y, asymptotic_distance);
		//limiting n.z to asymptotic distance is not straightforward
	}

	if (!g_vals_xy.assign(SZ3(n.x + 2, n.y + 2, 2 * n.z + 1), 0.0)) return false;
	if (!g_vals_xz.assign(SZ3(n.x + 2, 2 * n.z + 1, n.y + 2), 0.0)) return false;
	if (!g_vals_yz.assign(SZ3(n.y + 2, 2 * n.z + 1, n.x + 2), 0.0)) return false;

	//need additional del versions for irregular calculations -> these only need 1 elements along z direction
	if (!g_vals_xy_del.assign(SZ3(n.x + 2, n.y + 2, 1), 0.0)) return false;
	if (!g_vals_xz_del.assign(SZ3(n.x + 2, 1, n.y + 2), 0.0)) return false;
	if (!g_vals_yz_del.assign(SZ3(n.y + 2, 1, n.x + 2), 0.0)) return false;

	double del = h_src.z - h_dst.z;

#pragma omp parallel for
	for (int j = -1; j <= n.y; j++) {
		for (int i = -1; i <= n.x; i++) {

			//xy : +sz, -dz, +del
			//xz, yz : +sy, -dy, +del

			//k = -1
			g_vals_xy[INT3(i + 1, j + 1, 0)] = g(i * h_dst.x, j * h_dst.y, -h_dst.z + shift.z);
			g_vals_xz[INT3(i + 1, 0, j + 1)] = g(i * h_dst.x, -h_dst.z + shift.z, j * h_dst.y);
			g_vals_yz[INT3(j + 1, 0, i + 1)] = g(j * h_dst.y, -h_dst.z + shift.z, i * h_dst.x);

			//k = 0
			g_vals_xy[INT3(i + 1, j + 1, 1)] = g(i * h_dst.x, j * h_dst.y, shift.z);
			g_vals_xz[INT3(i + 1, 1, j + 1)] = g(i * h_dst.x, shift.z, j * h_dst.y);
			g_vals_yz[INT3(j + 1, 1, i + 1)] = g(j * h_dst.y, shift.z, i * h_dst.x);

			//k = 1
			g_vals_xy[INT3(i + 1, j + 1, 2)] = g(i * h_dst.x, j * h_dst.y, h_src.z + shift.z);
			g_vals_xz[INT3(i + 1, 2, j + 1)] = g(i * h_dst.x, h_src.z + shift.z, j * h_dst.y);
			g_vals_yz[INT3(j + 1, 2, i + 1)] = g(j * h_dst.y, h_src.z + shift.z, i * h_dst.x);

			//del versions
			g_vals_xy_del[INT3(i + 1, j + 1, 0)] = g(i * h_dst.x, j * h_dst.y, del + shift.z);
			g_vals_xz_del[INT3(i + 1, 0, j + 1)] = g(i * h_dst.x, del + shift.z, j * h_dst.y);
			g_vals_yz_del[INT3(j + 1, 0, i + 1)] = g(j * h_dst.y, del + shift.z, i * h_dst.x);
		}
	}

	return true;
}