#include "stdafx.h"
#include "DemagTFunc.h"

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
bool DemagTFunc::CalcDiagTens2D_Shifted_Irregular(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 s, DBL3 d, DBL3 shift, bool minus)
{
	//only use irregular version if you have to
	if (s == d) return CalcDiagTens2D_Shifted(Ddiag, n, N, s, shift, minus);

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted_irregular(n, s, d, shift)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted_irregular(n, s, d, shift)) return false;
	}

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val = DBL3(
					Ldia_zshifted_irregular_xx_yy(i, j, 0, n.x, n.y, n.z, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
					Ldia_zshifted_irregular_xx_yy(j, i, 0, n.y, n.x, n.z, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
					Ldia_zshifted_irregular_zz(0, j, i, n.z, n.y, n.x, d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;

				Ddiag[((i + I) % I) + ((j + J) % J)*I] = val;
				Ddiag[((-i + I) % I) + ((j + J) % J)*I] = val;
				Ddiag[((i + I) % I) + ((-j + J) % J)*I] = val;
				Ddiag[((-i + I) % I) + ((-j + J) % J)*I] = val;
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

#pragma omp parallel for
		for (int j = -n.y + 1; j < n.y; j++) {
			for (int i = -n.x + 1; i < n.x; i++) {

				DBL3 val = DBL3(
					Ldia_shifted_irregular_xx_yy(i, j, 0, n.x, n.y, n.z, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
					Ldia_shifted_irregular_xx_yy(j, i, 0, n.y, n.x, n.z, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
					Ldia_shifted_irregular_zz(0, j, i, n.z, n.y, n.x, d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;

				Ddiag[(i + I) % I + ((j + J) % J) * N.x] = val;
			}
		}
	}

	//free memory
	f_vals_xx.clear();
	f_vals_yy.clear();
	f_vals_zz.clear();
	f_vals_xx_del.clear();
	f_vals_yy_del.clear();
	f_vals_zz_del.clear();

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
bool DemagTFunc::CalcOffDiagTens2D_Shifted_Irregular(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 s, DBL3 d, DBL3 shift, bool minus)
{
	//only use irregular version if you have to
	if (s == d) return CalcOffDiagTens2D_Shifted(Dodiag, n, N, s, shift, minus);

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted_irregular(n, s, d, shift)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted_irregular(n, s, d, shift)) return false;
	}

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val = DBL3(
					Lodia_zshifted_irregular_xy(i, j, 0, n.x, n.y, n.z, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del),
					Lodia_zshifted_irregular_xz_yz(i, 0, j, n.x, n.z, n.y, d.x*d.y*d.z, g_vals_xz, g_vals_xz_del),
					Lodia_zshifted_irregular_xz_yz(j, 0, i, n.y, n.z, n.x, d.x*d.y*d.z, g_vals_yz, g_vals_yz_del)) * sign;

				Dodiag[((i + I) % I) + ((j + J) % J)*I] = val & DBL3(+1, +1, +1);
				Dodiag[((-i + I) % I) + ((j + J) % J)*I] = val & DBL3(-1, -1, +1);
				Dodiag[((i + I) % I) + ((-j + J) % J)*I] = val & DBL3(-1, +1, -1);
				Dodiag[((-i + I) % I) + ((-j + J) % J)*I] = val & DBL3(+1, -1, -1);
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

#pragma omp parallel for
		for (int j = -n.y + 1; j < n.y; j++) {
			for (int i = -n.x + 1; i < n.x; i++) {

				DBL3 val = DBL3(
					Lodia_shifted_irregular_xy(i, j, 0, n.x, n.y, n.z, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del),
					Lodia_shifted_irregular_xz_yz(i, 0, j, n.x, n.z, n.y, d.x*d.y*d.z, g_vals_xz, g_vals_xz_del),
					Lodia_shifted_irregular_xz_yz(j, 0, i, n.y, n.z, n.x, d.x*d.y*d.z, g_vals_yz, g_vals_yz_del)) * sign;

				Dodiag[(i + I) % I + ((j + J) % J) * N.x] = val;
			}
		}
	}

	//free memory
	g_vals_xy.clear();
	g_vals_xz.clear();
	g_vals_yz.clear();
	g_vals_xy_del.clear();
	g_vals_xz_del.clear();
	g_vals_yz_del.clear();

	return true;
}