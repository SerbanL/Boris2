#include "stdafx.h"
#include "DemagTFunc.h"

//---------------------SHIFTED VERSIONS (FOR STRAY FIELD)

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
bool DemagTFunc::CalcDiagTens3D_Shifted(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted(n, hRatios, shift)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted(n, hRatios, shift)) return false;
	}

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;
	int K = N.z;

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int k = -n.z + 1; k < n.z; k++) {
				for (int i = 0; i < n.x; i++) {

					DBL3 val = DBL3(
						Ldia_xx_yy_zshifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
						Ldia_xx_yy_zshifted(j, i, k, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
						Ldia_zz_zshifted(k, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;

					Ddiag[((i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val;
					Ddiag[((-i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val;
					Ddiag[((i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val;
					Ddiag[((-i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val;
				}
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

#pragma omp parallel for
		for (int j = -n.y + 1; j < n.y; j++) {
			for (int k = -n.z + 1; k < n.z; k++) {
				for (int i = -n.x + 1; i < n.x; i++) {

					DBL3 val = DBL3(
						Ldia_shifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
						Ldia_shifted(j, i, k, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
						Ldia_shifted(k, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;

					Ddiag[(i + I) % I + ((j + J) % J) * N.x + ((k + K) % K) * N.x * N.y] = val;
				}
			}
		}
	}

	//free memory
	f_vals_xx.clear();
	f_vals_yy.clear();
	f_vals_zz.clear();

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DemagTFunc::CalcOffDiagTens3D_Shifted(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus)
{
	//zero the tensor first
	Dodiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted(n, hRatios, shift)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted(n, hRatios, shift)) return false;
	}

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;
	int K = N.z;

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int k = -n.z + 1; k < n.z; k++) {
				for (int i = 0; i < n.x; i++) {

					DBL3 val = DBL3(
						Lodia_xy_zshifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
						Lodia_xz_yz_zshifted(i, k, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
						Lodia_xz_yz_zshifted(j, k, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;

					Dodiag[((i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(+1, +1, +1);
					Dodiag[((-i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(-1, -1, +1);
					Dodiag[((i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(-1, +1, -1);
					Dodiag[((-i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(+1, -1, -1);
				}
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with)

#pragma omp parallel for
		for (int j = -n.y + 1; j < n.y; j++) {
			for (int k = -n.z + 1; k < n.z; k++) {
				for (int i = -n.x + 1; i < n.x; i++) {

					DBL3 val = DBL3(
						Lodia_shifted(i, j, k, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
						Lodia_shifted(i, k, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
						Lodia_shifted(j, k, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;

					Dodiag[(i + I) % I + ((j + J) % J) * N.x + ((k + K) % K) * N.x * N.y] = val;
				}
			}
		}
	}

	//free memory
	g_vals_xy.clear();
	g_vals_xz.clear();
	g_vals_yz.clear();

	return true;
}

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DemagTFunc::CalcDiagTens2D_Shifted(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted(n, hRatios, shift)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted(n, hRatios, shift)) return false;
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
					Ldia_xx_yy_zshifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
					Ldia_xx_yy_zshifted(j, i, 0, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
					Ldia_zz_zshifted(0, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;

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
					Ldia_shifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
					Ldia_shifted(j, i, 0, n.y, n.x, n.z, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
					Ldia_shifted(0, j, i, n.z, n.y, n.x, hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;

				Ddiag[(i + I) % I + ((j + J) % J) * N.x] = val;
			}

		}
	}

	//free memory
	f_vals_xx.clear();
	f_vals_yy.clear();
	f_vals_zz.clear();

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy only) which has sizes given by N. 
bool DemagTFunc::CalcOffDiagTens2D_Shifted(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus)
{
	//zero the tensor first
	Dodiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted(n, hRatios, shift)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted(n, hRatios, shift)) return false;
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
					Lodia_xy_zshifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
					Lodia_xz_yz_zshifted(i, 0, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
					Lodia_xz_yz_zshifted(j, 0, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;

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
					Lodia_shifted(i, j, 0, n.x, n.y, n.z, hRatios.x, hRatios.y, hRatios.z, g_vals_xy),
					Lodia_shifted(i, 0, j, n.x, n.z, n.y, hRatios.x, hRatios.z, hRatios.y, g_vals_xz),
					Lodia_shifted(j, 0, i, n.y, n.z, n.x, hRatios.y, hRatios.z, hRatios.x, g_vals_yz)) * sign;

				Dodiag[(i + I) % I + ((j + J) % J) * N.x] = val;
			}
		}
	}

	//free memory
	g_vals_xy.clear();
	g_vals_xz.clear();
	g_vals_yz.clear();

	return true;
}