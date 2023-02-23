#include "stdafx.h"
#include "DemagTFunc.h"

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
bool DemagTFunc::CalcDiagTens2D_Shifted_Irregular(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 s, DBL3 d, DBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	//only use irregular version if you have to
	if (s * 1e-9 == d * 1e-9) return CalcDiagTens2D_Shifted(Ddiag, n, N, s, shift, minus, asymptotic_distance);

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted_irregular(n, s, d, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted_irregular(n, s, d, shift)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(d.x, d.y, d.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(d.y, d.x, d.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(d.z, d.y, d.x);

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val;

				//apply asymptotic equations?
				if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(mod(shift.z / d.z))) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

					//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
					//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
					val = DBL3(
						demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x, j * d.y, shift.z) * sign,
						demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y, i * d.x, shift.z) * sign,
						demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * d.y, i * d.x) * sign) * s.z;
				}
				else {

					val = DBL3(
						Ldia_zshifted_irregular_xx_yy(i, j, 0, n.x, n.y, n.z, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
						Ldia_zshifted_irregular_xx_yy(j, i, 0, n.y, n.x, n.z, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
						Ldia_zshifted_irregular_zz(0, j, i, n.z, n.y, n.x, d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;
				}

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

				DBL3 val;

				double is = mod(i + (shift.x / d.x));
				double js = mod(j + (shift.y / d.y));
				double ks = mod(shift.z / d.z);

				//apply asymptotic equations?
				if (asymptotic_distance > 0 && 
					(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
					int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {
				
					//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
					//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
					val = DBL3(
						demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x + shift.x, j * d.y + shift.y, shift.z) * sign,
						demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y + shift.y, i * d.x + shift.x, shift.z) * sign,
						demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * d.y + shift.y, i * d.x + shift.x) * sign) * s.z;
				}
				else {

					val = DBL3(
						Ldia_shifted_irregular_xx_yy(i, j, 0, n.x, n.y, n.z, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
						Ldia_shifted_irregular_xx_yy(j, i, 0, n.y, n.x, n.z, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
						Ldia_shifted_irregular_zz(0, j, i, n.z, n.y, n.x, d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;
				}

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
bool DemagTFunc::CalcOffDiagTens2D_Shifted_Irregular(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 s, DBL3 d, DBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	Dodiag.set(DBL3());

	//only use irregular version if you have to
	if (s * 1e-9 == d * 1e-9) return CalcOffDiagTens2D_Shifted(Dodiag, n, N, s, shift, minus, asymptotic_distance);

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted_irregular(n, s, d, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted_irregular(n, s, d, shift)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(d.x, d.y, d.z);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xz(d.x, d.z, d.y);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_yz(d.y, d.z, d.x);

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;

	if (IsZ(shift.x) && IsZ(shift.y)) {

		//z shift only, so use full xy plane symmetries

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val;

				//apply asymptotic equations?
				if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(mod(shift.z / d.z))) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

					//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
					//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
					//D12, D13, D23
					val = DBL3(
						demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x, j * d.y, shift.z) * sign,
						demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x, shift.z, j * d.y) * sign,
						demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y, shift.z, i * d.x) * sign) * s.z;
				}
				else {

					val = DBL3(
						Lodia_zshifted_irregular_xy(i, j, 0, n.x, n.y, n.z, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del),
						Lodia_zshifted_irregular_xz_yz(i, 0, j, n.x, n.z, n.y, d.x*d.y*d.z, g_vals_xz, g_vals_xz_del),
						Lodia_zshifted_irregular_xz_yz(j, 0, i, n.y, n.z, n.x, d.x*d.y*d.z, g_vals_yz, g_vals_yz_del)) * sign;
				}

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

				DBL3 val;

				double is = mod(i + (shift.x / d.x));
				double js = mod(j + (shift.y / d.y));
				double ks = mod(shift.z / d.z);

				//apply asymptotic equations?
				if (asymptotic_distance > 0 &&
					(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
					int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

					//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
					//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
					//D12, D13, D23
					val = DBL3(
						demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x + shift.x, j * d.y + shift.y, shift.z) * sign,
						demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x + shift.x, shift.z, j * d.y + shift.y) * sign,
						demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y + shift.y, shift.z, i * d.x + shift.x) * sign) * s.z;
				}
				else {

					val = DBL3(
						Lodia_shifted_irregular_xy(i, j, 0, n.x, n.y, n.z, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del),
						Lodia_shifted_irregular_xz_yz(i, 0, j, n.x, n.z, n.y, d.x*d.y*d.z, g_vals_xz, g_vals_xz_del),
						Lodia_shifted_irregular_xz_yz(j, 0, i, n.y, n.z, n.x, d.x*d.y*d.z, g_vals_yz, g_vals_yz_del)) * sign;
				}

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