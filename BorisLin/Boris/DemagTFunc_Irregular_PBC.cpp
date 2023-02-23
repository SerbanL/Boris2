#include "stdafx.h"
#include "DemagTFunc.h"

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
bool DemagTFunc::CalcDiagTens2D_Shifted_Irregular_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 s, DBL3 d, DBL3 shift,
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	//only use irregular version if you have to
	if (s * 1e-9 == d * 1e-9) return CalcDiagTens2D_Shifted_PBC(Ddiag, N, d, shift, minus, asymptotic_distance, x_images, y_images, z_images);

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift, and no pbc along z
		if (!fill_f_vals_zshifted_irregular(SZ3(
			(x_images ? asymptotic_distance : N.x / 2),
			(y_images ? asymptotic_distance : N.y / 2),
			1), s, d, shift, asymptotic_distance)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(d.x, d.y, d.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(d.y, d.x, d.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(d.z, d.y, d.x);

	int sign = 1;
	if (minus) sign = -1;

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift only, so use full xy plane symmetries. no pbc along z

#pragma omp parallel for
		for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

						//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
						//There is significant redundancy remaining, could be optimized further.
						int i = mod(i_img * N.x + i0);
						int j = mod(j_img * N.y + j0);

						//apply asymptotic equations?
						if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(mod(shift.z / d.z))) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

							//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
							//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
							val += DBL3(
								demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x, j * d.y, shift.z) * sign,
								demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y, i * d.x, shift.z) * sign,
								demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * d.y, i * d.x) * sign) * s.z;
						}
						else {

							val += DBL3(
								Ldia_zshifted_irregular_xx_yy(i, j, 0, (x_images ? N.x : N.x / 2), (y_images ? N.y : N.y / 2), 1, d.x * d.y * d.z, f_vals_xx, f_vals_xx_del),
								Ldia_zshifted_irregular_xx_yy(j, i, 0, (y_images ? N.y : N.y / 2), (x_images ? N.x : N.x / 2), 1, d.x * d.y * d.z, f_vals_yy, f_vals_yy_del),
								Ldia_zshifted_irregular_zz(0, j, i, 1, (y_images ? N.y : N.y / 2), (x_images ? N.x : N.x / 2), d.x * d.y * d.z, f_vals_zz, f_vals_zz_del)) * sign;
						}
					}
				}

				//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
				Ddiag[INT3(i0, j0, 0)] = val;

				if (!x_images) {

					if (i0) Ddiag[INT3((N.x - i0) % N.x, j0, 0)] = val;

					if (!y_images) {

						if (i0 && j0) Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, 0)] = val;
					}
				}

				if (!y_images) {

					if (j0) Ddiag[INT3(i0, (N.y - j0) % N.y, 0)] = val;
				}
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible,

#pragma omp parallel for
		for (int j0 = (y_images ? 0 : -N.y / 2 + 1); j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int i0 = (x_images ? 0 : -N.x / 2 + 1); i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
						for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

							int i = i_img * N.x + i0;
							int j = j_img * N.y + j0;
							int k = k_img;

							double is = mod(i + (shift.x / d.x));
							double js = mod(j + (shift.y / d.y));
							double ks = mod(k + (shift.z / d.z));

							//apply asymptotic equations?
							if (asymptotic_distance > 0 &&
								(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
									int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

								//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
								//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
								//D12, D13, D23
								val += DBL3(
									demagAsymptoticDiag_xx.AsymptoticLdia(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z) * sign,
									demagAsymptoticDiag_yy.AsymptoticLdia(j * d.y + shift.y, i * d.x + shift.x, k * d.z + shift.z) * sign,
									demagAsymptoticDiag_zz.AsymptoticLdia(k * d.z + shift.z, j * d.y + shift.y, i * d.x + shift.x) * sign) * s.z;
							}
							else {

								val += DBL3(
									Ldia_shifted_irregular_xx_yy_single(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z, d.x, d.y, s.z, d.z),
									Ldia_shifted_irregular_xx_yy_single(j * d.y + shift.y, i * d.x + shift.x, k * d.z + shift.z, d.y, d.x, s.z, d.z),
									Ldia_shifted_irregular_zz_single(k * d.z + shift.z, j * d.y + shift.y, i * d.x + shift.x, s.z, d.z, d.y, d.x)) * sign;
							}
						}
					}
				}

				//no symmetries used to speedup tensor calculation
				Ddiag[INT3((i0 + N.x) % N.x, (j0 + N.y) % N.y, 0)] = val;
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
bool DemagTFunc::CalcOffDiagTens2D_Shifted_Irregular_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 s, DBL3 d, DBL3 shift,
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Dodiag.set(DBL3());

	//only use irregular version if you have to
	if (s * 1e-9 == d * 1e-9) return CalcOffDiagTens2D_Shifted_PBC(Dodiag, N, s, shift, minus, asymptotic_distance, x_images, y_images, z_images);

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift, and no pbc along z
		if (!fill_g_vals_zshifted_irregular(SZ3(
			(x_images ? asymptotic_distance : N.x / 2),
			(y_images ? asymptotic_distance : N.y / 2),
			1), s, d, shift, asymptotic_distance)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(d.x, d.y, d.z);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xz(d.x, d.z, d.y);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_yz(d.y, d.z, d.x);

	int sign = 1;
	if (minus) sign = -1;

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift only, so use full xy plane symmetries. no pbc along z.

#pragma omp parallel for
		for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

						int i = i_img * N.x + i0;
						int j = j_img * N.y + j0;

						//use modulus of indexes and adjust for tensor element signs based on symmetries
						int sign_i = get_sign(i);
						int sign_j = get_sign(j);
						i = mod(i);
						j = mod(j);

						//apply asymptotic equations?
						if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(mod(shift.z / d.z))) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + (shift.z / d.z) * (shift.z / d.z))) >= asymptotic_distance * asymptotic_distance)) {

							//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
							//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
							//D12, D13, D23
							val += (DBL3(
								demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x, j * d.y, shift.z) * sign,
								demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x, shift.z, j * d.y) * sign,
								demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y, shift.z, i * d.x) * sign)
								& DBL3(sign_i * sign_j, sign_i, sign_j)) * s.z;
						}
						else {

							val += DBL3(
								Lodia_zshifted_irregular_xy(i, j, 0, (x_images ? N.x : N.x / 2), (y_images ? N.y : N.y / 2), 1, d.x*d.y*d.z, g_vals_xy, g_vals_xy_del) * sign,
								Lodia_zshifted_irregular_xz_yz(i, 0, j, (x_images ? N.x : N.x / 2), 1, (y_images ? N.y : N.y / 2), d.x*d.y*d.z, g_vals_xz, g_vals_xz_del) * sign,
								Lodia_zshifted_irregular_xz_yz(j, 0, i, (y_images ? N.y : N.y / 2), 1, (x_images ? N.x : N.x / 2), d.x*d.y*d.z, g_vals_yz, g_vals_yz_del) * sign)
								& DBL3(sign_i * sign_j, sign_i, sign_j);
						}
					}
				}

				//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
				Dodiag[INT3(i0, j0, 0)] = val;

				if (!x_images) {

					if (i0) Dodiag[INT3((N.x - i0) % N.x, j0, 0)] = val & DBL3(-1, -1, +1);

					if (!y_images) {

						if (i0 && j0) Dodiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, 0)] = val & DBL3(+1, -1, -1);
					}
				}

				if (!y_images) {

					if (j0) Dodiag[INT3(i0, (N.y - j0) % N.y, 0)] = val & DBL3(-1, +1, -1);
				}
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible,

#pragma omp parallel for
		for (int j0 = (y_images ? 0 : -N.y / 2 + 1); j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int i0 = (x_images ? 0 : -N.x / 2 + 1); i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
						for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

							int i = i_img * N.x + i0;
							int j = j_img * N.y + j0;
							int k = k_img;

							double is = mod(i + (shift.x / d.x));
							double js = mod(j + (shift.y / d.y));
							double ks = mod(k + (shift.z / d.z));

							//apply asymptotic equations?
							if (asymptotic_distance > 0 &&
								(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
									int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

								//asymptotic approximation : can apply equations for cells with same dimensions and simply multiply by source cell thickness to adjust (remember d and s can only differ in z component) - this is exact, not a further approximation
								//the demag field at the destination cell, in the asymptotic approximation regime, will scale with the source cell thickness, but not with the destination cell thickness
								//D12, D13, D23
								val += DBL3(
									demagAsymptoticOffDiag_xy.AsymptoticLodia(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z),
									demagAsymptoticOffDiag_xz.AsymptoticLodia(i * d.x + shift.x, k * d.z + shift.z, j * d.y + shift.y),
									demagAsymptoticOffDiag_yz.AsymptoticLodia(j * d.y + shift.y, k * d.z + shift.z, i * d.x + shift.x)) * sign * s.z;
							}
							else {

								val += DBL3(
									Lodia_shifted_irregular_xy_single(i * d.x + shift.x, j * d.y + shift.y, k * d.z + shift.z, d.x, d.y, s.z, d.z),
									Lodia_shifted_irregular_xz_yz_single(i * d.x + shift.x, k * d.z + shift.z, j * d.y + shift.y, d.x, s.z, d.z, d.y),
									Lodia_shifted_irregular_xz_yz_single(j * d.y + shift.y, k * d.z + shift.z, i * d.x + shift.x, d.y, s.z, d.z, d.x)) * sign;
							}
						}
					}
				}

				//no symmetries used to speedup tensor calculation
				Dodiag[INT3((i0 + N.x) % N.x, (j0 + N.y) % N.y, 0)] = val;
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