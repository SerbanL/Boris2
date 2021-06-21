#include "stdafx.h"
#include "DemagTFunc.h"

//---------------------SHIFTED VERSIONS (FOR STRAY FIELD) with PBC

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
bool DemagTFunc::CalcDiagTens3D_Shifted_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, DBL3 shift, 
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {
		
		//z shift, and no pbc along z
		if (!fill_f_vals_zshifted(SZ3(
			(x_images ? asymptotic_distance : N.x / 2),
			(y_images ? asymptotic_distance : N.y / 2),
			N.z / 2), hRatios, shift, asymptotic_distance)) return false;	
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(hRatios.y, hRatios.x, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(hRatios.z, hRatios.y, hRatios.x);

	int sign = 1;
	if (minus) sign = -1;

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift only, so use full xy plane symmetries. no pbc along z

#pragma omp parallel for
		for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int k0 = -N.z / 2 + 1; k0 < N.z / 2; k0++) {
				for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

					DBL3 val = DBL3();

					for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
						for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

							//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
							//There is significant redundancy remaining, could be optimized further.
							int i = mod(i_img * N.x + i0);
							int j = mod(j_img * N.y + j0);
							int k = k0;

							double ks = mod(k + shift.z / hRatios.z);

							//apply asymptotic equations?
							if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

								//D12, D13, D23
								val += DBL3(
									demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z) * sign,
									demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z + shift.z) * sign,
									demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y, i * hRatios.x) * sign);
							}
							else {
	
								val += DBL3(
									Ldia_zshifted_xx_yy(i, j, k, (x_images ? N.x : N.x / 2), (y_images ? N.y : N.y / 2), N.z / 2, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
									Ldia_zshifted_xx_yy(j, i, k, (y_images ? N.y : N.y / 2), (x_images ? N.x : N.x / 2), N.z / 2, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
									Ldia_zshifted_zz(k, j, i, N.z / 2, (y_images ? N.y : N.y / 2), (x_images ? N.x : N.x / 2), hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
							}
						}
					}

					//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
					Ddiag[INT3(i0, j0, (k0 + N.z) % N.z)] = val;

					if (!x_images) {

						if (i0) Ddiag[INT3((N.x - i0) % N.x, j0, (k0 + N.z) % N.z)] = val;

						if (!y_images) {

							if (i0 && j0) Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, (k0 + N.z) % N.z)] = val;
						}
					}

					if (!y_images) {

						if (j0) Ddiag[INT3(i0, (N.y - j0) % N.y, (k0 + N.z) % N.z)] = val;
					}
				}
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible,

#pragma omp parallel for
		for (int j0 = (y_images ? 0 : -N.y / 2 + 1); j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int k0 = (z_images ? 0 : -N.z / 2 + 1); k0 < (z_images ? N.z : N.z / 2); k0++) {
				for (int i0 = (x_images ? 0 : -N.x / 2 + 1); i0 < (x_images ? N.x : N.x / 2); i0++) {

					DBL3 val = DBL3();

					for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
						for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
							for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

								int i = i_img * N.x + i0;
								int j = j_img * N.y + j0;
								int k = k_img * N.z + k0;

								double is = mod(i + (shift.x / hRatios.x));
								double js = mod(j + (shift.y / hRatios.y));
								double ks = mod(k + (shift.z / hRatios.z));

								//apply asymptotic equations?
								if (asymptotic_distance > 0 &&
									(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
									int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

									//D12, D13, D23
									val += DBL3(
										demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z) * sign,
										demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z) * sign,
										demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x) * sign);
								}
								else {

									val += DBL3(
										Ldia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
										Ldia_single(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z, hRatios.y, hRatios.x, hRatios.z),
										Ldia_single(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x, hRatios.z, hRatios.y, hRatios.x)) * sign;
								}
							}
						}
					}

					//no symmetries used to speedup tensor calculation
					Ddiag[INT3((i0 + N.x) % N.x, (j0 + N.y) % N.y, (k0 + N.z) % N.z)] = val;
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
bool DemagTFunc::CalcOffDiagTens3D_Shifted_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 hRatios, DBL3 shift, 
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Dodiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift, and no pbc along z
		if (!fill_g_vals_zshifted(SZ3(
			(x_images ? asymptotic_distance : N.x / 2),
			(y_images ? asymptotic_distance : N.y / 2),
			N.z / 2), hRatios, shift, asymptotic_distance)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xz(hRatios.x, hRatios.z, hRatios.y);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_yz(hRatios.y, hRatios.z, hRatios.x);

	int sign = 1;
	if (minus) sign = -1;

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift only, so use full xy plane symmetries. no pbc along z.

#pragma omp parallel for
		for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int k0 = -N.z / 2 + 1; k0 < N.z / 2; k0++) {
				for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

					DBL3 val = DBL3();

					for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
						for (int j_img = -y_images; j_img < y_images + 1; j_img++) {

							int i = i_img * N.x + i0;
							int j = j_img * N.y + j0;
							int k = k0;

							//use modulus of indexes and adjust for tensor element signs based on symmetries
							int sign_i = get_sign(i);
							int sign_j = get_sign(j);
							i = mod(i);
							j = mod(j);

							double ks = mod(k + shift.z / hRatios.z);

							//apply asymptotic equations?
							if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

								//D12, D13, D23
								val += DBL3(
									demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z + shift.z) * sign,
									demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z + shift.z, j * hRatios.y) * sign,
									demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z + shift.z, i * hRatios.x) * sign)
									& DBL3(sign_i * sign_j, sign_i, sign_j);
							}
							else {

								val += DBL3(
									Lodia_xy_zshifted(i, j, k, (x_images ? N.x : N.x / 2), (y_images ? N.y : N.y / 2), N.z / 2, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
									Lodia_xz_yz_zshifted(i, k, j, (x_images ? N.x : N.x / 2), N.z / 2, (y_images ? N.y : N.y / 2), hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
									Lodia_xz_yz_zshifted(j, k, i, (y_images ? N.y : N.y / 2), N.z / 2, (x_images ? N.x : N.x / 2), hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign)
									& DBL3(sign_i * sign_j, sign_i, sign_j);
							}
						}
					}

					//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
					Dodiag[INT3(i0, j0, (k0 + N.z) % N.z)] = val;

					if (!x_images) {

						if (i0) Dodiag[INT3((N.x - i0) % N.x, j0, (k0 + N.z) % N.z)] = val & DBL3(-1, -1, +1);

						if (!y_images) {

							if (i0 && j0) Dodiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, (k0 + N.z) % N.z)] = val & DBL3(+1, -1, -1);
						}
					}

					if (!y_images) {

						if (j0) Dodiag[INT3(i0, (N.y - j0) % N.y, (k0 + N.z) % N.z)] = val & DBL3(-1, +1, -1);
					}
				}
			}
		}
	}
	else {

		//general case : no symmetries used (although it's still possible to use some restricted symmetries but not really worth bothering with). pbc along z possible,

#pragma omp parallel for
		for (int j0 = (y_images ? 0 : -N.y / 2 + 1); j0 < (y_images ? N.y : N.y / 2); j0++) {
			for (int k0 = (z_images ? 0 : -N.z / 2 + 1); k0 < (z_images ? N.z : N.z / 2); k0++) {
				for (int i0 = (x_images ? 0 : -N.x / 2 + 1); i0 < (x_images ? N.x : N.x / 2); i0++) {

					DBL3 val = DBL3();

					for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
						for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
							for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

								int i = i_img * N.x + i0;
								int j = j_img * N.y + j0;
								int k = k_img * N.z + k0;

								double is = mod(i + (shift.x / hRatios.x));
								double js = mod(j + (shift.y / hRatios.y));
								double ks = mod(k + (shift.z / hRatios.z));

								//apply asymptotic equations?
								if (asymptotic_distance > 0 &&
									(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
									int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

									//D12, D13, D23
									val += DBL3(
										demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z),
										demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y),
										demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x)) * sign;
								}
								else {

									val += DBL3(
										Lodia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
										Lodia_single(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y, hRatios.x, hRatios.z, hRatios.y),
										Lodia_single(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x, hRatios.y, hRatios.z, hRatios.x)) * sign;
								}
							}
						}
					}

					//no symmetries used to speedup tensor calculation
					Dodiag[INT3((i0 + N.x) % N.x, (j0 + N.y) % N.y, (k0 + N.z) % N.z)] = val;
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
bool DemagTFunc::CalcDiagTens2D_Shifted_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, DBL3 shift, 
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift, and no pbc along z
		if (!fill_f_vals_zshifted(SZ3(
			(x_images ? asymptotic_distance : N.x / 2),
			(y_images ? asymptotic_distance : N.y / 2),
			1), hRatios, shift, asymptotic_distance)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(hRatios.y, hRatios.x, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(hRatios.z, hRatios.y, hRatios.x);

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

						double ks = mod(shift.z / hRatios.z);

						//apply asymptotic equations?
						if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

							//D12, D13, D23
							val += DBL3(
								demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, shift.z) * sign,
								demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, shift.z) * sign,
								demagAsymptoticDiag_zz.AsymptoticLdia(shift.z, j * hRatios.y, i * hRatios.x) * sign);
						}
						else {

							val += DBL3(
								Ldia_zshifted_xx_yy(i, j, 0, (x_images ? N.x : N.x / 2), (y_images ? N.y : N.y / 2), 1, hRatios.x, hRatios.y, hRatios.z, f_vals_xx),
								Ldia_zshifted_xx_yy(j, i, 0, (y_images ? N.y : N.y / 2), (x_images ? N.x : N.x / 2), 1, hRatios.y, hRatios.x, hRatios.z, f_vals_yy),
								Ldia_zshifted_zz(0, j, i, 1, (y_images ? N.y : N.y / 2), (x_images ? N.x : N.x / 2), hRatios.z, hRatios.y, hRatios.x, f_vals_zz)) * sign;
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

							double is = mod(i + (shift.x / hRatios.x));
							double js = mod(j + (shift.y / hRatios.y));
							double ks = mod(k + (shift.z / hRatios.z));

							//apply asymptotic equations?
							if (asymptotic_distance > 0 &&
								(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
								int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

								//D12, D13, D23
								val += DBL3(
									demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z) * sign,
									demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z) * sign,
									demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x) * sign);
							}
							else {

								val += DBL3(
									Ldia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
									Ldia_single(j * hRatios.y + shift.y, i * hRatios.x + shift.x, k * hRatios.z + shift.z, hRatios.y, hRatios.x, hRatios.z),
									Ldia_single(k * hRatios.z + shift.z, j * hRatios.y + shift.y, i * hRatios.x + shift.x, hRatios.z, hRatios.y, hRatios.x)) * sign;
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

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy only) which has sizes given by N. 
bool DemagTFunc::CalcOffDiagTens2D_Shifted_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 hRatios, DBL3 shift, 
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Dodiag.set(DBL3());

	if (IsZ(shift.x) && IsZ(shift.y) && !z_images) {

		//z shift, and no pbc along z
		if (!fill_g_vals_zshifted(SZ3(
			(x_images ? asymptotic_distance : N.x / 2),
			(y_images ? asymptotic_distance : N.y / 2),
			1), hRatios, shift, asymptotic_distance)) return false;
	}

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xz(hRatios.x, hRatios.z, hRatios.y);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_yz(hRatios.y, hRatios.z, hRatios.x);

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

						double ks = mod(shift.z / hRatios.z);

						//apply asymptotic equations?
						if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance || int(floor_epsilon(i * i + j * j + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

							//D12, D13, D23
							val += DBL3(
								demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, shift.z) * sign,
								demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, shift.z, j * hRatios.y) * sign,
								demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, shift.z, i * hRatios.x) * sign)
								& DBL3(sign_i * sign_j, sign_i, sign_j);
						}
						else {

							val += DBL3(
								Lodia_xy_zshifted(i, j, 0, (x_images ? N.x : N.x / 2), (y_images ? N.y : N.y / 2), 1, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
								Lodia_xz_yz_zshifted(i, 0, j, (x_images ? N.x : N.x / 2), 1, (y_images ? N.y : N.y / 2), hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
								Lodia_xz_yz_zshifted(j, 0, i, (y_images ? N.y : N.y / 2), 1, (x_images ? N.x : N.x / 2), hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign)
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

							double is = mod(i + (shift.x / hRatios.x));
							double js = mod(j + (shift.y / hRatios.y));
							double ks = mod(k + (shift.z / hRatios.z));

							//apply asymptotic equations?
							if (asymptotic_distance > 0 &&
								(int(floor_epsilon(is)) >= asymptotic_distance || int(floor_epsilon(js)) >= asymptotic_distance || int(floor_epsilon(ks)) >= asymptotic_distance ||
								int(floor_epsilon(is * is + js * js + ks * ks)) >= asymptotic_distance * asymptotic_distance)) {

								//D12, D13, D23
								val += DBL3(
									demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z),
									demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y),
									demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x)) * sign;
							}
							else {

								val += DBL3(
									Lodia_single(i * hRatios.x + shift.x, j * hRatios.y + shift.y, k * hRatios.z + shift.z, hRatios.x, hRatios.y, hRatios.z),
									Lodia_single(i * hRatios.x + shift.x, k * hRatios.z + shift.z, j * hRatios.y + shift.y, hRatios.x, hRatios.z, hRatios.y),
									Lodia_single(j * hRatios.y + shift.y, k * hRatios.z + shift.z, i * hRatios.x + shift.x, hRatios.y, hRatios.z, hRatios.x)) * sign;
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

	return true;
}