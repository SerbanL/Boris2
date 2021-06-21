#include "stdafx.h"
#include "DemagTFunc.h"

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DemagTFunc::CalcDiagTens2D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, bool minus, int asymptotic_distance, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	// Set-up and calculate f values
	if (!fill_f_vals(SZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : 1)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(hRatios.y, hRatios.x, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(hRatios.z, hRatios.y, hRatios.x);

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

			DBL3 val = DBL3();

			for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
				for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
					for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

						//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
						//There is significant redundancy remaining, could be optimized further.
						int i = abs(i_img * N.x + i0);
						int j = abs(j_img * N.y + j0);
						int k = abs(k_img);

						//apply asymptotic equations?
						if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

							//D12, D13, D23
							val += DBL3(
								demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
								demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z) * sign,
								demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z, j * hRatios.y, i * hRatios.x) * sign);
						}
						else {

							val += DBL3(
								Ldia(i, j, k, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
								Ldia(j, i, k, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
								Ldia(k, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
						}
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

	//free memory
	f_vals_xx.clear();
	f_vals_yy.clear();
	f_vals_zz.clear();

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DemagTFunc::CalcOffDiagTens2D_PBC(std::vector<double> &Dodiag, INT3 N, DBL3 hRatios, bool minus, int asymptotic_distance, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	std::fill(Dodiag.begin(), Dodiag.end(), 0.0);

	// Set-up and calculate f values
	if (!fill_g2D_vals(SZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : 1)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(hRatios.x, hRatios.y, hRatios.z);

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

			double val = 0.0;

			for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
				for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
					for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

						int i = i_img * N.x + i0;
						int j = j_img * N.y + j0;
						int k = k_img;

						//use modulus of indexes and adjust for tensor element signs based on symmetries
						int sign_i = get_sign(i);
						int sign_j = get_sign(j);
						i = abs(i);
						j = abs(j);
						k = abs(k);

						//apply asymptotic equations?
						if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

							//D12, D13, D23
							val += demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign * sign_i * sign_j;
						}
						else {

							val += Lodia(i, j, k, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign * sign_i * sign_j;
						}
					}
				}
			}

			//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
			Dodiag[i0 + j0 * N.x] = val;

			if (!x_images) {

				if (i0) Dodiag[((N.x - i0) % N.x) + j0 * N.x] = -val;

				if (!y_images) {

					if (i0 && j0) Dodiag[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] = val;
				}
			}

			if (!y_images) {

				if (j0) Dodiag[i0 + ((N.y - j0) % N.y) * N.x] = -val;
			}
		}
	}

	//free memory
	g_vals_xy.clear();

	return true;
}

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DemagTFunc::CalcDiagTens3D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, bool minus, int asymptotic_distance, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	// Set-up and calculate f values
	if (!fill_f_vals(SZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : N.z / 2)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(hRatios.y, hRatios.x, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(hRatios.z, hRatios.y, hRatios.x);

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int k0 = 0; k0 < (z_images ? N.z : N.z / 2); k0++) {
			for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
						for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

							//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
							//There is significant redundancy remaining, could be optimized further.
							int i = abs(i_img * N.x + i0);
							int j = abs(j_img * N.y + j0);
							int k = abs(k_img * N.z + k0);

							//apply asymptotic equations?
							if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

								//D12, D13, D23
								val += DBL3(
									demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
									demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z) * sign,
									demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z, j * hRatios.y, i * hRatios.x) * sign);
							}
							else {

								val += DBL3(
									Ldia(i, j, k, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
									Ldia(j, i, k, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
									Ldia(k, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
							}
						}
					}
				}

				//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
				Ddiag[INT3(i0, j0, k0)] = val;

				if (!x_images) {

					if (i0) Ddiag[INT3((N.x - i0) % N.x, j0, k0)] = val;

					if (!y_images) {

						if (i0 && j0) Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, k0)] = val;

						if (!z_images) {

							if (i0 && j0 && k0) Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val;
						}
					}

					if (!z_images) {

						if (i0 && k0) Ddiag[INT3((N.x - i0) % N.x, j0, (N.z - k0) % N.z)] = val;
					}
				}

				if (!y_images) {

					if (j0) Ddiag[INT3(i0, (N.y - j0) % N.y, k0)] = val;

					if (!z_images) {

						if (j0 && k0) Ddiag[INT3(i0, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val;
					}
				}

				if (!z_images) {

					if (k0) Ddiag[INT3(i0, j0, (N.z - k0) % N.z)] = val;
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
bool DemagTFunc::CalcOffDiagTens3D_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 hRatios, bool minus, int asymptotic_distance, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Dodiag.set(DBL3());

	// Set-up and calculate f values
	if (!fill_g_vals(SZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : N.z / 2)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xz(hRatios.x, hRatios.z, hRatios.y);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_yz(hRatios.y, hRatios.z, hRatios.x);

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int k0 = 0; k0 < (z_images ? N.z : N.z / 2); k0++) {
			for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
						for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

							int i = i_img * N.x + i0;
							int j = j_img * N.y + j0;
							int k = k_img * N.z + k0;

							//use modulus of indexes and adjust for tensor element signs based on symmetries
							int sign_i = get_sign(i);
							int sign_j = get_sign(j);
							int sign_k = get_sign(k);
							i = abs(i);
							j = abs(j);
							k = abs(k);

							//apply asymptotic equations?
							if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

								//D12, D13, D23
								val += DBL3(
									demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
									demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z, j * hRatios.y) * sign,
									demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z, i * hRatios.x) * sign)
									& DBL3(sign_i * sign_j, sign_i * sign_k, sign_j * sign_k);
							}
							else {

								val += DBL3(
									Lodia(i, j, k, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
									Lodia(i, k, j, hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
									Lodia(j, k, i, hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign)
									& DBL3(sign_i * sign_j, sign_i * sign_k, sign_j * sign_k);
							}
						}
					}
				}

				//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
				Dodiag[INT3(i0, j0, k0)] = val;

				if (!x_images) {

					if (i0) Dodiag[INT3((N.x - i0) % N.x, j0, k0)] = val & DBL3(-1, -1, +1);

					if (!y_images) {

						if (i0 && j0) Dodiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, k0)] = val & DBL3(+1, -1, -1);

						if (!z_images) {

							if (i0 && j0 && k0) Dodiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val;
						}
					}

					if (!z_images) {

						if (i0 && k0) Dodiag[INT3((N.x - i0) % N.x, j0, (N.z - k0) % N.z)] = val & DBL3(-1, +1, -1);
					}
				}

				if (!y_images) {

					if (j0) Dodiag[INT3(i0, (N.y - j0) % N.y, k0)] = val & DBL3(-1, +1, -1);

					if (!z_images) {

						if (j0 && k0) Dodiag[INT3(i0, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val & DBL3(-1, -1, +1);
					}
				}

				if (!z_images) {

					if (k0) Dodiag[INT3(i0, j0, (N.z - k0) % N.z)] = val & DBL3(+1, -1, -1);
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