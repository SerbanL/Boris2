#include "stdafx.h"
#include "DemagTFunc.h"

DemagTFunc::DemagTFunc(void)
{
	int OmpThreads = omp_get_num_procs();

	std::vector<double> main_sum_sp(27, 0);
	std::vector<double> irregular_sum_sp(36, 0);
	std::vector<double> f_sum_sp(4, 0);
	std::vector<double> g_sum_sp(7, 0);

	main_sum.assign(OmpThreads, main_sum_sp);
	irregular_sum.assign(OmpThreads, irregular_sum_sp);
	f_sum.assign(OmpThreads, f_sum_sp);
	g_sum.assign(OmpThreads, g_sum_sp);
}

//---------------------ZERO SHIFT VERSIONS (FOR INTERNAL FIELD)

bool DemagTFunc::CalcDiagTens3D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, bool include_self_demag, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	// Set-up and calculate f values
	if (!fill_f_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(hRatios.y, hRatios.x, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(hRatios.z, hRatios.y, hRatios.x);

	//Calculate demag tensor diagonal elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;
	int K = N.z;

	// in the first octant only need to generate values up to n, not N / 2.
	//n may not be a power of 2, in which case any points between n and N/2 will be padded with zeroes, so we don't want the tensor coefficients from n to N/2 - keep them zero.
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val;

				//apply asymptotic equations?
				//Must include the separate i, j, k checks otherwise the i*i + j*j + k*k expression can overflow the 4 byte integer range
				if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

					val = DBL3(
						demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
						demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, k * hRatios.z) * sign,
						demagAsymptoticDiag_zz.AsymptoticLdia(k * hRatios.z, j * hRatios.y, i * hRatios.x) * sign);
				}
				else {

					val = DBL3(
						Ldia(i, j, k, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
						Ldia(j, i, k, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
						Ldia(k, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
				}

				Ddiag[((i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val;
				Ddiag[((i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val;
			}
		}
	}

	if (!include_self_demag) Ddiag[0] = DBL3();

	//free memory
	f_vals_xx.clear();
	f_vals_yy.clear();
	f_vals_zz.clear();

	return true;
}

bool DemagTFunc::CalcOffDiagTens3D(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	Dodiag.set(DBL3());

	// Set-up and calculate g values
	if (!fill_g_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xz(hRatios.x, hRatios.z, hRatios.y);
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_yz(hRatios.y, hRatios.z, hRatios.x);

	//Calculate demag tensor elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;
	int K = N.z;

	//in the first octant only need to generate values up to n, not N/2.
	//n may not be a power of 2, in which case any points between n and N/2 will be padded with zeroes, so we don't want the tensor coefficients from n to N/2 - keep them zero.
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val;

				//apply asymptotic equations?
				//Must include the separate i, j, k checks otherwise the i*i + j*j + k*k expression can overflow the 4 byte integer range
				if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || k >= asymptotic_distance || i * i + j * j + k * k >= asymptotic_distance * asymptotic_distance)) {

					//D12, D13, D23
					val = DBL3(
						demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, k * hRatios.z) * sign,
						demagAsymptoticOffDiag_xz.AsymptoticLodia(i * hRatios.x, k * hRatios.z, j * hRatios.y) * sign,
						demagAsymptoticOffDiag_yz.AsymptoticLodia(j * hRatios.y, k * hRatios.z, i * hRatios.x) * sign);
				}
				else {

					//D12, D13, D23
					val = DBL3(
						Lodia(i, j, k, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign,
						Lodia(i, k, j, hRatios.x, hRatios.z, hRatios.y, g_vals_xz) * sign,
						Lodia(j, k, i, hRatios.y, hRatios.z, hRatios.x, g_vals_yz) * sign);
				}

				Dodiag[((i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(+1, +1, +1);
				Dodiag[((-i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(-1, -1, +1);
				Dodiag[((i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(-1, +1, -1);
				Dodiag[((i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(+1, -1, -1);
				Dodiag[((-i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(+1, -1, -1);
				Dodiag[((-i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(-1, +1, -1);
				Dodiag[((i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(-1, -1, +1);
				Dodiag[((-i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(+1, +1, +1);
			}
		}
	}

	//free memory
	g_vals_xy.clear();
	g_vals_xz.clear();
	g_vals_yz.clear();

	return true;
}

bool DemagTFunc::CalcDiagTens2D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, bool include_self_demag, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	// Set-up and calculate f values
	if (!fill_f_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(hRatios.x, hRatios.y, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(hRatios.y, hRatios.x, hRatios.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(hRatios.z, hRatios.y, hRatios.x);

	//Calculate demag tensor diagonal elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			DBL3 val;

			//apply asymptotic equations?
			//Must include the separate i, j checks otherwise the i*i + j*j expression can overflow the 4 byte integer range
			if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || i * i + j * j >= asymptotic_distance * asymptotic_distance)) {

				val = DBL3(
					demagAsymptoticDiag_xx.AsymptoticLdia(i * hRatios.x, j * hRatios.y, 0) * sign,
					demagAsymptoticDiag_yy.AsymptoticLdia(j * hRatios.y, i * hRatios.x, 0) * sign,
					demagAsymptoticDiag_zz.AsymptoticLdia(0, j * hRatios.y, i * hRatios.x) * sign);
			}
			else {

				val = DBL3(
					Ldia(i, j, 0, hRatios.x, hRatios.y, hRatios.z, f_vals_xx) * sign,
					Ldia(j, i, 0, hRatios.y, hRatios.x, hRatios.z, f_vals_yy) * sign,
					Ldia(0, j, i, hRatios.z, hRatios.y, hRatios.x, f_vals_zz) * sign);
			}

			Ddiag[((i + I) % I) + ((j + J) % J)*I] = val;
			Ddiag[((-i + I) % I) + ((j + J) % J)*I] = val;
			Ddiag[((i + I) % I) + ((-j + J) % J)*I] = val;
			Ddiag[((-i + I) % I) + ((-j + J) % J)*I] = val;
		}
	}

	if (!include_self_demag) Ddiag[0] = DBL3();

	//free memory
	f_vals_xx.clear();
	f_vals_yy.clear();
	f_vals_zz.clear();

	return true;
}

bool DemagTFunc::CalcOffDiagTens2D(std::vector<double> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	std::fill(Dodiag.begin(), Dodiag.end(), 0.0);

	// Set-up and calculate g values
	if (!fill_g2D_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	DemagAsymptoticOffDiag demagAsymptoticOffDiag_xy(hRatios.x, hRatios.y, hRatios.z);

	//Calculate demag tensor elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

	int I = N.x;
	int J = N.y;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			double val;

			//apply asymptotic equations?
			//Must include the separate i, j checks otherwise the i*i + j*j expression can overflow the 4 byte integer range
			if (asymptotic_distance > 0 && (i >= asymptotic_distance || j >= asymptotic_distance || i * i + j * j >= asymptotic_distance * asymptotic_distance)) {

				val = demagAsymptoticOffDiag_xy.AsymptoticLodia(i * hRatios.x, j * hRatios.y, 0) * sign;
			}
			else {

				val = Lodia(i, j, 0, hRatios.x, hRatios.y, hRatios.z, g_vals_xy) * sign;
			}

			Dodiag[((i + I) % I) + ((j + J) % J)*I] = val;
			Dodiag[((-i + I) % I) + ((j + J) % J)*I] = -1 * val;
			Dodiag[((i + I) % I) + ((-j + J) % J)*I] = -1 * val;
			Dodiag[((-i + I) % I) + ((-j + J) % J)*I] = val;
		}
	}

	//free memory
	g_vals_xy.clear();

	return true;
}