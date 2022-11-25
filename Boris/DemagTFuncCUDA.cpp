#include "stdafx.h"
#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

DemagTFuncCUDA::DemagTFuncCUDA(void)
{
}

//---------------------ZERO SHIFT VERSIONS (FOR INTERNAL FIELD)

bool DemagTFuncCUDA::CalcDiagTens2D(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool include_self_demag, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	// Set-up and calculate f values
	if (!fill_f_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticDiag_yy()->setup(hRatios.y, hRatios.x, hRatios.z);
	demagAsymptoticDiag_zz()->setup(hRatios.z, hRatios.y, hRatios.x);

	//Calculate demag tensor diagonal elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Ldia(D11, D22, D33, n, N, hRatios, sign, asymptotic_distance);

	if (!include_self_demag) {

		D11.setvalue(0, 0.0);
		D22.setvalue(0, 0.0);
		D33.setvalue(0, 0.0);
	}

	return true;
}

bool DemagTFuncCUDA::CalcOffDiagTens2D(
	cu_arr<double>& D12,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D12.set(0.0);

	// Set-up and calculate g values
	if (!fill_g2D_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(hRatios.x, hRatios.y, hRatios.z);

	//Calculate demag tensor elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Lodia(D12, n, N, hRatios, sign, asymptotic_distance);

	return true;
}

//Compute diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
//Need mesh dimensions n, convolution mesh dimensions N (this must be a power of 2 and N/2 smallest integer >= n for all dimensions -> from n to N/2 pad with zeroes).
//Need cellsize hRatios -> these can be normalized (e.g. to have largest value 1)
//You can opt to set self-demag coefficients to zero (include_self_demag = false);
bool DemagTFuncCUDA::CalcDiagTens3D(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool include_self_demag, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	// Set-up and calculate f values
	if (!fill_f_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticDiag_yy()->setup(hRatios.y, hRatios.x, hRatios.z);
	demagAsymptoticDiag_zz()->setup(hRatios.z, hRatios.y, hRatios.x);

	//Calculate demag tensor diagonal elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens3D_Ldia(D11, D22, D33, n, N, hRatios, sign, asymptotic_distance);

	if (!include_self_demag) {

		D11.setvalue(0, 0.0);
		D22.setvalue(0, 0.0);
		D33.setvalue(0, 0.0);
	}

	return true;
}

//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DemagTFuncCUDA::CalcOffDiagTens3D(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D12.set(0.0);
	D13.set(0.0);
	D23.set(0.0);

	// Set-up and calculate g values
	if (!fill_g_vals(n, hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticOffDiag_xz()->setup(hRatios.x, hRatios.z, hRatios.y);
	demagAsymptoticOffDiag_yz()->setup(hRatios.y, hRatios.z, hRatios.x);

	//Calculate demag tensor elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens3D_Lodia(D12, D13, D23, n, N, hRatios, sign, asymptotic_distance);

	return true;
}

#endif