#include "stdafx.h"
#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

//---------------------ZERO SHIFT VERSIONS (FOR INTERNAL FIELD)

bool DemagTFuncCUDA::CalcDiagTens2D_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 hRatios, bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	// Set-up and calculate f values
	if (!fill_f_vals(cuSZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : 1)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticDiag_yy()->setup(hRatios.y, hRatios.x, hRatios.z);
	demagAsymptoticDiag_zz()->setup(hRatios.z, hRatios.y, hRatios.x);

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Ldia_PBC(D11, D22, D33, N, hRatios, sign, asymptotic_distance, x_images, y_images, z_images);

	return true;
}

bool DemagTFuncCUDA::CalcOffDiagTens2D_PBC(
	cu_arr<double>& D12,
	cuINT3 N, cuDBL3 hRatios,
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	D12.set(0.0);

	// Set-up and calculate f values
	if (!fill_g2D_vals(cuSZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : 1)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(hRatios.x, hRatios.y, hRatios.z);

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Lodia_PBC(D12, N, hRatios, sign, asymptotic_distance, x_images, y_images, z_images);

	return true;
}

//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DemagTFuncCUDA::CalcDiagTens3D_PBC(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 N, cuDBL3 hRatios,
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	// Set-up and calculate f values
	if (!fill_f_vals(cuSZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : N.z / 2)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticDiag_yy()->setup(hRatios.y, hRatios.x, hRatios.z);
	demagAsymptoticDiag_zz()->setup(hRatios.z, hRatios.y, hRatios.x);

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens3D_Ldia_PBC(D11, D22, D33, N, hRatios, sign, asymptotic_distance, x_images, y_images, z_images);

	return true;
}

//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DemagTFuncCUDA::CalcOffDiagTens3D_PBC(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 N, cuDBL3 hRatios,
	bool minus, int asymptotic_distance,
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but for demag we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	D12.set(0.0);
	D13.set(0.0);
	D23.set(0.0);

	// Set-up and calculate f values
	if (!fill_g_vals(cuSZ3(
		(x_images ? asymptotic_distance : N.x / 2),
		(y_images ? asymptotic_distance : N.y / 2),
		(z_images ? asymptotic_distance : N.z / 2)), hRatios, asymptotic_distance)) return false;

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticOffDiag_xz()->setup(hRatios.x, hRatios.z, hRatios.y);
	demagAsymptoticOffDiag_yz()->setup(hRatios.y, hRatios.z, hRatios.x);

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int sign = 1;
	if (minus) sign = -1;

	CalcTens3D_Lodia_PBC(D12, D13, D23, N, hRatios, sign, asymptotic_distance, x_images, y_images, z_images);

	return true;
}

#endif