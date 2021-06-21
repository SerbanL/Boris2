#include "stdafx.h"
#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

//---------------------SHIFTED VERSIONS (FOR STRAY FIELD)

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
bool DemagTFuncCUDA::CalcDiagTens3D_Shifted(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted(n, hRatios, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted(n, hRatios, shift)) return false;
	}

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticDiag_yy()->setup(hRatios.y, hRatios.x, hRatios.z);
	demagAsymptoticDiag_zz()->setup(hRatios.z, hRatios.y, hRatios.x);

	int sign = 1;
	if (minus) sign = -1;

	CalcTens3D_Shifted_Ldia(D11, D22, D33, n, N, hRatios, shift, sign, asymptotic_distance);

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DemagTFuncCUDA::CalcOffDiagTens3D_Shifted(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D12.set(0.0);
	D13.set(0.0);
	D23.set(0.0);

	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted(n, hRatios, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted(n, hRatios, shift)) return false;
	}

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticOffDiag_xz()->setup(hRatios.x, hRatios.z, hRatios.y);
	demagAsymptoticOffDiag_yz()->setup(hRatios.y, hRatios.z, hRatios.x);

	int sign = 1;
	if (minus) sign = -1;

	CalcTens3D_Shifted_Lodia(D12, D13, D23, n, N, hRatios, shift, sign, asymptotic_distance);

	return true;
}

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DemagTFuncCUDA::CalcDiagTens2D_Shifted(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted(n, hRatios, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted(n, hRatios, shift)) return false;
	}

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticDiag_yy()->setup(hRatios.y, hRatios.x, hRatios.z);
	demagAsymptoticDiag_zz()->setup(hRatios.z, hRatios.y, hRatios.x);

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Shifted_Ldia(D11, D22, D33, n, N, hRatios, shift, sign, asymptotic_distance);

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy only) which has sizes given by N. 
bool DemagTFuncCUDA::CalcOffDiagTens2D_Shifted(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
	cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D12.set(0.0);
	D13.set(0.0);
	D23.set(0.0);

	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted(n, hRatios, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted(n, hRatios, shift)) return false;
	}

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(hRatios.x, hRatios.y, hRatios.z);
	demagAsymptoticOffDiag_xz()->setup(hRatios.x, hRatios.z, hRatios.y);
	demagAsymptoticOffDiag_yz()->setup(hRatios.y, hRatios.z, hRatios.x);

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Shifted_Lodia(D12, D13, D23, n, N, hRatios, shift, sign, asymptotic_distance);

	return true;
}

#endif