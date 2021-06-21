#include "stdafx.h"
#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
bool DemagTFuncCUDA::CalcDiagTens2D_Shifted_Irregular(
	cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D11.set(0.0);
	D22.set(0.0);
	D33.set(0.0);

	//only use irregular version if you have to
	if (s * 1e-9 == d * 1e-9) return CalcDiagTens2D_Shifted(D11, D22, D33, n, N, s, shift, minus, asymptotic_distance);
	
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift
		if (!fill_f_vals_zshifted_irregular(n, s, d, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_f_vals_shifted_irregular(n, s, d, shift)) return false;
	}

	//Setup asymptotic approximation settings
	demagAsymptoticDiag_xx()->setup(d.x, d.y, d.z);
	demagAsymptoticDiag_yy()->setup(d.y, d.x, d.z);
	demagAsymptoticDiag_zz()->setup(d.z, d.y, d.x);

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Shifted_Irregular_Ldia(D11, D22, D33, n, N, s, d, shift, sign, asymptotic_distance);

	return true;
}

//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
bool DemagTFuncCUDA::CalcOffDiagTens2D_Shifted_Irregular(
	cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23, 
	cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, bool minus, int asymptotic_distance)
{
	//zero the tensor first
	D12.set(0.0);
	D13.set(0.0);
	D23.set(0.0);

	//only use irregular version if you have to
	if (s * 1e-9 == d * 1e-9) return CalcOffDiagTens2D_Shifted(D12, D13, D23, n, N, s, shift, minus, asymptotic_distance);
	
	if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

		//z shift
		if (!fill_g_vals_zshifted_irregular(n, s, d, shift, asymptotic_distance)) return false;
	}
	else {

		//not a z shift
		if (!fill_g_vals_shifted_irregular(n, s, d, shift)) return false;
	}

	//Setup asymptotic approximation settings
	demagAsymptoticOffDiag_xy()->setup(d.x, d.y, d.z);
	demagAsymptoticOffDiag_xz()->setup(d.x, d.z, d.y);
	demagAsymptoticOffDiag_yz()->setup(d.y, d.z, d.x);

	int sign = 1;
	if (minus) sign = -1;

	CalcTens2D_Shifted_Irregular_Lodia(D12, D13, D23, n, N, s, d, shift, sign, asymptotic_distance);

	return true;
}

#endif