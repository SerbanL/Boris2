#include "stdafx.h"
#include "DipoleDipoleTFunc.h"

double DipoleDipoleTFunc::SelfDemag(double hx, double hy, double hz)
{
	std::vector<double> main_sum(27, 0), f_sum(4, 0);

	//Newell f function
	auto f = [&](double x, double y, double z) -> double {

		double xsq = x * x;
		double ysq = y * y;
		double zsq = z * z;

		double Rxy_sq = xsq + ysq;
		double Rxz_sq = xsq + zsq;

		double R = sqrt(xsq + ysq + zsq);

		f_sum[0] = (2.0 * xsq - ysq - zsq) * R / 6.0;

		f_sum[1] = 0;
		if (IsNZ(Rxz_sq)) {

			double arg = 2.0 * y * (y + R) / Rxz_sq;
			if (arg > -1.0) f_sum[1] = y * (zsq - xsq) * log1p(arg) / 4.0;
		}

		f_sum[2] = 0;
		if (IsNZ(Rxy_sq)) {

			double arg = 2.0 * z * (z + R) / Rxy_sq;
			if (arg > -1.0) f_sum[2] = z * (ysq - xsq) * log1p(arg) / 4.0;
		}

		if (IsNZ(x)) f_sum[3] = -x * y*z*atan(y*z / (x*R));
		else f_sum[3] = 0;

		return sum_KahanNeumaier(f_sum);
	};

	main_sum[0] = +8 * f(0, 0, 0);

	main_sum[1] = -4 * f(hx, 0, 0);
	main_sum[2] = -4 * f(-hx, 0, 0);
	main_sum[3] = -4 * f(0, hy, 0);
	main_sum[4] = -4 * f(0, -hy, 0);
	main_sum[5] = -4 * f(0, 0, hz);
	main_sum[6] = -4 * f(0, 0, -hz);

	main_sum[7] = +2 * f(-hx, -hy, 0);
	main_sum[8] = +2 * f(-hx, hy, 0);
	main_sum[9] = +2 * f(hx, -hy, 0);
	main_sum[10] = +2 * f(hx, hy, 0);

	main_sum[11] = +2 * f(-hx, 0, -hz);
	main_sum[12] = +2 * f(-hx, 0, hz);
	main_sum[13] = +2 * f(hx, 0, -hz);
	main_sum[14] = +2 * f(hx, 0, hz);

	main_sum[15] = +2 * f(0, -hy, -hz);
	main_sum[16] = +2 * f(0, -hy, hz);
	main_sum[17] = +2 * f(0, hy, -hz);
	main_sum[18] = +2 * f(0, hy, hz);

	main_sum[19] = -1 * f(-hx, -hy, -hz);
	main_sum[20] = -1 * f(-hx, -hy, hz);
	main_sum[21] = -1 * f(-hx, hy, -hz);
	main_sum[22] = -1 * f(hx, -hy, -hz);
	main_sum[23] = -1 * f(-hx, hy, hz);
	main_sum[24] = -1 * f(hx, -hy, hz);
	main_sum[25] = -1 * f(hx, hy, -hz);
	main_sum[26] = -1 * f(hx, hy, hz);

	return sum_KahanNeumaier(main_sum) / (4 * PI * hx * hy * hz);
}

DBL3 DipoleDipoleTFunc::SelfDemag(DBL3 h)
{
	return DBL3(
		SelfDemag(h.x, h.y, h.z),
		SelfDemag(h.y, h.x, h.z),
		SelfDemag(h.z, h.y, h.x));
}