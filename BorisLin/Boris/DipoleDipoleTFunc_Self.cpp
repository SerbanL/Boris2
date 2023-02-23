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

DBL3 DipoleDipoleTFunc::SelfDemag_PBC(DBL3 h, DBL3 n, INT3 demag_pbc_images)
{
	DBL3 val = DBL3();

	for (int i_img = -demag_pbc_images.x; i_img < demag_pbc_images.x + 1; i_img++) {
		for (int j_img = -demag_pbc_images.y; j_img < demag_pbc_images.y + 1; j_img++) {
			for (int k_img = -demag_pbc_images.z; k_img < demag_pbc_images.z + 1; k_img++) {

				int i = i_img * n.x;
				int j = j_img * n.y;
				int k = k_img * n.z;

				if (i == 0 && j == 0 && k == 0) {

					val += DBL3(
						SelfDemag(h.x, h.y, h.z),
						SelfDemag(h.y, h.x, h.z),
						SelfDemag(h.z, h.y, h.x));
				}
				else {

					//displacement vector
					DBL3 r = DBL3(i * h.x, j * h.y, k * h.z);

					//length of displacement vector
					double r_norm = r.norm();

					//prefactor
					double c = MUB / (4 * PI * r_norm * r_norm * r_norm);

					//unit displacement vector
					DBL3 r_dir = r / r_norm;

					//D11, D22, D33
					val += DBL3(3 * r_dir.x*r_dir.x - 1.0, 3 * r_dir.y*r_dir.y - 1.0, 3 * r_dir.z*r_dir.z - 1.0) * c;
				}
			}
		}
	}

	return val;
}