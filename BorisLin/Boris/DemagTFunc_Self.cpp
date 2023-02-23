#include "stdafx.h"
#include "DemagTFunc.h"

//SELF DEMAG METHODS

double DemagTFunc::SelfDemag(double hx, double hy, double hz)
{
	int tn = omp_get_thread_num();

	main_sum[tn][0] = +8 * f(0, 0, 0);

	main_sum[tn][1] = -4 * f(hx, 0, 0);
	main_sum[tn][2] = -4 * f(-hx, 0, 0);
	main_sum[tn][3] = -4 * f(0, hy, 0);
	main_sum[tn][4] = -4 * f(0, -hy, 0);
	main_sum[tn][5] = -4 * f(0, 0, hz);
	main_sum[tn][6] = -4 * f(0, 0, -hz);

	main_sum[tn][7] = +2 * f(-hx, -hy, 0);
	main_sum[tn][8] = +2 * f(-hx, hy, 0);
	main_sum[tn][9] = +2 * f(hx, -hy, 0);
	main_sum[tn][10] = +2 * f(hx, hy, 0);

	main_sum[tn][11] = +2 * f(-hx, 0, -hz);
	main_sum[tn][12] = +2 * f(-hx, 0, hz);
	main_sum[tn][13] = +2 * f(hx, 0, -hz);
	main_sum[tn][14] = +2 * f(hx, 0, hz);

	main_sum[tn][15] = +2 * f(0, -hy, -hz);
	main_sum[tn][16] = +2 * f(0, -hy, hz);
	main_sum[tn][17] = +2 * f(0, hy, -hz);
	main_sum[tn][18] = +2 * f(0, hy, hz);

	main_sum[tn][19] = -1 * f(-hx, -hy, -hz);
	main_sum[tn][20] = -1 * f(-hx, -hy, hz);
	main_sum[tn][21] = -1 * f(-hx, hy, -hz);
	main_sum[tn][22] = -1 * f(hx, -hy, -hz);
	main_sum[tn][23] = -1 * f(-hx, hy, hz);
	main_sum[tn][24] = -1 * f(hx, -hy, hz);
	main_sum[tn][25] = -1 * f(hx, hy, -hz);
	main_sum[tn][26] = -1 * f(hx, hy, hz);

	return sum_KahanNeumaier(main_sum[tn]) / (4 * PI * hx * hy * hz);
}

DBL3 DemagTFunc::SelfDemag(DBL3 h, bool minus)
{
	int sign = 1;
	if (minus) sign = -1;

	return DBL3(
		SelfDemag(h.x, h.y, h.z),
		SelfDemag(h.y, h.x, h.z),
		SelfDemag(h.z, h.y, h.x)) * sign;
}

//As above but also add PBC contribution to self demag
DBL3 DemagTFunc::SelfDemag_PBC(DBL3 h, DBL3 n, INT3 demag_pbc_images, int asymptotic_distance, bool minus)
{
	int sign = 1;
	if (minus) sign = -1;

	//Setup asymptotic approximation settings
	DemagAsymptoticDiag demagAsymptoticDiag_xx(h.x, h.y, h.z);
	DemagAsymptoticDiag demagAsymptoticDiag_yy(h.x, h.y, h.z);
	DemagAsymptoticDiag demagAsymptoticDiag_zz(h.x, h.y, h.z);

	DBL3 val = DBL3();

	for (int i_img = -demag_pbc_images.x; i_img < demag_pbc_images.x + 1; i_img++) {
		for (int j_img = -demag_pbc_images.y; j_img < demag_pbc_images.y + 1; j_img++) {
			for (int k_img = -demag_pbc_images.z; k_img < demag_pbc_images.z + 1; k_img++) {

				//Diagonal elements are symmetric so take modulus of the indexes without any change in sign for the tensor elements
				//There is significant redundancy remaining, could be optimized further.
				int i = mod(i_img * n.x);
				int j = mod(j_img * n.y);
				int k = mod(k_img * n.z);

				if (i == 0 && j == 0 && k == 0) {

					val += DBL3(
						SelfDemag(h.x, h.y, h.z),
						SelfDemag(h.y, h.x, h.z),
						SelfDemag(h.z, h.y, h.x)) * sign;
				}
				else {

					val += DBL3(
						demagAsymptoticDiag_xx.AsymptoticLdia(i * h.x, j * h.y, k * h.z) * sign,
						demagAsymptoticDiag_yy.AsymptoticLdia(j * h.y, i * h.x, k * h.z) * sign,
						demagAsymptoticDiag_zz.AsymptoticLdia(k * h.z, j * h.y, i * h.x) * sign);
				}
			}
		}
	}

	return val;
}