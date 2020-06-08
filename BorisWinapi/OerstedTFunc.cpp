#include "stdafx.h"
#include "OerstedTFunc.h"

double OeFunc::G1(double x, double y, double z) 
{
	std::vector<double> sum_this(7, 0);

	double xsq = x * x;
	double ysq = y * y;
	double zsq = z * z;

	double R = sqrt(xsq + ysq + zsq);

	sum_this[0] = (3 * xsq + 3 * ysq - 2 * zsq) * z * R / 24;

	sum_this[1] = PI * z * fabs(x*y*z) / 4;

	if (z + R) sum_this[2] = (xsq*xsq - 6 * xsq*ysq + ysq * ysq) * log(z + R) / 24;
	else sum_this[2] = 0.0;

	if (y) sum_this[3] = x * y * (ysq - 3 * zsq) * atan(x*z / (y*R)) / 6;
	else sum_this[3] = 0.0;

	if (x) sum_this[4] = x * y * (xsq - 3 * zsq) * atan(y*z / (x*R)) / 6;
	else sum_this[4] = 0.0;

	if (x + R) sum_this[5] = z * x * (zsq - 3 * ysq) * log(x + R) / 6;
	else sum_this[5] = 0.0;

	if (y + R) sum_this[6] = z * y * (zsq - 3 * xsq) * log(y + R) / 6;
	else sum_this[6] = 0.0;

	return sum_KahanNeumaier(sum_this);
}

double OeFunc::I2(double x, double y, double z, double lx, double ly, double lz) 
{
	//x, y, z should be multiples of lx, ly, lz respectively. e.g. set lx = ly = lz = 1 for cubic cells then x, y, z are just the mesh cell numbers along x, y and z respectively.

	std::vector<double> sum_this(27, 0);

	int idx = 0;

	for (int n = -1; n <= 1; n++) {
		for (int m = -1; m <= 1; m++) {
			for (int k = -1; k <= 1; k++) {

				double c = double((2 - 3 * abs(n)) * (2 - 3 * abs(m)) * (2 - 3 * abs(k))) / (lx * ly * lz);

				sum_this[idx++] = c * G1(x + n * lx, y + m * ly, z + k * lz);
			}
		}
	}

	return sum_KahanNeumaier(sum_this);
}

//Kxy = -Kyx
double OeFunc::CalcKxy(DBL3 pos, DBL3 h) 
{
	return -I2(pos.x, pos.y, pos.z, h.x, h.y, h.z) / (4 * PI);
}

//Kyz = -Kzy
double OeFunc::CalcKyz(DBL3 pos, DBL3 h) 
{
	return -I2(pos.y, pos.z, pos.x, h.y, h.z, h.x) / (4 * PI);
}

//Kzx = -Kxz
double OeFunc::CalcKxz(DBL3 pos, DBL3 h) 
{
	return I2(pos.z, pos.x, pos.y, h.z, h.x, h.y) / (4 * PI);
}

void OeFunc::CalcOerstedTensors(VEC<DBL3>& DOe, INT3 n, INT3 N, DBL3 h)
{
	/*
	//Loop without use of symmetries to reduce computation time:

	#pragma omp parallel for
	for (int j = -(n.y - 1); j <= n.y - 1; j++) {
		for (int k = -(n.z - 1); k <= n.z - 1; k++) {
			for (int i = -(n.x - 1); i <= n.x - 1; i++) {

				DOe[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z) * N.x*N.y] = DBL3(CalcKxy(h & INT3(i, j, k), h), 
																									CalcKxz(h & INT3(i, j, k), h), 
																									CalcKyz(h & INT3(i, j, k), h));
			}
		}
	}
	*/

	//zero the tensor first
	DOe.set(DBL3());

	//use symmetries to reduce computation time:

	//Kxy : even in x, y, odd in z
	//Kxz : even in x, z, odd in y
	//Kyz : odd in x, even in y, z

	//NOTE : computation can be made more efficient by calculating the G1 values once then re-using them in the I2 function.

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				DBL3 val = DBL3(
					CalcKxy(h & INT3(i, j, k), h),
					CalcKxz(h & INT3(i, j, k), h),
					CalcKyz(h & INT3(i, j, k), h));

				DOe[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z) * N.x*N.y] = val & DBL3(1, 1, 1);
				DOe[(-i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((k + N.z) % N.z) * N.x*N.y] = val & DBL3(1, 1, -1);
				DOe[(i + N.x) % N.x + ((-j + N.y) % N.y) * N.x + ((k + N.z) % N.z) * N.x*N.y] = val & DBL3(1, -1, 1);
				DOe[(i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((-k + N.z) % N.z) * N.x*N.y] = val & DBL3(-1, 1, 1);

				DOe[(-i + N.x) % N.x + ((-j + N.y) % N.y) * N.x + ((k + N.z) % N.z) * N.x*N.y] = val & DBL3(1, -1, -1);
				DOe[(-i + N.x) % N.x + ((j + N.y) % N.y) * N.x + ((-k + N.z) % N.z) * N.x*N.y] = val & DBL3(-1, 1, -1);
				DOe[(i + N.x) % N.x + ((-j + N.y) % N.y) * N.x + ((-k + N.z) % N.z) * N.x*N.y] = val & DBL3(-1, -1, 1);
				DOe[(-i + N.x) % N.x + ((-j + N.y) % N.y) * N.x + ((-k + N.z) % N.z) * N.x*N.y] = val & DBL3(-1, -1, -1);
			}
		}
	}
}