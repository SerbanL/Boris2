#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include <cuda_runtime.h>

//---------------------f and g VECTORS COMPUTATION

//Basis f and g functions - these are correct
//The checks in them are absolutely essential, otherwise you'll run into NaNs sooner or later, especially in multi-layered convolution.
//Moreover don't check a double for zero, use the IsZ or IsNZ functions. If you do, guaranteed eventually the check will fail and result in NaNs.

inline __device__ double f(double x, double y, double z)
{
	double xsq = x * x;
	double ysq = y * y;
	double zsq = z * z;

	double Rxy_sq = xsq + ysq;
	double Rxz_sq = xsq + zsq;

	double R = sqrt(xsq + ysq + zsq);

	double f_sum = (2.0 * xsq - ysq - zsq) * R / 6.0;

	if (cuIsNZ(Rxz_sq)) {

		double arg = 2.0 * y * (y + R) / Rxz_sq;
		if (arg > -1.0) f_sum += y * (zsq - xsq) * log1p(arg) / 4.0;
	}

	if (cuIsNZ(Rxy_sq)) {

		double arg = 2.0 * z * (z + R) / Rxy_sq;
		if (arg > -1.0) f_sum += z * (ysq - xsq) * log1p(arg) / 4.0;
	}

	if (cuIsNZ(x)) f_sum += -x * y*z*atan(y*z / (x*R));

	return f_sum;
}

inline __device__ double g(double x, double y, double z)
{
	double xsq = x * x;
	double ysq = y * y;
	double zsq = z * z;

	double Rxy_sq = xsq + ysq;
	double Rxz_sq = xsq + zsq;
	double Ryz_sq = ysq + zsq;

	double R = sqrt(xsq + ysq + zsq);

	double g_sum = -x * y * R / 3.0;

	if (cuIsNZ(Rxy_sq)) {

		double arg = 2.0 * z * (z + R) / Rxy_sq;
		if (arg > -1.0) g_sum += x * y*z*log1p(arg) / 2.0;
	}

	if (cuIsNZ(Ryz_sq)) {

		double arg = 2.0 * x * (x + R) / Ryz_sq;
		if (arg > -1.0) g_sum += y * (3.0 * zsq - ysq) * log1p(arg) / 12.0;
	}

	if (cuIsNZ(Rxz_sq)) {

		double arg = 2.0 * y * (y + R) / Rxz_sq;
		if (arg > -1.0) g_sum += x * (3.0 * zsq - xsq) * log1p(arg) / 12.0;
	}

	if (cuIsNZ(z)) g_sum += -zsq * z * atan(x*y / (z*R)) / 6.0;

	if (cuIsNZ(y)) g_sum += -ysq * z * atan(x*z / (y*R)) / 2.0;

	if (cuIsNZ(x)) g_sum += -xsq * z * atan(y*z / (x*R)) / 2.0;

	return g_sum;
}

#endif
