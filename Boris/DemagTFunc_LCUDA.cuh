#include "DemagTFuncCUDA.h"

#if COMPILECUDA == 1

#include "DemagTFunc_fgCUDA.cuh"

//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

//REGULAR VERSIONS FOR INTERNAL FIELD

inline __device__ double Ldia(int i, int j, int k, double hx, double hy, double hz, cuVEC<double>& f_vals)
{
	//From Newell's paper, this is what the combination of F, F1, and F2 functions reduces to in terms of the f integrals
	//Note in Newell's paper we F2(x, y, z) = f(x, y, z) - f(x, 0, z) - f(x, y, 0) - f(x, 0, 0).
	//We can drop the f(x, 0, z), f(x, y, 0), f(x, 0, 0) since they don't affect the overall sum below, as they are cancelled out

	//this function should only be called with 0 <= i < n.x, 0 <= j < n.y, 0 <= k < n.z

	//to read values from f_vals increment i, j, k by 1
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (n.x + 2), (n.y + 2), (n.z + 2)

	i++; j++; k++;

	double main_sum = +8 * f_vals[cuINT3(i, j, k)];

	main_sum += -4 * f_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * f_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * f_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

inline __device__ double Lodia(int i, int j, int k, double hx, double hy, double hz, cuVEC<double>& g_vals)
{
	//From Newell's paper, this is what the combination of G, G1, and G2 functions reduces to in terms of the g integrals
	//Note in Newell's paper we G2(x, y, z) = g(x, y, z) - g(x, y, 0).
	//We can drop the g(x, y, 0) since it doesn't affect the overall sum below, as it's cancelled out (e.g. 8-4-4=0, -4+2+2=0, 2-1-1=0, etc)

	//this function should only be called with 0 <= i < n.x, 0 <= j < n.y, 0 <= k < n.z

	//to read values from g_vals increment i, j, k by 1
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (n.x + 2), (n.y + 2), (n.z + 2)

	i++; j++; k++;

	double main_sum = +8 * g_vals[cuINT3(i, j, k)];

	main_sum += -4 * g_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * g_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

//SHIFTED VERSIONS FOR STRAY FIELD

inline __device__ double Ldia_shifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, cuVEC<double>& f_vals)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	double main_sum = +8 * f_vals[cuINT3(i, j, k)];

	main_sum += -4 * f_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * f_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * f_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

inline __device__ double Lodia_shifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, cuVEC<double>& g_vals)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	double main_sum = +8 * g_vals[cuINT3(i, j, k)];

	main_sum += -4 * g_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * g_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

//special versions for z shift only
inline __device__ double Ldia_zshifted_xx_yy(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, cuVEC<double>& f_vals)
{
	i++; j++; k += nz;

	double main_sum = +8 * f_vals[cuINT3(i, j, k)];

	main_sum += -4 * f_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * f_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * f_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

inline __device__ double Ldia_zshifted_zz(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, cuVEC<double>& f_vals)
{
	i += nx; j++; k++;

	double main_sum = +8 * f_vals[cuINT3(i, j, k)];

	main_sum += -4 * f_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * f_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * f_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * f_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * f_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * f_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

inline __device__ double Lodia_xy_zshifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, cuVEC<double>& g_vals)
{
	i++; j++; k += nz;

	double main_sum = +8 * g_vals[cuINT3(i, j, k)];

	main_sum += -4 * g_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * g_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

inline __device__ double Lodia_xz_yz_zshifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, cuVEC<double>& g_vals)
{
	i++; j += ny; k++;

	double main_sum = +8 * g_vals[cuINT3(i, j, k)];

	main_sum += -4 * g_vals[cuINT3(i + 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i - 1, j, k)];
	main_sum += -4 * g_vals[cuINT3(i, j + 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j - 1, k)];
	main_sum += -4 * g_vals[cuINT3(i, j, k + 1)];
	main_sum += -4 * g_vals[cuINT3(i, j, k - 1)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j + 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j - 1, k)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j + 1, k)];

	main_sum += +2 * g_vals[cuINT3(i - 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i - 1, j, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i + 1, j, k + 1)];

	main_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	main_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];

	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	main_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return main_sum / (4 * PI * hx * hy * hz);
}

//IRREGULAR AND SHIFTED VERSIONS FOR STRAY FIELD

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
inline __device__ double Ldia_shifted_irregular_xx_yy(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& f_vals, cuVEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	double irregular_sum = +4 * f_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * f_vals_del[cuINT3(i, j, 0)];

	irregular_sum += -2 * f_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i + 1, j, 0)];
	irregular_sum += -2 * f_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i - 1, j, 0)];
	irregular_sum += -2 * f_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i, j + 1, 0)];
	irregular_sum += -2 * f_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i, j - 1, 0)];
	irregular_sum += -4 * f_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -4 * f_vals[cuINT3(i, j, k - 1)];

	irregular_sum += +f_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i + 1, j + 1, 0)];
	irregular_sum += +f_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i + 1, j - 1, 0)];
	irregular_sum += +f_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i - 1, j + 1, 0)];
	irregular_sum += +f_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i - 1, j - 1, 0)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j - 1, k - 1)];

	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
inline __device__ double Ldia_shifted_irregular_zz(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& f_vals, cuVEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	double irregular_sum = +4 * f_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * f_vals_del[cuINT3(0, j, k)];

	irregular_sum += -4 * f_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -4 * f_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * f_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j, k + 1)];
	irregular_sum += -2 * f_vals[cuINT3(i, j, k - 1)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j, k - 1)];
	irregular_sum += -2 * f_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j + 1, k)];
	irregular_sum += -2 * f_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j - 1, k)];

	irregular_sum += +2 * f_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +f_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j + 1, k + 1)];
	irregular_sum += +f_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j + 1, k - 1)];
	irregular_sum += +f_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j - 1, k + 1)];
	irregular_sum += +f_vals[cuINT3(i, j - 1, k - 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j - 1, k - 1)];

	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
inline __device__ double Lodia_shifted_irregular_xy(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& g_vals, cuVEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	double irregular_sum = +4 * g_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * g_vals_del[cuINT3(i, j, 0)];

	irregular_sum += -2 * g_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i + 1, j, 0)];
	irregular_sum += -2 * g_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i - 1, j, 0)];
	irregular_sum += -2 * g_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, j + 1, 0)];
	irregular_sum += -2 * g_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, j - 1, 0)];
	irregular_sum += -4 * g_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -4 * g_vals[cuINT3(i, j, k - 1)];

	irregular_sum += +g_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, j + 1, 0)];
	irregular_sum += +g_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, j - 1, 0)];
	irregular_sum += +g_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, j + 1, 0)];
	irregular_sum += +g_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, j - 1, 0)];
	irregular_sum += +2 * g_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];

	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
inline __device__ double Lodia_shifted_irregular_xz_yz(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& g_vals, cuVEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j += ny; k += nz;

	double irregular_sum = +4 * g_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * g_vals_del[cuINT3(i, 0, k)];

	irregular_sum += -2 * g_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i + 1, 0, k)];
	irregular_sum += -2 * g_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i - 1, 0, k)];
	irregular_sum += -4 * g_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -4 * g_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * g_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, 0, k + 1)];
	irregular_sum += -2 * g_vals[cuINT3(i, j, k - 1)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, 0, k - 1)];

	irregular_sum += +2 * g_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +2 * g_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +g_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, 0, k + 1)];
	irregular_sum += +g_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, 0, k - 1)];
	irregular_sum += +g_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, 0, k + 1)];
	irregular_sum += +g_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, 0, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];

	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
inline __device__ double Ldia_zshifted_irregular_xx_yy(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& f_vals, cuVEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i++; j++; k += nz;

	double irregular_sum = +4 * f_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * f_vals_del[cuINT3(i, j, 0)];

	irregular_sum += -2 * f_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i + 1, j, 0)];
	irregular_sum += -2 * f_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i - 1, j, 0)];
	irregular_sum += -2 * f_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i, j + 1, 0)];
	irregular_sum += -2 * f_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(i, j - 1, 0)];
	irregular_sum += -4 * f_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -4 * f_vals[cuINT3(i, j, k - 1)];

	irregular_sum += +f_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i + 1, j + 1, 0)];
	irregular_sum += +f_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i + 1, j - 1, 0)];
	irregular_sum += +f_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i - 1, j + 1, 0)];
	irregular_sum += +f_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +f_vals_del[cuINT3(i - 1, j - 1, 0)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i, j - 1, k - 1)];

	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
inline __device__ double Ldia_zshifted_irregular_zz(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& f_vals, cuVEC<double>& f_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from f_vals increment i, j, k by nx, ny, nz resp.
	//f_vals stores f values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i += nx; j++; k++;

	double irregular_sum = +4 * f_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * f_vals_del[cuINT3(0, j, k)];

	irregular_sum += -4 * f_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -4 * f_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * f_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j, k + 1)];
	irregular_sum += -2 * f_vals[cuINT3(i, j, k - 1)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j, k - 1)];
	irregular_sum += -2 * f_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j + 1, k)];
	irregular_sum += -2 * f_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * f_vals_del[cuINT3(0, j - 1, k)];

	irregular_sum += +2 * f_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +2 * f_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +f_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j + 1, k + 1)];
	irregular_sum += +f_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j + 1, k - 1)];
	irregular_sum += +f_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j - 1, k + 1)];
	irregular_sum += +f_vals[cuINT3(i, j - 1, k - 1)];
	irregular_sum += +f_vals_del[cuINT3(0, j - 1, k - 1)];

	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * f_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
inline __device__ double Lodia_zshifted_irregular_xy(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& g_vals, cuVEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i++; j++; k += nz;

	double irregular_sum = +4 * g_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * g_vals_del[cuINT3(i, j, 0)];

	irregular_sum += -2 * g_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i + 1, j, 0)];
	irregular_sum += -2 * g_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i - 1, j, 0)];
	irregular_sum += -2 * g_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, j + 1, 0)];
	irregular_sum += -2 * g_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, j - 1, 0)];
	irregular_sum += -4 * g_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -4 * g_vals[cuINT3(i, j, k - 1)];

	irregular_sum += +g_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, j + 1, 0)];
	irregular_sum += +g_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, j - 1, 0)];
	irregular_sum += +g_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, j + 1, 0)];
	irregular_sum += +g_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, j - 1, 0)];
	irregular_sum += +2 * g_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];

	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
inline __device__ double Lodia_zshifted_irregular_xz_yz(int i, int j, int k, int nx, int ny, int nz, double tau, cuVEC<double>& g_vals, cuVEC<double>& g_vals_del)
{
	//this function should only be called with -n.x + 1 <= i < n.x, -n.y + 1 <= j < n.y, -n.z + 1 <= k < n.z

	//to read values from g_vals increment i, j, k by nx, ny, nz resp.
	//g_vals stores g values on the mesh, including boundary points, thus it has dimensions (2*n.x + 1), (2*n.y + 1), (2*n.z + 1)

	i++; j += ny; k++;

	double irregular_sum = +4 * g_vals[cuINT3(i, j, k)];
	irregular_sum += +4 * g_vals_del[cuINT3(i, 0, k)];

	irregular_sum += -2 * g_vals[cuINT3(i + 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i + 1, 0, k)];
	irregular_sum += -2 * g_vals[cuINT3(i - 1, j, k)];
	irregular_sum += -2 * g_vals_del[cuINT3(i - 1, 0, k)];
	irregular_sum += -4 * g_vals[cuINT3(i, j + 1, k)];
	irregular_sum += -4 * g_vals[cuINT3(i, j - 1, k)];
	irregular_sum += -2 * g_vals[cuINT3(i, j, k + 1)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, 0, k + 1)];
	irregular_sum += -2 * g_vals[cuINT3(i, j, k - 1)];
	irregular_sum += -2 * g_vals_del[cuINT3(i, 0, k - 1)];

	irregular_sum += +2 * g_vals[cuINT3(i + 1, j + 1, k)];
	irregular_sum += +2 * g_vals[cuINT3(i + 1, j - 1, k)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j + 1, k)];
	irregular_sum += +2 * g_vals[cuINT3(i - 1, j - 1, k)];
	irregular_sum += +g_vals[cuINT3(i + 1, j, k + 1)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, 0, k + 1)];
	irregular_sum += +g_vals[cuINT3(i + 1, j, k - 1)];
	irregular_sum += +g_vals_del[cuINT3(i + 1, 0, k - 1)];
	irregular_sum += +g_vals[cuINT3(i - 1, j, k + 1)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, 0, k + 1)];
	irregular_sum += +g_vals[cuINT3(i - 1, j, k - 1)];
	irregular_sum += +g_vals_del[cuINT3(i - 1, 0, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j + 1, k - 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k + 1)];
	irregular_sum += +2 * g_vals[cuINT3(i, j - 1, k - 1)];

	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i - 1, j + 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j - 1, k + 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k - 1)];
	irregular_sum += -1 * g_vals[cuINT3(i + 1, j + 1, k + 1)];

	return irregular_sum / (4 * PI * tau);
}

// SINGLE VALUE COMPUTE (without pre-computed f, g values)

//get diagonal component at given position for given cellsizes
inline __device__ double Ldia_single(double x, double y, double z, double hx, double hy, double hz)
{
	double main_sum = +8 * f(x, y, z);

	main_sum += -4 * f(x + hx, y, z);
	main_sum += -4 * f(x - hx, y, z);
	main_sum += -4 * f(x, y + hy, z);
	main_sum += -4 * f(x, y - hy, z);
	main_sum += -4 * f(x, y, z + hz);
	main_sum += -4 * f(x, y, z - hz);

	main_sum += +2 * f(x - hx, y - hy, z);
	main_sum += +2 * f(x - hx, y + hy, z);
	main_sum += +2 * f(x + hx, y - hy, z);
	main_sum += +2 * f(x + hx, y + hy, z);

	main_sum += +2 * f(x - hx, y, z - hz);
	main_sum += +2 * f(x - hx, y, z + hz);
	main_sum += +2 * f(x + hx, y, z - hz);
	main_sum += +2 * f(x + hx, y, z + hz);

	main_sum += +2 * f(x, y - hy, z - hz);
	main_sum += +2 * f(x, y - hy, z + hz);
	main_sum += +2 * f(x, y + hy, z - hz);
	main_sum += +2 * f(x, y + hy, z + hz);

	main_sum += -1 * f(x - hx, y - hy, z - hz);
	main_sum += -1 * f(x - hx, y - hy, z + hz);
	main_sum += -1 * f(x - hx, y + hy, z - hz);
	main_sum += -1 * f(x + hx, y - hy, z - hz);
	main_sum += -1 * f(x - hx, y + hy, z + hz);
	main_sum += -1 * f(x + hx, y - hy, z + hz);
	main_sum += -1 * f(x + hx, y + hy, z - hz);
	main_sum += -1 * f(x + hx, y + hy, z + hz);

	return main_sum / (4 * PI * hx * hy * hz);
}

//get off-diagonal component at given position for given cellsizes
inline __device__ double Lodia_single(double x, double y, double z, double hx, double hy, double hz)
{
	double main_sum = +8 * g(x, y, z);

	main_sum += -4 * g(x + hx, y, z);
	main_sum += -4 * g(x - hx, y, z);
	main_sum += -4 * g(x, y + hy, z);
	main_sum += -4 * g(x, y - hy, z);
	main_sum += -4 * g(x, y, z + hz);
	main_sum += -4 * g(x, y, z - hz);

	main_sum += +2 * g(x - hx, y - hy, z);
	main_sum += +2 * g(x - hx, y + hy, z);
	main_sum += +2 * g(x + hx, y - hy, z);
	main_sum += +2 * g(x + hx, y + hy, z);

	main_sum += +2 * g(x - hx, y, z - hz);
	main_sum += +2 * g(x - hx, y, z + hz);
	main_sum += +2 * g(x + hx, y, z - hz);
	main_sum += +2 * g(x + hx, y, z + hz);

	main_sum += +2 * g(x, y - hy, z - hz);
	main_sum += +2 * g(x, y - hy, z + hz);
	main_sum += +2 * g(x, y + hy, z - hz);
	main_sum += +2 * g(x, y + hy, z + hz);

	main_sum += -1 * g(x - hx, y - hy, z - hz);
	main_sum += -1 * g(x - hx, y - hy, z + hz);
	main_sum += -1 * g(x - hx, y + hy, z - hz);
	main_sum += -1 * g(x + hx, y - hy, z - hz);
	main_sum += -1 * g(x - hx, y + hy, z + hz);
	main_sum += -1 * g(x + hx, y - hy, z + hz);
	main_sum += -1 * g(x + hx, y + hy, z - hz);
	main_sum += -1 * g(x + hx, y + hy, z + hz);

	return main_sum / (4 * PI * hx * hy * hz);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
inline __device__ double Ldia_shifted_irregular_xx_yy_single(double x, double y, double z, double hx, double hy, double sz, double dz)
{
	double del = sz - dz;

	double irregular_sum = +4 * f(x, y, z);
	irregular_sum += +4 * f(x, y, z - del);

	irregular_sum += -2 * f(x + hx, y, z);
	irregular_sum += -2 * f(x + hx, y, z - del);
	irregular_sum += -2 * f(x - hx, y, z);
	irregular_sum += -2 * f(x - hx, y, z - del);
	irregular_sum += -2 * f(x, y + hy, z);
	irregular_sum += -2 * f(x, y + hy, z - del);
	irregular_sum += -2 * f(x, y - hy, z);
	irregular_sum += -2 * f(x, y - hy, z - del);
	irregular_sum += -4 * f(x, y, z + dz);
	irregular_sum += -4 * f(x, y, z - sz);

	irregular_sum += +f(x + hx, y + hy, z);
	irregular_sum += +f(x + hx, y + hy, z - del);
	irregular_sum += +f(x + hx, y - hy, z);
	irregular_sum += +f(x + hx, y - hy, z - del);
	irregular_sum += +f(x - hx, y + hy, z);
	irregular_sum += +f(x - hx, y + hy, z - del);
	irregular_sum += +f(x - hx, y - hy, z);
	irregular_sum += +f(x - hx, y - hy, z - del);
	irregular_sum += +2 * f(x + hx, y, z + dz);
	irregular_sum += +2 * f(x + hx, y, z - sz);
	irregular_sum += +2 * f(x - hx, y, z + dz);
	irregular_sum += +2 * f(x - hx, y, z - sz);
	irregular_sum += +2 * f(x, y + hy, z + dz);
	irregular_sum += +2 * f(x, y + hy, z - sz);
	irregular_sum += +2 * f(x, y - hy, z + dz);
	irregular_sum += +2 * f(x, y - hy, z - sz);

	irregular_sum += -1 * f(x - hx, y - hy, z - sz);
	irregular_sum += -1 * f(x - hx, y - hy, z + dz);
	irregular_sum += -1 * f(x - hx, y + hy, z - sz);
	irregular_sum += -1 * f(x + hx, y - hy, z - sz);
	irregular_sum += -1 * f(x - hx, y + hy, z + dz);
	irregular_sum += -1 * f(x + hx, y - hy, z + dz);
	irregular_sum += -1 * f(x + hx, y + hy, z - sz);
	irregular_sum += -1 * f(x + hx, y + hy, z + dz);

	return irregular_sum / (4 * PI * hx * hy * dz);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
inline __device__ double Ldia_shifted_irregular_zz_single(double x, double y, double z, double sx, double dx, double hy, double hz)
{
	double del = sx - dx;

	double irregular_sum = +4 * f(x, y, z);
	irregular_sum += +4 * f(x + del, y, z);

	irregular_sum += -4 * f(x + sx, y, z);
	irregular_sum += -4 * f(x - dx, y, z);
	irregular_sum += -2 * f(x, y, z + hz);
	irregular_sum += -2 * f(x + del, y, z + hz);
	irregular_sum += -2 * f(x, y, z - hz);
	irregular_sum += -2 * f(x + del, y, z - hz);
	irregular_sum += -2 * f(x, y + hy, z);
	irregular_sum += -2 * f(x + del, y + hy, z);
	irregular_sum += -2 * f(x, y - hy, z);
	irregular_sum += -2 * f(x + del, y - hy, z);

	irregular_sum += +2 * f(x + sx, y + hy, z);
	irregular_sum += +2 * f(x + sx, y - hy, z);
	irregular_sum += +2 * f(x - dx, y + hy, z);
	irregular_sum += +2 * f(x - dx, y - hy, z);
	irregular_sum += +2 * f(x + sx, y, z + hz);
	irregular_sum += +2 * f(x + sx, y, z - hz);
	irregular_sum += +2 * f(x - dx, y, z + hz);
	irregular_sum += +2 * f(x - dx, y, z - hz);
	irregular_sum += +f(x, y + hy, z + hz);
	irregular_sum += +f(x + del, y + hy, z + hz);
	irregular_sum += +f(x, y + hy, z - hz);
	irregular_sum += +f(x + del, y + hy, z - hz);
	irregular_sum += +f(x, y - hy, z + hz);
	irregular_sum += +f(x + del, y - hy, z + hz);
	irregular_sum += +f(x, y - hy, z - hz);
	irregular_sum += +f(x + del, y - hy, z - hz);

	irregular_sum += -1 * f(x - dx, y - hy, z - hz);
	irregular_sum += -1 * f(x - dx, y - hy, z + hz);
	irregular_sum += -1 * f(x - dx, y + hy, z - hz);
	irregular_sum += -1 * f(x + sx, y - hy, z - hz);
	irregular_sum += -1 * f(x - dx, y + hy, z + hz);
	irregular_sum += -1 * f(x + sx, y - hy, z + hz);
	irregular_sum += -1 * f(x + sx, y + hy, z - hz);
	irregular_sum += -1 * f(x + sx, y + hy, z + hz);

	return irregular_sum / (4 * PI * dx * hy * hz);
}

//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
inline __device__ double Lodia_shifted_irregular_xy_single(double x, double y, double z, double hx, double hy, double sz, double dz)
{
	double del = sz - dz;

	double irregular_sum = +4 * g(x, y, z);
	irregular_sum += +4 * g(x, y, z + del);

	irregular_sum += -2 * g(x + hx, y, z);
	irregular_sum += -2 * g(x + hx, y, z + del);
	irregular_sum += -2 * g(x - hx, y, z);
	irregular_sum += -2 * g(x - hx, y, z + del);
	irregular_sum += -2 * g(x, y + hy, z);
	irregular_sum += -2 * g(x, y + hy, z + del);
	irregular_sum += -2 * g(x, y - hy, z);
	irregular_sum += -2 * g(x, y - hy, z + del);
	irregular_sum += -4 * g(x, y, z + sz);
	irregular_sum += -4 * g(x, y, z - dz);

	irregular_sum += +g(x + hx, y + hy, z);
	irregular_sum += +g(x + hx, y + hy, z + del);
	irregular_sum += +g(x + hx, y - hy, z);
	irregular_sum += +g(x + hx, y - hy, z + del);
	irregular_sum += +g(x - hx, y + hy, z);
	irregular_sum += +g(x - hx, y + hy, z + del);
	irregular_sum += +g(x - hx, y - hy, z);
	irregular_sum += +g(x - hx, y - hy, z + del);
	irregular_sum += +2 * g(x + hx, y, z + sz);
	irregular_sum += +2 * g(x + hx, y, z - dz);
	irregular_sum += +2 * g(x - hx, y, z + sz);
	irregular_sum += +2 * g(x - hx, y, z - dz);
	irregular_sum += +2 * g(x, y + hy, z + sz);
	irregular_sum += +2 * g(x, y + hy, z - dz);
	irregular_sum += +2 * g(x, y - hy, z + sz);
	irregular_sum += +2 * g(x, y - hy, z - dz);

	irregular_sum += -1 * g(x - hx, y - hy, z - dz);
	irregular_sum += -1 * g(x - hx, y - hy, z + sz);
	irregular_sum += -1 * g(x - hx, y + hy, z - dz);
	irregular_sum += -1 * g(x + hx, y - hy, z - dz);
	irregular_sum += -1 * g(x - hx, y + hy, z + sz);
	irregular_sum += -1 * g(x + hx, y - hy, z + sz);
	irregular_sum += -1 * g(x + hx, y + hy, z - dz);
	irregular_sum += -1 * g(x + hx, y + hy, z + sz);

	return irregular_sum / (4 * PI * hx * hy * dz);
}

//off-diagonal components for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
inline __device__ double Lodia_shifted_irregular_xz_yz_single(double x, double y, double z, double hx, double sy, double dy, double hz)
{
	double del = sy - dy;

	double irregular_sum = +4 * g(x, y, z);
	irregular_sum += +4 * g(x, y + del, z);

	irregular_sum += -2 * g(x + hx, y, z);
	irregular_sum += -2 * g(x + hx, y + del, z);
	irregular_sum += -2 * g(x - hx, y, z);
	irregular_sum += -2 * g(x - hx, y + del, z);
	irregular_sum += -4 * g(x, y + sy, z);
	irregular_sum += -4 * g(x, y - dy, z);
	irregular_sum += -2 * g(x, y, z + hz);
	irregular_sum += -2 * g(x, y + del, z + hz);
	irregular_sum += -2 * g(x, y, z - hz);
	irregular_sum += -2 * g(x, y + del, z - hz);

	irregular_sum += +2 * g(x + hx, y + sy, z);
	irregular_sum += +2 * g(x + hx, y - dy, z);
	irregular_sum += +2 * g(x - hx, y + sy, z);
	irregular_sum += +2 * g(x - hx, y - dy, z);
	irregular_sum += +g(x + hx, y, z + hz);
	irregular_sum += +g(x + hx, y + del, z + hz);
	irregular_sum += +g(x + hx, y, z - hz);
	irregular_sum += +g(x + hx, y + del, z - hz);
	irregular_sum += +g(x - hx, y, z + hz);
	irregular_sum += +g(x - hx, y + del, z + hz);
	irregular_sum += +g(x - hx, y, z - hz);
	irregular_sum += +g(x - hx, y + del, z - hz);
	irregular_sum += +2 * g(x, y + sy, z + hz);
	irregular_sum += +2 * g(x, y + sy, z - hz);
	irregular_sum += +2 * g(x, y - dy, z + hz);
	irregular_sum += +2 * g(x, y - dy, z - hz);

	irregular_sum += -1 * g(x - hx, y - dy, z - hz);
	irregular_sum += -1 * g(x - hx, y - dy, z + hz);
	irregular_sum += -1 * g(x - hx, y + sy, z - hz);
	irregular_sum += -1 * g(x + hx, y - dy, z - hz);
	irregular_sum += -1 * g(x - hx, y + sy, z + hz);
	irregular_sum += -1 * g(x + hx, y - dy, z + hz);
	irregular_sum += -1 * g(x + hx, y + sy, z - hz);
	irregular_sum += -1 * g(x + hx, y + sy, z + hz);

	return irregular_sum / (4 * PI * hx * dy * hz);
}

#endif