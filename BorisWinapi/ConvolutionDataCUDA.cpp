#include "stdafx.h"
#include "ConvolutionDataCUDA.h"

#if COMPILECUDA == 1

ConvolutionDataCUDA::~ConvolutionDataCUDA()
{
	//if n.z == 0 then no plan is currently allocated
	if (n.z == 1) {

		cufftDestroy(plan2D_fwd_x);
		cufftDestroy(plan2D_y);
		cufftDestroy(plan2D_inv_x);
	}
	else if (n.z > 1) {

		cufftDestroy(plan3D_fwd_x);
		cufftDestroy(plan3D_y);
		cufftDestroy(plan3D_z);
		cufftDestroy(plan3D_inv_x);
	}
}

BError ConvolutionDataCUDA::SetConvolutionDimensions(cuSZ3 n_, cuReal3 h_, bool embed_multiplication_)
{
	BError error(__FUNCTION__);

	//Turn off multiplication embedding by setting this to false - this will allocate full space memory for F and F2 scratch spaces. embed_multiplication = true by default.
	embed_multiplication = embed_multiplication_;

	//if n.z == 0 then no plan is currently allocated (since n has not been set yet)
	if (n.z == 1) {

		//currently using 2D plan : free resources as new ones will be remade
		cufftDestroy(plan2D_fwd_x);
		cufftDestroy(plan2D_y);
		cufftDestroy(plan2D_inv_x);
	}
	else if (n.z > 1) {

		//currently using 3D plan : free resources as new ones will be remade
		cufftDestroy(plan3D_fwd_x);
		cufftDestroy(plan3D_y);
		cufftDestroy(plan3D_z);
		cufftDestroy(plan3D_inv_x);
	}

	n = n_;
	h = h_;

	N = cuSZ3(1, 1, 1);

	//set N, M, K as smallest powers of 2 which will hold the demag kernel
	while (N.x < 2 * n.x - 1) { N.x = N.x * 2; };
	while (N.y < 2 * n.y - 1) { N.y = N.y * 2; };
	if (n.z > 1) while (N.z < 2 * n.z - 1) { N.z = N.z * 2; };

	//N in gpu memory
	cuN.from_cpu(N);

	//set size of complex arrays in gpu memory
	cuNc_xy.from_cpu(cuSZ3(N.x / 2 + 1, N.y, N.z));
	cuNc_yx.from_cpu(cuSZ3(N.y, N.x / 2 + 1, N.z));

	cuNc_xy_q2d.from_cpu(cuSZ3(N.x / 2 + 1, N.y, N.z / 2));
	cuNc_yx_q2d.from_cpu(cuSZ3(N.y, N.x / 2 + 1, N.z / 2));

	cuNcquart_xy.from_cpu(cuSZ3(N.x / 2 + 1, n.y, n.z));
	cuNcquart_yx.from_cpu(cuSZ3(n.y, N.x / 2 + 1, n.z));

	cuUpper_y_region.from_cpu(cuRect(cuReal3(0, n.y, 0), cuReal3(N.x / 2 + 1, N.y, n.z)));
	cuUpper_y_transposed_region.from_cpu(cuRect(cuReal3(n.y, 0, 0), cuReal3(N.y, N.x / 2 + 1, n.z)));
	cuUpper_z_region.from_cpu(cuRect(cuReal3(0, 0, n.z), cuReal3(N.x / 2 + 1, N.y, N.z)));

	cuUpper_z_region_q2d.from_cpu(cuRect(cuReal3(0, 0, n.z), cuReal3(N.x / 2 + 1, N.y, N.z / 2)));

	//Make cuda FFT plans
	cufftResult cuffterr = CUFFT_SUCCESS;

	if (n.z == 1) {

		//q2D mode not applicable for true 2D problems
		q2D_level = 0;
		
		//transpose xy before y fft / ifft?
		if (N.y >= TRANSPOSE_XY_YTHRESHOLD) transpose_xy = true;
		else transpose_xy = false;

#if SINGLEPRECISION == 1

		int embed[1] = { 0 };
		int ndims_x[1] = { (int)N.x };
		int ndims_y[1] = { (int)N.y };

		//Forward fft along x direction (out-of-place):
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_fwd_x, 1, ndims_x,
			embed, 1, N.x,
			embed, 1, (N.x / 2 + 1),
			CUFFT_R2C, (int)n.y);
		
		if (!transpose_xy) {

			//Forward fft along y direction
			if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_y, 1, ndims_y,
				embed, (N.x / 2 + 1), 1,
				embed, (N.x / 2 + 1), 1,
				CUFFT_C2C, (int)(N.x / 2 + 1));
		}
		else {

			if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_y, 1, ndims_y,
				embed, 1, N.y,
				embed, 1, N.y,
				CUFFT_C2C, (int)(N.x / 2 + 1));
		}

		//Inverse fft along x direction
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_inv_x, 1, ndims_x,
			embed, 1, (N.x / 2 + 1),
			embed, 1, N.x,
			CUFFT_C2R, (int)n.y);

#else

		int embed[1] = { 0 };
		int ndims_x[1] = { (int)N.x };
		int ndims_y[1] = { (int)N.y };

		//Forward fft along x direction (out-of-place):
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_fwd_x, 1, ndims_x,
			embed, 1, N.x,
			embed, 1, (N.x / 2 + 1),
			CUFFT_D2Z, (int)n.y);

		if (!transpose_xy) {

			//Forward fft along y direction
			if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_y, 1, ndims_y,
				embed, (N.x / 2 + 1), 1,
				embed, (N.x / 2 + 1), 1,
				CUFFT_Z2Z, (int)(N.x / 2 + 1));
		}
		else {

			if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_y, 1, ndims_y,
				embed, 1, N.y,
				embed, 1, N.y,
				CUFFT_Z2Z, (int)(N.x / 2 + 1));
		}

		//Inverse fft along x direction
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan2D_inv_x, 1, ndims_x,
			embed, 1, (N.x / 2 + 1),
			embed, 1, N.x,
			CUFFT_Z2D, (int)n.y);

#endif
	}
	else {

		//3D problem
		
		//quasi 2D mode? (not applicable if convolution not embedded)
		if (embed_multiplication && n.z == 2) {

			//N.z = 4
			q2D_level = 4;
		}
		else if (embed_multiplication && n.z <= 4) {

			//N.z = 8
			q2D_level = 8;
		}
		
		else if (embed_multiplication && n.z <= 8) {

			//N.z = 16
			q2D_level = 16;
		}
		
		else if (embed_multiplication && n.z <= 16) {

			//N.z = 32
			q2D_level = 32;
		}
		//above this level q2D mode is slower than full 3D mode due to inefficient use of gpu bandwidth
		else {

			//disable q2D mode
			q2D_level = 0;
		}

		//always use transpose_xy in 3D
		transpose_xy = true;

#if SINGLEPRECISION == 1

		int embed[1] = { 0 };		
		int ndims_x[1] = { (int)N.x };
		int ndims_y[1] = { (int)N.y };
		int ndims_z[1] = { (int)N.z };

		//batched x fft from cuIn to cuSquart
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_fwd_x, 1, ndims_x,
			embed, 1, N.x,
			embed, 1, (N.x / 2 + 1),
			CUFFT_R2C, (int)n.y * (int)n.z);

		//batched y fft from cuS (transposed) to cuS
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_y, 1, ndims_y,
			embed, 1, N.y,
			embed, 1, N.y,
			CUFFT_C2C, (int)(N.x / 2 + 1) * n.z);

		//batched inverse x fft from cuSquart to cuOut
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_inv_x, 1, ndims_x,
			embed, 1, (N.x / 2 + 1),
			embed, 1, N.x,
			CUFFT_C2R, (int)n.y * (int)n.z);
		
		if (!q2D_level) {

			//Forward and inverse fft along z direction all batched - not applicable for q2D mode

			if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_z, 1, ndims_z,
				embed, (N.x / 2 + 1)*N.y, 1,
				embed, (N.x / 2 + 1)*N.y, 1,
				CUFFT_C2C, (int)(N.x / 2 + 1) * N.y);
		}

#else
		
		int embed[1] = { 0 };
		int ndims_x[1] = { (int)N.x };
		int ndims_y[1] = { (int)N.y };
		int ndims_z[1] = { (int)N.z };

		//batched x fft from cuIn to cuSquart
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_fwd_x, 1, ndims_x,
			embed, 1, N.x,
			embed, 1, (N.x / 2 + 1),
			CUFFT_D2Z, (int)n.y * (int)n.z);

		//batched y fft from cuS (transposed) to cuS
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_y, 1, ndims_y,
			embed, 1, N.y,
			embed, 1, N.y,
			CUFFT_Z2Z, (int)(N.x / 2 + 1) * n.z);

		//batched inverse x fft from cuSquart to cuOut
		if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_inv_x, 1, ndims_x,
			embed, 1, (N.x / 2 + 1),
			embed, 1, N.x,
			CUFFT_Z2D, (int)n.y * (int)n.z);

		if (!q2D_level) {

			//Forward and inverse fft along z direction all batched - not applicable for q2D mode

			if (cuffterr == CUFFT_SUCCESS) cuffterr = cufftPlanMany(&plan3D_z, 1, ndims_z,
				embed, (N.x / 2 + 1)*N.y, 1,
				embed, (N.x / 2 + 1)*N.y, 1,
				CUFFT_Z2Z, (int)(N.x / 2 + 1) * N.y);
		}

#endif
	}

	if (cuffterr != CUFFT_SUCCESS) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//FFT spaces memory allocation

	//Input/Output space

	if (!cuIn_x.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIn_y.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIn_z.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	cuIn_x.set(cuReal(0.0));
	cuIn_y.set(cuReal(0.0));
	cuIn_z.set(cuReal(0.0));

	if (!cuOut_x.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuOut_y.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuOut_z.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	cuOut_x.set(cuReal(0.0));
	cuOut_y.set(cuReal(0.0));
	cuOut_z.set(cuReal(0.0));
	
	//full-size scratch space
	if (!q2D_level) {

		//applies for 3D problems as well as 2D problems (n.z = 1 and thus N.z = 1)

		if (!cuS_x.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_y.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_z.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		cuS_x.set(cuComplex());
		cuS_y.set(cuComplex());
		cuS_z.set(cuComplex());
	}
	else {

		//in q2D mode we can halve the space required for cuS scratch spaces.

		if (!cuS_x.resize((N.x / 2 + 1) * N.y * N.z / 2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_y.resize((N.x / 2 + 1) * N.y * N.z / 2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_z.resize((N.x / 2 + 1) * N.y * N.z / 2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		cuS_x.set(cuComplex());
		cuS_y.set(cuComplex());
		cuS_z.set(cuComplex());
	}

	if (!embed_multiplication) {

		//if not using embedded convolution then allocate cuS2 scratch spaces
		if (!cuS2_x.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS2_y.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS2_z.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		cuS2_x.set(cuComplex());
		cuS2_y.set(cuComplex());
		cuS2_z.set(cuComplex());
	}

	//quarter scratch space used only if transpose mode is on
	if (transpose_xy) {

		if (!cuSquart_x.resize((N.x / 2 + 1) * n.y * n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuSquart_y.resize((N.x / 2 + 1) * n.y * n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuSquart_z.resize((N.x / 2 + 1) * n.y * n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		cuSquart_x.set(cuComplex());
		cuSquart_y.set(cuComplex());
		cuSquart_z.set(cuComplex());
	}
	else {

		cuSquart_x.clear();
		cuSquart_y.clear();
		cuSquart_z.clear();
	}

	//Others

	transpose_xy_gpu.from_cpu(transpose_xy);

	return error;
}

#endif