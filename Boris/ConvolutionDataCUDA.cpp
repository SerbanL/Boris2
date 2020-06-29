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

BError ConvolutionDataCUDA::SetConvolutionDimensions(cuSZ3 n_, cuReal3 h_, bool embed_multiplication_, cuINT3 pbc_images)
{
	BError error(__FUNCTION__);

	this->pbc_images = pbc_images;

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
	
	//set N values for FFT dimensions
	if (pbc_images.x) {

		//pbc : can use wrap-around, but currently only even values of N are allowed. Thus if n is odd, N has an extra cell (user should be warned in this case to use only even values for n in pbc directions).
		N.x = n.x + (n.x % 2);
	}
	else {

		//no wrap-around thus double the size, with input to be zero-padded
		N.x = 2 * n.x;
	}

	if (pbc_images.y) {

		N.y = n.y + (n.y % 2);
	}
	else {

		N.y = 2 * n.y;
	}

	if (n.z > 1) {

		if (pbc_images.z) {

			N.z = n.z + (n.z % 2);

			//in z pbc mode we'll want to disable q2d mode as that only works for powers of 2 with N.z at least 2 times larger than n.z
			//if N.z is a power of 2 could adapt a q2D type mode for n.z = N.z but not worth the effort currently - is z pbc mode really that useful?
		}
		else {
			
			//if not in pbc mode we want to stick to powers of 2 here to take advantage of q2D mode
			while (N.z < 2 * n.z) { N.z = N.z * 2; };
		}
	}

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
		
		if (pbc_images.z) {

			//disable q2D mode if using z pbc
			q2D_level = 0;
		}
		else {

			//quasi 2D mode? (not applicable if convolution not embedded
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
	cuIn_x.set(cuBReal(0.0));
	cuIn_y.set(cuBReal(0.0));
	cuIn_z.set(cuBReal(0.0));

	if (!cuOut_x.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuOut_y.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuOut_z.resize(N.x*n.y*n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	cuOut_x.set(cuBReal(0.0));
	cuOut_y.set(cuBReal(0.0));
	cuOut_z.set(cuBReal(0.0));
	
	//full-size scratch space
	if (!q2D_level) {

		//applies for 3D problems as well as 2D problems (n.z = 1 and thus N.z = 1)

		if (!cuS_x.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_y.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_z.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		cuS_x.set(cuBComplex());
		cuS_y.set(cuBComplex());
		cuS_z.set(cuBComplex());
	}
	else {

		//in q2D mode we can halve the space required for cuS scratch spaces.

		if (!cuS_x.resize((N.x / 2 + 1) * N.y * N.z / 2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_y.resize((N.x / 2 + 1) * N.y * N.z / 2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS_z.resize((N.x / 2 + 1) * N.y * N.z / 2)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		cuS_x.set(cuBComplex());
		cuS_y.set(cuBComplex());
		cuS_z.set(cuBComplex());
	}

	if (!embed_multiplication) {

		//if not using embedded convolution then allocate cuS2 scratch spaces
		if (!cuS2_x.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS2_y.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuS2_z.resize((N.x / 2 + 1) * N.y * N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		cuS2_x.set(cuBComplex());
		cuS2_y.set(cuBComplex());
		cuS2_z.set(cuBComplex());
	}

	//quarter scratch space used only if transpose mode is on
	if (transpose_xy) {

		if (!cuSquart_x.resize((N.x / 2 + 1) * n.y * n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuSquart_y.resize((N.x / 2 + 1) * n.y * n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!cuSquart_z.resize((N.x / 2 + 1) * n.y * n.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		cuSquart_x.set(cuBComplex());
		cuSquart_y.set(cuBComplex());
		cuSquart_z.set(cuBComplex());
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