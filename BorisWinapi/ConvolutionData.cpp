#include "stdafx.h"
#include "ConvolutionData.h"

//-------------------------- CONSTRUCTORS

ConvolutionData::ConvolutionData(void)
{
	OmpThreads = omp_get_num_procs();
	
	plan_fwd_x.resize(OmpThreads);
	plan_fwd_y.resize(OmpThreads);
	plan_fwd_z.resize(OmpThreads);
	plan_inv_x.resize(OmpThreads);
	plan_inv_y.resize(OmpThreads);
	plan_inv_z.resize(OmpThreads);
	
	pline_zp_x.resize(OmpThreads);
	pline_zp_y.resize(OmpThreads);
	pline_zp_z.resize(OmpThreads);
	pline.resize(OmpThreads);
	pline_rev_x.resize(OmpThreads);
}

ConvolutionData::~ConvolutionData()
{
	//clean
	free_memory();
}

//-------------------------- HELPERS

void ConvolutionData::free_memory(void)
{
	if (fftw_plans_created) {

		//clean
		for (int idx = 0; idx < OmpThreads; idx++) {

			fftw_destroy_plan(plan_fwd_x[idx]);
			fftw_destroy_plan(plan_fwd_y[idx]);
			fftw_destroy_plan(plan_fwd_z[idx]);
			fftw_destroy_plan(plan_inv_x[idx]);
			fftw_destroy_plan(plan_inv_y[idx]);
			fftw_destroy_plan(plan_inv_z[idx]);

			fftw_free((double*)pline_zp_x[idx]);
			fftw_free((fftw_complex*)pline_zp_y[idx]);
			fftw_free((fftw_complex*)pline_zp_z[idx]);
			fftw_free((fftw_complex*)pline[idx]);
			fftw_free((double*)pline_rev_x[idx]);
		}
	}

	fftw_plans_created = false;
}

//Allocate memory for F and F2 (if needed) scratch spaces)
BError ConvolutionData::AllocateScratchSpaces(void)
{
	BError error(__FUNCTION__);

	if (embed_multiplication) {

		//if multiplication is embedded, we don't need the upper z-axis points (3D) or upper y-axis points (2D). We also don't need the F2 scratch space.

		if (n.z > 1) {

			//3D : don't need to allocate up to N.z as kernel multiplication is embedded with the z-direction fft / ifft
			if (!F.resize(SZ3(N.x / 2 + 1, N.y, n.z))) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else {

			//2D : don't need to allocate up to N.y as kernel multiplication is embedded with the y-direction fft / ifft
			if (!F.resize(SZ3(N.x / 2 + 1, n.y, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
		}

		F2.clear();
	}
	else {

		//if multiplication is not embedded, we need full spaces, and also the F2 scratch space.

		if (n.z > 1) {

			if (!F.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
			if (!F2.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else {

			if (!F.resize(SZ3(N.x / 2 + 1, N.y, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
			if (!F2.resize(SZ3(N.x / 2 + 1, N.y, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
		}
	}

	return error;
}

//-------------------------- CONFIGURATION

BError ConvolutionData::SetConvolutionDimensions(SZ3 n_, DBL3 h_, bool embed_multiplication_, INT3 pbc_images)
{
	BError error(__FUNCTION__);

	this->pbc_images = pbc_images;

	//Turn off multiplication embedding by setting this to false - this will allocate full space memory for F and F2 scratch spaces. embed_multiplication = true by default.
	embed_multiplication = embed_multiplication_;

	n = n_;
	h = h_;

	N = SZ3(1, 1, 1);

	//set N, M, K as smallest powers of 2 which will hold the demag kernel
	if (pbc_images.x) {

		//pbc : can use wrap-around thus half the size required
		while (N.x < n.x) { N.x *= 2; };
	}
	else {

		//no wrap-around thus double the size, with input to be zero-padded
		while (N.x < 2 * n.x) { N.x *= 2; };
	}

	if (pbc_images.y) {

		while (N.y < n.y) { N.y *= 2; };
	}
	else {

		while (N.y < 2 * n.y) { N.y *= 2; };
	}

	if (n.z > 1) {

		if (pbc_images.z) {

			while (N.z < n.z) { N.z *= 2; };
		}
		else {

			while (N.z < 2 * n.z) { N.z *= 2; };
		}
	}

	//now allocate memory

	//scratch spaces
	error = AllocateScratchSpaces();
	if (error) return error;

	//setup fftw for convolution
	
	//first clean any previously allocated memory
	free_memory();

	size_t maxN = maximum(N.x / 2 + 1, N.y, N.z);

	//allocate new fft lines
	for (int idx = 0; idx < OmpThreads; idx++) {

		pline_zp_x[idx] = fftw_alloc_real(N.x * 3);
		pline_rev_x[idx] = fftw_alloc_real(N.x * 3);

		pline_zp_y[idx] = fftw_alloc_complex(N.y * 3);
		pline_zp_z[idx] = fftw_alloc_complex(N.z * 3);

		pline[idx] = fftw_alloc_complex(maxN * 3);
	}

	//zero fft lines
	zero_fft_lines();

	//make fft plans
	int dims_x[1] = { (int)N.x };
	int dims_y[1] = { (int)N.y };
	int dims_z[1] = { (int)N.z };

	for (int idx = 0; idx < OmpThreads; idx++) {

		plan_fwd_x[idx] = fftw_plan_many_dft_r2c(1, dims_x, 3,
			pline_zp_x[idx], nullptr, 3, 1,
			pline[idx], nullptr, 3, 1,
			FFTW_PATIENT);

		plan_fwd_y[idx] = fftw_plan_many_dft(1, dims_y, 3,
			pline_zp_y[idx], nullptr, 3, 1,
			pline[idx], nullptr, 3, 1,
			FFTW_FORWARD, FFTW_PATIENT);

		plan_fwd_z[idx] = fftw_plan_many_dft(1, dims_z, 3,
			pline_zp_z[idx], nullptr, 3, 1,
			pline[idx], nullptr, 3, 1,
			FFTW_FORWARD, FFTW_PATIENT);

		plan_inv_z[idx] = fftw_plan_many_dft(1, dims_z, 3,
			pline[idx], nullptr, 3, 1,
			pline[idx], nullptr, 3, 1,
			FFTW_BACKWARD, FFTW_PATIENT);

		plan_inv_y[idx] = fftw_plan_many_dft(1, dims_y, 3,
			pline[idx], nullptr, 3, 1,
			pline[idx], nullptr, 3, 1,
			FFTW_BACKWARD, FFTW_PATIENT);

		plan_inv_x[idx] = fftw_plan_many_dft_c2r(1, dims_x, 3,
			pline[idx], nullptr, 3, 1,
			pline_rev_x[idx], nullptr, 3, 1,
			FFTW_PATIENT);
	}
	
	fftw_plans_created = true;

	return error;
}

//zero fftw memory
void ConvolutionData::zero_fft_lines(void)
{
	size_t maxN = maximum(N.x / 2 + 1, N.y, N.z);

	for (int idx = 0; idx < OmpThreads; idx++) {

		for (int i = 0; i < N.x; i++) {

			*reinterpret_cast<DBL3*>(pline_zp_x[idx] + i * 3) = DBL3();
			*reinterpret_cast<DBL3*>(pline_rev_x[idx] + i * 3) = DBL3();
		}

		for (int j = 0; j < N.y; j++) {

			*reinterpret_cast<ReIm3*>(pline_zp_y[idx] + j * 3) = ReIm3();
		}

		for (int k = 0; k < N.z; k++) {

			*reinterpret_cast<ReIm3*>(pline_zp_z[idx] + k * 3) = ReIm3();
		}

		for (int i = 0; i < maxN; i++) {

			*reinterpret_cast<ReIm3*>(pline[idx] + i * 3) = ReIm3();
		}
	}
}

//-------------------------- RUN-TIME METHODS
