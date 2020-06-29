#include "stdafx.h"
#include "OerstedKernelCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_OERSTED

#include "OerstedTFunc.h"

//-------------------------- MEMORY ALLOCATION

BError OerstedKernelCUDA::AllocateKernelMemory(void)
{
	BError error(__FUNCTION__);

	if (!transpose_xy) {

		if (!KOe()->resize(cuSZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}
	else {

		if (!KOe()->resize(cuSZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	return error;
}

//-------------------------- KERNEL CALCULATION

BError OerstedKernelCUDA::Calculate_Oersted_Kernels_2D(void)
{
	BError error(__FUNCTION__);

	//2D not used for Oersted module

	return error;
}

BError OerstedKernelCUDA::Calculate_Oersted_Kernels_3D(void)
{
	BError error(__FUNCTION__);

	//-------------- CALCULATE OERSTED TENSOR

	//Oersted tensor components
	//
	// 0    K12  K13
	// -K12 0    K23
	// -K13 -K23 0

	VEC<DBL3> DOe;

	if (!DOe.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//-------------- OERSTED KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<DBL3> KOe_cpu;

	if (!transpose_xy) {

		if (!KOe_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!KOe_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- SETUP FFT

	double* pline_real = fftw_alloc_real(maximum(N.x, N.y, N.z) * 3);
	fftw_complex* pline = fftw_alloc_complex(maximum(N.x / 2 + 1, N.y, N.z) * 3);

	//make fft plans
	int dims_x[1] = { (int)N.x };
	int dims_y[1] = { (int)N.y };
	int dims_z[1] = { (int)N.z };

	fftw_plan plan_fwd_x = fftw_plan_many_dft_r2c(1, dims_x, 3,
		pline_real, nullptr, 3, 1,
		pline, nullptr, 3, 1,
		FFTW_PATIENT);

	fftw_plan plan_fwd_y = fftw_plan_many_dft_r2c(1, dims_y, 3,
		pline_real, nullptr, 3, 1,
		pline, nullptr, 3, 1,
		FFTW_PATIENT);

	fftw_plan plan_fwd_z = fftw_plan_many_dft_r2c(1, dims_z, 3,
		pline_real, nullptr, 3, 1,
		pline, nullptr, 3, 1,
		FFTW_PATIENT);

	//lambda used to transform an input real tensor into an output real kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<DBL3>& kernel) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int k = 0; k < N.z; k++) {
			for (int j = 0; j < N.y; j++) {

				//write input into fft line (zero padding kept)
				for (int i = 0; i < N.x; i++) {

					int idx_in = i + j * N.x + k * N.x * N.y;

					*reinterpret_cast<DBL3*>(pline_real + i * 3) = tensor[idx_in];
				}

				//fft on line
				fftw_execute(plan_fwd_x);

				//pack into lower half of tensor row for next step (keep same row and plane strides)
				for (int i = 0; i < N.x / 2 + 1; i++) {

					ReIm3 value = *reinterpret_cast<ReIm3*>(pline + i * 3);

					//Dxy : even x, Dxz : even x, Dyz : odd x
					tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(value.x.Re, value.y.Re, value.z.Im);
				}
			}
		}

		for (int k = 0; k < N.z; k++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				//fetch line from fft array (zero padding kept)
				for (int j = 0; j < N.y; j++) {

					int idx_in = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

					*reinterpret_cast<DBL3*>(pline_real + j * 3) = tensor[idx_in];
				}

				//fft on line
				fftw_execute(plan_fwd_y);

				//pack into lower half of tensor column for next step (keep same row and plane strides)
				for (int j = 0; j < N.y / 2 + 1; j++) {

					ReIm3 value = *reinterpret_cast<ReIm3*>(pline + j * 3);

					//Dxy : even y, Dxz : odd y, Dyz : even y
					tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(value.x.Re, value.y.Im, value.z.Re);
				}
			}
		}

		for (int j = 0; j < N.y / 2 + 1; j++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				//fetch line from fft array (zero padding kept)
				for (int k = 0; k < N.z; k++) {

					int idx_in = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

					*reinterpret_cast<DBL3*>(pline_real + k * 3) = tensor[idx_in];
				}

				//fft on line
				fftw_execute(plan_fwd_z);

				//pack into output kernels with reduced strides
				for (int k = 0; k < N.z / 2 + 1; k++) {

					ReIm3 value = *reinterpret_cast<ReIm3*>(pline + k * 3);

					//Dxy : odd z, Dxz : even z, Dyz : even z
					//Note, all tensor elements where odd exactly once, thus the output will be purely real
					//This means when we use the kernel elements we need to multiply by i since we are saving the result as real values.

					if (!transpose_xy) {

						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(value.x.Im, value.y.Re, value.z.Re);
					}
					else {

						kernel[j + i * (N.y / 2 + 1) + k * (N.y / 2 + 1) * (N.x / 2 + 1)] = DBL3(value.x.Im, value.y.Re, value.z.Re);
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	OeFunc::CalcOerstedTensors(DOe, n, N, h);

	tensor_to_kernel(DOe, KOe_cpu);

	//transfer to GPU
	KOe()->copy_from_cpuvec(KOe_cpu);

	//-------------- CLEANUP

	fftw_destroy_plan(plan_fwd_x);
	fftw_destroy_plan(plan_fwd_y);
	fftw_destroy_plan(plan_fwd_z);

	fftw_free((double*)pline_real);
	fftw_free((fftw_complex*)pline);

	//Done
	return error;
}

#endif

#endif