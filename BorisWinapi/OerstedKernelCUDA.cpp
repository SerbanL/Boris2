#include "stdafx.h"
#include "OerstedKernelCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_OERSTED

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

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input real tensor into an output real kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<DBL3>& kernel) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int k = 0; k < N.z; k++) {
			for (int j = 0; j < N.y; j++) {

				fft.CopyRealShuffle(tensor.data() + j * N.x + k * N.x*N.y, fft_line.data(), N.x / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

				//pack into lower half of tensor row for next step (keep same row and plane strides)
				for (int i = 0; i < N.x / 2 + 1; i++) {

					//Dxy : even x, Dxz : even x, Dyz : odd x
					tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[i].x.Re, fft_line2[i].y.Re, fft_line2[i].z.Im);
				}
			}
		}

		for (int k = 0; k < N.z; k++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyRealShuffle(tensor.data() + i + k * (N.x / 2 + 1)*N.y, fft_line.data(), N.x / 2 + 1, N.y / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y) - 1, N.y / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.y / 2);

				//pack into lower half of tensor column for next step (keep same row and plane strides)
				for (int j = 0; j < N.y / 2 + 1; j++) {

					//Dxy : even y, Dxz : odd y, Dyz : even y
					tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Im, fft_line2[j].z.Re);
				}
			}
		}

		for (int j = 0; j < N.y / 2 + 1; j++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyRealShuffle(tensor.data() + i + j * (N.x / 2 + 1), fft_line.data(), (N.x / 2 + 1)*N.y, N.z / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.z) - 1, N.z / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.z / 2);

				//pack into output kernels with reduced strides
				for (int k = 0; k < N.z / 2 + 1; k++) {

					//Dxy : odd z, Dxz : even z, Dyz : even z
					//Note, all tensor elements where odd exactly once, thus the output will be purely real
					//This means when we use the kernel elements we need to multiply by i since we are saving the result as real values.

					if (!transpose_xy) {
						
						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Im, fft_line2[k].y.Re, fft_line2[k].z.Re);
					}
					else {

						kernel[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Im, fft_line2[k].y.Re, fft_line2[k].z.Re);
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

	//Done
	return error;
}

#endif

#endif