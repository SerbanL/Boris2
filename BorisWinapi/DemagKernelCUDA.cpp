#include "stdafx.h"
#include "DemagKernelCUDA.h"

#if COMPILECUDA == 1

#if defined MODULE_DEMAG || defined MODULE_SDEMAG

#include "DemagTFunc.h"

//-------------------------- MEMORY ALLOCATION

BError DemagKernelCUDA::AllocateKernelMemory(void)
{
	BError error(__FUNCTION__);

	//now allocate memory
	if (!transpose_xy) {

		if (n.z == 1) {

			Kodiag()->clear();

			if (!K2D_odiag()->resize(cuSZ3(N.x / 2 + 1, N.y / 2 + 1, 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Kdiag()->resize(cuSZ3(N.x / 2 + 1, N.y / 2 + 1, 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else {

			K2D_odiag()->clear();

			if (!Kodiag()->resize(cuSZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Kdiag()->resize(cuSZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
	}
	else {

		if (n.z == 1) {

			Kodiag()->clear();

			if (!K2D_odiag()->resize(cuSZ3(N.y / 2 + 1, N.x / 2 + 1, 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Kdiag()->resize(cuSZ3(N.y / 2 + 1, N.x / 2 + 1, 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else {

			K2D_odiag()->clear();
			
			if (!Kodiag()->resize(cuSZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
			if (!Kdiag()->resize(cuSZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
	}

	return error;
}

//-------------------------- KERNEL CALCULATION

BError DemagKernelCUDA::Calculate_Demag_Kernels_2D(bool include_self_demag)
{
	BError error(__FUNCTION__);

	//-------------- CALCULATE DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 0
	// D12 D22 0
	// 0   0   D33

	//Ddiag : D11, D22, D33 are the diagonal tensor elements
	VEC<DBL3> Ddiag;
	if (!Ddiag.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//off-diagonal tensor elements
	vector<double> Dodiag;
	if (!malloc_vector(Dodiag, N.x*N.y)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//use ratios instead of cellsizes directly - same result but better in terms of floating point errors
	DemagTFunc dtf;

	if (!dtf.CalcDiagTens2D(Ddiag, n, N, h / maximum(h.x, h.y, h.z), include_self_demag)) return error(BERROR_OUTOFMEMORY_NCRIT);
	if (!dtf.CalcOffDiagTens2D(Dodiag, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);

	//-------------- DEMAG KERNELS ON CPU-ADDRESSABLE MEMORY

	VEC<DBL3> Kdiag_cpu;
	vector<double> K2D_odiag_cpu;

	if (!transpose_xy) {

		if (!Kdiag_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!Kdiag_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	if (!malloc_vector(K2D_odiag_cpu, (N.x / 2 + 1) * (N.y / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);

	//-------------- SETUP FFT

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y / 2 + 1);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);
	vector<ReIm> fft_line2d(maxN);
	vector<ReIm> fft_line2d2(maxN);

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
	//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
	//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

	//FFT into Kernel forms ready for convolution multiplication
	for (int j = 0; j < N.y; j++) {

		fft.CopyRealShuffle(Ddiag.data() + j * N.x, fft_line.data(), N.x / 2);
		fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
		fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

		fft.CopyRealShuffle(Dodiag.data() + j * N.x, fft_line2d.data(), N.x / 2);
		fft.FFT_Radix4_DIT(fft_line2d.data(), log2(N.x) - 1, N.x / 2);
		fft.RealfromComplexFFT(fft_line2d.data(), fft_line2d2.data(), N.x / 2);

		//pack into Ddiag and Dodiag for next step
		for (int i = 0; i < N.x / 2 + 1; i++) {

			//even w.r.t. to x so output is purely real
			Ddiag[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[i].x.Re, fft_line2[i].y.Re, fft_line2[i].z.Re);

			//odd w.r.t. to x so output is purely imaginary
			Dodiag[i + j * (N.x / 2 + 1)] = fft_line2d2[i].Im;
		}
	}

	for (int i = 0; i < N.x / 2 + 1; i++) {

		fft.CopyRealShuffle(Ddiag.data() + i, fft_line.data(), N.x / 2 + 1, N.y / 2);
		fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y) - 1, N.y / 2);
		fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.y / 2);

		//the input sequence should actually be purely imaginary not purely real (but see below)
		fft.CopyRealShuffle(Dodiag.data() + i, fft_line2d.data(), N.x / 2 + 1, N.y / 2);
		fft.FFT_Radix4_DIT(fft_line2d.data(), log2(N.y) - 1, N.y / 2);
		fft.RealfromComplexFFT(fft_line2d.data(), fft_line2d2.data(), N.y / 2);

		//pack into output real kernels with reduced strides
		for (int j = 0; j < N.y / 2 + 1; j++) {

			if (!transpose_xy) {

				//even w.r.t. y so output is purely real
				Kdiag_cpu[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);

				//odd w.r.t. y so the purely imaginary input becomes purely real
				//however since we used CopyRealShuffle and RealfromComplexFFT, i.e. treating the input as purely real rather than purely imaginary we need to account for the i * i = -1 term, hence the - sign below
				K2D_odiag_cpu[i + j * (N.x / 2 + 1)] = -fft_line2d2[j].Im;
			}
			else {

				//even w.r.t. y so output is purely real
				Kdiag_cpu[j + i * (N.y / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);

				//odd w.r.t. y so the purely imaginary input becomes purely real
				//however since we used CopyRealShuffle and RealfromComplexFFT, i.e. treating the input as purely real rather than purely imaginary we need to account for the i * i = -1 term, hence the - sign below
				K2D_odiag_cpu[j + i * (N.y / 2 + 1)] = -fft_line2d2[j].Im;
			}
		}
	}

	//-------------- TRANSFER TO GPU KERNELS

	Kdiag()->copy_from_cpuvec(Kdiag_cpu);
	K2D_odiag()->copy_from_vector(K2D_odiag_cpu);

	//Done
	return error;
}

BError DemagKernelCUDA::Calculate_Demag_Kernels_3D(bool include_self_demag)
{
	BError error(__FUNCTION__);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	//D11, D22, D33 are the diagonal tensor elements
	//D12, D13, D23 are the off-diagonal tensor elements
	VEC<DBL3> D;

	if (!D.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//object used to compute tensor elements
	DemagTFunc dtf;

	//-------------- DEMAG KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<DBL3> K_cpu;

	if (!transpose_xy) {

		if (!K_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!K_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input real tensor into an output real kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<DBL3>& kernel, bool off_diagonal) -> void {

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

					if (!off_diagonal) {

						//even w.r.t. to x so output is purely real
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[i].x.Re, fft_line2[i].y.Re, fft_line2[i].z.Re);
					}
					else {

						//Dxy : odd x, Dxz : odd x, Dyz : even x
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[i].x.Im, fft_line2[i].y.Im, fft_line2[i].z.Re);
					}
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

					if (!off_diagonal) {

						//even w.r.t. to y so output is purely real
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);
					}
					else {

						//Dxy : odd y, Dxz : even y, Dyz : odd y
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[j].x.Im, fft_line2[j].y.Re, fft_line2[j].z.Im);
					}
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

					if (!off_diagonal) {

						//even w.r.t. to z so output is purely real

						if (!transpose_xy) {

							kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Re, fft_line2[k].y.Re, fft_line2[k].z.Re);
						}
						else {

							kernel[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Re, fft_line2[k].y.Re, fft_line2[k].z.Re);
						}
					}
					else {

						//Dxy : even z, Dxz : odd z, Dyz : odd z
						//Also multiply by -1 since all off-diagonal tensor elements have been odd twice
						//The final output is thus purely real but we always treated the input as purely real even when it should have been purely imaginary
						//This means we need to account for i * i = -1 at the end
						
						if (!transpose_xy) {

							kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(-fft_line2[k].x.Re, -fft_line2[k].y.Im, -fft_line2[k].z.Im);
						}
						else {

							kernel[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(-fft_line2[k].x.Re, -fft_line2[k].y.Im, -fft_line2[k].z.Im);
						}
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do
	if (!dtf.CalcDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z), include_self_demag)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, false);

	//transfer to GPU
	Kdiag()->copy_from_cpuvec(K_cpu);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (!dtf.CalcOffDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, true);

	//transfer to GPU
	Kodiag()->copy_from_cpuvec(K_cpu);

	//Done
	return error;
}

#endif

#endif