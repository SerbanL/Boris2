#include "stdafx.h"
#include "DemagKernel.h"

#if defined MODULE_DEMAG || defined MODULE_SDEMAG

//-------------------------- MEMORY ALLOCATION

BError DemagKernel::AllocateKernelMemory(void)
{
	BError error(__FUNCTION__);

	if (n.z == 1) {

		//2D
		Kodiag.clear();

		if (!Kdiag.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!malloc_vector(K2D_odiag, (N.x / 2 + 1) * (N.y / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		//3D
		K2D_odiag.clear();
		K2D_odiag.shrink_to_fit();

		if (!Kdiag.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!Kodiag.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	return error;
}

//-------------------------- RUN-TIME KERNEL MULTIPLICATION

//multiply kernels in line along y direction (so use stride of (N.x/2 + 1) to read from kernels), starting at given i index (this must be an index in the first x row)
void DemagKernel::KernelMultiplication_2D_line(ReIm3* pline, int i)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	//j = 0
	ReIm3 FM = pline[0];

	int idx_start = i;

	pline[0].x = (Kdiag[idx_start].x  * FM.x) + (K2D_odiag[idx_start] * FM.y);
	pline[0].y = (K2D_odiag[idx_start] * FM.x) + (Kdiag[idx_start].y  * FM.y);
	pline[0].z = (Kdiag[idx_start].z  * FM.z);

	//points between 1 and N.y / 2 - 1 inclusive
	for (int j = 1; j < N.y / 2; j++) {

		ReIm3 FM_l = pline[j];
		ReIm3 FM_h = pline[N.y - j];

		int ker_index = i + j * (N.x / 2 + 1);

		pline[j].x = (Kdiag[ker_index].x  * FM_l.x) + (K2D_odiag[ker_index] * FM_l.y);
		pline[j].y = (K2D_odiag[ker_index] * FM_l.x) + (Kdiag[ker_index].y  * FM_l.y);
		pline[j].z = (Kdiag[ker_index].z  * FM_l.z);

		pline[N.y - j].x = (Kdiag[ker_index].x  * FM_h.x) + (-K2D_odiag[ker_index] * FM_h.y);
		pline[N.y - j].y = (-K2D_odiag[ker_index] * FM_h.x) + (Kdiag[ker_index].y  * FM_h.y);
		pline[N.y - j].z = (Kdiag[ker_index].z  * FM_h.z);
	}

	//j = N.y / 2
	FM = pline[N.y / 2];

	int idx_mid = i + (N.y / 2) * (N.x / 2 + 1);

	pline[N.y / 2].x = (Kdiag[idx_mid].x  * FM.x) + (K2D_odiag[idx_mid] * FM.y);
	pline[N.y / 2].y = (K2D_odiag[idx_mid] * FM.x) + (Kdiag[idx_mid].y  * FM.y);
	pline[N.y / 2].z = (Kdiag[idx_mid].z  * FM.z);
}

//multiply kernels in line along z direction (so use stride of (N.x/2 + 1) * N.y to read from kernels), starting at given i and j indexes (these must be an indexes in the first xy plane)
void DemagKernel::KernelMultiplication_3D_line(ReIm3* pline, int i, int j)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.z/2 and N.y/2 points
	//Kxy is even about N.z/2 and odd about N.y/2
	//Kxz is odd about N.z/2 and even about N.y/2
	//Kyz is odd about N.z/2 and odd about N.y/2

	if (j <= N.y / 2) {

		//k = 0
		ReIm3 FM = pline[0];

		int idx_start = i + j * (N.x / 2 + 1);

		pline[0].x = (Kdiag[idx_start].x * FM.x) + (Kodiag[idx_start].x * FM.y) + (Kodiag[idx_start].y * FM.z);
		pline[0].y = (Kodiag[idx_start].x * FM.x) + (Kdiag[idx_start].y * FM.y) + (Kodiag[idx_start].z * FM.z);
		pline[0].z = (Kodiag[idx_start].y * FM.x) + (Kodiag[idx_start].z * FM.y) + (Kdiag[idx_start].z * FM.z);

		//points between 1 and N.z /2 - 1 inclusive
		for (int k = 1; k < N.z / 2; k++) {

			ReIm3 FM_l = pline[k];
			ReIm3 FM_h = pline[N.z - k];

			int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[k].x = (Kdiag[ker_index].x * FM_l.x) + (Kodiag[ker_index].x * FM_l.y) + (Kodiag[ker_index].y * FM_l.z);
			pline[k].y = (Kodiag[ker_index].x * FM_l.x) + (Kdiag[ker_index].y * FM_l.y) + (Kodiag[ker_index].z * FM_l.z);
			pline[k].z = (Kodiag[ker_index].y * FM_l.x) + (Kodiag[ker_index].z * FM_l.y) + (Kdiag[ker_index].z * FM_l.z);

			pline[N.z - k].x = (Kdiag[ker_index].x * FM_h.x) + (Kodiag[ker_index].x * FM_h.y) + (-Kodiag[ker_index].y * FM_h.z);
			pline[N.z - k].y = (Kodiag[ker_index].x * FM_h.x) + (Kdiag[ker_index].y * FM_h.y) + (-Kodiag[ker_index].z * FM_h.z);
			pline[N.z - k].z = (-Kodiag[ker_index].y * FM_h.x) + (-Kodiag[ker_index].z * FM_h.y) + (Kdiag[ker_index].z * FM_h.z);
		}

		//k = N.z / 2
		FM = pline[N.z / 2];

		int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

		pline[N.z / 2].x = (Kdiag[idx_mid].x * FM.x) + (Kodiag[idx_mid].x * FM.y) + (Kodiag[idx_mid].y * FM.z);
		pline[N.z / 2].y = (Kodiag[idx_mid].x * FM.x) + (Kdiag[idx_mid].y * FM.y) + (Kodiag[idx_mid].z * FM.z);
		pline[N.z / 2].z = (Kodiag[idx_mid].y * FM.x) + (Kodiag[idx_mid].z * FM.y) + (Kdiag[idx_mid].z * FM.z);
	}
	else {

		//k = 0
		ReIm3 FM = pline[0];

		int idx_start = i + (N.y - j) * (N.x / 2 + 1);

		pline[0].x = (Kdiag[idx_start].x * FM.x) + (-Kodiag[idx_start].x * FM.y) + (Kodiag[idx_start].y * FM.z);
		pline[0].y = (-Kodiag[idx_start].x * FM.x) + (Kdiag[idx_start].y * FM.y) + (-Kodiag[idx_start].z * FM.z);
		pline[0].z = (Kodiag[idx_start].y * FM.x) + (-Kodiag[idx_start].z * FM.y) + (Kdiag[idx_start].z * FM.z);

		//points between 1 and N.z /2 - 1 inclusive
		for (int k = 1; k < N.z / 2; k++) {

			ReIm3 FM_l = pline[k];
			ReIm3 FM_h = pline[N.z - k];

			int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[k].x = (Kdiag[ker_index].x * FM_l.x) + (-Kodiag[ker_index].x * FM_l.y) + (Kodiag[ker_index].y * FM_l.z);
			pline[k].y = (-Kodiag[ker_index].x * FM_l.x) + (Kdiag[ker_index].y * FM_l.y) + (-Kodiag[ker_index].z * FM_l.z);
			pline[k].z = (Kodiag[ker_index].y * FM_l.x) + (-Kodiag[ker_index].z * FM_l.y) + (Kdiag[ker_index].z * FM_l.z);

			pline[N.z - k].x = (Kdiag[ker_index].x * FM_h.x) + (-Kodiag[ker_index].x * FM_h.y) + (-Kodiag[ker_index].y * FM_h.z);
			pline[N.z - k].y = (-Kodiag[ker_index].x * FM_h.x) + (Kdiag[ker_index].y * FM_h.y) + (Kodiag[ker_index].z * FM_h.z);
			pline[N.z - k].z = (-Kodiag[ker_index].y * FM_h.x) + (Kodiag[ker_index].z * FM_h.y) + (Kdiag[ker_index].z * FM_h.z);
		}

		//k = N.z / 2
		FM = pline[N.z / 2];

		int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

		pline[N.z / 2].x = (Kdiag[idx_mid].x * FM.x) + (-Kodiag[idx_mid].x * FM.y) + (Kodiag[idx_mid].y * FM.z);
		pline[N.z / 2].y = (-Kodiag[idx_mid].x * FM.x) + (Kdiag[idx_mid].y * FM.y) + (-Kodiag[idx_mid].z * FM.z);
		pline[N.z / 2].z = (Kodiag[idx_mid].y * FM.x) + (-Kodiag[idx_mid].z * FM.y) + (Kdiag[idx_mid].z * FM.z);
	}
}

//-------------------------- KERNEL CALCULATION

BError DemagKernel::Calculate_Demag_Kernels_2D(bool include_self_demag)
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

			//even w.r.t. y so output is purely real
			Kdiag[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);
			
			//odd w.r.t. y so the purely imaginary input becomes purely real
			//however since we used CopyRealShuffle and RealfromComplexFFT, i.e. treating the input as purely real rather than purely imaginary we need to account for the i * i = -1 term, hence the - sign below
			K2D_odiag[i + j * (N.x / 2 + 1)] = -fft_line2d2[j].Im;
		}
	}

	return error;
}

BError DemagKernel::Calculate_Demag_Kernels_3D(bool include_self_demag)
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
						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Re, fft_line2[k].y.Re, fft_line2[k].z.Re);
					}
					else {

						//Dxy : even z, Dxz : odd z, Dyz : odd z
						//Also multiply by -1 since all off-diagonal tensor elements have been odd twice
						//The final output is thus purely real but we always treated the input as purely real even when it should have been purely imaginary
						//This means we need to account for i * i = -1 at the end
						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(-fft_line2[k].x.Re, -fft_line2[k].y.Im, -fft_line2[k].z.Im);
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do
	if (!dtf.CalcDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z), include_self_demag)) return error(BERROR_OUTOFMEMORY_NCRIT);
	
	tensor_to_kernel(D, Kdiag, false);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//important to zero D before calculating new tensor elements
	D.set(DBL3());

	if (!dtf.CalcOffDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, Kodiag, true);

	//Done
	return error;
}

#endif