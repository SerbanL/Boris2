#include "stdafx.h"
#include "RoughnessKernel.h"

#ifdef MODULE_ROUGHNESS

#include "Roughness.h"

template RoughnessKernel<Roughness>::RoughnessKernel(Roughness *pOwner_);

template <typename Owner>
RoughnessKernel<Owner>::RoughnessKernel(Owner* pOwner_)
{
	pOwner = pOwner_;
}

//-------------------------- MEMORY ALLOCATION

template BError RoughnessKernel<Roughness>::AllocateKernelMemory(void);

template <typename Owner>
BError RoughnessKernel<Owner>::AllocateKernelMemory(void)
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

template void RoughnessKernel<Roughness>::KernelMultiplication_2D(VEC<ReIm3>& In, VEC<ReIm3>& Out);

template <typename Owner>
void RoughnessKernel<Owner>::KernelMultiplication_2D(VEC<ReIm3>& In, VEC<ReIm3>& Out)
{
	//Full multiplication without use of kernel symmetries : testing only
	/*
	if (!pOwner->off_diagonal_convolution) {

#pragma omp parallel for
		for (int index = 0; index < (N.x / 2 + 1)*N.y; index++) {

			ReIm3 FM = In[index];

			Out[index].x = Kdiag[index].x  * FM.x;
			Out[index].y = Kdiag[index].y  * FM.y;
			Out[index].z = Kdiag[index].z  * FM.z;
		}
	}
	else {

#pragma omp parallel for
		for (int index = 0; index < (N.x / 2 + 1)*N.y; index++) {

			ReIm3 FM = In[index];

			Out[index].x = K2D_odiag[index] * FM.x;
			Out[index].y = ReIm();
			Out[index].z = ReIm();
		}
	}
	*/
}

template void RoughnessKernel<Roughness>::KernelMultiplication_3D(VEC<ReIm3>& In, VEC<ReIm3>& Out);

template <typename Owner>
void RoughnessKernel<Owner>::KernelMultiplication_3D(VEC<ReIm3>& In, VEC<ReIm3>& Out)
{
	//Full multiplication without use of kernel symmetries : testing only
	/*
	if (!pOwner->off_diagonal_convolution) {

#pragma omp parallel for
		for (int index = 0; index < (N.x / 2 + 1)*N.y*N.z; index++) {

			ReIm3 FM = In[index];

			Out[index].x = Kdiag[index].x * FM.x;
			Out[index].y = Kdiag[index].y * FM.y;
			Out[index].z = Kdiag[index].z * FM.z;
		}
	}
	else {

#pragma omp parallel for
		for (int index = 0; index < (N.x / 2 + 1)*N.y*N.z; index++) {

			ReIm3 FM = In[index];

			Out[index].x = Kodiag[index].x * FM.x;
			Out[index].y = Kodiag[index].y * FM.y;
			Out[index].z = Kodiag[index].z * FM.z;
		}
	}
	*/
}

template void RoughnessKernel<Roughness>::KernelMultiplication_2D_line(ReIm3* pline, int i);

//multiply kernels in line along y direction (so use stride of (N.x/2 + 1) to read from kernels), starting at given i index (this must be an index in the first x row)
template <typename Owner>
void RoughnessKernel<Owner>::KernelMultiplication_2D_line(ReIm3* pline, int i)
{
	//above N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.y/2 point
	//off-diagonal values are odd about the N.y/2 point

	if (!pOwner->off_diagonal_convolution) {

		//j = 0
		ReIm3 FM = pline[0];

		int idx_start = i;

		pline[0].x = Kdiag[idx_start].x  * FM.x;
		pline[0].y = Kdiag[idx_start].y  * FM.y;
		pline[0].z = Kdiag[idx_start].z  * FM.z;

		//points between 1 and N.y / 2 + 1 inclusive
		for (int j = 1; j <= N.y / 2; j++) {

			ReIm3 FM_l = pline[j];
			ReIm3 FM_h = pline[N.y - j];

			int ker_index = i + j * (N.x / 2 + 1);

			pline[j].x = Kdiag[ker_index].x  * FM_l.x;
			pline[j].y = Kdiag[ker_index].y  * FM_l.y;
			pline[j].z = Kdiag[ker_index].z  * FM_l.z;

			pline[N.y - j].x = Kdiag[ker_index].x  * FM_h.x;
			pline[N.y - j].y = Kdiag[ker_index].y  * FM_h.y;
			pline[N.y - j].z = Kdiag[ker_index].z  * FM_h.z;
		}

		//j = N.y / 2
		FM = pline[N.y / 2];

		int idx_mid = i + (N.y / 2) * (N.x / 2 + 1);

		pline[N.y / 2].x = Kdiag[idx_mid].x  * FM.x;
		pline[N.y / 2].y = Kdiag[idx_mid].y  * FM.y;
		pline[N.y / 2].z = Kdiag[idx_mid].z  * FM.z;
	}
	else {

		//j = 0
		ReIm3 FM = pline[0];

		int idx_start = i;

		pline[0].x = K2D_odiag[idx_start] * FM.x;
		pline[0].y = ReIm();
		pline[0].z = ReIm();

		//points between 1 and N.y / 2 + 1 inclusive
		for (int j = 1; j <= N.y / 2; j++) {

			ReIm3 FM_l = pline[j];
			ReIm3 FM_h = pline[N.y - j];

			int ker_index = i + j * (N.x / 2 + 1);

			pline[j].x = K2D_odiag[ker_index] * FM_l.x;
			pline[j].y = ReIm();
			pline[j].z = ReIm();

			pline[N.y - j].x = -K2D_odiag[ker_index] * FM_h.x;
			pline[N.y - j].y = ReIm();
			pline[N.y - j].z = ReIm();
		}

		//j = N.y / 2
		FM = pline[N.y / 2];

		int idx_mid = i + (N.y / 2) * (N.x / 2 + 1);

		pline[N.y / 2].x = K2D_odiag[idx_mid] * FM.x;
		pline[N.y / 2].y = ReIm();
		pline[N.y / 2].z = ReIm();
	}
}

template void RoughnessKernel<Roughness>::KernelMultiplication_3D_line(ReIm3* pline, int i, int j);

//multiply kernels in line along z direction (so use stride of (N.x/2 + 1) * N.y to read from kernels), starting at given i and j indexes (these must be an indexes in the first xy plane)
template <typename Owner>
void RoughnessKernel<Owner>::KernelMultiplication_3D_line(ReIm3* pline, int i, int j)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//diagonal components are even about the N.z/2 and N.y/2 points
	//Kxy is even about N.z/2 and odd about N.y/2
	//Kxz is odd about N.z/2 and even about N.y/2
	//Kyz is odd about N.z/2 and odd about N.y/2

	if (!pOwner->off_diagonal_convolution) {

		if (j <= N.y / 2) {

			//k = 0
			ReIm3 FM = pline[0];

			int idx_start = i + j * (N.x / 2 + 1);

			pline[0].x = Kdiag[idx_start].x * FM.x;
			pline[0].y = Kdiag[idx_start].y * FM.y;
			pline[0].z = Kdiag[idx_start].z * FM.z;

			//points between 1 and N.z /2 - 1 inclusive
			for (int k = 1; k < N.z / 2; k++) {

				ReIm3 FM_l = pline[k];
				ReIm3 FM_h = pline[N.z - k];

				int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				pline[k].x = Kdiag[ker_index].x * FM_l.x;
				pline[k].y = Kdiag[ker_index].y * FM_l.y;
				pline[k].z = Kdiag[ker_index].z * FM_l.z;

				pline[N.z - k].x = Kdiag[ker_index].x * FM_h.x;
				pline[N.z - k].y = Kdiag[ker_index].y * FM_h.y;
				pline[N.z - k].z = Kdiag[ker_index].z * FM_h.z;
			}

			//k = N.z / 2
			FM = pline[N.z / 2];

			int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[N.z / 2].x = Kdiag[idx_mid].x * FM.x;
			pline[N.z / 2].y = Kdiag[idx_mid].y * FM.y;
			pline[N.z / 2].z = Kdiag[idx_mid].z * FM.z;
		}
		else {

			//k = 0
			ReIm3 FM = pline[0];

			int idx_start = i + (N.y - j) * (N.x / 2 + 1);

			pline[0].x = Kdiag[idx_start].x * FM.x;
			pline[0].y = Kdiag[idx_start].y * FM.y;
			pline[0].z = Kdiag[idx_start].z * FM.z;

			//points between 1 and N.z /2 - 1 inclusive
			for (int k = 1; k < N.z / 2; k++) {

				ReIm3 FM_l = pline[k];
				ReIm3 FM_h = pline[N.z - k];

				int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				pline[k].x = Kdiag[ker_index].x * FM_l.x;
				pline[k].y = Kdiag[ker_index].y * FM_l.y;
				pline[k].z = Kdiag[ker_index].z * FM_l.z;

				pline[N.z - k].x = Kdiag[ker_index].x * FM_h.x;
				pline[N.z - k].y = Kdiag[ker_index].y * FM_h.y;
				pline[N.z - k].z = Kdiag[ker_index].z * FM_h.z;
			}

			//k = N.z / 2
			FM = pline[N.z / 2];

			int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[N.z / 2].x = Kdiag[idx_mid].x * FM.x;
			pline[N.z / 2].y = Kdiag[idx_mid].y * FM.y;
			pline[N.z / 2].z = Kdiag[idx_mid].z * FM.z;
		}
	}
	else {

		if (j <= N.y / 2) {

			//k = 0
			ReIm3 FM = pline[0];

			int idx_start = i + j * (N.x / 2 + 1);

			pline[0].x = Kodiag[idx_start].x * FM.x;
			pline[0].y = Kodiag[idx_start].y * FM.y;
			pline[0].z = Kodiag[idx_start].z * FM.z;

			//points between 1 and N.z /2 - 1 inclusive
			for (int k = 1; k < N.z / 2; k++) {

				ReIm3 FM_l = pline[k];
				ReIm3 FM_h = pline[N.z - k];

				int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				pline[k].x = Kodiag[ker_index].x * FM_l.x;
				pline[k].y = Kodiag[ker_index].y * FM_l.y;
				pline[k].z = Kodiag[ker_index].z * FM_l.z;

				pline[N.z - k].x = Kodiag[ker_index].x * FM_h.x;
				pline[N.z - k].y = -Kodiag[ker_index].y * FM_h.y;
				pline[N.z - k].z = -Kodiag[ker_index].z * FM_h.z;
			}

			//k = N.z / 2
			FM = pline[N.z / 2];

			int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[N.z / 2].x = Kodiag[idx_mid].x * FM.x;
			pline[N.z / 2].y = Kodiag[idx_mid].y * FM.y;
			pline[N.z / 2].z = Kodiag[idx_mid].z * FM.z;
		}
		else {

			//k = 0
			ReIm3 FM = pline[0];

			int idx_start = i + (N.y - j) * (N.x / 2 + 1);

			pline[0].x = -Kodiag[idx_start].x * FM.x;
			pline[0].y = Kodiag[idx_start].y * FM.y;
			pline[0].z = -Kodiag[idx_start].z * FM.z;

			//points between 1 and N.z /2 - 1 inclusive
			for (int k = 1; k < N.z / 2; k++) {

				ReIm3 FM_l = pline[k];
				ReIm3 FM_h = pline[N.z - k];

				int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				pline[k].x = -Kodiag[ker_index].x * FM_l.x;
				pline[k].y = Kodiag[ker_index].y * FM_l.y;
				pline[k].z = -Kodiag[ker_index].z * FM_l.z;

				pline[N.z - k].x = -Kodiag[ker_index].x * FM_h.x;
				pline[N.z - k].y = -Kodiag[ker_index].y * FM_h.y;
				pline[N.z - k].z = Kodiag[ker_index].z * FM_h.z;
			}

			//k = N.z / 2
			FM = pline[N.z / 2];

			int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[N.z / 2].x = -Kodiag[idx_mid].x * FM.x;
			pline[N.z / 2].y = Kodiag[idx_mid].y * FM.y;
			pline[N.z / 2].z = -Kodiag[idx_mid].z * FM.z;
		}
	}
}

//-------------------------- KERNEL CALCULATION

template BError RoughnessKernel<Roughness>::Calculate_Roughness_Kernels_2D(void);

template <typename Owner>
BError RoughnessKernel<Owner>::Calculate_Roughness_Kernels_2D(void)
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

	if (!dtf.CalcDiagTens2D(Ddiag, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);
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

template BError RoughnessKernel<Roughness>::Calculate_Roughness_Kernels_3D(void);

template <typename Owner>
BError RoughnessKernel<Owner>::Calculate_Roughness_Kernels_3D(void)
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
	if (!dtf.CalcDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);

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