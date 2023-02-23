#include "stdafx.h"
#include "RoughnessKernel.h"

#ifdef MODULE_COMPILATION_ROUGHNESS

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
	std::vector<double> Dodiag;
	if (!malloc_vector(Dodiag, N.x*N.y)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//use ratios instead of cellsizes directly - same result but better in terms of floating point errors
	DemagTFunc dtf;

	if (pbc_images.IsNull()) {

		//no pbcs in any dimension -> use these
		if (!dtf.CalcDiagTens2D(Ddiag, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);
		if (!dtf.CalcOffDiagTens2D(Dodiag, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);
	}
	else {

		//pbcs used in at least one dimension
		if (!dtf.CalcDiagTens2D_PBC(
			Ddiag, N, h / maximum(h.x, h.y, h.z), 
			true, ASYMPTOTIC_DISTANCE, 
			pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFMEMORY_NCRIT);

		if (!dtf.CalcOffDiagTens2D_PBC(
			Dodiag, N, h / maximum(h.x, h.y, h.z), 
			true, ASYMPTOTIC_DISTANCE, 
			pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	double* pline_real = fftw_alloc_real(maximum(N.x, N.y, N.z) * 3);
	double* pline_real_odiag = fftw_alloc_real(maximum(N.x, N.y, N.z));
	fftw_complex* pline_odiag = fftw_alloc_complex(maximum(N.x / 2 + 1, N.y, N.z));
	fftw_complex* pline = fftw_alloc_complex(maximum(N.x / 2 + 1, N.y, N.z) * 3);

	//make fft plans
	int dims_x[1] = { (int)N.x };
	int dims_y[1] = { (int)N.y };

	fftw_plan plan_fwd_x = fftw_plan_many_dft_r2c(1, dims_x, 3,
		pline_real, nullptr, 3, 1,
		pline, nullptr, 3, 1,
		FFTW_PATIENT);

	fftw_plan plan_fwd_x_odiag = fftw_plan_many_dft_r2c(1, dims_x, 1,
		pline_real_odiag, nullptr, 1, 1,
		pline_odiag, nullptr, 1, 1,
		FFTW_PATIENT);

	fftw_plan plan_fwd_y = fftw_plan_many_dft_r2c(1, dims_y, 3,
		pline_real, nullptr, 3, 1,
		pline, nullptr, 3, 1,
		FFTW_PATIENT);

	fftw_plan plan_fwd_y_odiag = fftw_plan_many_dft_r2c(1, dims_y, 1,
		pline_real_odiag, nullptr, 1, 1,
		pline_odiag, nullptr, 1, 1,
		FFTW_PATIENT);

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	//1. FFTs along x
	for (int j = 0; j < N.y; j++) {

		//write input into fft line
		for (int i = 0; i < N.x; i++) {

			int idx_in = i + j * N.x;

			*reinterpret_cast<DBL3*>(pline_real + i * 3) = Ddiag[idx_in];
			*reinterpret_cast<double*>(pline_real_odiag + i) = Dodiag[idx_in];
		}

		//fft on line
		fftw_execute(plan_fwd_x);
		fftw_execute(plan_fwd_x_odiag);

		//pack into tensor for next step
		for (int i = 0; i < N.x / 2 + 1; i++) {

			ReIm3 value = *reinterpret_cast<ReIm3*>(pline + i * 3);
			ReIm value_odiag = *reinterpret_cast<ReIm*>(pline_odiag + i);

			Ddiag[i + j * (N.x / 2 + 1)] = DBL3(value.x.Re, value.y.Re, value.z.Re);
			Dodiag[i + j * (N.x / 2 + 1)] = value_odiag.Im;
		}
	}

	//2. FFTs along y
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		//fetch line from array
		for (int j = 0; j < N.y; j++) {

			*reinterpret_cast<DBL3*>(pline_real + j * 3) = Ddiag[i + j * (N.x / 2 + 1)];
			*reinterpret_cast<double*>(pline_real_odiag + j) = Dodiag[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_fwd_y);
		fftw_execute(plan_fwd_y_odiag);

		//pack into output real kernels with reduced strides
		for (int j = 0; j < N.y / 2 + 1; j++) {

			ReIm3 value = *reinterpret_cast<ReIm3*>(pline + j * 3);
			ReIm value_odiag = *reinterpret_cast<ReIm*>(pline_odiag + j);

			//even w.r.t. y so output is purely real
			Kdiag[i + j * (N.x / 2 + 1)] = DBL3(value.x.Re, value.y.Re, value.z.Re);

			//odd w.r.t. y so the purely imaginary input becomes purely real
			//however since we used CopyRealShuffle and RealfromComplexFFT, i.e. treating the input as purely real rather than purely imaginary we need to account for the i * i = -1 term, hence the - sign below
			K2D_odiag[i + j * (N.x / 2 + 1)] = -value_odiag.Im;
		}
	}

	//-------------- CLEANUP

	fftw_destroy_plan(plan_fwd_x);
	fftw_destroy_plan(plan_fwd_x_odiag);
	fftw_destroy_plan(plan_fwd_y);
	fftw_destroy_plan(plan_fwd_y_odiag);

	fftw_free((double*)pline_real);
	fftw_free((double*)pline_real_odiag);
	fftw_free((fftw_complex*)pline);

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

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	//lambda used to transform an input tensor into an output kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<DBL3>& kernel, bool off_diagonal) -> void {

		//1. FFTs along x
		for (int k = 0; k < N.z; k++) {
			for (int j = 0; j < N.y; j++) {

				//write input into fft line (zero padding kept)
				for (int i = 0; i < N.x; i++) {

					int idx_in = i + j * N.x + k * N.x * N.y;

					*reinterpret_cast<DBL3*>(pline_real + i * 3) = tensor[idx_in];
				}

				//fft on line
				fftw_execute(plan_fwd_x);

				//pack into tensor for next step
				for (int i = 0; i < N.x / 2 + 1; i++) {

					ReIm3 value = *reinterpret_cast<ReIm3*>(pline + i * 3);

					if (!off_diagonal) {

						//even w.r.t. to x so output is purely real
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(value.x.Re, value.y.Re, value.z.Re);
					}
					else {

						//Dxy : odd x, Dxz : odd x, Dyz : even x
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(value.x.Im, value.y.Im, value.z.Re);
					}
				}
			}
		}

		//2. FFTs along y
		for (int k = 0; k < N.z; k++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

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

					if (!off_diagonal) {

						//even w.r.t. to y so output is purely real
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(value.x.Re, value.y.Re, value.z.Re);
					}
					else {

						//Dxy : odd y, Dxz : even y, Dyz : odd y
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(value.x.Im, value.y.Re, value.z.Im);
					}
				}
			}
		}

		//3. FFTs along z
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

					if (!off_diagonal) {

						//even w.r.t. to z so output is purely real
						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(value.x.Re, value.y.Re, value.z.Re);
					}
					else {

						//Dxy : even z, Dxz : odd z, Dyz : odd z
						//Also multiply by -1 since all off-diagonal tensor elements have been odd twice
						//The final output is thus purely real but we always treated the input as purely real even when it should have been purely imaginary
						//This means we need to account for i * i = -1 at the end
						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(-value.x.Re, -value.y.Im, -value.z.Im);
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do

	if (pbc_images.IsNull()) {

		//no pbcs in any dimension -> use these
		if (!dtf.CalcDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);
	}
	else {

		//pbcs used in at least one dimension
		if (!dtf.CalcDiagTens3D_PBC(D, N, h / maximum(h.x, h.y, h.z), true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFMEMORY_NCRIT);
	}

	tensor_to_kernel(D, Kdiag, false);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (pbc_images.IsNull()) {

		//no pbcs in any dimension -> use these
		if (!dtf.CalcOffDiagTens3D(D, n, N, h / maximum(h.x, h.y, h.z))) return error(BERROR_OUTOFMEMORY_NCRIT);
	}
	else {

		//pbcs used in at least one dimension
		if (!dtf.CalcOffDiagTens3D_PBC(D, N, h / maximum(h.x, h.y, h.z), true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFMEMORY_NCRIT);
	}

	tensor_to_kernel(D, Kodiag, true);

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