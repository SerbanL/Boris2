#include "stdafx.h"
#include "OerstedKernel.h"

#ifdef MODULE_COMPILATION_OERSTED

//-------------------------- MEMORY ALLOCATION

BError OerstedKernel::AllocateKernelMemory(void)
{
	BError error(__FUNCTION__);

	if (!KOe.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);

	return error;
}

//-------------------------- RUN-TIME KERNEL MULTIPLICATION

void OerstedKernel::KernelMultiplication_2D(VEC<ReIm3>& In, VEC<ReIm3>& Out)
{
	//2D not used for Oersted module
}

void OerstedKernel::KernelMultiplication_3D(VEC<ReIm3>& In, VEC<ReIm3>& Out)
{
	//Full multiplication without use of kernel symmetries : testing only
	/*
#pragma omp parallel for
	for (int index = 0; index < (N.x / 2 + 1) * N.y * N.z; index++) {

		ReIm3 FM = In[index];

		Out[index].x = !((KOe[index].x * FM.y) + (KOe[index].y * FM.z));
		Out[index].y = !((KOe[index].x * FM.x * (-1)) + (KOe[index].z * FM.z));
		Out[index].z = !((KOe[index].y * FM.x * (-1)) + (KOe[index].z * FM.y * (-1)));
	}
	*/
}

//multiply kernels in line along z direction (so use stride of (N.x/2 + 1) * N.y to read from kernels), starting at given idx_start (this must be an index in the first xy plane)
void OerstedKernel::KernelMultiplication_3D_line(ReIm3* pline, int i, int j)
{
	//above N.z/2 and N.y/2 use kernel symmetries to recover kernel values
	//Kxy is odd about N.z/2 and even about N.y/2
	//Kxz is even about N.z/2 and odd about N.y/2
	//Kyz is even about N.z/2 and even about N.y/2

	if (j <= N.y / 2) {

		//k = 0
		ReIm3 FM = pline[0];

		int idx_start = i + j * (N.x / 2 + 1);

		pline[0].x = !((KOe[idx_start].x * FM.y) + (KOe[idx_start].y * FM.z));
		pline[0].y = !((-KOe[idx_start].x * FM.x) + (KOe[idx_start].z * FM.z));
		pline[0].z = !((-KOe[idx_start].y * FM.x) + (-KOe[idx_start].z * FM.y));

		//points between 1 and N.z /2 - 1 inclusive
		for (int k = 1; k < N.z / 2; k++) {

			ReIm3 FM_l = pline[k];
			ReIm3 FM_h = pline[N.z - k];

			int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[k].x = !((KOe[ker_index].x * FM_l.y) + (KOe[ker_index].y * FM_l.z));
			pline[k].y = !((-KOe[ker_index].x * FM_l.x) + (KOe[ker_index].z * FM_l.z));
			pline[k].z = !((-KOe[ker_index].y * FM_l.x) + (-KOe[ker_index].z * FM_l.y));

			pline[N.z - k].x = !((-KOe[ker_index].x * FM_h.y) + (KOe[ker_index].y * FM_h.z));
			pline[N.z - k].y = !((KOe[ker_index].x * FM_h.x) + (KOe[ker_index].z * FM_h.z));
			pline[N.z - k].z = !((-KOe[ker_index].y * FM_h.x) + (-KOe[ker_index].z * FM_h.y));
		}

		//k = N.z / 2
		FM = pline[N.z / 2];

		int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

		pline[N.z / 2].x = !((KOe[idx_mid].x * FM.y) + (KOe[idx_mid].y * FM.z));
		pline[N.z / 2].y = !((-KOe[idx_mid].x * FM.x) + (KOe[idx_mid].z * FM.z));
		pline[N.z / 2].z = !((-KOe[idx_mid].y * FM.x) + (-KOe[idx_mid].z * FM.y));
	}
	else {

		//k = 0
		ReIm3 FM = pline[0];

		int idx_start = i + (N.y - j) * (N.x / 2 + 1);

		pline[0].x = !((KOe[idx_start].x * FM.y) + (-KOe[idx_start].y * FM.z));
		pline[0].y = !((-KOe[idx_start].x * FM.x) + (KOe[idx_start].z * FM.z));
		pline[0].z = !((KOe[idx_start].y * FM.x) + (-KOe[idx_start].z * FM.y));

		//points between 1 and N.z /2 - 1 inclusive
		for (int k = 1; k < N.z / 2; k++) {

			ReIm3 FM_l = pline[k];
			ReIm3 FM_h = pline[N.z - k];

			int ker_index = idx_start + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			pline[k].x = !((KOe[ker_index].x * FM_l.y) + (-KOe[ker_index].y * FM_l.z));
			pline[k].y = !((-KOe[ker_index].x * FM_l.x) + (KOe[ker_index].z * FM_l.z));
			pline[k].z = !((KOe[ker_index].y * FM_l.x) + (-KOe[ker_index].z * FM_l.y));

			pline[N.z - k].x = !((-KOe[ker_index].x * FM_h.y) + (-KOe[ker_index].y * FM_h.z));
			pline[N.z - k].y = !((KOe[ker_index].x * FM_h.x) + (KOe[ker_index].z * FM_h.z));
			pline[N.z - k].z = !((KOe[ker_index].y * FM_h.x) + (-KOe[ker_index].z * FM_h.y));
		}

		//k = N.z / 2
		FM = pline[N.z / 2];

		int idx_mid = idx_start + (N.z / 2) * (N.x / 2 + 1) * (N.y / 2 + 1);

		pline[N.z / 2].x = !((KOe[idx_mid].x * FM.y) + (-KOe[idx_mid].y * FM.z));
		pline[N.z / 2].y = !((-KOe[idx_mid].x * FM.x) + (KOe[idx_mid].z * FM.z));
		pline[N.z / 2].z = !((KOe[idx_mid].y * FM.x) + (-KOe[idx_mid].z * FM.y));
	}
}


//-------------------------- KERNEL CALCULATION

BError OerstedKernel::Calculate_Oersted_Kernels_2D(void)
{
	BError error(__FUNCTION__);

	//2D not used for Oersted module

	return error;
}

BError OerstedKernel::Calculate_Oersted_Kernels_3D(void)
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
					kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(value.x.Im, value.y.Re, value.z.Re);
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	OeFunc::CalcOerstedTensors(DOe, n, N, h);

	tensor_to_kernel(DOe, KOe);

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