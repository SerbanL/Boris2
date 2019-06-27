#include "stdafx.h"
#include "OerstedKernel.h"

#ifdef MODULE_OERSTED

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
			pline[N.z - k].y = !((KOe[ker_index].x * FM_h.x) + (-KOe[ker_index].z * FM_h.z));
			pline[N.z - k].z = !((-KOe[ker_index].y * FM_h.x) + (KOe[ker_index].z * FM_h.y));
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
					kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Im, fft_line2[k].y.Re, fft_line2[k].z.Re);
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	OeFunc::CalcOerstedTensors(DOe, n, N, h);

	tensor_to_kernel(DOe, KOe);

	//Done
	return error;
}

#endif