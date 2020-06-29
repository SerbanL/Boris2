#include "stdafx.h"
#include "DemagKernelCollection.h"

#ifdef MODULE_COMPILATION_SDEMAG

//-------------------------- RUN-TIME KERNEL MULTIPLICATION

//SELF VERSIONS

//These compute the self contributions using the real kernels with full use of symmetries
//These set the output, not add into it, so always call it first -> each kernel collection has exactly one self contribution
void DemagKernelCollection::KernelMultiplication_2D_Self(VEC<ReIm3>& In, VEC<ReIm3>& Out)
{
	//Full multiplication with use of kernel symmetries -> re-arranged for better cache use compared to the line versions

	VEC<DBL3>& Kdiag = kernels[self_contribution_index]->Kdiag_real;
	vector<double>& K2D_odiag = kernels[self_contribution_index]->K2D_odiag;

	//zero-th line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		ReIm3 FM = In[i];

		Out[i].x = (Kdiag[i].x  * FM.x) + (K2D_odiag[i] * FM.y);
		Out[i].y = (K2D_odiag[i] * FM.x) + (Kdiag[i].y  * FM.y);
		Out[i].z = (Kdiag[i].z  * FM.z);
	}

#pragma omp parallel for
	for (int j = 1; j < N.y / 2; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx_l = i + j * (N.x / 2 + 1);
			int idx_h = i + (N.y - j) * (N.x / 2 + 1);

			ReIm3 FM_l = In[idx_l];
			ReIm3 FM_h = In[idx_h];

			Out[idx_l].x = (Kdiag[idx_l].x  * FM_l.x) + (K2D_odiag[idx_l] * FM_l.y);
			Out[idx_l].y = (K2D_odiag[idx_l] * FM_l.x) + (Kdiag[idx_l].y  * FM_l.y);
			Out[idx_l].z = (Kdiag[idx_l].z  * FM_l.z);

			Out[idx_h].x = (Kdiag[idx_l].x  * FM_h.x) + (-K2D_odiag[idx_l] * FM_h.y);
			Out[idx_h].y = (-K2D_odiag[idx_l] * FM_h.x) + (Kdiag[idx_l].y  * FM_h.y);
			Out[idx_h].z = (Kdiag[idx_l].z  * FM_h.z);
		}
	}

	//half-way line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int idx = i + (N.y / 2) * (N.x / 2 + 1);

		ReIm3 FM = In[idx];

		Out[idx].x = (Kdiag[idx].x  * FM.x) + (K2D_odiag[idx] * FM.y);
		Out[idx].y = (K2D_odiag[idx] * FM.x) + (Kdiag[idx].y  * FM.y);
		Out[idx].z = (Kdiag[idx].z  * FM.z);
	}
}

void DemagKernelCollection::KernelMultiplication_3D_Self(VEC<ReIm3>& In, VEC<ReIm3>& Out)
{
	VEC<DBL3>& Kdiag = kernels[self_contribution_index]->Kdiag_real;
	VEC<DBL3>& Kodiag = kernels[self_contribution_index]->Kodiag_real;

	//Full multiplication with use of kernel symmetries -> re-arranged for better cache use compared to the line versions

	for (int k = 0; k <= N.z / 2; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x = (Kdiag[ker_index].x * FM.x) + (Kodiag[ker_index].x * FM.y) + (Kodiag[ker_index].y * FM.z);
			Out[idx].y = (Kodiag[ker_index].x * FM.x) + (Kdiag[ker_index].y * FM.y) + (Kodiag[ker_index].z * FM.z);
			Out[idx].z = (Kodiag[ker_index].y * FM.x) + (Kodiag[ker_index].z * FM.y) + (Kdiag[ker_index].z * FM.z);
		}

#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				ReIm3 FM_l = In[idx_l];
				ReIm3 FM_h = In[idx_h];

				int ker_index = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				Out[idx_l].x = (Kdiag[ker_index].x * FM_l.x) + (Kodiag[ker_index].x * FM_l.y) + (Kodiag[ker_index].y * FM_l.z);
				Out[idx_l].y = (Kodiag[ker_index].x * FM_l.x) + (Kdiag[ker_index].y * FM_l.y) + (Kodiag[ker_index].z * FM_l.z);
				Out[idx_l].z = (Kodiag[ker_index].y * FM_l.x) + (Kodiag[ker_index].z * FM_l.y) + (Kdiag[ker_index].z * FM_l.z);

				Out[idx_h].x = (Kdiag[ker_index].x * FM_h.x) + (-Kodiag[ker_index].x * FM_h.y) + (Kodiag[ker_index].y * FM_h.z);
				Out[idx_h].y = (-Kodiag[ker_index].x * FM_h.x) + (Kdiag[ker_index].y * FM_h.y) + (-Kodiag[ker_index].z * FM_h.z);
				Out[idx_h].z = (Kodiag[ker_index].y * FM_h.x) + (-Kodiag[ker_index].z * FM_h.y) + (Kdiag[ker_index].z * FM_h.z);
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x = (Kdiag[ker_index].x * FM.x) + (Kodiag[ker_index].x * FM.y) + (Kodiag[ker_index].y * FM.z);
			Out[idx].y = (Kodiag[ker_index].x * FM.x) + (Kdiag[ker_index].y * FM.y) + (Kodiag[ker_index].z * FM.z);
			Out[idx].z = (Kodiag[ker_index].y * FM.x) + (Kodiag[ker_index].z * FM.y) + (Kdiag[ker_index].z * FM.z);
		}
	}


	for (int k = N.z / 2 + 1; k < N.z; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x = (Kdiag[ker_index].x * FM.x) + (Kodiag[ker_index].x * FM.y) + (-Kodiag[ker_index].y * FM.z);
			Out[idx].y = (Kodiag[ker_index].x * FM.x) + (Kdiag[ker_index].y * FM.y) + (-Kodiag[ker_index].z * FM.z);
			Out[idx].z = (-Kodiag[ker_index].y * FM.x) + (-Kodiag[ker_index].z * FM.y) + (Kdiag[ker_index].z * FM.z);
		}

#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				ReIm3 FM_l = In[idx_l];
				ReIm3 FM_h = In[idx_h];

				int ker_index = i + j * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				Out[idx_l].x = (Kdiag[ker_index].x * FM_l.x) + (Kodiag[ker_index].x * FM_l.y) + (-Kodiag[ker_index].y * FM_l.z);
				Out[idx_l].y = (Kodiag[ker_index].x * FM_l.x) + (Kdiag[ker_index].y * FM_l.y) + (-Kodiag[ker_index].z * FM_l.z);
				Out[idx_l].z = (-Kodiag[ker_index].y * FM_l.x) + (-Kodiag[ker_index].z * FM_l.y) + (Kdiag[ker_index].z * FM_l.z);

				Out[idx_h].x = (Kdiag[ker_index].x * FM_h.x) + (-Kodiag[ker_index].x * FM_h.y) + (-Kodiag[ker_index].y * FM_h.z);
				Out[idx_h].y = (-Kodiag[ker_index].x * FM_h.x) + (Kdiag[ker_index].y * FM_h.y) + (Kodiag[ker_index].z * FM_h.z);
				Out[idx_h].z = (-Kodiag[ker_index].y * FM_h.x) + (Kodiag[ker_index].z * FM_h.y) + (Kdiag[ker_index].z * FM_h.z);
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x = (Kdiag[ker_index].x * FM.x) + (Kodiag[ker_index].x * FM.y) + (-Kodiag[ker_index].y * FM.z);
			Out[idx].y = (Kodiag[ker_index].x * FM.x) + (Kdiag[ker_index].y * FM.y) + (-Kodiag[ker_index].z * FM.z);
			Out[idx].z = (-Kodiag[ker_index].y * FM.x) + (-Kodiag[ker_index].z * FM.y) + (Kdiag[ker_index].z * FM.z);
		}
	}
}

//MULTIPLICATION FOR Z-SHIFTED DEMAG

void DemagKernelCollection::KernelMultiplication_2D_zShifted(VEC<ReIm3>& In, VEC<ReIm3>& Out, VEC<DBL3>& Kdiag, VEC<DBL3>& Kodiag)
{
	//Full multiplication with use of kernel symmetries, z-shifted version

	//zero-th line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		ReIm3 FM = In[i];

		Out[i].x += (Kdiag[i].x * FM.x) + (Kodiag[i].x * FM.y) + !(Kodiag[i].y * FM.z);
		Out[i].y += (Kodiag[i].x * FM.x) + (Kdiag[i].y * FM.y) + !(Kodiag[i].z * FM.z);
		Out[i].z += !(Kodiag[i].y * FM.x) + !(Kodiag[i].z * FM.y) + (Kdiag[i].z * FM.z);
	}

#pragma omp parallel for
	for (int j = 1; j < N.y / 2; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx_l = i + j * (N.x / 2 + 1);
			int idx_h = i + (N.y - j) * (N.x / 2 + 1);

			ReIm3 FM_l = In[idx_l];
			ReIm3 FM_h = In[idx_h];

			Out[idx_l].x += (Kdiag[idx_l].x * FM_l.x) + (Kodiag[idx_l].x * FM_l.y) + !(Kodiag[idx_l].y * FM_l.z);
			Out[idx_l].y += (Kodiag[idx_l].x * FM_l.x) + (Kdiag[idx_l].y * FM_l.y) + !(Kodiag[idx_l].z * FM_l.z);
			Out[idx_l].z += !(Kodiag[idx_l].y * FM_l.x) + !(Kodiag[idx_l].z * FM_l.y) + (Kdiag[idx_l].z * FM_l.z);

			Out[idx_h].x += (Kdiag[idx_l].x * FM_h.x) + (-Kodiag[idx_l].x * FM_h.y) + !(Kodiag[idx_l].y * FM_h.z);
			Out[idx_h].y += (-Kodiag[idx_l].x * FM_h.x) + (Kdiag[idx_l].y * FM_h.y) + !(-Kodiag[idx_l].z * FM_h.z);
			Out[idx_h].z += !(Kodiag[idx_l].y * FM_h.x) + !(-Kodiag[idx_l].z * FM_h.y) + (Kdiag[idx_l].z * FM_h.z);
		}
	}

	//half-way line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int idx = i + (N.y / 2) * (N.x / 2 + 1);

		ReIm3 FM = In[idx];

		Out[idx].x += (Kdiag[idx].x * FM.x) + (Kodiag[idx].x * FM.y) + !(Kodiag[idx].y * FM.z);
		Out[idx].y += (Kodiag[idx].x * FM.x) + (Kdiag[idx].y * FM.y) + !(Kodiag[idx].z * FM.z);
		Out[idx].z += !(Kodiag[idx].y * FM.x) + !(Kodiag[idx].z * FM.y) + (Kdiag[idx].z * FM.z);
	}
}

//z shifted but with kernel calculated for the other direction shift, so adjust multiplications
void DemagKernelCollection::KernelMultiplication_2D_inversezShifted(VEC<ReIm3>& In, VEC<ReIm3>& Out, VEC<DBL3>& Kdiag, VEC<DBL3>& Kodiag)
{
	//Full multiplication with use of kernel symmetries, inverse z-shifted version
	
	//zero-th line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		ReIm3 FM = In[i];

		Out[i].x += (Kdiag[i].x * FM.x) + (Kodiag[i].x * FM.y) + !(-Kodiag[i].y * FM.z);
		Out[i].y += (Kodiag[i].x * FM.x) + (Kdiag[i].y * FM.y) + !(-Kodiag[i].z * FM.z);
		Out[i].z += !(-Kodiag[i].y * FM.x) + !(-Kodiag[i].z * FM.y) + (Kdiag[i].z * FM.z);
	}

#pragma omp parallel for
	for (int j = 1; j < N.y / 2; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx_l = i + j * (N.x / 2 + 1);
			int idx_h = i + (N.y - j) * (N.x / 2 + 1);

			ReIm3 FM_l = In[idx_l];
			ReIm3 FM_h = In[idx_h];

			Out[idx_l].x += (Kdiag[idx_l].x * FM_l.x) + (Kodiag[idx_l].x * FM_l.y) + !(-Kodiag[idx_l].y * FM_l.z);
			Out[idx_l].y += (Kodiag[idx_l].x * FM_l.x) + (Kdiag[idx_l].y * FM_l.y) + !(-Kodiag[idx_l].z * FM_l.z);
			Out[idx_l].z += !(-Kodiag[idx_l].y * FM_l.x) + !(-Kodiag[idx_l].z * FM_l.y) + (Kdiag[idx_l].z * FM_l.z);

			Out[idx_h].x += (Kdiag[idx_l].x * FM_h.x) + (-Kodiag[idx_l].x * FM_h.y) + !(-Kodiag[idx_l].y * FM_h.z);
			Out[idx_h].y += (-Kodiag[idx_l].x * FM_h.x) + (Kdiag[idx_l].y * FM_h.y) + !(Kodiag[idx_l].z * FM_h.z);
			Out[idx_h].z += !(-Kodiag[idx_l].y * FM_h.x) + !(Kodiag[idx_l].z * FM_h.y) + (Kdiag[idx_l].z * FM_h.z);
		}
	}

	//half-way line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int idx = i + (N.y / 2) * (N.x / 2 + 1);

		ReIm3 FM = In[idx];

		Out[idx].x += (Kdiag[idx].x * FM.x) + (Kodiag[idx].x * FM.y) + !(-Kodiag[idx].y * FM.z);
		Out[idx].y += (Kodiag[idx].x * FM.x) + (Kdiag[idx].y * FM.y) + !(-Kodiag[idx].z * FM.z);
		Out[idx].z += !(-Kodiag[idx].y * FM.x) + !(-Kodiag[idx].z * FM.y) + (Kdiag[idx].z * FM.z);
	}
}

//z shifted for 3D : complex kernels, but use kernel symmetries
void DemagKernelCollection::KernelMultiplication_3D_zShifted(VEC<ReIm3>& In, VEC<ReIm3>& Out, VEC<ReIm3>& Kdiag, VEC<ReIm3>& Kodiag)
{
	//z shifted for 3D : can use kernels of reduced dimensions but must be complex
	//
	//Kxx : y - symmetrical (+), z - Re part symmetrical (+), Im part inv. symmetric (-)
	//Kyy : y - symmetrical (+), z - Re part symmetrical (+), Im part inv. symmetric (-)
	//Kzz : y - symmetrical (+), z - Re part symmetrical (+), Im part inv. symmetric (-)
	//
	//Kxy : y - inv. symmetric (-), z - Re part symmetrical  (+), Im part inv. symmetric (-)
	//Kxz : y - symmetrical  (+), z - Re part inv. symmetric (-), Im part symmetrical  (+)
	//Kyz : y - inv. symmetric (-), z - Re part inv. symmetric (-), Im part symmetrical  (+)

	for (int k = 0; k <= N.z / 2; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x += (Kdiag[ker_index].x * FM.x) + (Kodiag[ker_index].x * FM.y) + (Kodiag[ker_index].y * FM.z);
			Out[idx].y += (Kodiag[ker_index].x * FM.x) + (Kdiag[ker_index].y * FM.y) + (Kodiag[ker_index].z * FM.z);
			Out[idx].z += (Kodiag[ker_index].y * FM.x) + (Kodiag[ker_index].z * FM.y) + (Kdiag[ker_index].z * FM.z);
		}

		//in between 0 and middle
#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				ReIm3 FM_l = In[idx_l];
				ReIm3 FM_h = In[idx_h];

				int ker_index = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				//lower z, lower y
				Out[idx_l].x += (Kdiag[ker_index].x * FM_l.x) + (Kodiag[ker_index].x * FM_l.y) + (Kodiag[ker_index].y * FM_l.z);
				Out[idx_l].y += (Kodiag[ker_index].x * FM_l.x) + (Kdiag[ker_index].y * FM_l.y) + (Kodiag[ker_index].z * FM_l.z);
				Out[idx_l].z += (Kodiag[ker_index].y * FM_l.x) + (Kodiag[ker_index].z * FM_l.y) + (Kdiag[ker_index].z * FM_l.z);

				//lower z, upper y
				Out[idx_h].x += (Kdiag[ker_index].x * FM_h.x) - (Kodiag[ker_index].x * FM_h.y) + (Kodiag[ker_index].y * FM_h.z);
				Out[idx_h].y += -1.0 * (Kodiag[ker_index].x * FM_h.x) + (Kdiag[ker_index].y * FM_h.y) - (Kodiag[ker_index].z * FM_h.z);
				Out[idx_h].z += (Kodiag[ker_index].y * FM_h.x) - (Kodiag[ker_index].z * FM_h.y) + (Kdiag[ker_index].z * FM_h.z);
			}
		}

		//mid line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x += (Kdiag[ker_index].x * FM.x) + (Kodiag[ker_index].x * FM.y) + (Kodiag[ker_index].y * FM.z);
			Out[idx].y += (Kodiag[ker_index].x * FM.x) + (Kdiag[ker_index].y * FM.y) + (Kodiag[ker_index].z * FM.z);
			Out[idx].z += (Kodiag[ker_index].y * FM.x) + (Kodiag[ker_index].z * FM.y) + (Kdiag[ker_index].z * FM.z);
		}
	}

	//upper z
	for (int k = N.z / 2 + 1; k < N.z; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			//upper z, lower y
			Out[idx].x += ((~Kdiag[ker_index].x) * FM.x) + ((~Kodiag[ker_index].x) * FM.y) - ((~Kodiag[ker_index].y) * FM.z);
			Out[idx].y += ((~Kodiag[ker_index].x) * FM.x) + ((~Kdiag[ker_index].y) * FM.y) - ((~Kodiag[ker_index].z) * FM.z);
			Out[idx].z += -1.0 * ((~Kodiag[ker_index].y) * FM.x) - ((~Kodiag[ker_index].z) * FM.y) + ((~Kdiag[ker_index].z) * FM.z);
		}

#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				ReIm3 FM_l = In[idx_l];
				ReIm3 FM_h = In[idx_h];

				int ker_index = i + j * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				Out[idx_l].x += ((~Kdiag[ker_index].x) * FM_l.x) + ((~Kodiag[ker_index].x) * FM_l.y) - ((~Kodiag[ker_index].y) * FM_l.z);
				Out[idx_l].y += ((~Kodiag[ker_index].x) * FM_l.x) + ((~Kdiag[ker_index].y) * FM_l.y) - ((~Kodiag[ker_index].z) * FM_l.z);
				Out[idx_l].z += -1.0 * ((~Kodiag[ker_index].y) * FM_l.x) - ((~Kodiag[ker_index].z) * FM_l.y) + ((~Kdiag[ker_index].z) * FM_l.z);

				//upper z, upper y
				Out[idx_h].x += ((~Kdiag[ker_index].x) * FM_h.x) - ((~Kodiag[ker_index].x) * FM_h.y) - ((~Kodiag[ker_index].y) * FM_h.z);
				Out[idx_h].y += -1.0 * ((~Kodiag[ker_index].x) * FM_h.x) + ((~Kdiag[ker_index].y) * FM_h.y) + ((~Kodiag[ker_index].z) * FM_h.z);
				Out[idx_h].z += -1.0 * ((~Kodiag[ker_index].y) * FM_h.x) + ((~Kodiag[ker_index].z) * FM_h.y) + ((~Kdiag[ker_index].z) * FM_h.z);
			}
		}

		//mid line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

			ReIm3 FM = In[idx];

			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			Out[idx].x += ((~Kdiag[ker_index].x) * FM.x) + ((~Kodiag[ker_index].x) * FM.y) - ((~Kodiag[ker_index].y) * FM.z);
			Out[idx].y += ((~Kodiag[ker_index].x) * FM.x) + ((~Kdiag[ker_index].y) * FM.y) - ((~Kodiag[ker_index].z) * FM.z);
			Out[idx].z += -1.0 * ((~Kodiag[ker_index].y) * FM.x) - ((~Kodiag[ker_index].z) * FM.y) + ((~Kdiag[ker_index].z) * FM.z);
		}
	}
}

//MULTIPLE INPUTS

//multiple input spaces version, used when not embedding the kenel multiplication
void DemagKernelCollection::KernelMultiplication_2D(std::vector<VEC<ReIm3>*>& Incol, VEC<ReIm3>& Out)
{
	//first compute the self contribution -> this sets Out
	KernelMultiplication_2D_Self(*Incol[self_contribution_index], Out);

	for (int mesh_index = 0; mesh_index < Incol.size(); mesh_index++) {

		//z-shifted : use symmetries
		if (kernels[mesh_index]->zshifted) {

			//inverse : adjust signs
			if (inverse_shifted[mesh_index]) {

				KernelMultiplication_2D_inversezShifted(*Incol[mesh_index], Out, kernels[mesh_index]->Kdiag_real, kernels[mesh_index]->Kodiag_real);
			}
			//regular : as is
			else {

				KernelMultiplication_2D_zShifted(*Incol[mesh_index], Out, kernels[mesh_index]->Kdiag_real, kernels[mesh_index]->Kodiag_real);
			}
		}

		//now compute the other contributions by adding to Out : general kernel multiplication without any symmetries used
		else if (!kernels[mesh_index]->internal_demag) {

#pragma omp parallel for
			for (int index = 0; index < (N.x / 2 + 1)*N.y; index++) {

				ReIm3 FM = (*Incol[mesh_index])[index];

				Out[index].x += (kernels[mesh_index]->Kdiag_cmpl[index].x * FM.x) + (kernels[mesh_index]->Kodiag_cmpl[index].x * FM.y) + (kernels[mesh_index]->Kodiag_cmpl[index].y * FM.z);
				Out[index].y += (kernels[mesh_index]->Kodiag_cmpl[index].x * FM.x) + (kernels[mesh_index]->Kdiag_cmpl[index].y * FM.y) + (kernels[mesh_index]->Kodiag_cmpl[index].z * FM.z);
				Out[index].z += (kernels[mesh_index]->Kodiag_cmpl[index].y * FM.x) + (kernels[mesh_index]->Kodiag_cmpl[index].z * FM.y) + (kernels[mesh_index]->Kdiag_cmpl[index].z * FM.z);
			}
		}
	}
}

void DemagKernelCollection::KernelMultiplication_3D(std::vector<VEC<ReIm3>*>& Incol, VEC<ReIm3>& Out)
{
	//first compute the self contribution -> this sets Out
	KernelMultiplication_3D_Self(*Incol[self_contribution_index], Out);

	//now compute the other contribution by adding to Out
	for (int mesh_index = 0; mesh_index < Incol.size(); mesh_index++) {

		//z-shifted : use symmetries
		if (kernels[mesh_index]->zshifted) {

			KernelMultiplication_3D_zShifted(*Incol[mesh_index], Out, kernels[mesh_index]->Kdiag_cmpl, kernels[mesh_index]->Kodiag_cmpl);
		}

		//now compute the other contributions by adding to Out : general kernel multiplication without any symmetries used
		else if (!kernels[mesh_index]->internal_demag) {

#pragma omp parallel for
			for (int index = 0; index < (N.x / 2 + 1)*N.y*N.z; index++) {

				ReIm3 FM = (*Incol[mesh_index])[index];

				Out[index].x += (kernels[mesh_index]->Kdiag_cmpl[index].x * FM.x) + (kernels[mesh_index]->Kodiag_cmpl[index].x * FM.y) + (kernels[mesh_index]->Kodiag_cmpl[index].y * FM.z);
				Out[index].y += (kernels[mesh_index]->Kodiag_cmpl[index].x * FM.x) + (kernels[mesh_index]->Kdiag_cmpl[index].y * FM.y) + (kernels[mesh_index]->Kodiag_cmpl[index].z * FM.z);
				Out[index].z += (kernels[mesh_index]->Kodiag_cmpl[index].y * FM.x) + (kernels[mesh_index]->Kodiag_cmpl[index].z * FM.y) + (kernels[mesh_index]->Kdiag_cmpl[index].z * FM.z);
			}
		}
	}
}

#endif