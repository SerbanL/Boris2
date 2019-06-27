#include "stdafx.h"
#include "SDemag_KernelCollection.h"

#ifdef MODULE_SDEMAG

//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE

void KerTypeCollection::KernelMultiplication_2D_Self(void)
{
	//Full multiplication with use of kernel symmetries -> re-arranged for better cache use compared to the line versions

	VEC<DBL3>& Kdiag = kernel->Kdiag_real;
	vector<double>& K2D_odiag = kernel->K2D_odiag;

	SZ3& N = kernel->N;
	int size = In.size();

	//zero-th line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		DBL3 Kdiag_val = Kdiag[i];
		double Kodiag_val = K2D_odiag[i];

		for (int vidx = 0; vidx < size; vidx++) {

			ReIm3 FM = (*In[vidx])[i];

			(*Out[vidx])[i].x = (Kdiag_val.x  * FM.x) + (Kodiag_val * FM.y);
			(*Out[vidx])[i].y = (Kodiag_val * FM.x) + (Kdiag_val.y  * FM.y);
			(*Out[vidx])[i].z = (Kdiag_val.z  * FM.z);
		}
	}

#pragma omp parallel for
	for (int j = 1; j < N.y / 2; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx_l = i + j * (N.x / 2 + 1);
			int idx_h = i + (N.y - j) * (N.x / 2 + 1);

			DBL3 Kdiag_val = Kdiag[idx_l];
			double Kodiag_val = K2D_odiag[idx_l];

			//putting this loop on the outside actually makes the code a little bit faster in testing (but still alot slower than the multiple inputs version)
			//a bit surprising since the point of having this loop here is to minimise kernel reading operations, but reality differs. It must be In reading and Out writing trump kernel reading operations.
			//Same goes for all other places like this in this file
			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM_l = (*In[vidx])[idx_l];
				ReIm3 FM_h = (*In[vidx])[idx_h];

				(*Out[vidx])[idx_l].x = (Kdiag_val.x  * FM_l.x) + (Kodiag_val * FM_l.y);
				(*Out[vidx])[idx_l].y = (Kodiag_val * FM_l.x) + (Kdiag_val.y  * FM_l.y);
				(*Out[vidx])[idx_l].z = (Kdiag_val.z  * FM_l.z);

				(*Out[vidx])[idx_h].x = (Kdiag_val.x  * FM_h.x) + (-Kodiag_val * FM_h.y);
				(*Out[vidx])[idx_h].y = (-Kodiag_val * FM_h.x) + (Kdiag_val.y  * FM_h.y);
				(*Out[vidx])[idx_h].z = (Kdiag_val.z  * FM_h.z);
			}
		}
	}

	//half-way line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int idx = i + (N.y / 2) * (N.x / 2 + 1);

		DBL3 Kdiag_val = Kdiag[idx];
		double Kodiag_val = K2D_odiag[idx];

		for (int vidx = 0; vidx < size; vidx++) {

			ReIm3 FM = (*In[vidx])[idx];

			(*Out[vidx])[idx].x = (Kdiag_val.x  * FM.x) + (Kodiag_val * FM.y);
			(*Out[vidx])[idx].y = (Kodiag_val * FM.x) + (Kdiag_val.y  * FM.y);
			(*Out[vidx])[idx].z = (Kdiag_val.z  * FM.z);
		}
	}
}

void KerTypeCollection::KernelMultiplication_3D_Self(void)
{
	VEC<DBL3>& Kdiag = kernel->Kdiag_real;
	VEC<DBL3>& Kodiag = kernel->Kodiag_real;

	SZ3& N = kernel->N;
	int size = In.size();

	//Full multiplication with use of kernel symmetries -> re-arranged for better cache use compared to the line versions

	for (int k = 0; k <= N.z / 2; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			int ker_index = i + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			DBL3 Kdiag_val = Kdiag[ker_index];
			DBL3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x = (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y = (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z = (Kodiag_val.y * FM.x) + (Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}

#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				int ker_index = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				DBL3 Kdiag_val = Kdiag[ker_index];
				DBL3 Kodiag_val = Kodiag[ker_index];

				for (int vidx = 0; vidx < size; vidx++) {

					ReIm3 FM_l = (*In[vidx])[idx_l];
					ReIm3 FM_h = (*In[vidx])[idx_h];

					(*Out[vidx])[idx_l].x = (Kdiag_val.x * FM_l.x) + (Kodiag_val.x * FM_l.y) + (Kodiag_val.y * FM_l.z);
					(*Out[vidx])[idx_l].y = (Kodiag_val.x * FM_l.x) + (Kdiag_val.y * FM_l.y) + (Kodiag_val.z * FM_l.z);
					(*Out[vidx])[idx_l].z = (Kodiag_val.y * FM_l.x) + (Kodiag_val.z * FM_l.y) + (Kdiag_val.z * FM_l.z);

					(*Out[vidx])[idx_h].x = (Kdiag_val.x * FM_h.x) + (-Kodiag_val.x * FM_h.y) + (Kodiag_val.y * FM_h.z);
					(*Out[vidx])[idx_h].y = (-Kodiag_val.x * FM_h.x) + (Kdiag_val.y * FM_h.y) + (-Kodiag_val.z * FM_h.z);
					(*Out[vidx])[idx_h].z = (Kodiag_val.y * FM_h.x) + (-Kodiag_val.z * FM_h.y) + (Kdiag_val.z * FM_h.z);
				}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
			
			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			DBL3 Kdiag_val = Kdiag[ker_index];
			DBL3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x = (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y = (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z = (Kodiag_val.y * FM.x) + (Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}
	}


	for (int k = N.z / 2 + 1; k < N.z; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			int ker_index = i + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			DBL3 Kdiag_val = Kdiag[ker_index];
			DBL3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x = (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (-Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y = (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (-Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z = (-Kodiag_val.y * FM.x) + (-Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}

#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				int ker_index = i + j * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				DBL3 Kdiag_val = Kdiag[ker_index];
				DBL3 Kodiag_val = Kodiag[ker_index];

				for (int vidx = 0; vidx < size; vidx++) {

					ReIm3 FM_l = (*In[vidx])[idx_l];
					ReIm3 FM_h = (*In[vidx])[idx_h];

					(*Out[vidx])[idx_l].x = (Kdiag_val.x * FM_l.x) + (Kodiag_val.x * FM_l.y) + (-Kodiag_val.y * FM_l.z);
					(*Out[vidx])[idx_l].y = (Kodiag_val.x * FM_l.x) + (Kdiag_val.y * FM_l.y) + (-Kodiag_val.z * FM_l.z);
					(*Out[vidx])[idx_l].z = (-Kodiag_val.y * FM_l.x) + (-Kodiag_val.z * FM_l.y) + (Kdiag_val.z * FM_l.z);

					(*Out[vidx])[idx_h].x = (Kdiag_val.x * FM_h.x) + (-Kodiag_val.x * FM_h.y) + (-Kodiag_val.y * FM_h.z);
					(*Out[vidx])[idx_h].y = (-Kodiag_val.x * FM_h.x) + (Kdiag_val.y * FM_h.y) + (Kodiag_val.z * FM_h.z);
					(*Out[vidx])[idx_h].z = (-Kodiag_val.y * FM_h.x) + (Kodiag_val.z * FM_h.y) + (Kdiag_val.z * FM_h.z);
				}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
			
			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			DBL3 Kdiag_val = Kdiag[ker_index];
			DBL3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x = (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (-Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y = (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (-Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z = (-Kodiag_val.y * FM.x) + (-Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}
	}
}

void KerTypeCollection::KernelMultiplication_2D_zShifted(void)
{
	//Full multiplication with use of kernel symmetries, z-shifted version

	VEC<DBL3>& Kdiag = kernel->Kdiag_real;
	VEC<DBL3>& Kodiag = kernel->Kodiag_real;

	SZ3& N = kernel->N;
	int size = In.size();

	//zero-th line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		DBL3 Kdiag_val = Kdiag[i];
		DBL3 Kodiag_val = Kodiag[i];

		for (int vidx = 0; vidx < size; vidx++) {

			ReIm3 FM = (*In[vidx])[i];

			//you can improve this by sorting the In and Out vector by the inverse_shifted vector.
			//thus you wouldn't need this check here, but can separate the vidx loop into two - will improve speed slightly but not enough to get anywhere near the multiple inputs method so no point putting more work here
			if (inverse_shifted[vidx]) {

				(*Out[vidx])[i].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + !(-Kodiag_val.y * FM.z);
				(*Out[vidx])[i].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + !(-Kodiag_val.z * FM.z);
				(*Out[vidx])[i].z += !(-Kodiag_val.y * FM.x) + !(-Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
			else {

				(*Out[vidx])[i].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + !(Kodiag_val.y * FM.z);
				(*Out[vidx])[i].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + !(Kodiag_val.z * FM.z);
				(*Out[vidx])[i].z += !(Kodiag_val.y * FM.x) + !(Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}
	}

#pragma omp parallel for
	for (int j = 1; j < N.y / 2; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx_l = i + j * (N.x / 2 + 1);
			int idx_h = i + (N.y - j) * (N.x / 2 + 1);

			DBL3 Kdiag_val = Kdiag[idx_l];
			DBL3 Kodiag_val = Kodiag[idx_l];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM_l = (*In[vidx])[idx_l];
				ReIm3 FM_h = (*In[vidx])[idx_h];

				if (inverse_shifted[vidx]) {

					(*Out[vidx])[idx_l].x += (Kdiag_val.x * FM_l.x) + (Kodiag_val.x * FM_l.y) + !(-Kodiag_val.y * FM_l.z);
					(*Out[vidx])[idx_l].y += (Kodiag_val.x * FM_l.x) + (Kdiag_val.y * FM_l.y) + !(-Kodiag_val.z * FM_l.z);
					(*Out[vidx])[idx_l].z += !(-Kodiag_val.y * FM_l.x) + !(-Kodiag_val.z * FM_l.y) + (Kdiag_val.z * FM_l.z);

					(*Out[vidx])[idx_h].x += (Kdiag_val.x * FM_h.x) + (-Kodiag_val.x * FM_h.y) + !(-Kodiag_val.y * FM_h.z);
					(*Out[vidx])[idx_h].y += (-Kodiag_val.x * FM_h.x) + (Kdiag_val.y * FM_h.y) + !(Kodiag_val.z * FM_h.z);
					(*Out[vidx])[idx_h].z += !(-Kodiag_val.y * FM_h.x) + !(Kodiag_val.z * FM_h.y) + (Kdiag_val.z * FM_h.z);
				}
				else {

					(*Out[vidx])[idx_l].x += (Kdiag_val.x * FM_l.x) + (Kodiag_val.x * FM_l.y) + !(Kodiag_val.y * FM_l.z);
					(*Out[vidx])[idx_l].y += (Kodiag_val.x * FM_l.x) + (Kdiag_val.y * FM_l.y) + !(Kodiag_val.z * FM_l.z);
					(*Out[vidx])[idx_l].z += !(Kodiag_val.y * FM_l.x) + !(Kodiag_val.z * FM_l.y) + (Kdiag_val.z * FM_l.z);

					(*Out[vidx])[idx_h].x += (Kdiag_val.x * FM_h.x) + (-Kodiag_val.x * FM_h.y) + !(Kodiag_val.y * FM_h.z);
					(*Out[vidx])[idx_h].y += (-Kodiag_val.x * FM_h.x) + (Kdiag_val.y * FM_h.y) + !(-Kodiag_val.z * FM_h.z);
					(*Out[vidx])[idx_h].z += !(Kodiag_val.y * FM_h.x) + !(-Kodiag_val.z * FM_h.y) + (Kdiag_val.z * FM_h.z);
				}
			}
		}
	}

	//half-way line
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int idx = i + (N.y / 2) * (N.x / 2 + 1);

		DBL3 Kdiag_val = Kdiag[idx];
		DBL3 Kodiag_val = Kodiag[idx];

		for (int vidx = 0; vidx < size; vidx++) {

			ReIm3 FM = (*In[vidx])[idx];

			if (inverse_shifted[vidx]) {

				(*Out[vidx])[idx].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + !(-Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + !(-Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z += !(-Kodiag_val.y * FM.x) + !(-Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
			else {

				(*Out[vidx])[idx].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + !(Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + !(Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z += !(Kodiag_val.y * FM.x) + !(Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}
	}
}

void KerTypeCollection::KernelMultiplication_3D_zShifted(void)
{
	//z shifted for 3D : can use kernels of reduced dimensions but must be complex
	//
	//Kxx : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
	//Kyy : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
	//Kzz : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
	//
	//Kxy : y - asymmetrical (-), z - Re part symmetrical  (+), Im part asymmetrical (-)
	//Kxz : y - symmetrical  (+), z - Re part asymmetrical (-), Im part symmetrical  (+)
	//Kyz : y - asymmetrical (-), z - Re part asymmetrical (-), Im part symmetrical  (+)

	VEC<ReIm3>& Kdiag = kernel->Kdiag_cmpl;
	VEC<ReIm3>& Kodiag = kernel->Kodiag_cmpl;

	SZ3& N = kernel->N;
	int size = In.size();

	for (int k = 0; k <= N.z / 2; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			int ker_index = i + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			ReIm3 Kdiag_val = Kdiag[ker_index];
			ReIm3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z += (Kodiag_val.y * FM.x) + (Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}

		//in between 0 and middle
#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				int ker_index = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

				ReIm3 Kdiag_val = Kdiag[ker_index];
				ReIm3 Kodiag_val = Kodiag[ker_index];

				for (int vidx = 0; vidx < size; vidx++) {

					ReIm3 FM_l = (*In[vidx])[idx_l];
					ReIm3 FM_h = (*In[vidx])[idx_h];

					//lower z, lower y
					(*Out[vidx])[idx_l].x += (Kdiag_val.x * FM_l.x) + (Kodiag_val.x * FM_l.y) + (Kodiag_val.y * FM_l.z);
					(*Out[vidx])[idx_l].y += (Kodiag_val.x * FM_l.x) + (Kdiag_val.y * FM_l.y) + (Kodiag_val.z * FM_l.z);
					(*Out[vidx])[idx_l].z += (Kodiag_val.y * FM_l.x) + (Kodiag_val.z * FM_l.y) + (Kdiag_val.z * FM_l.z);

					//lower z, upper y
					(*Out[vidx])[idx_h].x += (Kdiag_val.x * FM_h.x) - (Kodiag_val.x * FM_h.y) + (Kodiag_val.y * FM_h.z);
					(*Out[vidx])[idx_h].y += -1.0 * (Kodiag_val.x * FM_h.x) + (Kdiag_val.y * FM_h.y) - (Kodiag_val.z * FM_h.z);
					(*Out[vidx])[idx_h].z += (Kodiag_val.y * FM_h.x) - (Kodiag_val.z * FM_h.y) + (Kdiag_val.z * FM_h.z);
				}
			}
		}

		//mid line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1);

			ReIm3 Kdiag_val = Kdiag[ker_index];
			ReIm3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (Kodiag_val.y * FM.z);
				(*Out[vidx])[idx].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (Kodiag_val.z * FM.z);
				(*Out[vidx])[idx].z += (Kodiag_val.y * FM.x) + (Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
			}
		}
	}

	//upper z
	for (int k = N.z / 2 + 1; k < N.z; k++) {

		//zero-th line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + k * (N.x / 2 + 1) * N.y;

			int ker_index = i + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			ReIm3 Kdiag_val = Kdiag[ker_index];
			ReIm3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];
				
				//upper z, lower y
				(*Out[vidx])[idx].x += ((~Kdiag_val.x) * FM.x) + ((~Kodiag_val.x) * FM.y) - ((~Kodiag_val.y) * FM.z);
				(*Out[vidx])[idx].y += ((~Kodiag_val.x) * FM.x) + ((~Kdiag_val.y) * FM.y) - ((~Kodiag_val.z) * FM.z);
				(*Out[vidx])[idx].z += -1.0 * ((~Kodiag_val.y) * FM.x) - ((~Kodiag_val.z) * FM.y) + ((~Kdiag_val.z) * FM.z);
			}
		}

#pragma omp parallel for
		for (int j = 1; j < N.y / 2; j++) {
			for (int i = 0; i < (N.x / 2 + 1); i++) {

				int idx_l = i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;
				int idx_h = i + (N.y - j) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

				int ker_index = i + j * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

				ReIm3 Kdiag_val = Kdiag[ker_index];
				ReIm3 Kodiag_val = Kodiag[ker_index];

				for (int vidx = 0; vidx < size; vidx++) {

					ReIm3 FM_l = (*In[vidx])[idx_l];
					ReIm3 FM_h = (*In[vidx])[idx_h];

					(*Out[vidx])[idx_l].x += ((~Kdiag_val.x) * FM_l.x) + ((~Kodiag_val.x) * FM_l.y) - ((~Kodiag_val.y) * FM_l.z);
					(*Out[vidx])[idx_l].y += ((~Kodiag_val.x) * FM_l.x) + ((~Kdiag_val.y) * FM_l.y) - ((~Kodiag_val.z) * FM_l.z);
					(*Out[vidx])[idx_l].z += -1.0 * ((~Kodiag_val.y) * FM_l.x) - ((~Kodiag_val.z) * FM_l.y) + ((~Kdiag_val.z) * FM_l.z);

					//upper z, upper y
					(*Out[vidx])[idx_h].x += ((~Kdiag_val.x) * FM_h.x) - ((~Kodiag_val.x) * FM_h.y) - ((~Kodiag_val.y) * FM_h.z);
					(*Out[vidx])[idx_h].y += -1.0 * ((~Kodiag_val.x) * FM_h.x) + ((~Kdiag_val.y) * FM_h.y) + ((~Kodiag_val.z) * FM_h.z);
					(*Out[vidx])[idx_h].z += -1.0 * ((~Kodiag_val.y) * FM_h.x) + ((~Kodiag_val.z) * FM_h.y) + ((~Kdiag_val.z) * FM_h.z);
				}
			}
		}

		//mid line
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int idx = i + (N.y / 2) * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y;

			int ker_index = i + (N.y / 2) * (N.x / 2 + 1) + (N.z - k) * (N.x / 2 + 1) * (N.y / 2 + 1);

			ReIm3 Kdiag_val = Kdiag[ker_index];
			ReIm3 Kodiag_val = Kodiag[ker_index];

			for (int vidx = 0; vidx < size; vidx++) {

				ReIm3 FM = (*In[vidx])[idx];

				(*Out[vidx])[idx].x += ((~Kdiag_val.x) * FM.x) + ((~Kodiag_val.x) * FM.y) - ((~Kodiag_val.y) * FM.z);
				(*Out[vidx])[idx].y += ((~Kodiag_val.x) * FM.x) + ((~Kdiag_val.y) * FM.y) - ((~Kodiag_val.z) * FM.z);
				(*Out[vidx])[idx].z += -1.0 * ((~Kodiag_val.y) * FM.x) - ((~Kodiag_val.z) * FM.y) + ((~Kdiag_val.z) * FM.z);
			}
		}
	}
}

void KerTypeCollection::KernelMultiplication_2D_Regular(void)
{
	VEC<ReIm3>& Kdiag = kernel->Kdiag_cmpl;
	VEC<ReIm3>& Kodiag = kernel->Kodiag_cmpl;

	SZ3& N = kernel->N;
	int size = In.size();

#pragma omp parallel for
	for (int idx = 0; idx < (N.x / 2 + 1)*N.y; idx++) {

		ReIm3 Kdiag_val = Kdiag[idx];
		ReIm3 Kodiag_val = Kodiag[idx];

		for (int vidx = 0; vidx < size; vidx++) {

			ReIm3 FM = (*In[vidx])[idx];

			(*Out[vidx])[idx].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (Kodiag_val.y * FM.z);
			(*Out[vidx])[idx].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (Kodiag_val.z * FM.z);
			(*Out[vidx])[idx].z += (Kodiag_val.y * FM.x) + (Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
		}
	}
}

void KerTypeCollection::KernelMultiplication_3D_Regular(void)
{
	VEC<ReIm3>& Kdiag = kernel->Kdiag_cmpl;
	VEC<ReIm3>& Kodiag = kernel->Kodiag_cmpl;

	SZ3& N = kernel->N;
	int size = In.size();

#pragma omp parallel for
	for (int idx = 0; idx < (N.x / 2 + 1)*N.y*N.z; idx++) {

		ReIm3 Kdiag_val = Kdiag[idx];
		ReIm3 Kodiag_val = Kodiag[idx];

		for (int vidx = 0; vidx < size; vidx++) {

			ReIm3 FM = (*In[vidx])[idx];

			(*Out[vidx])[idx].x += (Kdiag_val.x * FM.x) + (Kodiag_val.x * FM.y) + (Kodiag_val.y * FM.z);
			(*Out[vidx])[idx].y += (Kodiag_val.x * FM.x) + (Kdiag_val.y * FM.y) + (Kodiag_val.z * FM.z);
			(*Out[vidx])[idx].z += (Kodiag_val.y * FM.x) + (Kodiag_val.z * FM.y) + (Kdiag_val.z * FM.z);
		}
	}
}

void KerTypeCollection::Kernel_Multiplication_2D(void)
{
	if (kernel->internal_demag) {

		KernelMultiplication_2D_Self();
	}
	else {

		if (kernel->zshifted) {

			KernelMultiplication_2D_zShifted();
		}
		else {

			KernelMultiplication_2D_Regular();
		}
	}
}

void KerTypeCollection::Kernel_Multiplication_3D(void)
{
	if (kernel->internal_demag) {

		KernelMultiplication_3D_Self();
	}
	else {

		if (kernel->zshifted) {

			KernelMultiplication_3D_zShifted();
		}
		else {

			KernelMultiplication_3D_Regular();
		}
	}
}

#endif