#pragma once

#include "BorisLib.h"
#include "ErrorHandler.h"

#include "ConvolutionData.h"

using namespace std;

//Convolution based on convolution theorem : 
//Out is the convolution of In with some tensor. 
//FFT of the tensor is called a kernel.
//Take FFT on In, multiply pointwise with the kernels, take IFFT of the result to obtain Out.

//The method for doing this is contained in the Convolution class.
//1. Pre-set dimensions of In and Out VECs so memory can be allocated correctly : SZ3 n for number of cells in each dimension, DBL3 h for cellsize (h needed for Kernel calculation since h may not be cubic)

//2. Construct this with the required Kernel class.
//This defines exactly what kernels are used, how they are calculated - Calculate_Demag_Kernels() - and how the pointwise multiplication is done (see e.g. DemagKernel.h for an example)

//3. Before starting a computation call Calculate_Demag_Kernels() - provided by the Kernel class - to pre-compute the Kernels

//4. Call Convolute with In and Out to compute the convolution during run-time

//template <typename Owner, typename Kernel>
template <typename Owner, typename Kernel>
class Convolution :
	public virtual ConvolutionData,		//virtual to solve the diamond problem as it's also inherited by Kernel
	public Kernel
{

private:

	//if the object couldn't be created properly in the constructor an error is set here
	BError convolution_error_on_create;

private:

	//Embedded (default)

	//convolute In with kernels, set output in Out. 2D is for n.z == 1.
	//If clearOut flag is true then Out is set, otherwise Out is added into.
	//SINGLE INPUT, SINGLE OUTPUT
	double Convolute_2D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut);
	double Convolute_3D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut);

	//AVERAGED INPUTS, SINGLE OUTPUT
	//Same as Convolution with (In1 + In2) / 2 as input.
	double Convolute_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut);
	double Convolute_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut);

	//AVERAGED INPUTS, DUPLICATED OUTPUTS
	//Same as Convolution with (In1 + In2) / 2 as input and output copied to both Out1 and Out2.
	double Convolute_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut);
	double Convolute_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut);

	//Not embedded

	//SINGLE INPUT

	void ForwardFFT_2D(VEC<DBL3> &In);
	void ForwardFFT_3D(VEC<DBL3> &In);

	//AVERAGED INPUTS

	void ForwardFFT_2D(VEC<DBL3> &In1, VEC<DBL3> &In2);
	void ForwardFFT_3D(VEC<DBL3> &In1, VEC<DBL3> &In2);

	//SINGLE OUTPUT

	double InverseFFT_2D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut);
	double InverseFFT_3D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut);

	//AVERAGED INPUTS, SINGLE OUTPUT

	double InverseFFT_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut);
	double InverseFFT_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut);

	//AVERAGED INPUTS, DUPLICATED OUTPUTS

	double InverseFFT_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut);
	double InverseFFT_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut);

protected:

	//-------------------------- CONSTRUCTORS

	//empty constructor only used when loading from file : FFT Dimensions will be set when calling RepairObjectState shortly after
	Convolution(void) :
		ConvolutionData(),
		Kernel()
	{}

	Convolution(INT3 n_, DBL3 h_) :
		ConvolutionData(),
		Kernel()
	{
		convolution_error_on_create = SetDimensions(n_, h_);
	}

	//in some cases the Kernel needs a pointer to its owner class so pass this in addition to dimensions
	//template <typename Owner>
	Convolution(Owner* pOwner) :
		ConvolutionData(),
		Kernel(pOwner)
	{
	}

	//in some cases the Kernel needs a pointer to its owner class so pass this in addition to dimensions
	//template <typename Owner>
	Convolution(INT3 n_, DBL3 h_, Owner* pOwner) :
		ConvolutionData(),
		Kernel(pOwner)
	{
		convolution_error_on_create = SetDimensions(n_, h_);
	}

	virtual ~Convolution() {}

	//-------------------------- ERROR CHECKER

	BError Convolution_Error_on_Create(void) { return convolution_error_on_create; }

	//-------------------------- CONFIGURATION

	//This methods sets all values from and h, including allocating memory - call this before initializing kernels or doing any convolutions
	BError SetDimensions(SZ3 n_, DBL3 h_, bool multiplication_embedding_ = true, INT3 pbc_images_ = INT3());

	//-------------------------- CHECK

	//return true only if both n_ and h_ match the current FFT dimensions (n and h); also number of pbc images must match
	bool CheckDimensions(SZ3 n_, DBL3 h_, INT3 pbc_images_) { return (n == n_ && h == h_ && pbc_images == pbc_images_); }

	//-------------------------- RUN-TIME CONVOLUTION

	//Embedded (default)

	//SINGLE INPUT, SINGLE OUTPUT

	//Run a 2D or 3D convolution depending on set n (n.z == 1 for 2D).
	//Return dot product of In with Out
	double Convolute(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut)
	{
		if (n.z == 1) return Convolute_2D(In, Out, clearOut);
		else return Convolute_3D(In, Out, clearOut);
	}

	//AVERAGED INPUTS, SINGLE OUTPUT

	//Run a 2D or 3D convolution depending on set n (n.z == 1 for 2D).
	//Return dot product of In with Out
	//Same as Convolution with (In1 + In2) / 2 as input.
	double Convolute_AveragedInputs(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut)
	{
		if (n.z == 1) return Convolute_2D(In1, In2, Out, clearOut);
		else return Convolute_3D(In1, In2, Out, clearOut);
	}

	//AVERAGED INPUTS, DUPLICATED OUTPUTS

	//Run a 2D or 3D convolution depending on set n (n.z == 1 for 2D).
	//Return dot product of In with Out
	//Same as Convolution with (In1 + In2) / 2 as input and output copied to both Out1 and Out2.
	double Convolute_AveragedInputs_DuplicatedOutputs(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut)
	{
		if (n.z == 1) return Convolute_2D(In1, In2, Out1, Out2, clearOut);
		else return Convolute_3D(In1, In2, Out1, Out2, clearOut);
	}

	//Not embedded (must set multiplication_embedding = false before using these)

	//Rather than using a single call to the Convolute method above, you can also break down the convolution into its 3 main parts:
	//Forward FFT, Kernel Multiplcation, Inverse FFT
	//This must be used only after setting multiplcation embedding off : multiplication_embedding = false
	//This is useful if you need several forward FFTs from different sources before performing a set of kernel multiplications
	
	//SINGLE INPUT, SINGLE OUTPUT

	//1. Forward (In -> F)
	void ForwardFFT(VEC<DBL3> &In)
	{
		if (n.z == 1) ForwardFFT_2D(In);
		else ForwardFFT_3D(In);
	}

	//AVERAGED INPUTS

	//1. Forward (In -> F)
	//Same as ForwardFFT with (In1 + In2) / 2 as input
	void ForwardFFT_AveragedInputs(VEC<DBL3> &In1, VEC<DBL3> &In2)
	{
		if (n.z == 1) ForwardFFT_2D(In1, In2);
		else ForwardFFT_3D(In1, In2);
	}

	//MULTIPLICATION

	//2. (F -> F2) -> single space version
	void KernelMultiplication(void)
	{
		
		if (n.z == 1) static_cast<Owner*>(this)->KernelMultiplication_2D(F, F2);
		else static_cast<Owner*>(this)->KernelMultiplication_3D(F, F2);
	}

	//2. (F -> F2) -> multiple input spaces version using a collection of FFT spaces (Kernel must be configured for this)
	void KernelMultiplication_MultipleInputs(std::vector<VEC<ReIm3>*>& Fcol)
	{
		if (n.z == 1) static_cast<Owner*>(this)->KernelMultiplication_2D(Fcol, F2);
		else static_cast<Owner*>(this)->KernelMultiplication_3D(Fcol, F2);
	}

	//SINGLE OUTPUT

	//3. Inverse. Return dot product of In with Out. (F2 -> Out)
	double InverseFFT(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut)
	{
		if (n.z == 1) return InverseFFT_2D(In, Out, clearOut);
		else return InverseFFT_3D(In, Out, clearOut);
	}

	//AVERAGED INPUTS, SINGLE OUTPUT

	//3. Inverse. Return dot product of In with Out. (F2 -> Out)
	//Same as InverseFFT with (In1 + In2) / 2 as input and output copied to both Out1 and Out2.
	double InverseFFT_AveragedInputs(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut)
	{
		if (n.z == 1) return InverseFFT_2D(In1, In2, Out, clearOut);
		else return InverseFFT_3D(In1, In2, Out, clearOut);
	}

	//AVERAGED INPUTS, DUPLICATED OUTPUTS

	//3. Inverse. Return dot product of In with Out. (F2 -> Out)
	//Same as InverseFFT with (In1 + In2) / 2 as input and output copied to both Out1 and Out2.
	double InverseFFT_AveragedInputs_DuplicatedOutputs(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut)
	{
		if (n.z == 1) return InverseFFT_2D(In1, In2, Out1, Out2, clearOut);
		else return InverseFFT_3D(In1, In2, Out1, Out2, clearOut);
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------- CONFIGURATION

template <typename Owner, typename Kernel>
BError Convolution<Owner, Kernel>::SetDimensions(SZ3 n_, DBL3 h_, bool multiplication_embedding_, INT3 pbc_images_)
{
	BError error(__FUNCTION__);

	error = SetConvolutionDimensions(n_, h_, multiplication_embedding_, pbc_images_);
	if (!error) error = static_cast<Owner*>(this)->AllocateKernelMemory();

	return error;
}

//-------------------------- RUN-TIME CONVOLUTION : 2D (multiplication embedded)

//SINGLE INPUT, SINGLE OUTPUT
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::Convolute_2D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut)
{
	//2D

	//1. FFTs along x
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {

		int tn = omp_get_thread_num();

		//write input into fft line (zero padding kept)
		for (int i = 0; i < n.x; i++) {

			int idx_in = i + j * n.x;

			*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = In[idx_in];
		}

		//fft on line
		fftw_execute(plan_fwd_x[tn]);

		//write line to fft array
		for (int i = 0; i < N.x / 2 + 1; i++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
		}
	}

	//2. FFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array (zero padding kept)
		for (int j = 0; j < n.y; j++) {

			*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_fwd_y[tn]);

		//3. kernel multiplication on line
		static_cast<Owner*>(this)->KernelMultiplication_2D_line(reinterpret_cast<ReIm3*>(pline[tn]), i);

		//4. ifft on line
		fftw_execute(plan_inv_y[tn]);

		//write line to fft array, truncating upper part (from n.y to N.y if different)
		for (int j = 0; j < n.y; j++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//5. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//write line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = In[i + j * n.x];

				Out[i + j * n.x] = Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}
	else {

		//5. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);


			//add line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = In[i + j * n.x];

				Out[i + j * n.x] += Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS, SINGLE OUTPUT
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::Convolute_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut)
{
	//2D

	//1. FFTs along x
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {

		int tn = omp_get_thread_num();

		//write input into fft line (zero padding kept)
		for (int i = 0; i < n.x; i++) {

			int idx_in = i + j * n.x;

			*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = (In1[idx_in] + In2[idx_in]) / 2;
		}

		//fft on line
		fftw_execute(plan_fwd_x[tn]);

		//write line to fft array
		for (int i = 0; i < N.x / 2 + 1; i++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
		}
	}

	//2. FFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array (zero padding kept)
		for (int j = 0; j < n.y; j++) {

			*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_fwd_y[tn]);

		//3. kernel multiplication on line
		static_cast<Owner*>(this)->KernelMultiplication_2D_line(reinterpret_cast<ReIm3*>(pline[tn]), i);

		//4. ifft on line
		fftw_execute(plan_inv_y[tn]);

		//write line to fft array, truncating upper part (from n.y to N.y if different)
		for (int j = 0; j < n.y; j++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//5. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//write line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out[i + j * n.x] = Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}
	else {

		//5. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);


			//add line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out[i + j * n.x] += Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS, DUPLICATED OUTPUTS
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::Convolute_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut)
{
	//2D

	//1. FFTs along x
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {

		int tn = omp_get_thread_num();

		//write input into fft line (zero padding kept)
		for (int i = 0; i < n.x; i++) {

			int idx_in = i + j * n.x;

			*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = (In1[idx_in] + In2[idx_in]) / 2;
		}

		//fft on line
		fftw_execute(plan_fwd_x[tn]);

		//write line to fft array
		for (int i = 0; i < N.x / 2 + 1; i++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
		}
	}

	//2. FFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array (zero padding kept)
		for (int j = 0; j < n.y; j++) {

			*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_fwd_y[tn]);

		//3. kernel multiplication on line
		static_cast<Owner*>(this)->KernelMultiplication_2D_line(reinterpret_cast<ReIm3*>(pline[tn]), i);

		//4. ifft on line
		fftw_execute(plan_inv_y[tn]);

		//write line to fft array, truncating upper part (from n.y to N.y if different)
		for (int j = 0; j < n.y; j++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//5. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//write line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x] ) / 2;

				Out1[i + j * n.x] = Out_val;
				Out2[i + j * n.x] = Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}
	else {

		//5. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);


			//add line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out1[i + j * n.x] += Out_val;
				Out2[i + j * n.x] += Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}

	return dot_product;
}

//-------------------------- RUN-TIME CONVOLUTION : 3D (multiplication embedded)

//SINGLE INPUT, SINGLE OUTPUT
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::Convolute_3D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut)
{	
	//1. FFTs along x
	for (int k = 0; k < n.z; k++) {
	#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line (zero padding kept)
			for (int i = 0; i < n.x; i++) {

				int idx_in = i + j * n.x + k * n.x * n.y;

				*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = In[idx_in];
			}

			//fft on line
			fftw_execute(plan_fwd_x[tn]);

			//write line to fft array
			for (int i = 0; i < N.x / 2 + 1; i++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
			}
		}
	}

	//2. FFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int j = 0; j < n.y; j++) {

				*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_y[tn]);

			//write line to fft array
			for (int j = 0; j < N.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	//3. FFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int k = 0; k < n.z; k++) {

				*reinterpret_cast<ReIm3*>(pline_zp_z[tn] + k * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_z[tn]);

			//4. kernel multiplication on line
			static_cast<Owner*>(this)->KernelMultiplication_3D_line(reinterpret_cast<ReIm3*>(pline[tn]), i, j);

			//5. ifft on line
			fftw_execute(plan_inv_z[tn]);

			//write line to fft array, truncating upper half
			for (int k = 0; k < n.z; k++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}

	//6. IFFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int j = 0; j < N.y; j++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_inv_y[tn]);

			//write line to fft array, truncating upper half
			for (int j = 0; j < n.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//7. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//write line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = In[i + j * n.x + k * n.x * n.y];

					Out[i + j * n.x + k * n.x * n.y] = Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}
	else {

		//7. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//add line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = In[i + j * n.x + k * n.x * n.y];

					Out[i + j * n.x + k * n.x * n.y] += Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS, SINGLE OUTPUT
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::Convolute_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut)
{
	//1. FFTs along x
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line (zero padding kept)
			for (int i = 0; i < n.x; i++) {

				int idx_in = i + j * n.x + k * n.x * n.y;

				*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = (In1[idx_in] + In2[idx_in]) / 2;
			}

			//fft on line
			fftw_execute(plan_fwd_x[tn]);

			//write line to fft array
			for (int i = 0; i < N.x / 2 + 1; i++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
			}
		}
	}

	//2. FFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int j = 0; j < n.y; j++) {

				*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_y[tn]);

			//write line to fft array
			for (int j = 0; j < N.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	//3. FFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int k = 0; k < n.z; k++) {

				*reinterpret_cast<ReIm3*>(pline_zp_z[tn] + k * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_z[tn]);

			//4. kernel multiplication on line
			static_cast<Owner*>(this)->KernelMultiplication_3D_line(reinterpret_cast<ReIm3*>(pline[tn]), i, j);

			//5. ifft on line
			fftw_execute(plan_inv_z[tn]);

			//write line to fft array, truncating upper half
			for (int k = 0; k < n.z; k++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}

	//6. IFFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int j = 0; j < N.y; j++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_inv_y[tn]);

			//write line to fft array, truncating upper half
			for (int j = 0; j < n.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//7. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//write line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out[i + j * n.x + k * n.x * n.y] = Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}
	else {

		//7. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//add line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out[i + j * n.x + k * n.x * n.y] += Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS, DUPLICATED OUTPUTS
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::Convolute_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut)
{
	//1. FFTs along x
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line (zero padding kept)
			for (int i = 0; i < n.x; i++) {

				int idx_in = i + j * n.x + k * n.x * n.y;

				*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = (In1[idx_in] + In2[idx_in]) / 2;
			}

			//fft on line
			fftw_execute(plan_fwd_x[tn]);

			//write line to fft array
			for (int i = 0; i < N.x / 2 + 1; i++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
			}
		}
	}

	//2. FFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int j = 0; j < n.y; j++) {

				*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_y[tn]);

			//write line to fft array
			for (int j = 0; j < N.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	//3. FFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int k = 0; k < n.z; k++) {

				*reinterpret_cast<ReIm3*>(pline_zp_z[tn] + k * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_z[tn]);

			//4. kernel multiplication on line
			static_cast<Owner*>(this)->KernelMultiplication_3D_line(reinterpret_cast<ReIm3*>(pline[tn]), i, j);

			//5. ifft on line
			fftw_execute(plan_inv_z[tn]);

			//write line to fft array, truncating upper half
			for (int k = 0; k < n.z; k++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}

	//6. IFFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int j = 0; j < N.y; j++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_inv_y[tn]);

			//write line to fft array, truncating upper half
			for (int j = 0; j < n.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//7. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//write line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out1[i + j * n.x + k * n.x * n.y] = Out_val;
					Out2[i + j * n.x + k * n.x * n.y] = Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}
	else {

		//7. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//add line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out1[i + j * n.x + k * n.x * n.y] += Out_val;
					Out2[i + j * n.x + k * n.x * n.y] += Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}

	return dot_product;
}

//-------------------------- RUN-TIME CONVOLUTION : 2D (multiplication not embedded)

//SINGLE INPUT

//1. Forward (In -> F)
template <typename Owner, typename Kernel>
void Convolution<Owner, Kernel>::ForwardFFT_2D(VEC<DBL3> &In)
{
	//2D

	//1. FFTs along x
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {

		int tn = omp_get_thread_num();

		//write input into fft line (zero padding kept)
		for (int i = 0; i < n.x; i++) {

			int idx_in = i + j * n.x;

			*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = In[idx_in];
		}

		//fft on line
		fftw_execute(plan_fwd_x[tn]);

		//write line to fft array
		for (int i = 0; i < N.x / 2 + 1; i++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
		}
	}

	//2. FFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array (zero padding kept)
		for (int j = 0; j < n.y; j++) {

			*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_fwd_y[tn]);

		//write line to fft array
		for (int j = 0; j < N.y; j++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}
}

//SINGLE OUTPUT

//3. Inverse. Return dot product of In with Out. (F2 -> Out)
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::InverseFFT_2D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut)
{
	//1. IFFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array
		for (int j = 0; j < N.y; j++) {

			*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F2[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_inv_y[tn]);

		//write line to fft array, truncating upper half
		for (int j = 0; j < n.y; j++) {

			F2[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//2. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//write line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = In[i + j * n.x];

				Out[i + j * n.x] = Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}
	else {

		//2. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//add line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = In[i + j * n.x];

				Out[i + j * n.x] += Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS

//1. Forward (In -> F)
template <typename Owner, typename Kernel>
void Convolution<Owner, Kernel>::ForwardFFT_2D(VEC<DBL3> &In1, VEC<DBL3> &In2)
{
	//2D

	//1. FFTs along x
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {

		int tn = omp_get_thread_num();

		//write input into fft line (zero padding kept)
		for (int i = 0; i < n.x; i++) {

			int idx_in = i + j * n.x;

			*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = (In1[idx_in] + In2[idx_in]) / 2;
		}

		//fft on line
		fftw_execute(plan_fwd_x[tn]);

		//write line to fft array
		for (int i = 0; i < N.x / 2 + 1; i++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
		}
	}

	//2. FFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array (zero padding kept)
		for (int j = 0; j < n.y; j++) {

			*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_fwd_y[tn]);

		//write line to fft array
		for (int j = 0; j < N.y; j++) {

			F[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}
}

//AVERAGED INPUTS, SINGLE OUTPUT

//3. Inverse. Return dot product of In with Out. (F2 -> Out)
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::InverseFFT_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut)
{
	//1. IFFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array
		for (int j = 0; j < N.y; j++) {

			*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F2[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_inv_y[tn]);

		//write line to fft array, truncating upper half
		for (int j = 0; j < n.y; j++) {

			F2[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//2. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//write line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out[i + j * n.x] = Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}
	else {

		//2. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//add line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out[i + j * n.x] += Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS, DUPLICATED OUTPUTS

//3. Inverse. Return dot product of In with Out. (F2 -> Out)
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::InverseFFT_2D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut)
{
	//1. IFFTs along y
#pragma omp parallel for
	for (int i = 0; i < (N.x / 2 + 1); i++) {

		int tn = omp_get_thread_num();

		//fetch line from fft array
		for (int j = 0; j < N.y; j++) {

			*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F2[i + j * (N.x / 2 + 1)];
		}

		//fft on line
		fftw_execute(plan_inv_y[tn]);

		//write line to fft array, truncating upper half
		for (int j = 0; j < n.y; j++) {

			F2[i + j * (N.x / 2 + 1)] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//2. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//write line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out1[i + j * n.x] = Out_val;
				Out2[i + j * n.x] = Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}
	else {

		//2. IFFTs along x
#pragma omp parallel for reduction(+:dot_product)
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line
			for (int i = 0; i < N.x / 2 + 1; i++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1)];
			}

			//fft on line
			fftw_execute(plan_inv_x[tn]);

			//add line to output
			for (int i = 0; i < n.x; i++) {

				DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
				DBL3 In_val = (In1[i + j * n.x] + In2[i + j * n.x]) / 2;

				Out1[i + j * n.x] += Out_val;
				Out2[i + j * n.x] += Out_val;

				dot_product += In_val * Out_val;
			}
		}
	}

	return dot_product;
}

//-------------------------- RUN-TIME CONVOLUTION : 3D (multiplication not embedded)

//SINGLE INPUT

//1. Forward (In -> F)
template <typename Owner, typename Kernel>
void Convolution<Owner, Kernel>::ForwardFFT_3D(VEC<DBL3> &In)
{
	//1. FFTs along x
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line (zero padding kept)
			for (int i = 0; i < n.x; i++) {

				int idx_in = i + j * n.x + k * n.x * n.y;

				*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = In[idx_in];
			}

			//fft on line
			fftw_execute(plan_fwd_x[tn]);

			//write line to fft array
			for (int i = 0; i < N.x / 2 + 1; i++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
			}
		}
	}

	//2. FFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int j = 0; j < n.y; j++) {

				*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_y[tn]);

			//write line to fft array
			for (int j = 0; j < N.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	//3. FFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int k = 0; k < n.z; k++) {

				*reinterpret_cast<ReIm3*>(pline_zp_z[tn] + k * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_z[tn]);

			//write line to fft array
			for (int k = 0; k < N.z; k++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}
}

//SINGLE OUTPUT

//3. Inverse. Return dot product of In with Out. (F2 -> Out)
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::InverseFFT_3D(VEC<DBL3> &In, VEC<DBL3> &Out, bool clearOut)
{
	//1. IFFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int k = 0; k < N.z; k++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + k * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//ifft on line
			fftw_execute(plan_inv_z[tn]);

			//write line to fft array, truncating upper half
			for (int k = 0; k < n.z; k++) {

				F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}

	//2. IFFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int j = 0; j < N.y; j++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_inv_y[tn]);

			//write line to fft array, truncating upper half
			for (int j = 0; j < n.y; j++) {

				F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//3. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//write line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = In[i + j * n.x + k * n.x * n.y];

					Out[i + j * n.x + k * n.x * n.y] = Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}
	else {

		//3. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//add line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = In[i + j * n.x + k * n.x * n.y];

					Out[i + j * n.x + k * n.x * n.y] += Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS

//1. Forward (In -> F)
template <typename Owner, typename Kernel>
void Convolution<Owner, Kernel>::ForwardFFT_3D(VEC<DBL3> &In1, VEC<DBL3> &In2)
{
	//1. FFTs along x
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {

			int tn = omp_get_thread_num();

			//write input into fft line (zero padding kept)
			for (int i = 0; i < n.x; i++) {

				int idx_in = i + j * n.x + k * n.x * n.y;

				*reinterpret_cast<DBL3*>(pline_zp_x[tn] + i * 3) = (In1[idx_in] + In2[idx_in]) / 2;
			}

			//fft on line
			fftw_execute(plan_fwd_x[tn]);

			//write line to fft array
			for (int i = 0; i < N.x / 2 + 1; i++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + i * 3);
			}
		}
	}

	//2. FFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int j = 0; j < n.y; j++) {

				*reinterpret_cast<ReIm3*>(pline_zp_y[tn] + j * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_y[tn]);

			//write line to fft array
			for (int j = 0; j < N.y; j++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	//3. FFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array (zero padding kept)
			for (int k = 0; k < n.z; k++) {

				*reinterpret_cast<ReIm3*>(pline_zp_z[tn] + k * 3) = F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_fwd_z[tn]);

			//write line to fft array
			for (int k = 0; k < N.z; k++) {

				F[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}
}

//AVERAGED INPUTS, SINGLE OUTPUT

//3. Inverse. Return dot product of In with Out. (F2 -> Out)
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::InverseFFT_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out, bool clearOut)
{
	//1. IFFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int k = 0; k < N.z; k++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + k * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//ifft on line
			fftw_execute(plan_inv_z[tn]);

			//write line to fft array, truncating upper half
			for (int k = 0; k < n.z; k++) {

				F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}

	//2. IFFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int j = 0; j < N.y; j++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_inv_y[tn]);

			//write line to fft array, truncating upper half
			for (int j = 0; j < n.y; j++) {

				F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//3. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//write line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out[i + j * n.x + k * n.x * n.y] = Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}
	else {

		//3. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//add line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out[i + j * n.x + k * n.x * n.y] += Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}

	return dot_product;
}

//AVERAGED INPUTS, DUPLICATED OUTPUTS

//3. Inverse. Return dot product of In with Out. (F2 -> Out)
template <typename Owner, typename Kernel>
double Convolution<Owner, Kernel>::InverseFFT_3D(VEC<DBL3> &In1, VEC<DBL3> &In2, VEC<DBL3> &Out1, VEC<DBL3> &Out2, bool clearOut)
{
	//1. IFFTs along z
#pragma omp parallel for
	for (int j = 0; j < N.y; j++) {
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int k = 0; k < N.z; k++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + k * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//ifft on line
			fftw_execute(plan_inv_z[tn]);

			//write line to fft array, truncating upper half
			for (int k = 0; k < n.z; k++) {

				F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + k * 3);
			}
		}
	}

	//2. IFFTs along y
	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int i = 0; i < (N.x / 2 + 1); i++) {

			int tn = omp_get_thread_num();

			//fetch line from fft array
			for (int j = 0; j < N.y; j++) {

				*reinterpret_cast<ReIm3*>(pline[tn] + j * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
			}

			//fft on line
			fftw_execute(plan_inv_y[tn]);

			//write line to fft array, truncating upper half
			for (int j = 0; j < n.y; j++) {

				F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = *reinterpret_cast<ReIm3*>(pline[tn] + j * 3);
			}
		}
	}

	double dot_product = 0;

	if (clearOut) {

		//3. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//write line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out1[i + j * n.x + k * n.x * n.y] = Out_val;
					Out2[i + j * n.x + k * n.x * n.y] = Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}
	else {

		//3. IFFTs along x
		for (int k = 0; k < n.z; k++) {
#pragma omp parallel for reduction(+:dot_product)
			for (int j = 0; j < n.y; j++) {

				int tn = omp_get_thread_num();

				//write input into fft line
				for (int i = 0; i < N.x / 2 + 1; i++) {

					*reinterpret_cast<ReIm3*>(pline[tn] + i * 3) = F2[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y];
				}

				//fft on line
				fftw_execute(plan_inv_x[tn]);

				//add line to output
				for (int i = 0; i < n.x; i++) {

					DBL3 Out_val = *reinterpret_cast<DBL3*>(pline_rev_x[tn] + i * 3) / N.dim();
					DBL3 In_val = (In1[i + j * n.x + k * n.x * n.y] + In2[i + j * n.x + k * n.x * n.y]) / 2;

					Out1[i + j * n.x + k * n.x * n.y] += Out_val;
					Out2[i + j * n.x + k * n.x * n.y] += Out_val;

					dot_product += In_val * Out_val;
				}
			}
		}
	}

	return dot_product;
}