#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined MODULE_COMPILATION_DEMAG || defined MODULE_COMPILATION_SDEMAG

#include "BorisCUDALib.h"

#include "ConvolutionDataCUDA.h"

#include "ErrorHandler.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag Kernel calculated from demag tensor

//This must be used as a template parameter in Convolution class.

class DemagKernelCUDA :
	public virtual ConvolutionDataCUDA		//virtual to solve the diamond problem as it's also inherited by the Convolution class which inherits from this Kernel
{

private:

	//off-diagonal Kernel used for 2D only (real parts only, imaginary parts are zero)
	cu_obj<cuVEC<cuBReal>> K2D_odiag;

	//Kernels for 3D : GPU memory
	//Kdiag : Kx, Ky, Kz; Kodiag : Kxy, Kxz, Kyz; (real parts only, imaginary parts are zero)
	cu_obj<cuVEC<cuReal3>> Kdiag, Kodiag;

private:

	//-------------------------- KERNEL CALCULATION

	BError Calculate_Demag_Kernels_2D(bool include_self_demag, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_3D(bool include_self_demag, bool initialize_on_gpu);

	//these versions compute the kernels entirely on the GPU, but in double precision
	BError Calculate_Demag_Kernels_2D_onGPU(bool include_self_demag);
	BError Calculate_Demag_Kernels_3D_onGPU(bool include_self_demag);

	//Auxiliary for kernel computations on the GPU

	//copy Re or Im parts of cuOut to cuIn
	void cuOut_to_cuIn_Re(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut);
	void cuOut_to_cuIn_Im(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut);

	//Copy Re parts of cuOut to Kdiag component (1: Kx, 2: Ky, 3: Kz). Takes into account transpose_xy flag.
	void cuOut_to_Kdiagcomponent(cu_arr<cufftDoubleComplex>& cuOut, int component);

	//Copy -(Re, Im, Im) parts of cuOut to Kodiag component (1: Kxy, 2: Kxz, 3: Kyz). Takes into account transpose_xy flag.
	void cuOut_to_Kodiagcomponent(cu_arr<cufftDoubleComplex>& cuOut, int component);
	
	//Copy -Im parts of cuOut to K2D_odiag. Takes into account transpose_xy flag.
	void cuOut_to_K2D_odiag(cu_arr<cufftDoubleComplex>& cuOut);

protected:

	//-------------------------- CONSTRUCTOR

	DemagKernelCUDA(void) {}

	virtual ~DemagKernelCUDA() {}

	//-------------------------- MEMORY ALLOCATION

	//Called by SetDimensions in ConvolutionCUDA class
	BError AllocateKernelMemory(void);

	//-------------------------- KERNEL CALCULATION

	//this initializes the convolution kernels for the given mesh dimensions. 2D is for n.z == 1.
	BError Calculate_Demag_Kernels(bool include_self_demag = true, bool initialize_on_gpu = true)
	{
		if (n.z == 1) return Calculate_Demag_Kernels_2D(include_self_demag, initialize_on_gpu);
		else return Calculate_Demag_Kernels_3D(include_self_demag, initialize_on_gpu);
	}

	//-------------------------- RUN-TIME KERNEL MULTIPLICATION

	//Called by Convolute_2D/Convolute_3D methods in ConvolutionCUDA class : define pointwise multiplication with Kernels (using the cuSx, cuSy and cuSz arrays)
	void KernelMultiplication_2D(void);
	void KernelMultiplication_3D(void);

	//Kernel multiplication in quasi-2D mode : z-axis fft / kernel multiplication / z-axis ifft rolled into one (but do not divide by N for the ifft)
	void KernelMultiplication_q2D(int q2D_level);
};

#endif

#endif

