#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_OERSTED

#include "BorisCUDALib.h"

#include "ConvolutionDataCUDA.h"

#include "ErrorHandler.h"



////////////////////////////////////////////////////////////////////////////////////////////////
//
// Oersted Kernel calculated from tensors

//This must be used as a template parameter in Convolution class.

class OerstedKernelCUDA :
	public virtual ConvolutionDataCUDA		//virtual to solve the diamond problem as it's also inherited by the Convolution class which inherits from this Kernel
{

private:

	// Kernels for 3D Oersted field : Kxy, Kxz, Kyz; (imaginary parts only, real parts are zero)
	cu_obj<cuVEC<cuReal3>> KOe;

private:

	//-------------------------- KERNEL CALCULATION

	BError Calculate_Oersted_Kernels_2D(void);
	BError Calculate_Oersted_Kernels_3D(void);

protected:

	//-------------------------- CONSTRUCTOR

	OerstedKernelCUDA(void) {}

	virtual ~OerstedKernelCUDA() {}

	//-------------------------- MEMORY ALLOCATION

	//Called by SetDimensions in ConvolutionCUDA class
	BError AllocateKernelMemory(void);

	//-------------------------- KERNEL CALCULATION

	//this initializes the convolution kernels for the given mesh dimensions. 2D is for n.z == 1.
	BError Calculate_Oersted_Kernels(void)
	{
		if (n.z == 1) return Calculate_Oersted_Kernels_2D();
		else return Calculate_Oersted_Kernels_3D();
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
