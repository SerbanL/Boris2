#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include <memory>

#include "BorisCUDALib.h"

#include "ConvolutionDataCUDA.h"
#include "DemagKernelCollectionCUDA_KerType.h"
#include "ErrorHandler.h"



////////////////////////////////////////////////////////////////////////////////////////////////
//
// Collection of demag kernels calculated from demag tensor, which can also use shifts. Intended for multi-layered demag calculations.

//This must be used as a template parameter in Convolution class.

class DemagKernelCollectionCUDA :
	public virtual ConvolutionDataCUDA		//virtual to solve the diamond problem as it's also inherited by the Convolution class which inherits from this Kernel
{

private:

	//For the vectors below the order must match, i.e. the nth rect in Rect_collection corresponds to the nth kernel
	//These are the kernels used to obtain values in this mesh using contributions from those respective meshes.

	//kernels used to multiply fft-ed magnetization from other layers
	std::vector<std::shared_ptr<cu_obj<cuKerType>>> kernels;

	//kernels collected in a cu_arr so we can pass it to a __global__ function
	cu_arr<cuKerType> kernels_gpu;

	//number of kernels in gpu memory
	cu_obj<int> num_kernels;

	//if z-shifted, re-use already calculated kernel if the only difference is just the sign of the shift : use kernel symmetries to recover correct multiplication
	//this vector has 1-2-1 correspondence with kernels
	std::vector<bool> inverse_shifted;
	
	//this also needs to be passed to a __global__ so build it from inverse_shifted
	cu_arr<bool> inverse_shifted_gpu;

	//also keep a copy of this in cpu memory for convenience (these flags also found in kernels)
	std::vector<bool> zshifted;

	//collection of all mesh rectangles participating in convolution
	//All meshes must have same number of cells, thus you can determine the cellsize used in each mesh
	std::vector<cuRect> Rect_collection;

	//also need to know what the rectangle for this mesh is (i.e. the mesh that this DemagKernelCollection object applies to)
	cuRect this_rect;

	//each kernel collection has exactly one self contribution : store the index here (index in kernels and hence also Rect_collection)
	int self_contribution_index;

	//maximum common cellsize dimension used to normalize dimensions
	cuBReal h_max;

	//configure shifts and source and destination cellsizes so kernels apply for input to multiple outputs - true - (rather than multiple inputs to single output)
	//Always use  false here - it's faster! Included for testing only.
	bool use_multiple_outputs = false;

private:

	//-------------------------- KERNEL CALCULATION

	//2D layers, real kernels for self demag (Kdiag_real, and K2D_odiag, with full use of kernel symmetries)
	BError Calculate_Demag_Kernels_2D_Self(int index, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_2D_Self_onGPU(int index);

	//2D layers, z shift only : Kernels can be stored as real with use of kernel symmetries. Kxx, Kyy, Kzz, Kxy real, Kxz, Kyz imaginary
	BError Calculate_Demag_Kernels_2D_zShifted(int index, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_2D_zShifted_onGPU(int index);

	//2D layers, complex kernels most general case (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
	BError Calculate_Demag_Kernels_2D_Complex_Full(int index, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_2D_Complex_Full_onGPU(int index);

	//3D layers, real kernels for self demag (Kdiag_real, and Kodiag_real, with full use of kernel symmetries)
	BError Calculate_Demag_Kernels_3D_Self(int index, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_3D_Self_onGPU(int index);

	//3D layers, z shift only : Kernels can be stored with use of kernel symmetries (but still complex).
	BError Calculate_Demag_Kernels_3D_zShifted(int index, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_3D_zShifted_onGPU(int index);

	//3D layers, complex kernels most general case (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
	BError Calculate_Demag_Kernels_3D_Complex_Full(int index, bool initialize_on_gpu);
	BError Calculate_Demag_Kernels_3D_Complex_Full_onGPU(int index);

	//search to find a matching kernel that has already been computed and return pointer to it -> kernel can be identified from shift, source and destination discretisation
	std::shared_ptr<cu_obj<cuKerType>> KernelAlreadyComputed(cuReal3 shift, cuReal3 h_src, cuReal3 h_dst);

	//Auxiliary for kernel computations on the GPU

	//copy Re or Im parts of cuOut to cuIn
	void cuOut_to_cuIn_Re(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut);
	void cuOut_to_cuIn_Im(size_t size, cu_arr<cufftDoubleReal>& cuIn, cu_arr<cufftDoubleComplex>& cuOut);

protected:

	//-------------------------- CONSTRUCTOR

	DemagKernelCollectionCUDA(void) {}

	virtual ~DemagKernelCollectionCUDA() {}

	//-------------------------- MEMORY ALLOCATION

	//Called by SetDimensions in ConvolutionCUDA class
	BError AllocateKernelMemory(void);

	//-------------------------- SETTERS

	//Set all the rectangles participating in convolution. This determines the number of kernels needed : one for each mesh.
	BError Set_Rect_Collection(std::vector<cuRect>& Rect_collection_, cuRect this_rect_, cuBReal h_max_);

	//-------------------------- GETTERS

	std::shared_ptr<cu_obj<cuKerType>> Get_Kernel(int index) { return kernels[index]; }

	bool is_inverse_shifted(int index) { return inverse_shifted[index]; }

	//-------------------------- KERNEL CALCULATION

	//this initializes all the convolution kernels for the given mesh dimensions.
	//use information for other DemagKernelCollection objects in the set so we re-use kernels as much as possible
	BError Calculate_Demag_Kernels(std::vector<DemagKernelCollectionCUDA*>& kernelCollection, bool initialize_on_gpu = true);

	//-------------------------- RUN-TIME KERNEL MULTIPLICATION

	//Called by Convolute_2D/Convolute_3D methods in ConvolutionCUDA class : define pointwise multiplication with Kernels (using the cuSx, cuSy and cuSz arrays)
	//single kernel multiplication versions not used here
	void KernelMultiplication_2D(void) {}
	void KernelMultiplication_3D(void) {}

	//Kernel multiplication in quasi-2D mode : z-axis fft / kernel multiplication / z-axis ifft rolled into one (but do not divide by N for the ifft)
	//not use here
	void KernelMultiplication_q2D(int q2D_level) {}

	//Multiple inputs to single output - use stand-alone __global__ functions for multiplication
	//These are faster than the single input to multiple output versions below, with speedup factor increasing the larger the number of layers
	//These versions involve n^2 kernel launches, whilst the multiple output versions involve n kernel launches (but each n times larger and with poor load balancing)

	void KernelMultiplication_2D(
		std::vector<cu_arr<cuBComplex>*>& Incol_x, std::vector<cu_arr<cuBComplex>*>& Incol_y, std::vector<cu_arr<cuBComplex>*>& Incol_z,
		cu_arr<cuBComplex>& Out_x, cu_arr<cuBComplex>& Out_y, cu_arr<cuBComplex>& Out_z);

	void KernelMultiplication_3D(
		std::vector<cu_arr<cuBComplex>*>& Incol_x, std::vector<cu_arr<cuBComplex>*>& Incol_y, std::vector<cu_arr<cuBComplex>*>& Incol_z,
		cu_arr<cuBComplex>& Out_x, cu_arr<cuBComplex>& Out_y, cu_arr<cuBComplex>& Out_z);
};

#endif

#endif

