#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

//size of N.y required to switch to transpose_xy mode
//This is also the threshold of N.x/2 required to switch to interleaving transposing operations for the 3 vector components, rather than launching a single kernel call
#define TRANSPOSE_XY_YTHRESHOLD	2048

class ConvolutionDataCUDA
{

protected:

	//mesh dimensions
	cuSZ3 n;

	//mesh cellsize
	cuReal3 h;

	//dimensions of FFT spaces
	cuSZ3 N;

	//dimensions of FFT space in GPU memory
	cu_obj<cuSZ3> cuN;
	
	//dimensions of complex fft spaces in GPU memory (N.x/2 + 1)*N.y*N.z, and N.y*(N.x/2 + 1)*N.z respectively
	cu_obj<cuSZ3> cuNc_xy, cuNc_yx;
	//dimensions of complex fft spaces in GPU memory for q2d mode (N.x/2 + 1)*N.y*N.z/2, and N.y*(N.x/2 + 1)*N.z/2 respectively
	cu_obj<cuSZ3> cuNc_xy_q2d, cuNc_yx_q2d;
	//quarter size scratch space sizes as: (N.x/2 + 1)*n.y*n.z for xy transpose and n*y*(N.x/2 + 1)*n.z for yx transpose respectively
	cu_obj<cuSZ3> cuNcquart_xy, cuNcquart_yx;

	//regions for zero padding
	cu_obj<cuRect> cuUpper_y_region, cuUpper_y_transposed_region, cuUpper_z_region, cuUpper_z_region_q2d;

	//cuFFT plans : 2D
	cufftHandle plan2D_fwd_x, plan2D_y;
	cufftHandle plan2D_inv_x;
	
	//cuFFT plans : 3D
	cufftHandle plan3D_fwd_x, plan3D_y, plan3D_z;
	cufftHandle plan3D_inv_x;

	//FFT arrays in GPU memory

	//Input / output space of size N.x * n.y * n.z; zero padding for input from n.x to N.x must be kept by computations
	cu_arr<cuReal> cuIn_x, cuIn_y, cuIn_z;
	cu_arr<cuReal> cuOut_x, cuOut_y, cuOut_z;

	//Scratch space (N.x / 2 + 1) * N.y * N.z (except in q2D mode where last dimension is n.z, same as in 2D mode -> q2D mode is disabled if not using the embedded convolution pipeline)
	cu_arr<cuComplex> cuS_x, cuS_y, cuS_z;

	//by default the embedded convolution pipeline is used (i.e. FFT -> Mult -> iFFT are done one after another)
	//if you want to break down the convolution into separate calls then set embed_multiplication = false.
	//This disables q2D mode and allocates additional scratch space cuS2
	bool embed_multiplication = true;

	//additiona scratch space used if not using embedded convolution
	cu_arr<cuComplex> cuS2_x, cuS2_y, cuS2_z;

	//Quarter-size scratch space used for x fft/ifft when the transpose_xy mode is used.
	//After the x ffts this space is transposed into the main scratch space.
	//Similarly before the x iffts, the required part from the main scratch space is transposed into the quarter scratch space.
	cu_arr<cuComplex> cuSquart_x, cuSquart_y, cuSquart_z;

	//transpose xy planes before doing the y direction fft/ifft?
	//for 2D mode we have a choice : transpose_xy mode triggered if N.y is above a set threshold
	//for 3D (and q2D) mode always use transpose_xy as it turns out it's a very good catch-all approach
	bool transpose_xy = true;

	//transpose_xy but in gpu memory so we can pass it to a __global__ quickly
	cu_obj<bool> transpose_xy_gpu;

	//quasi 2D mode : 3D mode but with the z-axis fft / kernel multiplication / z-axis ifft rolled into one step
	//currently handled values:
	//0 : q2D disabled
	//4 : N.z = 4 (n.z = 2)
	//8 : N.z = 8 (n.z = 3, 4)
	//16 : N.z = 16 (n.z = 5, 6, 7, 8)
	//32 : N.z = 16 (n.z = 9, 10, ..., 16)
	int q2D_level = 0;

protected:

	//-------------------------- CONSTRUCTORS

	ConvolutionDataCUDA(void) {}

	virtual ~ConvolutionDataCUDA();

	//-------------------------- CONFIGURATION

	BError SetConvolutionDimensions(cuSZ3 n_, cuReal3 h_, bool embed_multiplication_ = true);

	//-------------------------- RUN-TIME METHODS

	//Input

	//Copy data to cuSx, cuSy, cuSz arrays at start of convolution iteration
	template <typename cuVECIn>
	void CopyInputData(cu_obj<cuVECIn>& In);

	//FFTs

	void forward_fft_2D(void);
	void inverse_fft_2D(void);
	void inverse_fft_2D_2(void);

	void forward_fft_3D(void);
	void inverse_fft_3D(void);
	void inverse_fft_3D_2(void);

	void forward_fft_q2D(void);
	void inverse_fft_q2D(void);

	//Output

	//Copy convolution result (in cuS arrays) to output and obtain energy value : product of In with Out times -MU0 / (2 * non_empty_points), where non_empty_points = In.get_nonempty_points();
	template <typename cuVECIn, typename cuVECOut>
	void FinishConvolution_Set(cu_obj<cuVECIn>& In, cu_obj<cuVECOut>& Out, cu_obj<cuReal>& energy, bool get_energy);
	
	template <typename cuVECIn, typename cuVECOut>
	void FinishConvolution_Add(cu_obj<cuVECIn>& In, cu_obj<cuVECOut>& Out, cu_obj<cuReal>& energy, bool get_energy);

	//-------------------------- GETTERS

	//Get pointer to the S scratch space
	cu_arr<cuComplex>* Get_Input_Scratch_Space_x(void) { return &cuS_x; }
	cu_arr<cuComplex>* Get_Input_Scratch_Space_y(void) { return &cuS_y; }
	cu_arr<cuComplex>* Get_Input_Scratch_Space_z(void) { return &cuS_z; }

	//Get pointer to the S2 scratch space
	cu_arr<cuComplex>* Get_Output_Scratch_Space_x(void) { return &cuS2_x; }
	cu_arr<cuComplex>* Get_Output_Scratch_Space_y(void) { return &cuS2_y; }
	cu_arr<cuComplex>* Get_Output_Scratch_Space_z(void) { return &cuS2_z; }
};

#endif