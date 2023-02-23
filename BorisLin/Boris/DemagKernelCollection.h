#pragma once

#include "Boris_Enums_Defs.h"
#ifdef MODULE_COMPILATION_SDEMAG

#include "BorisLib.h"
#include "DemagTFunc.h"

#include "ConvolutionData.h"

//NOTES:

//Multiple options available for kernel multiplications:
//
//1. Kernels calculated to collect multiple sources (inputs) into single destination (this) - this is the currently used method.
//
//2. Kernel calculated to contribute to multiple destinations from single source (this) - tested this method but comes out slightly slower than 1. and not as elegant to implement.
//Note to calculate these kernels all you need to do is inverse the shift direction, and swap source and destination cellsizes (in 2D irregular mode they can differ).
//Code deleted to reduce clutter but really quick to redo.
//
//3. Partially embedding kernel multiplication into forward FFT. 
//After taking forward FFT, at the last step rather than writing FFT result to F scratch space, perform kernel multiplications for multiple destinations and add result in F2 scratch spaces.
//Thus kernel multiplications embedded with forward FFT, but inverse IFFT completely separate. This means the F scratch space can be half-sized.
//However, this is much slower than 1. and 2. in all scenarios tested. Probably due to constantly flushing cache lines used for forward FFT.
//Abandon this entirely (slight memory saving not worth it); also ugly to implement.
//Code deleted to reduce clutter but really quick to redo.
//
//4. Prioritise kernel cache use. This means kernel multiplications are done so as to minimise kernel reading operations.
//Thus identify kernels which are used multiple times, making a list of sources and destinations for them, then perform multiplications in this order.
//First instinct this should be the best method since it minimises the read operations from kernels; it's also fairly elegant to implement.
//I've done this in SDemag_KernelCollection.h (deleted now - dead code policing!), but it comes out a fair bit slower than method 1. so not in current use.

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Collection of demag kernels calculated from demag tensor, which can also use shifts. Intended for multi-layered demag calculations.

//collect various kernel types in this struct used in demag kernel collections
struct KerType {

	//2D off-diagonal for self-demag (only need real version)
	std::vector<double> K2D_odiag;

	//3D Complex
	VEC<ReIm3> Kdiag_cmpl, Kodiag_cmpl;

	//3D Real
	VEC<DBL3> Kdiag_real, Kodiag_real;

	//if z shifted only then kernel multiplications and storage can be made more efficient for 2D. For 3D we can make storage more efficient.
	bool zshifted = false;

	//has this kernel been calculated and ready to use?
	bool kernel_calculated = false;

	//shift, source and destination cellsizes for which this kernel applies (not normalized)
	DBL3 shift, h_src, h_dst;

	//The dimensions of the FFT space (In and Out in the kernel multiplication methods)
	SZ3 N;

	//-----------------------

	KerType(void) {}

	BError AllocateKernels(Rect from_rect, Rect this_rect, SZ3 N_);

	void FreeKernels(void);
};

//This must be used as a template parameter in Convolution class.

class DemagKernelCollection :
	public virtual ConvolutionData		//virtual to solve the diamond problem as it's also inherited by the Convolution class which inherits from this Kernel
{

private:

	//For the vectors below the order must match, i.e. the nth rect in Rect_collection corresponds to the nth kernel
	//These are the kernels used to obtain values in this mesh using contributions from those respective meshes.

	//kernels used to multiply fft-ed magnetization from other layers
	std::vector<std::shared_ptr<KerType>> kernels;

	//if z-shifted, re-use already calculated kernel if the only difference is just the sign of the shift : use kernel symmetries to recover correct multiplication
	//this vector has 1-2-1 correspondence with kernels
	std::vector<bool> inverse_shifted;

	//collection of all mesh rectangles participating in convolution
	//All meshes must have same number of cells, thus you can determine the cellsize used in each mesh
	std::vector<Rect> Rect_collection;

	//also need to know what the rectangle for this mesh is (i.e. the mesh that this DemagKernelCollection object applies to)
	Rect this_rect;

	//each kernel collection has exactly one self contribution : store the index here (index in kernels and hence also Rect_collection)
	int self_contribution_index;

	//maximum common cellsize dimension used to normalize dimensions
	double h_max;

private:

	//-------------------------- KERNEL CALCULATION

	//2D layers, real kernels for self demag (Kdiag_real, and K2D_odiag, with full use of kernel symmetries)
	BError Calculate_Demag_Kernels_2D_Self(int index);
	
	//2D layers, z shift only : Kernels can be stored as real with use of kernel symmetries. Kxx, Kyy, Kzz, Kxy real, Kxz, Kyz imaginary
	BError Calculate_Demag_Kernels_2D_zShifted(int index);

	//2D layers, complex kernels most general case (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
	BError Calculate_Demag_Kernels_2D_Complex_Full(int index);

	//3D layers, real kernels for self demag (Kdiag_real, and Kodiag_real, with full use of kernel symmetries)
	BError Calculate_Demag_Kernels_3D_Self(int index);

	//3D layers, z shift only : Kernels can be stored with use of kernel symmetries (but still complex).
	BError Calculate_Demag_Kernels_3D_zShifted(int index);

	//3D layers, complex kernels most general case (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
	BError Calculate_Demag_Kernels_3D_Complex_Full(int index);

	//search to find a matching kernel that has already been computed and return pointer to it -> kernel can be identified from shift, source and destination discretisation
	std::shared_ptr<KerType> KernelAlreadyComputed(DBL3 shift, DBL3 h_src, DBL3 h_dst);

	//-------------------------- RUN-TIME KERNEL MULTIPLICATION

	//These compute the self contributions using the real kernels with full use of symmetries
	//These set the output (unless specified otherwise), not add into it, so always call it first -> each kernel collection has exactly one self contribution
	void KernelMultiplication_2D_Self(VEC<ReIm3>& In, VEC<ReIm3>& Out, bool set_output = true);
	void KernelMultiplication_3D_Self(VEC<ReIm3>& In, VEC<ReIm3>& Out, bool set_output = true);

	//z shifted regular version
	void KernelMultiplication_2D_zShifted(VEC<ReIm3>& In, VEC<ReIm3>& Out, VEC<DBL3>& Kdiag, VEC<DBL3>& Kodiag);
	//z shifted but with kernel calculated for the other direction shift, so adjust multiplications
	void KernelMultiplication_2D_inversezShifted(VEC<ReIm3>& In, VEC<ReIm3>& Out, VEC<DBL3>& Kdiag, VEC<DBL3>& Kodiag);

	//z shifted for 3D : complex kernels, but use kernel symmetries
	void KernelMultiplication_3D_zShifted(VEC<ReIm3>& In, VEC<ReIm3>& Out, VEC<ReIm3>& Kdiag, VEC<ReIm3>& Kodiag);

protected:

	//-------------------------- CONSTRUCTOR

	DemagKernelCollection(void) {}

	virtual ~DemagKernelCollection() {}

	//-------------------------- MEMORY ALLOCATION

	//Called by SetDimensions in Convolution class
	BError AllocateKernelMemory(void);

	//-------------------------- SETTERS

	//Set all the rectangles participating in convolution. This determines the number of kernels needed : one for each mesh.
	BError Set_Rect_Collection(std::vector<Rect>& Rect_collection_, Rect this_rect_, double h_max_, int self_contribution_index_);

	//-------------------------- GETTERS

	std::shared_ptr<KerType> Get_Kernel(int index) { return kernels[index]; }

	bool is_inverse_shifted(int index) { return inverse_shifted[index]; }

	//-------------------------- KERNEL CALCULATION

	//this initializes all the convolution kernels for the given mesh dimensions.
	//use information for other DemagKernelCollection objects in the set so we re-use kernels as much as possible
	BError Calculate_Demag_Kernels(std::vector<DemagKernelCollection*>& kernelCollection);

	//-------------------------- RUN-TIME KERNEL MULTIPLICATION

	//Called by Convolute_2D/Convolute_3D methods in Convolution class : define pointwise multiplication of In with Kernels and set result in Out

	//These versions should be used when not embedding the kernel multiplication
	//not used here
	void KernelMultiplication_2D(VEC<ReIm3>& In, VEC<ReIm3>& Out) {}
	void KernelMultiplication_3D(VEC<ReIm3>& In, VEC<ReIm3>& Out) {}

	//multiple input spaces version, used without embedding the kernel multiplication
	void KernelMultiplication_2D(std::vector<VEC<ReIm3>*>& Incol, VEC<ReIm3>& Out);
	void KernelMultiplication_3D(std::vector<VEC<ReIm3>*>& Incol, VEC<ReIm3>& Out);

	//kernel multiplications for a single line : used for embedding. Not used by DemagKernelCollection.
	void KernelMultiplication_2D_line(ReIm3* pline, int i) {}
	void KernelMultiplication_3D_line(ReIm3* pline, int i, int j) {}
};

#endif