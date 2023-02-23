#pragma once

#include "Boris_Enums_Defs.h"
#ifdef MODULE_COMPILATION_OERSTED

#include "BorisLib.h"
#include "OerstedTFunc.h"

#include "ConvolutionData.h"



////////////////////////////////////////////////////////////////////////////////////////////////
//
// Oersted Kernel calculated from tensors

//This must be used as a template parameter in Convolution class.

class OerstedKernel :
	public virtual ConvolutionData		//virtual to solve the diamond problem as it's also inherited by the Convolution class which inherits from this Kernel
{

private:

	//Kernels for 3D Oersted field : Kxy, Kxz, Kyz (imaginary parts only, real parts are zero)
	VEC<DBL3> KOe;

private:

	//-------------------------- KERNEL CALCULATION

	BError Calculate_Oersted_Kernels_2D(void);
	BError Calculate_Oersted_Kernels_3D(void);

protected:

	//-------------------------- CONSTRUCTOR

	OerstedKernel(void) {}

	virtual ~OerstedKernel() {}

	//-------------------------- MEMORY ALLOCATION

	//Called by SetDimensions in Convolution class
	BError AllocateKernelMemory(void);

	//-------------------------- KERNEL CALCULATION

	//this initializes the convolution kernels for the given mesh dimensions. 2D is for n.z == 1.
	BError Calculate_Oersted_Kernels(void)
	{
		if (n.z == 1) return Calculate_Oersted_Kernels_2D();
		else return Calculate_Oersted_Kernels_3D();
	}

	//-------------------------- RUN-TIME KERNEL MULTIPLICATION

	//Called by Convolute_2D/Convolute_3D methods in Convolution class : define pointwise multiplication of In with Kernels and set result in Out
	void KernelMultiplication_2D(VEC<ReIm3>& In, VEC<ReIm3>& Out);
	void KernelMultiplication_3D(VEC<ReIm3>& In, VEC<ReIm3>& Out);

	//multiple input spaces version, used when not embedding the kenel multiplication
	void KernelMultiplication_2D(std::vector<VEC<ReIm3>*>& Incol, VEC<ReIm3>& Out) {}
	void KernelMultiplication_3D(std::vector<VEC<ReIm3>*>& Incol, VEC<ReIm3>& Out) {}

	//kernel multiplications for a single line : used for embedding

	//multiply kernels in line along y direction (so use stride of (N.x/2 + 1) to read from kernels), starting at given i index (this must be an index in the first x row)
	void KernelMultiplication_2D_line(ReIm3* pline, int i) {}

	//multiply kernels in line along z direction (so use stride of (N.x/2 + 1) * N.y to read from kernels), starting at given i and j indexes (these must be an indexes in the first xy plane)
	void KernelMultiplication_3D_line(ReIm3* pline, int i, int j);
};

#endif

