#pragma once

#include "Boris_Enums_Defs.h"
#ifdef MODULE_SDEMAG

#include "DemagKernelCollection.h"

//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE

//NOTE, some improvements could be made but not enough to catch up to the multiple inputs version

//gather kernels, input and output spaces in this object from DemagKernelCollection.h
//the idea is to re-arrange the kernel multiplications to optimise kernel cache use
//Thus for each kernel type make a list of input and output spaces, so we can perform the multiplications by kernel type first
//SDemag will hold a vector of KerTypeCollection

class KerTypeCollection {

private:

	//the unique kernel
	std::shared_ptr<KerType> kernel;

	//input and output spaces vectors
	//Note, for each input there is exactly one output and conversely (otherwise the kernel would differ)
	std::vector<VEC<ReIm3>*> In, Out;

	//do we need to correct for an inverse z shift? in the In to Out multiplication?
	//Note, the first entry in the In, Out vectors cannot be z shifted, since the kernels are gathered in order
	std::vector<bool> inverse_shifted;

private:

	//2D

	void KernelMultiplication_2D_Self(void);

	void KernelMultiplication_2D_zShifted(void);

	void KernelMultiplication_2D_Regular(void);

	//3D

	void KernelMultiplication_3D_Self(void);

	void KernelMultiplication_3D_zShifted(void);

	void KernelMultiplication_3D_Regular(void);

public:

	KerTypeCollection(std::shared_ptr<KerType> kernel_, VEC<ReIm3>* pIn, VEC<ReIm3>* pOut, bool inverse_shifted_)
	{
		kernel = kernel_;

		In.push_back(pIn);
		Out.push_back(pOut);
		inverse_shifted.push_back(inverse_shifted_);
	}

	//---------------------------------------------------------------------

	//if the kernel matches the one already stored here then add new entry in In, Out, and inverse_shifted vectors; return true
	//if it doesn't match return false
	bool add_entry_if_kernel_matches(std::shared_ptr<KerType> kernel_, VEC<ReIm3>* pIn, VEC<ReIm3>* pOut, bool inverse_shifted_)
	{
		if (kernel == kernel_) {

			In.push_back(pIn);
			Out.push_back(pOut);
			inverse_shifted.push_back(inverse_shifted_);

			return true;
		}
		else return false;
	}

	//---------------------------------------------------------------------

	int get_size(void) { return In.size(); }

	//---------------------------------------------------------------------

	void Kernel_Multiplication_2D(void);

	void Kernel_Multiplication_3D(void);
};

#endif
