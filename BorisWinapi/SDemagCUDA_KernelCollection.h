#pragma once

//TESTING ONLY - SLOWER THAN MULTIPLE INPUTS VERSION SO NOT IN CURRENT USE

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "BorisCUDALib.h"
#include "DemagKernelCollectionCUDA.h"

class KerTypeCollectionCUDA {

private:

	//the unique kernel
	std::shared_ptr<cu_obj<cuKerType>> kernel;

	//flags in cpu memory so we can quickly launch the required __global__ from Kernel_Multiplication_2D / Kernel_Multiplication_3D
	bool internal_demag, zshifted;

	//do we need to correct for an inverse z shift? in the In to Out multiplication?
	//Note, the first entry in the In, Out vectors cannot be z shifted, since the kernels are gathered in order
	cu_arr<bool> inverse_shifted;

	//input and output spaces vectors
	//Note, for each input there is exactly one output and conversely (otherwise the kernel would differ)
	cu_arr<cuComplex*> InCol_x, InCol_y, InCol_z, OutCol_x, OutCol_y, OutCol_z;

	//number of arrays in In, Out
	cu_obj<int> size;

	//cpu memory version of size
	int size_cpu;

	//need this in cpu memory
	SZ3 N;

public:

	KerTypeCollectionCUDA(
		std::shared_ptr<cu_obj<cuKerType>> kernel_, 
		cu_arr<cuComplex>* pIn_x, cu_arr<cuComplex>* pIn_y, cu_arr<cuComplex>* pIn_z,
		cu_arr<cuComplex>* pOut_x, cu_arr<cuComplex>* pOut_y, cu_arr<cuComplex>* pOut_z,
		bool inverse_shifted_)
	{
		kernel = kernel_;

		internal_demag = (*kernel)()->GetFlag_InternalDemag();
		zshifted = (*kernel)()->GetFlag_zShifted();
		N = (*kernel)()->Get_N();

		InCol_x.push_back(pIn_x->get_managed_array());
		InCol_y.push_back(pIn_y->get_managed_array());
		InCol_z.push_back(pIn_z->get_managed_array());

		OutCol_x.push_back(pOut_x->get_managed_array());
		OutCol_y.push_back(pOut_y->get_managed_array());
		OutCol_z.push_back(pOut_z->get_managed_array());

		size_cpu = 1;
		size.from_cpu(size_cpu);

		cu_obj<bool> gpu_value;
		gpu_value.from_cpu(inverse_shifted_);
		inverse_shifted.push_back((bool*&)gpu_value.get_managed_object());
	}

	//if the kernel matches the one already stored here then add new entry in In, Out, and inverse_shifted vectors; return true
	//if it doesn't match return false
	bool add_entry_if_kernel_matches(
		std::shared_ptr<cu_obj<cuKerType>> kernel_,
		cu_arr<cuComplex>* pIn_x, cu_arr<cuComplex>* pIn_y, cu_arr<cuComplex>* pIn_z,
		cu_arr<cuComplex>* pOut_x, cu_arr<cuComplex>* pOut_y, cu_arr<cuComplex>* pOut_z,
		bool inverse_shifted_)
	{
		if (kernel == kernel_) {

			InCol_x.push_back(pIn_x->get_managed_array());
			InCol_y.push_back(pIn_y->get_managed_array());
			InCol_z.push_back(pIn_z->get_managed_array());

			OutCol_x.push_back(pOut_x->get_managed_array());
			OutCol_y.push_back(pOut_y->get_managed_array());
			OutCol_z.push_back(pOut_z->get_managed_array());

			size_cpu++;
			size.from_cpu(size_cpu);

			cu_obj<bool> gpu_value;
			gpu_value.from_cpu(inverse_shifted_);
			inverse_shifted.push_back((bool*&)gpu_value.get_managed_object());

			return true;
		}
		else return false;
	}

	//---------------------------------------------------------------------

	void Kernel_Multiplication_2D(bool transpose_xy);

	void Kernel_Multiplication_3D(void);
};

#endif

#endif
