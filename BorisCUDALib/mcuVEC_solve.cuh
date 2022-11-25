#pragma once

#include "mcuVEC.h"

//--------------------------------------------LAPLACE / POISSON EQUATION : mcuVEC_solve.h

//NOTE : these are placed in mcuVEC_solve.cuh, since they can only be called from cu files, not from h or cpp.
//Reason for this, IteratePoisson_SOR_red and IteratePoisson_SOR_black, etc., are not declared with explicit template parameters (it would be very restrictive/cumbersome to specify all possible Class_Poisson_RHS)
//Then, without including cuVEC_VC_solve.cuh a linker error would be issued. The latter can only be included in cu files since it contains cuda global kernels.

//POISSON with homogeneous Neumann boundaries

template <typename VType, typename MType>
template <typename Class_Poisson_RHS, typename Policy_Class_Poisson_RHS, typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::IteratePoisson_SOR(mcu_obj<Class_Poisson_RHS, Policy_Class_Poisson_RHS>& obj, mcu_val<cuBReal>& damping, mcu_val<cuBReal>& max_error, mcu_val<cuBReal>& max_val)
{
	//for multiple devices we need halos before each red / black passes.
	exchange_halos();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->IteratePoisson_SOR_red(pn_d[mGPU].dim(), obj.get_deviceobject(mGPU), damping(mGPU));
	}

	exchange_halos();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->IteratePoisson_SOR_black(pn_d[mGPU].dim(), obj.get_deviceobject(mGPU), damping(mGPU), max_error(mGPU), max_val(mGPU));
	}
}

//POISSON with non-homogeneous Neumann boundaries

template <typename VType, typename MType>
template <typename Class_Poisson_NNeu, typename Policy_Class_Poisson_NNeu, typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::IteratePoisson_NNeu_SOR(mcu_obj<Class_Poisson_NNeu, Policy_Class_Poisson_NNeu>& obj, mcu_val<cuBReal>& damping, mcu_val<cuBReal>& max_error, mcu_val<cuBReal>& max_val)
{
	//for multiple devices we need halos before each red / black passes.
	exchange_halos();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->IteratePoisson_NNeu_SOR_red(pn_d[mGPU].dim(), obj.get_deviceobject(mGPU), damping(mGPU));
	}

	exchange_halos();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->IteratePoisson_NNeu_SOR_black(pn_d[mGPU].dim(), obj.get_deviceobject(mGPU), damping(mGPU), max_error(mGPU), max_val(mGPU));
	}
}