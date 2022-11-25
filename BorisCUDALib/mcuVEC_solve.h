#pragma once

#include "mcuVEC.h"

//--------------------------------------------LAPLACE / POISSON EQUATION : mcuVEC_solve.h

//LAPLACE

//NOTE : this is placed in mcuVEC_solve.h, since it can be called from h or cpp files (could put it in mcuVEC_solve.cuh, but keeping to convention)

//Take one SOR iteration for Laplace equation on this VEC. Return error (maximum change in quantity from one iteration to the next) by reference.
//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
//Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::IterateLaplace_SOR(mcu_val<cuBReal>& damping, mcu_val<cuBReal>& max_error, mcu_val<cuBReal>& max_val)
{
	//for multiple devices we need halos before each red / black passes.
	exchange_halos();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->IterateLaplace_SOR_red(pn_d[mGPU].dim(), damping(mGPU));
	}

	exchange_halos();
	
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->IterateLaplace_SOR_black(pn_d[mGPU].dim(), damping(mGPU), max_error(mGPU), max_val(mGPU));
	}
}