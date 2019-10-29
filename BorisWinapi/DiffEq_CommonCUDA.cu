#include "DiffEq_CommonCUDA.h"

#if COMPILECUDA == 1

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- AUXILIARY

__global__ void Zerovalues_kernel(cuBReal& mxh, cuReal3& mxh_av, size_t& avpoints, cuBReal& dmdt, cuReal3& dmdt_av, size_t& avpoints2, cuBReal& lte)
{
	if (threadIdx.x == 0) mxh = 0.0;
	else if (threadIdx.x == 1) mxh_av = cuReal3(0.0);
	else if (threadIdx.x == 2) avpoints = 0;
	else if (threadIdx.x == 3) dmdt = 0.0;
	else if (threadIdx.x == 4) dmdt_av = cuReal3(0.0);
	else if (threadIdx.x == 5) avpoints2 = 0;
	else if (threadIdx.x == 6) lte = 0.0;
}

void ODECommonCUDA::Zero_reduction_values(void)
{
	Zerovalues_kernel << < 1, CUDATHREADS >> > (*pmxh, *pmxh_av, *pavpoints, *pdmdt, *pdmdt_av, *pavpoints2, *plte);
}

//-----------------------------------------

__global__ void Zeromxh_kernel(cuBReal& mxh, cuReal3& mxh_av, size_t& avpoints, cuBReal& lte)
{
	if (threadIdx.x == 0) mxh = 0.0;
	else if (threadIdx.x == 1) mxh_av = cuReal3(0.0);
	else if (threadIdx.x == 2) avpoints = 0;
	else if (threadIdx.x == 3) lte = 0.0;
}

void ODECommonCUDA::Zero_mxh_lte_values(void)
{
	Zeromxh_kernel << < 1, CUDATHREADS >> > (*pmxh, *pmxh_av, *pavpoints, *plte);
}

//-----------------------------------------

__global__ void Zerodmdt_kernel(cuBReal& dmdt, cuReal3& dmdt_av, size_t& avpoints2, cuBReal& lte)
{
	if (threadIdx.x == 0) dmdt = 0.0;
	else if (threadIdx.x == 1) dmdt_av = cuReal3(0.0);
	else if (threadIdx.x == 2) avpoints2 = 0;
	else if (threadIdx.x == 3) lte = 0.0;
}

void ODECommonCUDA::Zero_dmdt_lte_values(void)
{
	Zerodmdt_kernel << < 1, CUDATHREADS >> > (*pdmdt, *pdmdt_av, *pavpoints2, *plte);
}

//-----------------------------------------

__global__ void Zerolte_kernel(cuBReal& lte)
{
	if (threadIdx.x == 0) lte = 0.0;
}

void ODECommonCUDA::Zero_lte_value(void)
{
	Zerolte_kernel << < 1, CUDATHREADS >> > (*plte);
}

//-----------------------------------------

__global__ void mxhav_to_mxh_kernel(cuBReal& mxh, cuReal3& mxh_av, size_t& avpoints)
{
	if (threadIdx.x == 0) {

		if (avpoints) {

			mxh = cu_GetMagnitude(mxh_av) / avpoints;
		}
		else {

			mxh = 0.0;
		}
	}
}

void ODECommonCUDA::mxhav_to_mxh(void)
{
	mxhav_to_mxh_kernel << < 1, CUDATHREADS >> > (*pmxh, *pmxh_av, *pavpoints);
}

//-----------------------------------------

__global__ void dmdtav_to_dmdt_kernel(cuBReal& dmdt, cuReal3& dmdt_av, size_t& avpoints2)
{
	if (threadIdx.x == 0) {

		if (avpoints2) {

			dmdt = cu_GetMagnitude(dmdt_av) / avpoints2;
		}
		else {

			dmdt = 0.0;
		}
	}
}

void ODECommonCUDA::dmdtav_to_dmdt(void)
{
	dmdtav_to_dmdt_kernel << < 1, CUDATHREADS >> > (*pdmdt, *pdmdt_av, *pavpoints2);
}

//-----------------------------------------

__global__ void Zero_SD_Solver_BB_Values_kernel(cuBReal& delta_M_sq, cuBReal& delta_G_sq, cuBReal& delta_M_dot_delta_G, cuBReal& delta_M2_sq, cuBReal& delta_G2_sq, cuBReal& delta_M2_dot_delta_G2)
{
	if (threadIdx.x == 0) delta_M_sq = 0.0;
	else if (threadIdx.x == 1) delta_G_sq = 0.0;
	else if (threadIdx.x == 2) delta_M_dot_delta_G = 0.0;
	else if (threadIdx.x == 3) delta_M2_sq = 0.0;
	else if (threadIdx.x == 4) delta_G2_sq = 0.0;
	else if (threadIdx.x == 5) delta_M2_dot_delta_G2 = 0.0;
}

void ODECommonCUDA::Zero_SD_Solver_BB_Values(void)
{
	Zero_SD_Solver_BB_Values_kernel << < 1, CUDATHREADS >> > (*pdelta_M_sq, *pdelta_G_sq, *pdelta_M_dot_delta_G, *pdelta_M2_sq, *pdelta_G2_sq, *pdelta_M2_dot_delta_G2);
}

#endif