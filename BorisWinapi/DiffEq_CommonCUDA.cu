#include "DiffEq_CommonCUDA.h"

#if COMPILECUDA == 1

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

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
	Zero_SD_Solver_BB_Values_kernel <<< 1, CUDATHREADS >>> (*pdelta_M_sq, *pdelta_G_sq, *pdelta_M_dot_delta_G, *pdelta_M2_sq, *pdelta_G2_sq, *pdelta_M2_dot_delta_G2);
}

#endif