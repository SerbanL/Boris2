#include "DiffEqAFMCUDA.h"
#include "DiffEqAFM_EquationsCUDA.h"
#include "DiffEqAFM_SEquationsCUDA.h"
#include "DiffEq_Defs.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

__global__ void set_ODEFunc_pointers(ManagedDiffEqAFMCUDA& cuDiffEq)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		switch (*cuDiffEq.psetODE) {

		case ODE_LLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLG;
			break;

		case ODE_LLGSTATIC:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLGStatic;
			break;

		case ODE_LLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLGSTT;
			break;

		case ODE_LLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLB;
			break;

		case ODE_LLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLBSTT;
			break;

		case ODE_SLLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::SLLG;
			break;

		case ODE_SLLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::SLLGSTT;
			break;

		case ODE_SLLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::SLLB;
			break;

		case ODE_SLLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::SLLBSTT;
			break;

		case ODE_LLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLG;
			break;

		case ODE_LLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLB;
			break;

		case ODE_SLLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::SLLG;
			break;

		case ODE_SLLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::SLLB;
			break;

		default:
			cuDiffEq.pODEFunc = &ManagedDiffEqAFMCUDA::LLG;
			break;
		}
	}
}

void DifferentialEquationAFMCUDA::SetODEMethodPointers(void)
{
	set_ODEFunc_pointers <<<1, 1>>> (cuDiffEq);
}

#endif