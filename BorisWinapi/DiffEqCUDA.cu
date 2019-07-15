#include "DiffEqCUDA.h"
#include "DiffEq_Defs.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

__global__ void set_ODEFunc_pointers(ManagedDiffEqCUDA& cuDiffEq)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		switch (*cuDiffEq.psetODE) {

		case ODE_LLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLG;
			break;

		case ODE_LLGSTATIC:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLGStatic;
			break;

		case ODE_LLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLGSTT;
			break;

		case ODE_LLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLB;
			break;

		case ODE_LLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLBSTT;
			break;

		case ODE_SLLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::SLLG;
			break;

		case ODE_SLLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::SLLGSTT;
			break;

		case ODE_SLLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::SLLB;
			break;

		case ODE_SLLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::SLLBSTT;
			break;

		case ODE_LLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLG;
			break;

		case ODE_LLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLB;
			break;

		case ODE_SLLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::SLLG;
			break;

		case ODE_SLLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::SLLB;
			break;

		default:
			cuDiffEq.pODEFunc = &ManagedDiffEqCUDA::LLG;
			break;
		}
	}
}

void DifferentialEquationCUDA::SetODEMethodPointers(void)
{
	set_ODEFunc_pointers <<<1, 1>>> (cuDiffEq);
}

#endif