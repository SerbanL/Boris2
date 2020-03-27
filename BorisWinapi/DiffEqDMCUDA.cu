#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "DiffEqDM_EquationsCUDA.h"
#include "DiffEqDM_SEquationsCUDA.h"
#include "DiffEq_Defs.h"

#include "BorisCUDALib.cuh"

__global__ void set_ODEFunc_pointers(ManagedDiffEqDMCUDA& cuDiffEq)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		switch (*cuDiffEq.psetODE) {

		case ODE_LLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLG;
			break;

		case ODE_LLGSTATIC:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLGStatic;
			break;

		case ODE_LLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLGSTT;
			break;

		case ODE_LLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLB;
			break;

		case ODE_LLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLBSTT;
			break;

		case ODE_SLLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::SLLG;
			break;

		case ODE_SLLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::SLLGSTT;
			break;

		case ODE_SLLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::SLLB;
			break;

		case ODE_SLLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::SLLBSTT;
			break;

		case ODE_LLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLG;
			break;

		case ODE_LLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLB;
			break;

		case ODE_SLLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::SLLG;
			break;

		case ODE_SLLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::SLLB;
			break;

		default:
			cuDiffEq.pODEFunc = &ManagedDiffEqDMCUDA::LLG;
			break;
		}
	}
}

void DifferentialEquationDMCUDA::SetODEMethodPointers(void)
{
	set_ODEFunc_pointers <<<1, 1>>> (cuDiffEq);
}

#endif
#endif