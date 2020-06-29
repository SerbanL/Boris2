#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "DiffEqFM_EquationsCUDA.h"
#include "DiffEqFM_SEquationsCUDA.h"
#include "DiffEq_Defs.h"

#include "BorisCUDALib.cuh"

__global__ void set_ODEFunc_pointers(ManagedDiffEqFMCUDA& cuDiffEq)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		switch (*cuDiffEq.psetODE) {

		case ODE_LLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLG;
			break;

		case ODE_LLGSTATIC:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLGStatic;
			break;

		case ODE_LLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLGSTT;
			break;

		case ODE_LLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLB;
			break;

		case ODE_LLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLBSTT;
			break;

		case ODE_SLLG:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::SLLG;
			break;

		case ODE_SLLGSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::SLLGSTT;
			break;

		case ODE_SLLB:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::SLLB;
			break;

		case ODE_SLLBSTT:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::SLLBSTT;
			break;

		case ODE_LLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLG;
			break;

		case ODE_LLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLB;
			break;

		case ODE_SLLGSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::SLLG;
			break;

		case ODE_SLLBSA:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::SLLB;
			break;

		default:
			cuDiffEq.pODEFunc = &ManagedDiffEqFMCUDA::LLG;
			break;
		}
	}
}

void DifferentialEquationFMCUDA::SetODEMethodPointers(void)
{
	set_ODEFunc_pointers <<<1, 1>>> (cuDiffEq);
}

#endif
#endif