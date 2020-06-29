#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_DiffEqCubic_EquationsCUDA.h"
#include "Atom_DiffEqCubic_SEquationsCUDA.h"
#include "DiffEq_Defs.h"

#include "BorisCUDALib.cuh"

__global__ void set_ODEFunc_pointers(ManagedAtom_DiffEqCubicCUDA& cuaDiffEq)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		switch (*cuaDiffEq.psetODE) {

		case ODE_LLG:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::LLG;
			break;

		case ODE_LLGSTATIC:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::LLGStatic;
			break;

		case ODE_LLGSTT:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::LLGSTT;
			break;

		case ODE_SLLG:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::SLLG;
			break;

		case ODE_SLLGSTT:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::SLLGSTT;
			break;

		case ODE_LLGSA:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::LLG;
			break;

		case ODE_SLLGSA:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::SLLG;
			break;

		default:
			cuaDiffEq.pODEFunc = &ManagedAtom_DiffEqCubicCUDA::LLG;
			break;
		}
	}
}

void Atom_DifferentialEquationCubicCUDA::SetODEMethodPointers(void)
{
	set_ODEFunc_pointers <<<1, 1>>> (cuaDiffEq);
}

#endif
#endif