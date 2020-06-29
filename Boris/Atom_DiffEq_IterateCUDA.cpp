#include "stdafx.h"
#include "Atom_DiffEq.h"

#if COMPILECUDA == 1

//---------------------------------------- ITERATE METHODS : CUDA

//restore ODE solvers so iteration can be redone (call if AdvanceIteration has failed)
void Atom_ODECommon::RestoreCUDA(void)
{
	for (int idx = 0; idx < (int)pODE.size(); idx++) {

		pODE[idx]->pameshODECUDA->RestoreMoments();
	}
}

#endif