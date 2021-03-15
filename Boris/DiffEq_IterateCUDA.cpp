#include "stdafx.h"
#include "DiffEq.h"

#if COMPILECUDA == 1

//---------------------------------------- ITERATE METHODS : CUDA

//restore ODE solvers so iteration can be redone (call if AdvanceIteration has failed)
void ODECommon::RestoreCUDA(void)
{
	for (int idx = 0; idx < (int)pODE.size(); idx++) {

		pODE[idx]->pmeshODECUDA->RestoreMagnetization();
	}
}

#endif