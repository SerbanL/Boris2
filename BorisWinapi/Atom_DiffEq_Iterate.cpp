#include "stdafx.h"
#include "Atom_DiffEq.h"

//---------------------------------------- ITERATE METHODS

//restore ODE solvers so iteration can be redone (call if AdvanceIteration has failed)
void Atom_ODECommon::Restore(void)
{
	for (int idx = 0; idx < (int)pODE.size(); idx++) {

		pODE[idx]->RestoreMoments();
	}
}