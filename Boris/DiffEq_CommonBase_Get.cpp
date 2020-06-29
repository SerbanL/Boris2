#include "stdafx.h"
#include "DiffEq_CommonBase.h"

#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

double ODECommon_Base::Get_mxh(void)
{
	calculate_mxh = true;

#if COMPILECUDA == 1
	//mxh same for gpu memory also across micromagnetic and atomistic diffeq, so only need to use one pODECUDA
	if (podeSolver->pODECUDA) return podeSolver->pODECUDA->Get_mxh();
#endif

	return mxh;
}

double ODECommon_Base::Get_dmdt(void)
{
	calculate_dmdt = true;

#if COMPILECUDA == 1
	//mxh same for gpu memory also across micromagnetic and atomistic diffeq, so only need to use one pODECUDA
	if (podeSolver->pODECUDA) return podeSolver->pODECUDA->Get_dmdt();
#endif

	return dmdt;
}
