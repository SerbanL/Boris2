#include "HeatBaseCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.cuh"

//-------------------Setters

BError HeatBaseCUDA::SetQEquation(const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	BError error(CLASS_STR(HeatBaseCUDA));

	if (!Q_equation.make_scalar(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif

#endif