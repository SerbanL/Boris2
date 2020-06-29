#include "stdafx.h"
#include "ManagedDiffEq_CommonCUDA.h"
#include "DiffEq_CommonCUDA.h"

#if COMPILECUDA == 1

BError ManagedDiffEq_CommonCUDA::set_pointers(ODECommonCUDA* pDiffEqCUDA)
{
	BError error(__FUNCTION__);

	//Pointers to data in ODECommonCUDA

	if (set_gpu_value(ptime, pDiffEqCUDA->ptime->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pstagetime, pDiffEqCUDA->pstagetime->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdT, pDiffEqCUDA->pdT->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdT_last, pDiffEqCUDA->pdT_last->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pmxh, pDiffEqCUDA->pmxh->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pmxh_av, pDiffEqCUDA->pmxh_av->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pavpoints, pDiffEqCUDA->pavpoints->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdmdt, pDiffEqCUDA->pdmdt->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdmdt_av, pDiffEqCUDA->pdmdt_av->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pavpoints2, pDiffEqCUDA->pavpoints2->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(plte, pDiffEqCUDA->plte->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(prenormalize, pDiffEqCUDA->prenormalize->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psolve_spin_current, pDiffEqCUDA->psolve_spin_current->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psetODE, pDiffEqCUDA->psetODE->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(palternator, pDiffEqCUDA->palternator->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdelta_M_sq, pDiffEqCUDA->pdelta_M_sq->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelta_G_sq, pDiffEqCUDA->pdelta_G_sq->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelta_M_dot_delta_G, pDiffEqCUDA->pdelta_M_dot_delta_G->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdelta_M2_sq, pDiffEqCUDA->pdelta_M2_sq->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelta_G2_sq, pDiffEqCUDA->pdelta_G2_sq->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelta_M2_dot_delta_G2, pDiffEqCUDA->pdelta_M2_dot_delta_G2->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

