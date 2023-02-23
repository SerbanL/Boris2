#include "stdafx.h"
#include "ManagedAtom_DiffEq_CommonCUDA.h"

#if COMPILECUDA == 1

#include "Atom_DiffEq_CommonCUDA.h"

BError ManagedAtom_DiffEq_CommonCUDA::set_pointers(Atom_ODECommonCUDA* paDiffEqCUDA)
{
	BError error(__FUNCTION__);

	//Pointers to data in Atom_ODECommonCUDA

	if (set_gpu_value(ptime, paDiffEqCUDA->ptime->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pstagetime, paDiffEqCUDA->pstagetime->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdT, paDiffEqCUDA->pdT->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdT_last, paDiffEqCUDA->pdT_last->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pmxh, paDiffEqCUDA->pmxh->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pmxh_av, paDiffEqCUDA->pmxh_av->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pavpoints, paDiffEqCUDA->pavpoints->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdmdt, paDiffEqCUDA->pdmdt->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdmdt_av, paDiffEqCUDA->pdmdt_av->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pavpoints2, paDiffEqCUDA->pavpoints2->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(plte, paDiffEqCUDA->plte->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(prenormalize, paDiffEqCUDA->prenormalize->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psolve_spin_current, paDiffEqCUDA->psolve_spin_current->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psetODE, paDiffEqCUDA->psetODE->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(palternator, paDiffEqCUDA->palternator->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdelta_M_sq, paDiffEqCUDA->pdelta_M_sq->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelta_G_sq, paDiffEqCUDA->pdelta_G_sq->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelta_M_dot_delta_G, paDiffEqCUDA->pdelta_M_dot_delta_G->get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

