#include "stdafx.h"
#include "ManagedAtom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_DiffEqCubicCUDA.h"

BError ManagedAtom_DiffEqCubicCUDA::set_pointers(Atom_DifferentialEquationCubicCUDA* paDiffEqCUDA)
{
	BError error(__FUNCTION__);
	
	//Pointers to data in Atom_ODECommonCUDA

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

	//Pointers to data in DifferentialEquationCUDA

	if (set_gpu_value(psM1, paDiffEqCUDA->sM1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psEval0, paDiffEqCUDA->sEval0.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval1, paDiffEqCUDA->sEval1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval2, paDiffEqCUDA->sEval2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval3, paDiffEqCUDA->sEval3.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval4, paDiffEqCUDA->sEval4.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval5, paDiffEqCUDA->sEval5.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pH_Thermal, paDiffEqCUDA->H_Thermal.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	//Managed cuda mesh pointer so all mesh data can be accessed in device code

	if (set_gpu_value(pcuaMesh, paDiffEqCUDA->paMeshCUDA->cuaMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	return error;
}

#endif
#endif

