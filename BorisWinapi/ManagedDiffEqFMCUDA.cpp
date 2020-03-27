#include "stdafx.h"
#include "ManagedDiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "DiffEqFMCUDA.h"

BError ManagedDiffEqFMCUDA::set_pointers(DifferentialEquationFMCUDA* pDiffEqCUDA)
{
	BError error(__FUNCTION__);
	
	//Pointers to data in ODECommonCUDA

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

	//Pointers to data in DifferentialEquationCUDA

	if (set_gpu_value(psM1, pDiffEqCUDA->sM1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psEval0, pDiffEqCUDA->sEval0.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval1, pDiffEqCUDA->sEval1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval2, pDiffEqCUDA->sEval2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval3, pDiffEqCUDA->sEval3.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval4, pDiffEqCUDA->sEval4.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psEval5, pDiffEqCUDA->sEval5.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pH_Thermal, pDiffEqCUDA->H_Thermal.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pTorque_Thermal, pDiffEqCUDA->Torque_Thermal.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	//Managed cuda mesh pointer so all mesh data can be accessed in device code

	if (set_gpu_value(pcuMesh, pDiffEqCUDA->pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	return error;
}

#endif
#endif

