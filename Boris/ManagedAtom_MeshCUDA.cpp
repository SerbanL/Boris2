#include "stdafx.h"
#include "ManagedAtom_MeshCUDA.h"
#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "SuperMesh.h"

BError ManagedAtom_MeshCUDA::set_pointers(Atom_MeshCUDA* paMeshCUDA)
{
	BError error(__FUNCTION__);

	//Material Parameters

	//-----------SIMPLE CUBIC

	if (set_gpu_value(pgrel, paMeshCUDA->grel.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(palpha, paMeshCUDA->alpha.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pmu_s, paMeshCUDA->mu_s.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pJ, paMeshCUDA->J.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pD, paMeshCUDA->D.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pD_dir, paMeshCUDA->D_dir.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pJs, paMeshCUDA->Js.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pK1, paMeshCUDA->K1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pK2, paMeshCUDA->K2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pK3, paMeshCUDA->K3.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pmcanis_ea1, paMeshCUDA->mcanis_ea1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pmcanis_ea2, paMeshCUDA->mcanis_ea2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pmcanis_ea3, paMeshCUDA->mcanis_ea3.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pKt, paMeshCUDA->Kt.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	//-----------Others

	if (set_gpu_value(pNxy, paMeshCUDA->Nxy.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pcHA, paMeshCUDA->cHA.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pcHmo, paMeshCUDA->cHmo.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pelecCond, paMeshCUDA->elecCond.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pbase_temperature, paMeshCUDA->base_temperature.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pthermCond, paMeshCUDA->thermCond.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pdensity, paMeshCUDA->density.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pshc, paMeshCUDA->shc.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pshc_e, paMeshCUDA->shc_e.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pG_e, paMeshCUDA->G_e.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pcT, paMeshCUDA->cT.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pQ, paMeshCUDA->Q.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	//Mesh quantities

	if (set_gpu_value(pM1, paMeshCUDA->M1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pHeff1, paMeshCUDA->Heff1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pV, paMeshCUDA->V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pelC, paMeshCUDA->elC.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pE, paMeshCUDA->E.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pS, paMeshCUDA->S.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pTemp, paMeshCUDA->Temp.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pTemp_l, paMeshCUDA->Temp_l.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pu_disp, paMeshCUDA->u_disp.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pstrain_diag, paMeshCUDA->strain_diag.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pstrain_odiag, paMeshCUDA->strain_odiag.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	//Managed Atom_DiffEq_CommonCUDA pointer so all common diffeq data can be accessed in device code
	if (set_gpu_value(pcuaDiffEq, paMeshCUDA->Get_ManagedAtom_DiffEq_CommonCUDA().get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif