#include "stdafx.h"
#include "ManagedMeshCUDA.h"
#include "MeshCUDA.h"

#if COMPILECUDA == 1

#include "SuperMesh.h"

BError ManagedMeshCUDA::set_pointers(MeshCUDA* pMeshCUDA)
{
	BError error(__FUNCTION__);
	
	//Material Parameters
	
	if (set_gpu_value(pgrel, pMeshCUDA->grel.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pgrel_AFM, pMeshCUDA->grel_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(palpha, pMeshCUDA->alpha.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(palpha_AFM, pMeshCUDA->alpha_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pMs, pMeshCUDA->Ms.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pMs_AFM, pMeshCUDA->Ms_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pNxy, pMeshCUDA->Nxy.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pA, pMeshCUDA->A.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pA_AFM, pMeshCUDA->A_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pAh, pMeshCUDA->Ah.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pAnh, pMeshCUDA->Anh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pD, pMeshCUDA->D.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pD_AFM, pMeshCUDA->D_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(ptau_ii, pMeshCUDA->tau_ii.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(ptau_ij, pMeshCUDA->tau_ij.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pJ1, pMeshCUDA->J1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pJ2, pMeshCUDA->J2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pneta_dia, pMeshCUDA->neta_dia.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pK1, pMeshCUDA->K1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pK2, pMeshCUDA->K2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pmcanis_ea1, pMeshCUDA->mcanis_ea1.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pmcanis_ea2, pMeshCUDA->mcanis_ea2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pK1_AFM, pMeshCUDA->K1_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pK2_AFM, pMeshCUDA->K2_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(psusrel, pMeshCUDA->susrel.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psusrel_AFM, pMeshCUDA->susrel_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(psusprel, pMeshCUDA->susprel.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pcHA, pMeshCUDA->cHA.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pcHmo, pMeshCUDA->cHmo.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pelecCond, pMeshCUDA->elecCond.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pamrPercentage, pMeshCUDA->amrPercentage.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pP, pMeshCUDA->P.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pbeta, pMeshCUDA->beta.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pDe, pMeshCUDA->De.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pn_density, pMeshCUDA->n_density.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pbetaD, pMeshCUDA->betaD.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pSHA, pMeshCUDA->SHA.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(piSHA, pMeshCUDA->iSHA.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pflSOT, pMeshCUDA->flSOT.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pl_sf, pMeshCUDA->l_sf.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pl_ex, pMeshCUDA->l_ex.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pl_ph, pMeshCUDA->l_ph.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pGi, pMeshCUDA->Gi.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pGmix, pMeshCUDA->Gmix.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pts_eff, pMeshCUDA->ts_eff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(ptsi_eff, pMeshCUDA->tsi_eff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(ppump_eff, pMeshCUDA->pump_eff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pcpump_eff, pMeshCUDA->cpump_eff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pthe_eff, pMeshCUDA->the_eff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pbase_temperature, pMeshCUDA->base_temperature.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pT_Curie, pMeshCUDA->T_Curie.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(patomic_moment, pMeshCUDA->atomic_moment.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(patomic_moment_AFM, pMeshCUDA->atomic_moment_AFM.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pthermCond, pMeshCUDA->thermCond.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdensity, pMeshCUDA->density.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pMEc, pMeshCUDA->MEc.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pYm, pMeshCUDA->Ym.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pPr, pMeshCUDA->Pr.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pshc, pMeshCUDA->shc.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pshc_e, pMeshCUDA->shc_e.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pG_e, pMeshCUDA->G_e.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pcT, pMeshCUDA->cT.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pQ, pMeshCUDA->Q.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	//Mesh quantities
	
	if (set_gpu_value(pM, pMeshCUDA->M.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pM2, pMeshCUDA->M2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pHeff, pMeshCUDA->Heff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pHeff2, pMeshCUDA->Heff2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pV, pMeshCUDA->V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pelC, pMeshCUDA->elC.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pE, pMeshCUDA->E.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pS, pMeshCUDA->S.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	if (set_gpu_value(pTemp, pMeshCUDA->Temp.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pTemp_l, pMeshCUDA->Temp_l.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pu_disp, pMeshCUDA->u_disp.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pstrain_diag, pMeshCUDA->strain_diag.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pstrain_odiag, pMeshCUDA->strain_odiag.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	//Managed DiffEq_CommonCUDA pointer so all common diffeq data can be accessed in device code
	if (set_gpu_value(pcuDiffEq, pMeshCUDA->Get_ManagedDiffEq_CommonCUDA().get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif