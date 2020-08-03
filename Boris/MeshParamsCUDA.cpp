#include "stdafx.h"
#include "MeshParamsCUDA.h"
#include "MeshParams.h"

#if COMPILECUDA == 1

MeshParamsCUDA::MeshParamsCUDA(MeshParams *pmeshParams)
{
	this->pmeshParams = pmeshParams;
	
	grel()->set_from_cpu(pmeshParams->grel);
	pmeshParams->grel.set_p_cu_obj_mpcuda(&grel);
	grel_AFM()->set_from_cpu(pmeshParams->grel_AFM);
	pmeshParams->grel_AFM.set_p_cu_obj_mpcuda(&grel_AFM);

	alpha()->set_from_cpu(pmeshParams->alpha);
	pmeshParams->alpha.set_p_cu_obj_mpcuda(&alpha);
	alpha_AFM()->set_from_cpu(pmeshParams->alpha_AFM);
	pmeshParams->alpha_AFM.set_p_cu_obj_mpcuda(&alpha_AFM);

	Ms()->set_from_cpu(pmeshParams->Ms);
	pmeshParams->Ms.set_p_cu_obj_mpcuda(&Ms);
	Ms_AFM()->set_from_cpu(pmeshParams->Ms_AFM);
	pmeshParams->Ms_AFM.set_p_cu_obj_mpcuda(&Ms_AFM);

	Nxy()->set_from_cpu(pmeshParams->Nxy);
	pmeshParams->Nxy.set_p_cu_obj_mpcuda(&Nxy);

	A()->set_from_cpu(pmeshParams->A);
	pmeshParams->A.set_p_cu_obj_mpcuda(&A);
	A_AFM()->set_from_cpu(pmeshParams->A_AFM);
	pmeshParams->A_AFM.set_p_cu_obj_mpcuda(&A_AFM);
	
	Ah()->set_from_cpu(pmeshParams->Ah);
	pmeshParams->Ah.set_p_cu_obj_mpcuda(&Ah);
	Anh()->set_from_cpu(pmeshParams->Anh);
	pmeshParams->Anh.set_p_cu_obj_mpcuda(&Anh);

	D()->set_from_cpu(pmeshParams->D);
	pmeshParams->D.set_p_cu_obj_mpcuda(&D);
	D_AFM()->set_from_cpu(pmeshParams->D_AFM);
	pmeshParams->D_AFM.set_p_cu_obj_mpcuda(&D_AFM);

	tau_ii()->set_from_cpu(pmeshParams->tau_ii);
	pmeshParams->tau_ii.set_p_cu_obj_mpcuda(&tau_ii);
	tau_ij()->set_from_cpu(pmeshParams->tau_ij);
	pmeshParams->tau_ij.set_p_cu_obj_mpcuda(&tau_ij);

	J1()->set_from_cpu(pmeshParams->J1);
	pmeshParams->J1.set_p_cu_obj_mpcuda(&J1);
	J2()->set_from_cpu(pmeshParams->J2);
	pmeshParams->J2.set_p_cu_obj_mpcuda(&J2);
	neta_dia()->set_from_cpu(pmeshParams->neta_dia);
	pmeshParams->neta_dia.set_p_cu_obj_mpcuda(&neta_dia);

	K1()->set_from_cpu(pmeshParams->K1);
	pmeshParams->K1.set_p_cu_obj_mpcuda(&K1);
	K2()->set_from_cpu(pmeshParams->K2);
	pmeshParams->K2.set_p_cu_obj_mpcuda(&K2);
	mcanis_ea1()->set_from_cpu(pmeshParams->mcanis_ea1);
	pmeshParams->mcanis_ea1.set_p_cu_obj_mpcuda(&mcanis_ea1);
	mcanis_ea2()->set_from_cpu(pmeshParams->mcanis_ea2);
	pmeshParams->mcanis_ea2.set_p_cu_obj_mpcuda(&mcanis_ea2);

	K1_AFM()->set_from_cpu(pmeshParams->K1_AFM);
	pmeshParams->K1_AFM.set_p_cu_obj_mpcuda(&K1_AFM);
	K2_AFM()->set_from_cpu(pmeshParams->K2_AFM);
	pmeshParams->K2_AFM.set_p_cu_obj_mpcuda(&K2_AFM);

	susrel()->set_from_cpu(pmeshParams->susrel);
	pmeshParams->susrel.set_p_cu_obj_mpcuda(&susrel);
	susrel_AFM()->set_from_cpu(pmeshParams->susrel_AFM);
	pmeshParams->susrel_AFM.set_p_cu_obj_mpcuda(&susrel_AFM);
	susprel()->set_from_cpu(pmeshParams->susprel);
	pmeshParams->susprel.set_p_cu_obj_mpcuda(&susprel);

	cHA()->set_from_cpu(pmeshParams->cHA);
	pmeshParams->cHA.set_p_cu_obj_mpcuda(&cHA);

	cHmo()->set_from_cpu(pmeshParams->cHmo);
	pmeshParams->cHmo.set_p_cu_obj_mpcuda(&cHmo);

	elecCond()->set_from_cpu(pmeshParams->elecCond);
	pmeshParams->elecCond.set_p_cu_obj_mpcuda(&elecCond);
	amrPercentage()->set_from_cpu(pmeshParams->amrPercentage);
	pmeshParams->amrPercentage.set_p_cu_obj_mpcuda(&amrPercentage);

	P()->set_from_cpu(pmeshParams->P);
	pmeshParams->P.set_p_cu_obj_mpcuda(&P);
	beta()->set_from_cpu(pmeshParams->beta);
	pmeshParams->beta.set_p_cu_obj_mpcuda(&beta);

	De()->set_from_cpu(pmeshParams->De);
	pmeshParams->De.set_p_cu_obj_mpcuda(&De);

	n_density()->set_from_cpu(pmeshParams->n_density);
	pmeshParams->n_density.set_p_cu_obj_mpcuda(&n_density);
	
	betaD()->set_from_cpu(pmeshParams->betaD);
	pmeshParams->betaD.set_p_cu_obj_mpcuda(&betaD);
	
	SHA()->set_from_cpu(pmeshParams->SHA);
	pmeshParams->SHA.set_p_cu_obj_mpcuda(&SHA);
	iSHA()->set_from_cpu(pmeshParams->iSHA);
	pmeshParams->iSHA.set_p_cu_obj_mpcuda(&iSHA);
	flSOT()->set_from_cpu(pmeshParams->flSOT);
	pmeshParams->flSOT.set_p_cu_obj_mpcuda(&flSOT);

	STq()->set_from_cpu(pmeshParams->STq);
	pmeshParams->STq.set_p_cu_obj_mpcuda(&STq);
	STa()->set_from_cpu(pmeshParams->STa);
	pmeshParams->STa.set_p_cu_obj_mpcuda(&STa);
	STp()->set_from_cpu(pmeshParams->STp);
	pmeshParams->STp.set_p_cu_obj_mpcuda(&STp);
	
	l_sf()->set_from_cpu(pmeshParams->l_sf);
	pmeshParams->l_sf.set_p_cu_obj_mpcuda(&l_sf);
	l_ex()->set_from_cpu(pmeshParams->l_ex);
	pmeshParams->l_ex.set_p_cu_obj_mpcuda(&l_ex);
	l_ph()->set_from_cpu(pmeshParams->l_ph);
	pmeshParams->l_ph.set_p_cu_obj_mpcuda(&l_ph);
	
	Gi()->set_from_cpu(pmeshParams->Gi);
	pmeshParams->Gi.set_p_cu_obj_mpcuda(&Gi);
	Gmix()->set_from_cpu(pmeshParams->Gmix);
	pmeshParams->Gmix.set_p_cu_obj_mpcuda(&Gmix);
	
	ts_eff()->set_from_cpu(pmeshParams->ts_eff);
	pmeshParams->ts_eff.set_p_cu_obj_mpcuda(&ts_eff);
	tsi_eff()->set_from_cpu(pmeshParams->tsi_eff);
	pmeshParams->tsi_eff.set_p_cu_obj_mpcuda(&tsi_eff);
	
	pump_eff()->set_from_cpu(pmeshParams->pump_eff);
	pmeshParams->pump_eff.set_p_cu_obj_mpcuda(&pump_eff);

	cpump_eff()->set_from_cpu(pmeshParams->cpump_eff);
	pmeshParams->cpump_eff.set_p_cu_obj_mpcuda(&cpump_eff);

	the_eff()->set_from_cpu(pmeshParams->the_eff);
	pmeshParams->the_eff.set_p_cu_obj_mpcuda(&the_eff);

	base_temperature.from_cpu(pmeshParams->base_temperature);
	T_Curie.from_cpu(pmeshParams->T_Curie);

	atomic_moment()->set_from_cpu(pmeshParams->atomic_moment);
	pmeshParams->atomic_moment.set_p_cu_obj_mpcuda(&atomic_moment);
	atomic_moment_AFM()->set_from_cpu(pmeshParams->atomic_moment_AFM);
	pmeshParams->atomic_moment_AFM.set_p_cu_obj_mpcuda(&atomic_moment_AFM);

	thermCond()->set_from_cpu(pmeshParams->thermCond);
	pmeshParams->thermCond.set_p_cu_obj_mpcuda(&thermCond);
	density()->set_from_cpu(pmeshParams->density);
	pmeshParams->density.set_p_cu_obj_mpcuda(&density);
	
	shc()->set_from_cpu(pmeshParams->shc);
	pmeshParams->shc.set_p_cu_obj_mpcuda(&shc);
	shc_e()->set_from_cpu(pmeshParams->shc_e);
	pmeshParams->shc_e.set_p_cu_obj_mpcuda(&shc_e);
	G_e()->set_from_cpu(pmeshParams->G_e);
	pmeshParams->G_e.set_p_cu_obj_mpcuda(&G_e);

	cT()->set_from_cpu(pmeshParams->cT);
	pmeshParams->cT.set_p_cu_obj_mpcuda(&cT);
	Q()->set_from_cpu(pmeshParams->Q);
	pmeshParams->Q.set_p_cu_obj_mpcuda(&Q);

	MEc()->set_from_cpu(pmeshParams->MEc);
	pmeshParams->MEc.set_p_cu_obj_mpcuda(&MEc);
	Ym()->set_from_cpu(pmeshParams->Ym);
	pmeshParams->Ym.set_p_cu_obj_mpcuda(&Ym);
	Pr()->set_from_cpu(pmeshParams->Pr);
	pmeshParams->Pr.set_p_cu_obj_mpcuda(&Pr);

	//setup CUDA special functions to corresponding data held in the cpu objects
	set_special_functions_data();

	//make sure special functions are set by default for all material parameters text equations
	set_special_functions();
}

//set pre-calculated Funcs_Special objects in material parameters
void MeshParamsCUDA::set_special_functions(PARAM_ paramID)
{
	auto set_param_special_functions = [&](PARAM_ update_paramID) {

		auto code = [&](auto& MatP_object) -> void {

			MatP_object.set_t_scaling_special_functions_CUDA(&CurieWeiss_CUDA, &LongRelSus_CUDA, &CurieWeiss1_CUDA, &CurieWeiss2_CUDA, &LongRelSus1_CUDA, &LongRelSus2_CUDA, &Alpha1_CUDA, &Alpha2_CUDA);
		};

		pmeshParams->run_on_param<void>(update_paramID, code);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < pmeshParams->meshParams.size(); index++) {

			set_param_special_functions((PARAM_)pmeshParams->meshParams.get_ID_from_index(index));
		}
	}
	else set_param_special_functions(paramID);
}

void MeshParamsCUDA::set_special_functions_data(void)
{
	//setup CUDA special functions to corresponding data held in the cpu objects
	CurieWeiss_CUDA()->set_data(pmeshParams->pCurieWeiss->get_data(), pmeshParams->pCurieWeiss->get_start(), pmeshParams->pCurieWeiss->get_resolution());
	LongRelSus_CUDA()->set_data(pmeshParams->pLongRelSus->get_data(), pmeshParams->pLongRelSus->get_start(), pmeshParams->pLongRelSus->get_resolution());

	CurieWeiss1_CUDA()->set_data(pmeshParams->pCurieWeiss1->get_data(), pmeshParams->pCurieWeiss1->get_start(), pmeshParams->pCurieWeiss1->get_resolution());
	CurieWeiss2_CUDA()->set_data(pmeshParams->pCurieWeiss2->get_data(), pmeshParams->pCurieWeiss2->get_start(), pmeshParams->pCurieWeiss2->get_resolution());
	LongRelSus1_CUDA()->set_data(pmeshParams->pLongRelSus1->get_data(), pmeshParams->pLongRelSus1->get_start(), pmeshParams->pLongRelSus1->get_resolution());
	LongRelSus2_CUDA()->set_data(pmeshParams->pLongRelSus2->get_data(), pmeshParams->pLongRelSus2->get_start(), pmeshParams->pLongRelSus2->get_resolution());

	Alpha1_CUDA()->set_data(pmeshParams->pAlpha1->get_data(), pmeshParams->pAlpha1->get_start(), pmeshParams->pAlpha1->get_resolution());
	Alpha2_CUDA()->set_data(pmeshParams->pAlpha2->get_data(), pmeshParams->pAlpha2->get_start(), pmeshParams->pAlpha2->get_resolution());
}

MeshParamsCUDA::~MeshParamsCUDA()
{
	//fine to access data in MeshParams here : Mesh inherits from MeshParams, so in the destruction process Mesh gets destroyed first. 
	//Mesh destructor then calls for MeshCUDA implementation to be destroyed, then we get here since MeshCUDA inherits from MeshParamsCUDA. After this we return back to Mesh destructor to continue destruction down the list.
	
	pmeshParams->grel.null_p_cu_obj_mpcuda();
	pmeshParams->grel_AFM.null_p_cu_obj_mpcuda();

	pmeshParams->alpha.null_p_cu_obj_mpcuda();
	pmeshParams->alpha_AFM.null_p_cu_obj_mpcuda();

	pmeshParams->Ms.null_p_cu_obj_mpcuda();
	pmeshParams->Ms_AFM.null_p_cu_obj_mpcuda();

	pmeshParams->Nxy.null_p_cu_obj_mpcuda();

	pmeshParams->A.null_p_cu_obj_mpcuda();
	pmeshParams->A_AFM.null_p_cu_obj_mpcuda();
	
	pmeshParams->Ah.null_p_cu_obj_mpcuda();
	pmeshParams->Anh.null_p_cu_obj_mpcuda();

	pmeshParams->D.null_p_cu_obj_mpcuda();
	pmeshParams->D_AFM.null_p_cu_obj_mpcuda();

	pmeshParams->tau_ii.null_p_cu_obj_mpcuda();
	pmeshParams->tau_ij.null_p_cu_obj_mpcuda();

	pmeshParams->J1.null_p_cu_obj_mpcuda();
	pmeshParams->J2.null_p_cu_obj_mpcuda();
	pmeshParams->neta_dia.null_p_cu_obj_mpcuda();

	pmeshParams->K1.null_p_cu_obj_mpcuda();
	pmeshParams->K2.null_p_cu_obj_mpcuda();
	pmeshParams->mcanis_ea1.null_p_cu_obj_mpcuda();
	pmeshParams->mcanis_ea2.null_p_cu_obj_mpcuda();

	pmeshParams->K1_AFM.null_p_cu_obj_mpcuda();
	pmeshParams->K2_AFM.null_p_cu_obj_mpcuda();

	pmeshParams->susrel.null_p_cu_obj_mpcuda();
	pmeshParams->susrel_AFM.null_p_cu_obj_mpcuda();
	pmeshParams->susprel.null_p_cu_obj_mpcuda();

	pmeshParams->cHA.null_p_cu_obj_mpcuda();
	pmeshParams->cHmo.null_p_cu_obj_mpcuda();

	pmeshParams->elecCond.null_p_cu_obj_mpcuda();
	pmeshParams->amrPercentage.null_p_cu_obj_mpcuda();

	pmeshParams->P.null_p_cu_obj_mpcuda();
	pmeshParams->beta.null_p_cu_obj_mpcuda();

	pmeshParams->De.null_p_cu_obj_mpcuda();
	pmeshParams->n_density.null_p_cu_obj_mpcuda();
	pmeshParams->betaD.null_p_cu_obj_mpcuda();
	
	pmeshParams->SHA.null_p_cu_obj_mpcuda();
	pmeshParams->iSHA.null_p_cu_obj_mpcuda();
	pmeshParams->flSOT.null_p_cu_obj_mpcuda();

	pmeshParams->STq.null_p_cu_obj_mpcuda();
	pmeshParams->STa.null_p_cu_obj_mpcuda();
	pmeshParams->STp.null_p_cu_obj_mpcuda();
	
	pmeshParams->l_sf.null_p_cu_obj_mpcuda();
	pmeshParams->l_ex.null_p_cu_obj_mpcuda();
	pmeshParams->l_ph.null_p_cu_obj_mpcuda();
	
	pmeshParams->Gi.null_p_cu_obj_mpcuda();
	pmeshParams->Gmix.null_p_cu_obj_mpcuda();
	
	pmeshParams->ts_eff.null_p_cu_obj_mpcuda();
	pmeshParams->tsi_eff.null_p_cu_obj_mpcuda();
	
	pmeshParams->pump_eff.null_p_cu_obj_mpcuda();
	pmeshParams->cpump_eff.null_p_cu_obj_mpcuda();
	pmeshParams->the_eff.null_p_cu_obj_mpcuda();

	pmeshParams->atomic_moment.null_p_cu_obj_mpcuda();
	pmeshParams->atomic_moment_AFM.null_p_cu_obj_mpcuda();

	pmeshParams->thermCond.null_p_cu_obj_mpcuda();
	pmeshParams->density.null_p_cu_obj_mpcuda();
	
	pmeshParams->shc.null_p_cu_obj_mpcuda();
	pmeshParams->shc_e.null_p_cu_obj_mpcuda();
	pmeshParams->G_e.null_p_cu_obj_mpcuda();

	pmeshParams->cT.null_p_cu_obj_mpcuda();
	pmeshParams->Q.null_p_cu_obj_mpcuda();

	pmeshParams->MEc.null_p_cu_obj_mpcuda();
	pmeshParams->Ym.null_p_cu_obj_mpcuda();
	pmeshParams->Pr.null_p_cu_obj_mpcuda();
}

#endif