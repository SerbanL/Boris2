#include "stdafx.h"
#include "MeshParamsCUDA.h"
#include "MeshParams.h"

#if COMPILECUDA == 1

MeshParamsCUDA::MeshParamsCUDA(MeshParams *pmeshParams)
{
	this->pmeshParams = pmeshParams;
	
	grel()->set_from_cpu(pmeshParams->grel);
	pmeshParams->grel.set_p_cu_obj_mpcuda(&grel);

	alpha()->set_from_cpu(pmeshParams->alpha);
	pmeshParams->alpha.set_p_cu_obj_mpcuda(&alpha);

	Ms()->set_from_cpu(pmeshParams->Ms);
	pmeshParams->Ms.set_p_cu_obj_mpcuda(&Ms);

	Nxy()->set_from_cpu(pmeshParams->Nxy);
	pmeshParams->Nxy.set_p_cu_obj_mpcuda(&Nxy);

	A()->set_from_cpu(pmeshParams->A);
	pmeshParams->A.set_p_cu_obj_mpcuda(&A);
	D()->set_from_cpu(pmeshParams->D);
	pmeshParams->D.set_p_cu_obj_mpcuda(&D);

	J1()->set_from_cpu(pmeshParams->J1);
	pmeshParams->J1.set_p_cu_obj_mpcuda(&J1);
	J2()->set_from_cpu(pmeshParams->J2);
	pmeshParams->J2.set_p_cu_obj_mpcuda(&J2);

	K1()->set_from_cpu(pmeshParams->K1);
	pmeshParams->K1.set_p_cu_obj_mpcuda(&K1);
	K2()->set_from_cpu(pmeshParams->K2);
	pmeshParams->K2.set_p_cu_obj_mpcuda(&K2);
	mcanis_ea1()->set_from_cpu(pmeshParams->mcanis_ea1);
	pmeshParams->mcanis_ea1.set_p_cu_obj_mpcuda(&mcanis_ea1);
	mcanis_ea2()->set_from_cpu(pmeshParams->mcanis_ea2);
	pmeshParams->mcanis_ea2.set_p_cu_obj_mpcuda(&mcanis_ea2);

	susrel()->set_from_cpu(pmeshParams->susrel);
	pmeshParams->susrel.set_p_cu_obj_mpcuda(&susrel);
	susprel()->set_from_cpu(pmeshParams->susprel);
	pmeshParams->susprel.set_p_cu_obj_mpcuda(&susprel);

	cHA()->set_from_cpu(pmeshParams->cHA);
	pmeshParams->cHA.set_p_cu_obj_mpcuda(&cHA);

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
	
	betaD()->set_from_cpu(pmeshParams->betaD);
	pmeshParams->betaD.set_p_cu_obj_mpcuda(&betaD);
	
	SHA()->set_from_cpu(pmeshParams->SHA);
	pmeshParams->SHA.set_p_cu_obj_mpcuda(&SHA);
	iSHA()->set_from_cpu(pmeshParams->iSHA);
	pmeshParams->iSHA.set_p_cu_obj_mpcuda(&iSHA);
	flSOT()->set_from_cpu(pmeshParams->flSOT);
	pmeshParams->flSOT.set_p_cu_obj_mpcuda(&flSOT);
	
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

	base_temperature.from_cpu(pmeshParams->base_temperature);
	T_Curie.from_cpu(pmeshParams->T_Curie);

	thermCond()->set_from_cpu(pmeshParams->thermCond);
	pmeshParams->thermCond.set_p_cu_obj_mpcuda(&thermCond);
	density()->set_from_cpu(pmeshParams->density);
	pmeshParams->density.set_p_cu_obj_mpcuda(&density);
	shc()->set_from_cpu(pmeshParams->shc);
	pmeshParams->shc.set_p_cu_obj_mpcuda(&shc);

	cT()->set_from_cpu(pmeshParams->cT);
	pmeshParams->cT.set_p_cu_obj_mpcuda(&cT);
	Q()->set_from_cpu(pmeshParams->Q);
	pmeshParams->Q.set_p_cu_obj_mpcuda(&Q);
}

MeshParamsCUDA::~MeshParamsCUDA()
{
	//fine to access data in MeshParams here : Mesh inherits from MeshParams, so in the destruction process Mesh gets destroyed first. 
	//Mesh destructor then calls for MeshCUDA implementation to be destroyed, then we get here since MeshCUDA inherits from MeshParamsCUDA. After this we return back to Mesh destructor to continue destruction down the list.
	
	pmeshParams->grel.null_p_cu_obj_mpcuda();
	pmeshParams->alpha.null_p_cu_obj_mpcuda();

	pmeshParams->Ms.null_p_cu_obj_mpcuda();

	pmeshParams->Nxy.null_p_cu_obj_mpcuda();

	pmeshParams->A.null_p_cu_obj_mpcuda();
	pmeshParams->D.null_p_cu_obj_mpcuda();

	pmeshParams->J1.null_p_cu_obj_mpcuda();
	pmeshParams->J2.null_p_cu_obj_mpcuda();

	pmeshParams->K1.null_p_cu_obj_mpcuda();
	pmeshParams->K2.null_p_cu_obj_mpcuda();
	pmeshParams->mcanis_ea1.null_p_cu_obj_mpcuda();
	pmeshParams->mcanis_ea2.null_p_cu_obj_mpcuda();

	pmeshParams->susrel.null_p_cu_obj_mpcuda();
	pmeshParams->susprel.null_p_cu_obj_mpcuda();

	pmeshParams->cHA.null_p_cu_obj_mpcuda();

	pmeshParams->elecCond.null_p_cu_obj_mpcuda();
	pmeshParams->amrPercentage.null_p_cu_obj_mpcuda();

	pmeshParams->P.null_p_cu_obj_mpcuda();
	pmeshParams->beta.null_p_cu_obj_mpcuda();

	pmeshParams->De.null_p_cu_obj_mpcuda();
	pmeshParams->betaD.null_p_cu_obj_mpcuda();
	
	pmeshParams->SHA.null_p_cu_obj_mpcuda();
	pmeshParams->iSHA.null_p_cu_obj_mpcuda();
	pmeshParams->flSOT.null_p_cu_obj_mpcuda();
	
	pmeshParams->l_sf.null_p_cu_obj_mpcuda();
	pmeshParams->l_ex.null_p_cu_obj_mpcuda();
	pmeshParams->l_ph.null_p_cu_obj_mpcuda();
	
	pmeshParams->Gi.null_p_cu_obj_mpcuda();
	pmeshParams->Gmix.null_p_cu_obj_mpcuda();
	
	pmeshParams->ts_eff.null_p_cu_obj_mpcuda();
	pmeshParams->tsi_eff.null_p_cu_obj_mpcuda();
	
	pmeshParams->pump_eff.null_p_cu_obj_mpcuda();

	pmeshParams->thermCond.null_p_cu_obj_mpcuda();
	pmeshParams->density.null_p_cu_obj_mpcuda();
	pmeshParams->shc.null_p_cu_obj_mpcuda();

	pmeshParams->cT.null_p_cu_obj_mpcuda();
	pmeshParams->Q.null_p_cu_obj_mpcuda();
}

#endif