#include "stdafx.h"
#include "Atom_MeshParamsCUDA.h"
#include "Atom_MeshParams.h"

#if COMPILECUDA == 1

Atom_MeshParamsCUDA::Atom_MeshParamsCUDA(Atom_MeshParams *pameshParams)
{
	this->pameshParams = pameshParams;

	//-----------SIMPLE CUBIC

	grel()->set_from_cpu(pameshParams->grel);
	pameshParams->grel.set_p_cu_obj_mpcuda(&grel);

	alpha()->set_from_cpu(pameshParams->alpha);
	pameshParams->alpha.set_p_cu_obj_mpcuda(&alpha);

	mu_s()->set_from_cpu(pameshParams->mu_s);
	pameshParams->mu_s.set_p_cu_obj_mpcuda(&mu_s);

	J()->set_from_cpu(pameshParams->J);
	pameshParams->J.set_p_cu_obj_mpcuda(&J);
	
	D()->set_from_cpu(pameshParams->D);
	pameshParams->D.set_p_cu_obj_mpcuda(&D);

	D_dir()->set_from_cpu(pameshParams->D_dir);
	pameshParams->D_dir.set_p_cu_obj_mpcuda(&D_dir);

	Js()->set_from_cpu(pameshParams->Js);
	pameshParams->Js.set_p_cu_obj_mpcuda(&Js);
	Js2()->set_from_cpu(pameshParams->Js2);
	pameshParams->Js2.set_p_cu_obj_mpcuda(&Js2);

	K1()->set_from_cpu(pameshParams->K1);
	pameshParams->K1.set_p_cu_obj_mpcuda(&K1);
	K2()->set_from_cpu(pameshParams->K2);
	pameshParams->K2.set_p_cu_obj_mpcuda(&K2);
	K3()->set_from_cpu(pameshParams->K3);
	pameshParams->K3.set_p_cu_obj_mpcuda(&K3);

	mcanis_ea1()->set_from_cpu(pameshParams->mcanis_ea1);
	pameshParams->mcanis_ea1.set_p_cu_obj_mpcuda(&mcanis_ea1);
	mcanis_ea2()->set_from_cpu(pameshParams->mcanis_ea2);
	pameshParams->mcanis_ea2.set_p_cu_obj_mpcuda(&mcanis_ea2);
	mcanis_ea3()->set_from_cpu(pameshParams->mcanis_ea3);
	pameshParams->mcanis_ea3.set_p_cu_obj_mpcuda(&mcanis_ea3);

	//-----------Others

	Nxy()->set_from_cpu(pameshParams->Nxy);
	pameshParams->Nxy.set_p_cu_obj_mpcuda(&Nxy);

	cHA()->set_from_cpu(pameshParams->cHA);
	pameshParams->cHA.set_p_cu_obj_mpcuda(&cHA);
	cHmo()->set_from_cpu(pameshParams->cHmo);
	pameshParams->cHmo.set_p_cu_obj_mpcuda(&cHmo);

	s_eff()->set_from_cpu(pameshParams->s_eff);
	pameshParams->s_eff.set_p_cu_obj_mpcuda(&s_eff);

	elecCond()->set_from_cpu(pameshParams->elecCond);
	pameshParams->elecCond.set_p_cu_obj_mpcuda(&elecCond);
	amrPercentage()->set_from_cpu(pameshParams->amrPercentage);
	pameshParams->amrPercentage.set_p_cu_obj_mpcuda(&amrPercentage);
	RAtmr_p()->set_from_cpu(pameshParams->RAtmr_p);
	pameshParams->RAtmr_p.set_p_cu_obj_mpcuda(&RAtmr_p);
	RAtmr_ap()->set_from_cpu(pameshParams->RAtmr_ap);
	pameshParams->RAtmr_ap.set_p_cu_obj_mpcuda(&RAtmr_ap);

	P()->set_from_cpu(pameshParams->P);
	pameshParams->P.set_p_cu_obj_mpcuda(&P);
	beta()->set_from_cpu(pameshParams->beta);
	pameshParams->beta.set_p_cu_obj_mpcuda(&beta);

	De()->set_from_cpu(pameshParams->De);
	pameshParams->De.set_p_cu_obj_mpcuda(&De);

	n_density()->set_from_cpu(pameshParams->n_density);
	pameshParams->n_density.set_p_cu_obj_mpcuda(&n_density);

	betaD()->set_from_cpu(pameshParams->betaD);
	pameshParams->betaD.set_p_cu_obj_mpcuda(&betaD);

	SHA()->set_from_cpu(pameshParams->SHA);
	pameshParams->SHA.set_p_cu_obj_mpcuda(&SHA);
	flSOT()->set_from_cpu(pameshParams->flSOT);
	pameshParams->flSOT.set_p_cu_obj_mpcuda(&flSOT);

	STq()->set_from_cpu(pameshParams->STq);
	pameshParams->STq.set_p_cu_obj_mpcuda(&STq);
	STa()->set_from_cpu(pameshParams->STa);
	pameshParams->STa.set_p_cu_obj_mpcuda(&STa);
	STp()->set_from_cpu(pameshParams->STp);
	pameshParams->STp.set_p_cu_obj_mpcuda(&STp);

	l_sf()->set_from_cpu(pameshParams->l_sf);
	pameshParams->l_sf.set_p_cu_obj_mpcuda(&l_sf);
	l_ex()->set_from_cpu(pameshParams->l_ex);
	pameshParams->l_ex.set_p_cu_obj_mpcuda(&l_ex);
	l_ph()->set_from_cpu(pameshParams->l_ph);
	pameshParams->l_ph.set_p_cu_obj_mpcuda(&l_ph);

	Gi()->set_from_cpu(pameshParams->Gi);
	pameshParams->Gi.set_p_cu_obj_mpcuda(&Gi);
	Gmix()->set_from_cpu(pameshParams->Gmix);
	pameshParams->Gmix.set_p_cu_obj_mpcuda(&Gmix);

	ts_eff()->set_from_cpu(pameshParams->ts_eff);
	pameshParams->ts_eff.set_p_cu_obj_mpcuda(&ts_eff);
	tsi_eff()->set_from_cpu(pameshParams->tsi_eff);
	pameshParams->tsi_eff.set_p_cu_obj_mpcuda(&tsi_eff);

	pump_eff()->set_from_cpu(pameshParams->pump_eff);
	pameshParams->pump_eff.set_p_cu_obj_mpcuda(&pump_eff);

	cpump_eff()->set_from_cpu(pameshParams->cpump_eff);
	pameshParams->cpump_eff.set_p_cu_obj_mpcuda(&cpump_eff);

	the_eff()->set_from_cpu(pameshParams->the_eff);
	pameshParams->the_eff.set_p_cu_obj_mpcuda(&the_eff);

	thermCond()->set_from_cpu(pameshParams->thermCond);
	pameshParams->thermCond.set_p_cu_obj_mpcuda(&thermCond);

	base_temperature.from_cpu(pameshParams->base_temperature);
	
	density()->set_from_cpu(pameshParams->density);
	pameshParams->density.set_p_cu_obj_mpcuda(&density);

	shc()->set_from_cpu(pameshParams->shc);
	pameshParams->shc.set_p_cu_obj_mpcuda(&shc);
	shc_e()->set_from_cpu(pameshParams->shc_e);
	pameshParams->shc_e.set_p_cu_obj_mpcuda(&shc_e);
	G_e()->set_from_cpu(pameshParams->G_e);
	pameshParams->G_e.set_p_cu_obj_mpcuda(&G_e);

	cT()->set_from_cpu(pameshParams->cT);
	pameshParams->cT.set_p_cu_obj_mpcuda(&cT);
	
	Q()->set_from_cpu(pameshParams->Q);
	pameshParams->Q.set_p_cu_obj_mpcuda(&Q);
}

Atom_MeshParamsCUDA::~Atom_MeshParamsCUDA()
{
	//fine to access data in MeshParams here : Mesh inherits from MeshParams, so in the destruction process Mesh gets destroyed first. 
	//Mesh destructor then calls for MeshCUDA implementation to be destroyed, then we get here since MeshCUDA inherits from MeshParamsCUDA. After this we return back to Mesh destructor to continue destruction down the list.

	//-----------SIMPLE CUBIC

	pameshParams->grel.null_p_cu_obj_mpcuda();

	pameshParams->alpha.null_p_cu_obj_mpcuda();

	pameshParams->mu_s.null_p_cu_obj_mpcuda();

	pameshParams->J.null_p_cu_obj_mpcuda();
	
	pameshParams->D.null_p_cu_obj_mpcuda();
	pameshParams->D_dir.null_p_cu_obj_mpcuda();

	pameshParams->Js.null_p_cu_obj_mpcuda();
	pameshParams->Js2.null_p_cu_obj_mpcuda();

	pameshParams->K1.null_p_cu_obj_mpcuda();
	pameshParams->K2.null_p_cu_obj_mpcuda();
	pameshParams->K3.null_p_cu_obj_mpcuda();

	pameshParams->mcanis_ea1.null_p_cu_obj_mpcuda();
	pameshParams->mcanis_ea2.null_p_cu_obj_mpcuda();
	pameshParams->mcanis_ea3.null_p_cu_obj_mpcuda();

	//-----------Others

	pameshParams->Nxy.null_p_cu_obj_mpcuda();

	pameshParams->cHA.null_p_cu_obj_mpcuda();
	pameshParams->cHmo.null_p_cu_obj_mpcuda();

	pameshParams->s_eff.null_p_cu_obj_mpcuda();

	pameshParams->elecCond.null_p_cu_obj_mpcuda();
	pameshParams->amrPercentage.null_p_cu_obj_mpcuda();
	pameshParams->RAtmr_p.null_p_cu_obj_mpcuda();
	pameshParams->RAtmr_ap.null_p_cu_obj_mpcuda();

	pameshParams->P.null_p_cu_obj_mpcuda();
	pameshParams->beta.null_p_cu_obj_mpcuda();

	pameshParams->De.null_p_cu_obj_mpcuda();
	pameshParams->n_density.null_p_cu_obj_mpcuda();
	pameshParams->betaD.null_p_cu_obj_mpcuda();

	pameshParams->SHA.null_p_cu_obj_mpcuda();
	pameshParams->flSOT.null_p_cu_obj_mpcuda();

	pameshParams->STq.null_p_cu_obj_mpcuda();
	pameshParams->STa.null_p_cu_obj_mpcuda();
	pameshParams->STp.null_p_cu_obj_mpcuda();

	pameshParams->l_sf.null_p_cu_obj_mpcuda();
	pameshParams->l_ex.null_p_cu_obj_mpcuda();
	pameshParams->l_ph.null_p_cu_obj_mpcuda();

	pameshParams->Gi.null_p_cu_obj_mpcuda();
	pameshParams->Gmix.null_p_cu_obj_mpcuda();

	pameshParams->ts_eff.null_p_cu_obj_mpcuda();
	pameshParams->tsi_eff.null_p_cu_obj_mpcuda();

	pameshParams->pump_eff.null_p_cu_obj_mpcuda();
	pameshParams->cpump_eff.null_p_cu_obj_mpcuda();
	pameshParams->the_eff.null_p_cu_obj_mpcuda();

	pameshParams->thermCond.null_p_cu_obj_mpcuda();
	
	pameshParams->density.null_p_cu_obj_mpcuda();

	pameshParams->shc.null_p_cu_obj_mpcuda();
	pameshParams->shc_e.null_p_cu_obj_mpcuda();
	pameshParams->G_e.null_p_cu_obj_mpcuda();

	pameshParams->cT.null_p_cu_obj_mpcuda();
	
	pameshParams->Q.null_p_cu_obj_mpcuda();
}

#endif