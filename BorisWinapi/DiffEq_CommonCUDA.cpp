#include "stdafx.h"
#include "DiffEq_CommonCUDA.h"
#include "DiffEq_Common.h"

#if COMPILECUDA == 1

cu_obj<cuBReal>* ODECommonCUDA::ptime = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pstagetime = nullptr;

cu_obj<cuBReal>* ODECommonCUDA::pdT = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pdT_last = nullptr;

cu_obj<cuBReal>* ODECommonCUDA::pmxh = nullptr;
cu_obj<cuReal3>* ODECommonCUDA::pmxh_av = nullptr;
cu_obj<size_t>* ODECommonCUDA::pavpoints = nullptr;

cu_obj<cuBReal>* ODECommonCUDA::pdmdt = nullptr;
cu_obj<cuReal3>* ODECommonCUDA::pdmdt_av = nullptr;
cu_obj<size_t>* ODECommonCUDA::pavpoints2 = nullptr;

cu_obj<cuBReal>* ODECommonCUDA::plte = nullptr;

cu_obj<bool>* ODECommonCUDA::prenormalize = nullptr;

cu_obj<bool>* ODECommonCUDA::psolve_spin_current = nullptr;

cu_obj<int>* ODECommonCUDA::psetODE = nullptr;

cu_obj<bool>* ODECommonCUDA::palternator = nullptr;

cu_obj<cuBReal>* ODECommonCUDA::pdelta_M_sq = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pdelta_G_sq = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pdelta_M_dot_delta_G = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pdelta_M2_sq = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pdelta_G2_sq = nullptr;
cu_obj<cuBReal>* ODECommonCUDA::pdelta_M2_dot_delta_G2 = nullptr;

ODECommon* ODECommonCUDA::pODE = nullptr;

ODECommonCUDA::ODECommonCUDA(ODECommon *pODE_)
{
	pODE = pODE_;

	if (!ptime) ptime = new cu_obj<cuBReal>();
	if (!pstagetime) pstagetime = new cu_obj<cuBReal>();

	if (!pdT) pdT = new cu_obj<cuBReal>();
	if (!pdT_last) pdT_last = new cu_obj<cuBReal>();

	if (!pmxh) pmxh = new cu_obj<cuBReal>();
	if (!pmxh_av) pmxh_av = new cu_obj<cuReal3>();
	if (!pavpoints) pavpoints = new cu_obj<size_t>();

	if (!pdmdt) pdmdt = new cu_obj<cuBReal>();
	if (!pdmdt_av) pdmdt_av = new cu_obj<cuReal3>();
	if (!pavpoints2) pavpoints2 = new cu_obj<size_t>();

	if (!plte) plte = new cu_obj<cuBReal>();

	if (!prenormalize) prenormalize = new cu_obj<bool>();

	if (!psolve_spin_current) psolve_spin_current = new cu_obj<bool>();

	if (!psetODE) psetODE = new cu_obj<int>();

	if (!palternator) palternator = new cu_obj<bool>();

	if (!pdelta_M_sq) pdelta_M_sq = new cu_obj<cuBReal>();
	if (!pdelta_G_sq) pdelta_G_sq = new cu_obj<cuBReal>();
	if (!pdelta_M_dot_delta_G) pdelta_M_dot_delta_G = new cu_obj<cuBReal>();

	if (!pdelta_M2_sq) pdelta_M2_sq = new cu_obj<cuBReal>();
	if (!pdelta_G2_sq) pdelta_G2_sq = new cu_obj<cuBReal>();
	if (!pdelta_M2_dot_delta_G2) pdelta_M2_dot_delta_G2 = new cu_obj<cuBReal>();

	SyncODEValues();
}

ODECommonCUDA::~ODECommonCUDA()
{
	if (pODE) {

		if (ptime) {

			pODE->time = ptime->to_cpu();
			delete ptime;
			ptime = nullptr;
		}

		if (pstagetime) {

			pODE->stagetime = pstagetime->to_cpu();
			delete pstagetime;
			pstagetime = nullptr;
		}

		if (pdT) {

			pODE->dT = pdT->to_cpu();
			delete pdT;
			pdT = nullptr;
		}

		if (pdT_last) {

			pODE->dT_last = pdT_last->to_cpu();
			delete pdT_last;
			pdT_last = nullptr;
		}

		if (pmxh) {

			pODE->mxh = pmxh->to_cpu();
			delete pmxh;
			pmxh = nullptr;
		}

		if (pdmdt) {

			pODE->dmdt = pdmdt->to_cpu();
			delete pdmdt;
			pdmdt = nullptr;
		}

		if (pmxh_av) { delete pmxh_av; pmxh_av = nullptr; }
		if (pavpoints) { delete pavpoints; pavpoints = nullptr; }

		if (pdmdt_av) { delete pdmdt_av; pdmdt_av = nullptr; }
		if (pavpoints2) { delete pavpoints2; pavpoints2 = nullptr; }
		if (plte) { delete plte; plte = nullptr; }

		if (prenormalize) { delete  prenormalize; prenormalize = nullptr; }

		if (psolve_spin_current) { delete psolve_spin_current; psolve_spin_current = nullptr; }

		if (psetODE) { delete psetODE; psetODE = nullptr; }

		if (palternator) { delete palternator; palternator = nullptr; }

		if (pdelta_M_sq) { delete pdelta_M_sq; pdelta_M_sq = nullptr; }
		if (pdelta_G_sq) { delete pdelta_G_sq; pdelta_G_sq = nullptr; }
		if (pdelta_M_dot_delta_G) { delete pdelta_M_dot_delta_G; pdelta_M_dot_delta_G = nullptr; }

		if (pdelta_M2_sq) { delete pdelta_M2_sq; pdelta_M2_sq = nullptr; }
		if (pdelta_G2_sq) { delete pdelta_G2_sq; pdelta_G2_sq = nullptr; }
		if (pdelta_M2_dot_delta_G2) { delete pdelta_M2_dot_delta_G2; pdelta_M2_dot_delta_G2 = nullptr; }
	}
}

BError ODECommonCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(__FUNCTION__);

	return error;
}

void ODECommonCUDA::SyncODEValues(void)
{
	ptime->from_cpu(pODE->time);
	pstagetime->from_cpu(pODE->stagetime);

	pdT->from_cpu(pODE->dT);
	pdT_last->from_cpu(pODE->dT_last);

	pmxh->from_cpu(pODE->mxh);
	pdmdt->from_cpu(pODE->dmdt);

	prenormalize->from_cpu(pODE->renormalize);

	psolve_spin_current->from_cpu(pODE->solve_spin_current);

	psetODE->from_cpu(pODE->setODE);

	palternator->from_cpu(pODE->alternator);
}

//set specific cuda values (used often)
void ODECommonCUDA::Sync_time(void)
{
	ptime->from_cpu(pODE->time);
	pstagetime->from_cpu(pODE->stagetime);
}

void ODECommonCUDA::Sync_dT(void)
{
	pdT->from_cpu(pODE->dT);
}

void ODECommonCUDA::Sync_dT_last(void)
{
	pdT_last->from_cpu(pODE->dT_last);
}

void ODECommonCUDA::Sync_alternator(void)
{
	palternator->from_cpu(pODE->alternator);
}

#endif