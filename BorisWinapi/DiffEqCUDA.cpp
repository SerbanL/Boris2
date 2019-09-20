#include "stdafx.h"
#include "DiffEqCUDA.h"

#include "DiffEq.h"

#include "Mesh_Ferromagnetic.h"
#include "Mesh_FerromagneticCUDA.h"

#if COMPILECUDA == 1

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

ODECommon* ODECommonCUDA::pODE = nullptr;

ODECommonCUDA::ODECommonCUDA(ODECommon *pODE_)
{
	pODE = pODE_;
	
	if(!pdT) pdT = new cu_obj<cuBReal>();
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

	SyncODEValues();
}

ODECommonCUDA::~ODECommonCUDA()
{
	if (pODE) {

		pODE->dT = pdT->to_cpu();
		pODE->dT_last = pdT_last->to_cpu();
		pODE->mxh = pmxh->to_cpu();
		pODE->dmdt = pdmdt->to_cpu();
	}
}

BError ODECommonCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(__FUNCTION__);

	return error;
}

void ODECommonCUDA::SyncODEValues(void)
{
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

//---------------------------------------------------------------

DifferentialEquationCUDA::DifferentialEquationCUDA(DifferentialEquation *pmeshODE_) :
	ODECommonCUDA()
{
	pmeshODE = pmeshODE_;
	
	pMesh = pmeshODE->pMesh;
	pMeshCUDA = dynamic_cast<FMeshCUDA*>(pmeshODE->pMesh->pMeshCUDA);

	AllocateMemory();
}

BError DifferentialEquationCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquationCUDA));
	
	error = AllocateMemory();

	if (!error) {

		//set skip cells flags for moving mesh if enabled
		if (pmeshODE->moving_mesh) {

			Rect mesh_rect = pMesh->GetMeshRect();

			DBL3 end_size = mesh_rect.size() & DBL3(MOVEMESH_ENDRATIO, 1, 1);

			Rect end_rect_left = Rect(mesh_rect.s, mesh_rect.s + end_size);
			Rect end_rect_right = Rect(mesh_rect.e - end_size, mesh_rect.e);

			pMeshCUDA->M()->set_skipcells((cuRect)end_rect_left);
			pMeshCUDA->M()->set_skipcells((cuRect)end_rect_right);
		}
		else {

			pMeshCUDA->M()->clear_skipcells();
		}
	}
	
	return error;
}

//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
BError DifferentialEquationCUDA::AllocateMemory(void)
{
	BError error(CLASS_STR(DifferentialEquationCUDA));
	
	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	//copy values from cpu : it's possible the user switches to CUDA during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	else sM1()->copy_from_cpuvec(pmeshODE->sM1);

	switch (pmeshODE->evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_AHEUN:
	case EVAL_TEULER:
		break;

	case EVAL_RK4:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2()->copy_from_cpuvec(pmeshODE->sEval2);
		break;

	case EVAL_ABM:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);
		break;

	case EVAL_RK23:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2()->copy_from_cpuvec(pmeshODE->sEval2);
		break;

	case EVAL_RKF:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2()->copy_from_cpuvec(pmeshODE->sEval2);

		if (!sEval3()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval3()->copy_from_cpuvec(pmeshODE->sEval3);

		if (!sEval4()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval4()->copy_from_cpuvec(pmeshODE->sEval4);
		break;

	case EVAL_RKCK:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2()->copy_from_cpuvec(pmeshODE->sEval2);

		if (!sEval3()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval3()->copy_from_cpuvec(pmeshODE->sEval3);

		if (!sEval4()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval4()->copy_from_cpuvec(pmeshODE->sEval4);
		break;

	case EVAL_RKDP:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2()->copy_from_cpuvec(pmeshODE->sEval2);

		if (!sEval3()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval3()->copy_from_cpuvec(pmeshODE->sEval3);

		if (!sEval4()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval4()->copy_from_cpuvec(pmeshODE->sEval4);

		if (!sEval5()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval5()->copy_from_cpuvec(pmeshODE->sEval5);
		break;

	case EVAL_SD:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	bool prng_used = false;

	switch (pmeshODE->setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if (!H_Thermal()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else H_Thermal()->copy_from_cpuvec(pmeshODE->H_Thermal);
		prng_used = true;
		break;

	case ODE_SLLB:
	case ODE_SLLBSTT:
	case ODE_SLLBSA:
		if (!H_Thermal()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else H_Thermal()->copy_from_cpuvec(pmeshODE->H_Thermal);

		if (!Torque_Thermal()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else Torque_Thermal()->copy_from_cpuvec(pmeshODE->Torque_Thermal);
		prng_used = true;
		break;
	}

	if (prng_used) {

		//initialize the pseudo-random number generator with a seed and memory size - recommended use kernel size divided by 128
		if (prng()->initialize(GetTickCount(), pMesh->n.dim() / 128) != cudaSuccess) error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	return error;
}

//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
void DifferentialEquationCUDA::CleanupMemory(void)
{
	//copy values to cpu before erasing : it's possible the user switches CUDA off during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	//CleanupMemory will be called by destructor in this case, so before destroying gpu data copy it over to cpu if possible
	//CleanupMemory may also be called in other circumstances, in particular from the cpu version of CleanupMemory, after having cleaned cpu vecs, thus in this case the copy methods will not run.

	//Only clear vectors not used for current evaluation method
	if (sM1()->size_cpu() == pmeshODE->sM1.size()) sM1()->copy_to_cpuvec(pmeshODE->sM1);
	sM1()->clear();

	if (
		pmeshODE->evalMethod != EVAL_RK4 && 
		pmeshODE->evalMethod != EVAL_ABM && 
		pmeshODE->evalMethod != EVAL_RK23 && 
		pmeshODE->evalMethod != EVAL_RKF && 
		pmeshODE->evalMethod != EVAL_RKCK && 
		pmeshODE->evalMethod != EVAL_RKDP && 
		pmeshODE->evalMethod != EVAL_SD) {

		if (sEval0()->size_cpu() == pmeshODE->sEval0.size()) sEval0()->copy_to_cpuvec(pmeshODE->sEval0);
		sEval0()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 && 
		pmeshODE->evalMethod != EVAL_ABM && 
		pmeshODE->evalMethod != EVAL_RK23 && 
		pmeshODE->evalMethod != EVAL_RKF &&
		pmeshODE->evalMethod != EVAL_RKCK &&
		pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval1()->size_cpu() == pmeshODE->sEval1.size()) sEval1()->copy_to_cpuvec(pmeshODE->sEval1);
		sEval1()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 && 
		pmeshODE->evalMethod != EVAL_RK23 && 
		pmeshODE->evalMethod != EVAL_RKF &&
		pmeshODE->evalMethod != EVAL_RKCK &&
		pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval2()->size_cpu() == pmeshODE->sEval2.size()) sEval2()->copy_to_cpuvec(pmeshODE->sEval2);
		sEval2()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 && 
		pmeshODE->evalMethod != EVAL_RKF &&
		pmeshODE->evalMethod != EVAL_RKCK &&
		pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval3()->size_cpu() == pmeshODE->sEval3.size()) sEval3()->copy_to_cpuvec(pmeshODE->sEval3);
		sEval3()->clear();

		if (sEval4()->size_cpu() == pmeshODE->sEval4.size()) sEval4()->copy_to_cpuvec(pmeshODE->sEval4);
		sEval4()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval5()->size_cpu() == pmeshODE->sEval5.size()) sEval5()->copy_to_cpuvec(pmeshODE->sEval5);
		sEval5()->clear();
	}

	//For thermal vecs only clear if not used for current set ODE
	if (pmeshODE->setODE != ODE_SLLG && 
		pmeshODE->setODE != ODE_SLLGSTT && 
		pmeshODE->setODE != ODE_SLLB && 
		pmeshODE->setODE != ODE_SLLBSTT && 
		pmeshODE->setODE != ODE_SLLGSA && 
		pmeshODE->setODE != ODE_SLLBSA) {

		if (H_Thermal()->size_cpu() == pmeshODE->H_Thermal.size()) H_Thermal()->copy_to_cpuvec(pmeshODE->H_Thermal);
		H_Thermal()->clear();
	}

	if (pmeshODE->setODE != ODE_SLLB && 
		pmeshODE->setODE != ODE_SLLBSTT && 
		pmeshODE->setODE != ODE_SLLBSA) {

		if (Torque_Thermal()->size_cpu() == pmeshODE->Torque_Thermal.size()) Torque_Thermal()->copy_to_cpuvec(pmeshODE->Torque_Thermal);
		Torque_Thermal()->clear();
	}

	prng()->clear();
}

#endif