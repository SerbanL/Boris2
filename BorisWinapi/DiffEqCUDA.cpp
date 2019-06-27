#include "stdafx.h"
#include "DiffEqCUDA.h"

#include "DiffEq.h"

#include "Mesh_Ferromagnetic.h"
#include "Mesh_FerromagneticCUDA.h"

#if COMPILECUDA == 1

cu_obj<cuReal>* ODECommonCUDA::pdT = nullptr;
cu_obj<cuReal>* ODECommonCUDA::pdT_last = nullptr;

cu_obj<cuReal>* ODECommonCUDA::pmxh = nullptr;
cu_obj<cuReal3>* ODECommonCUDA::pmxh_av = nullptr;
cu_obj<size_t>* ODECommonCUDA::pavpoints = nullptr;

cu_obj<cuReal>* ODECommonCUDA::plte = nullptr;

cu_obj<bool>* ODECommonCUDA::prenormalize = nullptr;

cu_obj<bool>* ODECommonCUDA::psolve_spin_current = nullptr;

cu_obj<int>* ODECommonCUDA::psetODE = nullptr;

cu_obj<bool>* ODECommonCUDA::palternator = nullptr;

ODECommon* ODECommonCUDA::pODE = nullptr;

ODECommonCUDA::ODECommonCUDA(ODECommon *pODE_)
{
	pODE = pODE_;
	
	if(!pdT) pdT = new cu_obj<cuReal>();
	if (!pdT_last) pdT_last = new cu_obj<cuReal>();
	
	if (!pmxh) pmxh = new cu_obj<cuReal>();
	if (!pmxh_av) pmxh_av = new cu_obj<cuReal3>();
	if (!pavpoints) pavpoints = new cu_obj<size_t>();
	
	if (!plte) plte = new cu_obj<cuReal>();
	
	if (!prenormalize) prenormalize = new cu_obj<bool>();
	
	if (!psolve_spin_current) psolve_spin_current = new cu_obj<bool>();
	
	if (!psetODE) psetODE = new cu_obj<int>();
	
	if (!palternator) palternator = new cu_obj<bool>();

	SyncODEValues();
}

ODECommonCUDA::~ODECommonCUDA()
{
	if (pODE) {

		pODE->dT = pdT->to_cpu();
		pODE->dT_last = pdT_last->to_cpu();
		pODE->mxh = pmxh->to_cpu();
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

	switch (pmeshODE->evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_TEULER:
		break;

	case EVAL_RK4:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		break;

	case EVAL_ABM:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		break;

	case EVAL_RKF:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval3()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!sEval4()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	bool prng_used = false;

	switch (pmeshODE->setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if (!H_Thermal()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		prng_used = true;
		break;

	case ODE_SLLB:
	case ODE_SLLBSTT:
	case ODE_SLLBSA:
		if (!H_Thermal()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!Torque_Thermal()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
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
	//Only clear vectors not used for current evaluation method
	sM1()->clear();

	if (pmeshODE->evalMethod != EVAL_RK4 && pmeshODE->evalMethod != EVAL_ABM && pmeshODE->evalMethod != EVAL_RKF) {

		sEval0()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 && pmeshODE->evalMethod != EVAL_ABM && pmeshODE->evalMethod != EVAL_RKF) {

		sEval1()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 && pmeshODE->evalMethod != EVAL_RKF) {

		sEval2()->clear();
		sEval3()->clear();
		sEval4()->clear();
	}

	//For thermal vecs only clear if not used for current set ODE
	if (pmeshODE->setODE != ODE_SLLG && pmeshODE->setODE != ODE_SLLGSTT && pmeshODE->setODE != ODE_SLLB && pmeshODE->setODE != ODE_SLLBSTT && pmeshODE->setODE != ODE_SLLGSA && pmeshODE->setODE != ODE_SLLBSA) {

		H_Thermal()->clear();
	}

	if (pmeshODE->setODE != ODE_SLLB && pmeshODE->setODE != ODE_SLLBSTT && pmeshODE->setODE != ODE_SLLBSA) {

		Torque_Thermal()->clear();
	}

	prng()->clear();
}

#endif