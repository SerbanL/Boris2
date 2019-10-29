#include "stdafx.h"
#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1

#include "DiffEqAFM.h"

#include "Mesh_AntiFerromagnetic.h"
#include "Mesh_AntiFerromagneticCUDA.h"

DifferentialEquationAFMCUDA::DifferentialEquationAFMCUDA(DifferentialEquation *pmeshODE) :
	DifferentialEquationCUDA(pmeshODE)
{
	error_on_create = AllocateMemory();
}

DifferentialEquationAFMCUDA::~DifferentialEquationAFMCUDA()
{
	CleanupMemory();
}

//---------------------------------------- SET-UP METHODS  : DiffEqCUDA.cpp and DiffEqCUDA.cu

//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
BError DifferentialEquationAFMCUDA::AllocateMemory(void)
{
	BError error(CLASS_STR(DifferentialEquationAFMCUDA));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	//copy values from cpu : it's possible the user switches to CUDA during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	else sM1()->copy_from_cpuvec(pmeshODE->sM1);

	if (!sM1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	//copy values from cpu : it's possible the user switches to CUDA during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	else sM1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sM1_2);

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

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);

		if (!sEval1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);

		if (!sEval2_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2);
		break;

	case EVAL_ABM:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);

		if (!sEval1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);
		break;

	case EVAL_RK23:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1()->copy_from_cpuvec(pmeshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2()->copy_from_cpuvec(pmeshODE->sEval2);

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);

		if (!sEval1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);

		if (!sEval2_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2);
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

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);

		if (!sEval1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);

		if (!sEval2_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2);

		if (!sEval3_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval3_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval3_2);

		if (!sEval4_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval4_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval4_2);
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

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);

		if (!sEval1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);

		if (!sEval2_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2);

		if (!sEval3_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval3_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval3_2);

		if (!sEval4_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval4_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval4_2);
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

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);

		if (!sEval1_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval1_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);

		if (!sEval2_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval2_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2);

		if (!sEval3_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval3_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval3_2);

		if (!sEval4_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval4_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval4_2);

		if (!sEval5_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval5_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval5_2);
		break;

	case EVAL_SD:
		if (!sEval0()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0()->copy_from_cpuvec(pmeshODE->sEval0);

		if (!sEval0_2()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else sEval0_2()->copy_from_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);
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
void DifferentialEquationAFMCUDA::CleanupMemory(void)
{
	//copy values to cpu before erasing : it's possible the user switches CUDA off during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	//CleanupMemory will be called by destructor in this case, so before destroying gpu data copy it over to cpu if possible
	//CleanupMemory may also be called in other circumstances, in particular from the cpu version of CleanupMemory, after having cleaned cpu vecs, thus in this case the copy methods will not run.

	//Only clear vectors not used for current evaluation method
	if (sM1()->size_cpu() == pmeshODE->sM1.size()) sM1()->copy_to_cpuvec(pmeshODE->sM1);
	sM1()->clear();

	if (sM1_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sM1_2.size()) sM1_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sM1_2);
	sM1_2()->clear();

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

		if (sEval0_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2.size()) sEval0_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval0_2);
		sEval0_2()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 &&
		pmeshODE->evalMethod != EVAL_ABM &&
		pmeshODE->evalMethod != EVAL_RK23 &&
		pmeshODE->evalMethod != EVAL_RKF &&
		pmeshODE->evalMethod != EVAL_RKCK &&
		pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval1()->size_cpu() == pmeshODE->sEval1.size()) sEval1()->copy_to_cpuvec(pmeshODE->sEval1);
		sEval1()->clear();

		if (sEval1_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2.size()) sEval1_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval1_2);
		sEval1_2()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 &&
		pmeshODE->evalMethod != EVAL_RK23 &&
		pmeshODE->evalMethod != EVAL_RKF &&
		pmeshODE->evalMethod != EVAL_RKCK &&
		pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval2()->size_cpu() == pmeshODE->sEval2.size()) sEval2()->copy_to_cpuvec(pmeshODE->sEval2);
		sEval2()->clear();

		if (sEval2_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2.size()) sEval2_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval2_2);
		sEval2_2()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RK4 &&
		pmeshODE->evalMethod != EVAL_RKF &&
		pmeshODE->evalMethod != EVAL_RKCK &&
		pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval3()->size_cpu() == pmeshODE->sEval3.size()) sEval3()->copy_to_cpuvec(pmeshODE->sEval3);
		sEval3()->clear();

		if (sEval4()->size_cpu() == pmeshODE->sEval4.size()) sEval4()->copy_to_cpuvec(pmeshODE->sEval4);
		sEval4()->clear();

		if (sEval3_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval3_2.size()) sEval3_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval3_2);
		sEval3_2()->clear();

		if (sEval4_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval4_2.size()) sEval4_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval4_2);
		sEval4_2()->clear();
	}

	if (pmeshODE->evalMethod != EVAL_RKDP) {

		if (sEval5()->size_cpu() == pmeshODE->sEval5.size()) sEval5()->copy_to_cpuvec(pmeshODE->sEval5);
		sEval5()->clear();

		if (sEval5_2()->size_cpu() == reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval5_2.size()) sEval5_2()->copy_to_cpuvec(reinterpret_cast<DifferentialEquationAFM*>(pmeshODE)->sEval5_2);
		sEval5_2()->clear();
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


BError DifferentialEquationAFMCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquationFMCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		error = AllocateMemory();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_ODE_MOVEMESH)) {

		if (!error) {

			//set skip cells flags for moving mesh if enabled
			if (pmeshODE->moving_mesh) {

				Rect mesh_rect = pMesh->GetMeshRect();

				DBL3 end_size = mesh_rect.size() & DBL3(MOVEMESH_ENDRATIO, 1, 1);

				Rect end_rect_left = Rect(mesh_rect.s, mesh_rect.s + end_size);
				Rect end_rect_right = Rect(mesh_rect.e - end_size, mesh_rect.e);

				pMeshCUDA->M()->set_skipcells((cuRect)end_rect_left);
				pMeshCUDA->M()->set_skipcells((cuRect)end_rect_right);

				pMeshCUDA->M2()->set_skipcells((cuRect)end_rect_left);
				pMeshCUDA->M2()->set_skipcells((cuRect)end_rect_right);
			}
			else {

				pMeshCUDA->M()->clear_skipcells();
				pMeshCUDA->M2()->clear_skipcells();
			}
		}
	}

	return error;
}

#endif