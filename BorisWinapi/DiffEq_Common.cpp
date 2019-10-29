#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh.h"
#include "SuperMesh.h"

///////////////////////////////////////////////////////////////////////
// Static variables in ODECommon

int ODECommon::setODE = ODE_LLG;
int ODECommon::evalMethod = EVAL_EULER;

int ODECommon::iteration = 0;
int ODECommon::stageiteration = 0;

double ODECommon::dT = 0.0;
double ODECommon::dT_last = 0.0;
double ODECommon::time = 0.0;
double ODECommon::stagetime = 0.0;
double ODECommon::mxh = 1.0;
double ODECommon::dmdt = 1.0;

bool ODECommon::calculate_mxh = true;
bool ODECommon::calculate_dmdt = true;

double ODECommon::err_high_fail = RKF_RELERRFAIL;
double ODECommon::err_high = RKF_RELERRMAX;
double ODECommon::err_low = RKF_RELERRMIN;
double ODECommon::dT_increase = RKF_DTINCREASE;
double ODECommon::dT_max = RKF_MAXDT;
double ODECommon::dT_min = RKF_MINDT;

double ODECommon::delta_M_sq = 0.0;
double ODECommon::delta_G_sq = 0.0;
double ODECommon::delta_M_dot_delta_G = 0.0;

double ODECommon::delta_M2_sq = 0.0;
double ODECommon::delta_G2_sq = 0.0;
double ODECommon::delta_M2_dot_delta_G2 = 0.0;

int ODECommon::evalStep = 0;

bool ODECommon::available = true;

bool ODECommon::primed = false;

bool ODECommon::alternator = false;

vector_lut<DifferentialEquation*> ODECommon::pODE;

Equation ODECommon::equation;

bool ODECommon::renormalize = true;

bool ODECommon::solve_spin_current = false;

int ODECommon::use_evaluation_speedup = (int)EVALSPEEDUP_NONE;

bool ODECommon::moving_mesh = false;
bool ODECommon::moving_mesh_antisymmetric = true;
double ODECommon::moving_mesh_threshold = MOVEMESH_ANTISYMMETRIC_THRESHOLD;
double ODECommon::moving_mesh_dwshift = 0.0;

#if COMPILECUDA == 1
ODECommonCUDA* ODECommon::pODECUDA = nullptr;
#endif

///////////////////////////////////////////////////////////////////////

ODECommon::ODECommon(bool called_from_derived) :
	ProgramStateNames(this,
		{
			VINFO(iteration), VINFO(stageiteration),
			VINFO(time), VINFO(stagetime),
			VINFO(mxh), VINFO(dmdt),
			VINFO(setODE), VINFO(evalMethod), VINFO(dT),
			VINFO(err_high_fail), VINFO(err_high), VINFO(err_low), VINFO(dT_increase), VINFO(dT_min), VINFO(dT_max),
			VINFO(use_evaluation_speedup),
			VINFO(moving_mesh), VINFO(moving_mesh_antisymmetric), VINFO(moving_mesh_threshold), VINFO(moving_mesh_dwshift)
		}, {})
{
	//when a new ferromagnetic mesh is added this constructor is called with called_from_derived = true
	//we only need to set ODE when ODECommon is first created as these are common to all ferromagnetic meshes (ODE can be changed separately of course later)
	if (!called_from_derived) {

		SetODE((ODE_)setODE, (EVAL_)evalMethod);
	}

	//ODECommon is held in SuperMesh. No need to make pODECUDA here as cuda must be switched off when ODECommon is made (start of program)
}

ODECommon::~ODECommon()
{
#if COMPILECUDA == 1
	if (pODECUDA) delete pODECUDA;
	pODECUDA = nullptr;
#endif
}

void ODECommon::RepairObjectState(void)
{
	//Here need to make sure everything is correctly conigured from primary data (which was just loaded)

	//must remake equation: do not set eval method yet. As meshes are loaded later, they'll each make their own settings for the current evaluation method
	SetODE((ODE_)setODE, (EVAL_)evalMethod, false);
}

//---------------------------------------- SET-UP METHODS

BError ODECommon::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(ODECommon));

	//make sure moving_mesh flag is set correctly
	moving_mesh = false;

	for (int idx = 0; idx < pODE.size(); idx++) {

		if (pODE[idx]->pMesh->GetMoveMeshTrigger()) {

			moving_mesh = true;
			break;
		}
	}

	//reset evaluation method priming and step
	evalStep = 0;
	alternator = false;
	primed = false;

#if COMPILECUDA == 1
	if (pODECUDA) pODECUDA->UpdateConfiguration(cfgMessage);
#endif

	return error;
}

//switch CUDA state on/off
BError ODECommon::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(ODECommon));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		///delete cuda module object and null (just in case)
		if (pODECUDA) delete pODECUDA;
		pODECUDA = nullptr;

		pODECUDA = new ODECommonCUDA(this);
	}
	else {

		//cuda switched off so delete cuda module object
		if (pODECUDA) delete pODECUDA;
		pODECUDA = nullptr;
	}

	primed = false;

#endif

	return error;
}

BError ODECommon::SetODE(ODE_ setODE_, EVAL_ evalMethod_, bool set_eval_method)
{
	BError error(__FUNCTION__);

	//use a function pointer to assign equation to solve
	//this approach saves on having to write out the evaluation methods for every equation, with possible impact on performance due to equation evaluation not being manually inlined
	//tests for typical problems show this has virtually no impact on performance! Maybe the compiler has inlined the code for all the different equations used? (could check this)
	//Or maybe the overhead from function calls is just insignificant for typical problems

	setODE = setODE_;

	switch (setODE) {

	case ODE_LLG:
		equation = &DifferentialEquation::LLG;
		renormalize = true;
		solve_spin_current = false;
		break;

	case ODE_LLGSTATIC:
		equation = &DifferentialEquation::LLGStatic;
		renormalize = true;
		solve_spin_current = false;
		break;

	case ODE_LLGSTT:
		equation = &DifferentialEquation::LLGSTT;
		renormalize = true;
		solve_spin_current = false;
		break;

	case ODE_LLB:
		equation = &DifferentialEquation::LLB;
		renormalize = false;
		solve_spin_current = false;
		break;

	case ODE_LLBSTT:
		equation = &DifferentialEquation::LLBSTT;
		renormalize = false;
		solve_spin_current = false;
		break;

	case ODE_SLLG:
		equation = &DifferentialEquation::SLLG;
		renormalize = true;
		solve_spin_current = false;
		break;

	case ODE_SLLGSTT:
		equation = &DifferentialEquation::SLLGSTT;
		renormalize = true;
		solve_spin_current = false;
		break;

	case ODE_SLLB:
		equation = &DifferentialEquation::SLLB;
		renormalize = false;
		solve_spin_current = false;
		break;

	case ODE_SLLBSTT:
		equation = &DifferentialEquation::SLLBSTT;
		renormalize = false;
		solve_spin_current = false;
		break;

	case ODE_LLGSA:
		equation = &DifferentialEquation::LLG;
		renormalize = true;
		solve_spin_current = true;
		break;

	case ODE_SLLGSA:
		equation = &DifferentialEquation::SLLG;
		renormalize = true;
		solve_spin_current = true;
		break;

	case ODE_LLBSA:
		equation = &DifferentialEquation::LLB;
		renormalize = false;
		solve_spin_current = true;
		break;

	case ODE_SLLBSA:
		equation = &DifferentialEquation::SLLB;
		renormalize = false;
		solve_spin_current = true;
		break;

	default:
		equation = &DifferentialEquation::LLG;
		renormalize = true;
		solve_spin_current = false;
		break;
	}

	if (set_eval_method)
		error = SetEvaluationMethod(evalMethod_);

#if COMPILECUDA == 1
	if (pODECUDA) {

		pODECUDA->SyncODEValues();

		//based on setODE value, setup device method pointers in cuDiffEq held in each DifferentialEquationCUDA object
		for (int idx = 0; idx < (int)pODE.size(); idx++) {

			pODE[idx]->pmeshODECUDA->SetODEMethodPointers();
		}
	}
#endif

	return error;
}

BError ODECommon::SetEvaluationMethod(EVAL_ evalMethod_)
{
	BError error(__FUNCTION__);

	evalMethod = evalMethod_;

	//allocate memory for this evaluation method in all currently set ODEs.
	for (int idx = 0; idx < (int)pODE.size(); idx++)
		if (!error) error = pODE[idx]->AllocateMemory();

	//set default parameters for given evaluation method
	switch (evalMethod) {

	case EVAL_EULER:
	{
		dT = EULER_DEFAULT_DT;
	}
	break;

	case EVAL_TEULER:
	{
		dT = TEULER_DEFAULT_DT;
	}
	break;

	case EVAL_AHEUN:
	{
		dT = AHEUN_DEFAULT_DT;

		err_high_fail = AHEUN_RELERRFAIL;
		err_high = AHEUN_RELERRMAX;
		err_low = AHEUN_RELERRMIN;
		dT_increase = AHEUN_DTINCREASE;
		dT_max = AHEUN_MAXDT;
		dT_min = AHEUN_MINDT;
	}
	break;

	case EVAL_RK4:
	{
		dT = RK4_DEFAULT_DT;
	}
	break;

	case EVAL_ABM:
	{
		dT = ABM_DEFAULT_DT;

		err_high_fail = ABM_RELERRFAIL;
		err_high = ABM_RELERRMAX;
		err_low = ABM_RELERRMIN;
		dT_increase = ABM_DTINCREASE;
		dT_max = ABM_MAXDT;
		dT_min = ABM_MINDT;
	}
	break;

	case EVAL_RK23:
	{
		dT = RK23_DEFAULT_DT;

		err_high_fail = RK23_RELERRFAIL;
		err_high = RK23_RELERRMAX;
		err_low = RK23_RELERRMIN;
		dT_increase = RK23_DTINCREASE;
		dT_max = RK23_MAXDT;
		dT_min = RK23_MINDT;
	}
	break;

	case EVAL_SD:
	{
		//set starting dT - set a very conservative initial stepsize otherwise solver priming can be bad.
		dT = SD_DEFAULT_DT;
	}
	break;

	default:
	case EVAL_RKF:
	{
		dT = RKF_DEFAULT_DT;

		err_high_fail = RKF_RELERRFAIL;
		err_high = RKF_RELERRMAX;
		err_low = RKF_RELERRMIN;
		dT_increase = RKF_DTINCREASE;
		dT_max = RKF_MAXDT;
		dT_min = RKF_MINDT;
	}
	break;

	case EVAL_RKCK:
	{
		dT = RKCK_DEFAULT_DT;

		err_high_fail = RKCK_RELERRFAIL;
		err_high = RKCK_RELERRMAX;
		err_low = RKCK_RELERRMIN;
		dT_increase = RKCK_DTINCREASE;
		dT_max = RKCK_MAXDT;
		dT_min = RKCK_MINDT;
	}
	break;

	case EVAL_RKDP:
	{
		dT = RKDP_DEFAULT_DT;

		err_high_fail = RKDP_RELERRFAIL;
		err_high = RKDP_RELERRMAX;
		err_low = RKDP_RELERRMIN;
		dT_increase = RKDP_DTINCREASE;
		dT_max = RKDP_MAXDT;
		dT_min = RKDP_MINDT;
	}
	break;
	}

	//initial settings
	dT_last = dT;

	mxh = 1;
	dmdt = 1;

	available = true;
	evalStep = 0;

	alternator = false;
	primed = false;

	calculate_mxh = true;
	calculate_dmdt = true;

#if COMPILECUDA == 1
	if (pODECUDA) pODECUDA->SyncODEValues();
#endif

	return error;
}

void ODECommon::Reset(void)
{
	iteration = 0;
	stageiteration = 0;
	time = 0;
	stagetime = 0;

	mxh = 1;
	dmdt = 1;

	available = true;
	evalStep = 0;

	alternator = false;
	primed = false;

	calculate_mxh = true;
	calculate_dmdt = true;

	moving_mesh_dwshift = 0.0;

#if COMPILECUDA == 1
	if (pODECUDA) pODECUDA->SyncODEValues();
#endif
}

void ODECommon::NewStage(void)
{
	stagetime = 0;
	stageiteration = 0;

	available = true;
	evalStep = 0;

	mxh = 1;
	dmdt = 1;

	alternator = false;
	primed = false;

	calculate_mxh = true;
	calculate_dmdt = true;

#if COMPILECUDA == 1
	if (pODECUDA) pODECUDA->SyncODEValues();
#endif
}

void ODECommon::SetdT(double dT) 
{
	this->dT = dT;

#if COMPILECUDA == 1
	if (pODECUDA) pODECUDA->SyncODEValues();
#endif
}

void ODECommon::SetAdaptiveTimeStepCtrl(double err_high_fail, double err_high, double err_low, double dT_increase, double dT_min, double dT_max)
{
	this->err_high_fail = err_high_fail;
	//decrease time step above this relative error - decrease factor depends on error ratio.
	this->err_high = err_high;
	//increase time step above this relative error
	this->err_low = err_low;
	//when increasing dT multiply by this (> 1.0)
	this->dT_increase = dT_increase;
	//maximum and minimum dT values allowed
	this->dT_min = dT_min;
	this->dT_max = dT_max;
}

void ODECommon::SetMoveMeshTrigger(bool status, int meshId)
{
	//first turn it off
	for (int idx = 0; idx < pODE.size(); idx++) {

		pODE[idx]->pMesh->SetMoveMeshTrigger(false);
	}

	moving_mesh = false;

	//set to trigger on a specified ferromagnetic mesh
	if (meshId >= 0 && status) {

		for (int idx = 0; idx < pODE.size(); idx++) {

			if (pODE[idx]->pMesh->get_id() == meshId) {

				pODE[idx]->pMesh->SetMoveMeshTrigger(true);
				moving_mesh = true;
				break;
			}
		}
	}
	//set to trigger on first ferromagnetic mesh
	else if (status) {

		if (pODE.size()) {

			pODE[0]->pMesh->SetMoveMeshTrigger(true);
			moving_mesh = true;
		}
	}
}

//---------------------------------------- GET / SET METHODS

void ODECommon::Set_mxh(void)
{
	//set mxh as the maximum values from all the set meshes
	if (pODE.size()) mxh = pODE[0]->mxh_reduction.max;

	for (int idx = 1; idx < pODE.size(); idx++) {

		if (pODE[idx]->mxh_reduction.max > mxh) mxh = pODE[idx]->mxh_reduction.max;
	}
}

void ODECommon::Set_dmdt(void)
{
	//set dmdt as the maximum values from all the set meshes
	if (pODE.size()) dmdt = pODE[0]->dmdt_reduction.max;

	for (int idx = 1; idx < pODE.size(); idx++) {

		if (pODE[idx]->dmdt_reduction.max > dmdt) dmdt = pODE[idx]->dmdt_reduction.max;
	}
}

double ODECommon::Get_lte(void)
{
	double lte = 0.0;

	//set lte as the maximum values from all the set meshes
	if (pODE.size()) lte = pODE[0]->lte_reduction.max;

	for (int idx = 1; idx < pODE.size(); idx++) {

		if (pODE[idx]->lte_reduction.max > lte) lte = pODE[idx]->lte_reduction.max;
	}

	return lte;
}

int ODECommon::GetId_of_MoveMeshTrigger(void)
{
	if (!moving_mesh) return -1;

	for (int idx = 0; idx < pODE.size(); idx++) {

		if (pODE[idx]->pMesh->GetMoveMeshTrigger())
			return pODE[idx]->pMesh->get_id();
	}

	return -1;
}

double ODECommon::Get_mxh(void)
{
	calculate_mxh = true;

#if COMPILECUDA == 1
	if (pODECUDA) return pODECUDA->Get_mxh();
#endif

	return mxh;
}

double ODECommon::Get_dmdt(void)
{
	calculate_dmdt = true;

#if COMPILECUDA == 1
	if (pODECUDA) return pODECUDA->Get_dmdt();
#endif

	return dmdt;
}