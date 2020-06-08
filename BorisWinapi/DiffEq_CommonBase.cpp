#include "stdafx.h"
#include "DiffEq_CommonBase.h"

#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

#include "Mesh.h"
#include "Atom_Mesh.h"

///////////////////////////////////////////////////////////////////////
// Static variables in ODECommon_Base

//-----------------------------------Pointers

ODECommon* ODECommon_Base::podeSolver = nullptr;
Atom_ODECommon* ODECommon_Base::patom_odeSolver = nullptr;

//-----------------------------------Primary Data

int ODECommon_Base::iteration = 0;
int ODECommon_Base::stageiteration = 0;

double ODECommon_Base::time = 0.0;
double ODECommon_Base::stagetime = 0.0;

//-----------------------------------Time step

double ODECommon_Base::dT = 0.0;
double ODECommon_Base::dT_last = 0.0;

double ODECommon_Base::dTstoch = 0.0;
double ODECommon_Base::time_stoch = 0.0;
bool ODECommon_Base::link_dTstoch = true;

//-----------------------------------Evaluation Method Data

bool ODECommon_Base::available = true;

int ODECommon_Base::evalMethod = EVAL_EULER;

int ODECommon_Base::evalStep = 0;

int ODECommon_Base::use_evaluation_speedup = (int)EVALSPEEDUP_NONE;

//-----------------------------------mxh and dmdt

double ODECommon_Base::mxh = 1.0;
double ODECommon_Base::dmdt = 1.0;

bool ODECommon_Base::calculate_mxh = true;
bool ODECommon_Base::calculate_dmdt = true;

//-----------------------------------Adaptive time step control

double ODECommon_Base::lte = 0.0;

double ODECommon_Base::err_high_fail = RKF_RELERRFAIL;
double ODECommon_Base::err_high = RKF_RELERRMAX;
double ODECommon_Base::err_low = RKF_RELERRMIN;
double ODECommon_Base::dT_increase = RKF_DTINCREASE;
double ODECommon_Base::dT_max = RKF_MAXDT;
double ODECommon_Base::dT_min = RKF_MINDT;

//-----------------------------------Special values

bool ODECommon_Base::alternator = false;

bool ODECommon_Base::primed = false;

//-----------------------------------Moving mesh data

bool ODECommon_Base::moving_mesh = false;
bool ODECommon_Base::moving_mesh_antisymmetric = true;
double ODECommon_Base::moving_mesh_threshold = MOVEMESH_ANTISYMMETRIC_THRESHOLD;
double ODECommon_Base::moving_mesh_dwshift = 0.0;

//-----------------------------------Special Properties

bool ODECommon_Base::solve_spin_current = false;

//----------------------------------- Setup

//once ODECommon and Atom_ODECommon have been constructed, this function should be called to set pointers here (in SuperMesh constructor).
void ODECommon_Base::set_pointers(ODECommon& odeSolver, Atom_ODECommon& atom_odeSolver)
{
	podeSolver = &odeSolver;
	patom_odeSolver = &atom_odeSolver;
}

//----------------------------------- Important Control Methods

BError ODECommon_Base::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(ODECommon_Base));

	//make sure moving_mesh flag is set correctly
	moving_mesh = false;

	for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

		if (podeSolver->pODE[idx]->pMesh->GetMoveMeshTrigger()) {

			moving_mesh = true;
			break;
		}
	}

	if (!moving_mesh) {

		for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

			if (patom_odeSolver->pODE[idx]->paMesh->GetMoveMeshTrigger()) {

				moving_mesh = true;
				break;
			}
		}
	}

	//reset evaluation method priming and step
	evalStep = 0;
	alternator = false;
	primed = false;

#if COMPILECUDA == 1
	if (podeSolver->pODECUDA) podeSolver->pODECUDA->UpdateConfiguration(cfgMessage);
	if (patom_odeSolver->pODECUDA) patom_odeSolver->pODECUDA->UpdateConfiguration(cfgMessage);
#endif

	return error;
}

//switch CUDA state on/off
BError ODECommon_Base::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(ODECommon_Base));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		if (!podeSolver->pODECUDA) {

			podeSolver->pODECUDA = new ODECommonCUDA(podeSolver);
		}

		if (!patom_odeSolver->pODECUDA) {

			patom_odeSolver->pODECUDA = new Atom_ODECommonCUDA(patom_odeSolver);
		}
	}
	else {

		//cuda switched off so delete cuda module object
		if (podeSolver->pODECUDA) delete podeSolver->pODECUDA;
		podeSolver->pODECUDA = nullptr;

		if (patom_odeSolver->pODECUDA) delete patom_odeSolver->pODECUDA;
		patom_odeSolver->pODECUDA = nullptr;
	}

	primed = false;

#endif

	return error;
}