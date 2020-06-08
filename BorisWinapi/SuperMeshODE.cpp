#include "stdafx.h"
#include "SuperMesh.h"

//----------------------------------- ODE SOLVER CONTROL

//Reset all ODE solvers in meshes with on ODE
void SuperMesh::ResetODE(void) 
{ 
	odeSolver.Reset();
}

//set new stage for ODE solvers
void SuperMesh::NewStageODE(void)
{
	//if setting new stage first clear all text equation objects - if this stage uses a text equation then it will be set elsewhere after this call
	UpdateConfiguration_Values(UPDATECONFIG_TEQUATION_CLEAR);
	
	//new stage in ode
	odeSolver.NewStage();
}

//set the ode and evaluation method. Any new ODE in a magnetic mesh will use these settings. Currently Micromagnetic and Atomistic ODEs use the same evaluation method.
BError SuperMesh::SetODE(ODE_ setOde, EVAL_ evalMethod)
{
	BError error(__FUNCTION__);

	if (setOde <= ODE_ERROR || evalMethod <= EVAL_ERROR) return error(BERROR_INCORRECTNAME);

	error = odeSolver.SetODE(setOde, evalMethod);

	//now get currently set atomistic equation, so we can change its evaluation method to the one used by the micromagnetics ODE
	atom_odeSolver.QueryODE(setOde);
	error = atom_odeSolver.SetODE(setOde, evalMethod);

	error = UpdateConfiguration(UPDATECONFIG_ODE_SOLVER);

	return error;
}

//same for the atomistic ODE. Currently Micromagnetic and Atomistic ODEs use the same evaluation method.
BError SuperMesh::SetAtomisticODE(ODE_ setOde, EVAL_ evalMethod)
{
	BError error(__FUNCTION__);

	if (setOde <= ODE_ERROR || evalMethod <= EVAL_ERROR) return error(BERROR_INCORRECTNAME);

	error = atom_odeSolver.SetODE(setOde, evalMethod);

	//now get currently set micromagnetics equation, so we can change its evaluation method to the one used by the atomistic ODE
	odeSolver.QueryODE(setOde);
	error = odeSolver.SetODE(setOde, evalMethod);

	error = UpdateConfiguration(UPDATECONFIG_ODE_SOLVER);

	return error;
}

void SuperMesh::SetEvaluationSpeedup(int status) 
{ 
	odeSolver.SetEvaluationSpeedup(status); 
	
	UpdateConfiguration(UPDATECONFIG_ODE_SOLVER); 
}

//set the time step for the magnetisation solver
void SuperMesh::SetTimeStep(double dT) 
{ 
	odeSolver.SetdT(dT); 
}

//set parameters for adaptive time step control
void SuperMesh::SetAdaptiveTimeStepCtrl(double err_fail, double err_high, double err_low, double dT_incr, double dT_min, double dT_max) 
{ 
	odeSolver.SetAdaptiveTimeStepCtrl(err_fail, err_high, err_low, dT_incr, dT_min, dT_max); 
}

void SuperMesh::SetStochTimeStep(double dTstoch) 
{ 
	odeSolver.SetdTstoch(dTstoch); 
}

double SuperMesh::GetStochTimeStep(void) 
{ 
	return odeSolver.GetStochTimeStep(); 
}

void SuperMesh::SetLinkdTStochastic(bool flag) 
{ 
	odeSolver.SetLink_dTstoch(flag);;
}

bool SuperMesh::GetLink_dTstoch(void) 
{ 
	return odeSolver.GetLink_dTstoch(); 
}

//is the current time step fully finished? - most evaluation schemes need multiple sub-steps
bool SuperMesh::CurrentTimeStepSolved(void) 
{ 
	return odeSolver.TimeStepSolved(); 
}

//check evaluation speedup settings in ode solver
int SuperMesh::EvaluationSpeedup(void) 
{ 
	//assumed to be the same for the atomistic one
	return odeSolver.EvaluationSpeedup();
}

//check in ODECommon the type of field update we need to do depending on the ODE evaluation step
int SuperMesh::Check_Step_Update(void) 
{ 
	return odeSolver.Check_Step_Update();
}

//check if ODE solver needs spin accumulation solved
bool SuperMesh::SolveSpinCurrent(void)
{ 
	return odeSolver.SolveSpinCurrent();
}

void SuperMesh::SetMoveMeshAntisymmetric(bool antisymmetric) 
{ 
	odeSolver.SetMoveMeshAntisymmetric(antisymmetric); 
}

void SuperMesh::SetMoveMeshThreshold(double threshold) 
{ 
	odeSolver.SetMoveMeshThreshold(threshold); 
}

void SuperMesh::SetMoveMeshTrigger(bool status, string meshName)
{
	//if meshName not contained then we want meshId to be -1 : this means we'll be setting trigger on first magnetic mesh
	if (contains(meshName)) {

		int meshId = pMesh[meshName]->get_id();

		odeSolver.SetMoveMeshTrigger(status, meshId);
	}
	else {

		//if calling with meshId = -1 then first magnetic mesh will be used
		odeSolver.SetMoveMeshTrigger(status, -1);
	}

	UpdateConfiguration(UPDATECONFIG_ODE_MOVEMESH);
}

//---Other ODE Getters

//get set ODE and evaluation method
void SuperMesh::QueryODE(ODE_ &setODE, EVAL_ &evalMethod) 
{ 
	odeSolver.QueryODE(setODE, evalMethod); 
}

void SuperMesh::QueryODE(ODE_ &setODE) 
{ 
	odeSolver.QueryODE(setODE); 
}

//get set atomistic ODE and evaluation method
void SuperMesh::QueryAtomODE(ODE_ &setODE, EVAL_ &evalMethod)
{
	atom_odeSolver.QueryODE(setODE, evalMethod);
}

void SuperMesh::QueryAtomODE(ODE_ &setODE)
{
	atom_odeSolver.QueryODE(setODE);
}

int SuperMesh::GetIteration(void) 
{ 
	return odeSolver.GetIteration();
}

int SuperMesh::GetStageIteration(void) 
{ 
	return odeSolver.GetStageIteration();
}

double SuperMesh::GetTime(void) 
{ 
	return odeSolver.GetTime();
}

double SuperMesh::GetStageTime(void) 
{ 
	return odeSolver.GetStageTime();
}

double SuperMesh::GetTimeStep(void) 
{ 
	return odeSolver.GetTimeStep();
}

double SuperMesh::Get_mxh(void) 
{ 
	return odeSolver.Get_mxh();
}

double SuperMesh::Get_dmdt(void) 
{ 
	return odeSolver.Get_dmdt(); 
}

DBL3 SuperMesh::Get_AStepRelErrCtrl(void) 
{ 
	return odeSolver.Get_AStepRelErrCtrl();
}

DBL3 SuperMesh::Get_AStepdTCtrl(void) 
{ 
	return odeSolver.Get_AStepdTCtrl();
}

bool SuperMesh::IsMovingMeshSet(void) 
{ 
	return odeSolver.IsMovingMeshSet();
}

int SuperMesh::GetId_of_MoveMeshTrigger(void) 
{ 
	return odeSolver.GetId_of_MoveMeshTrigger(); 
}

double SuperMesh::Get_dwshift(void) 
{ 
	return odeSolver.Get_dwshift();
}

bool SuperMesh::MoveMeshAntisymmetric(void) 
{ 
	return odeSolver.MoveMeshAntisymmetric(); 
}

double SuperMesh::MoveMeshThreshold(void) 
{ 
	return odeSolver.MoveMeshThreshold();
}