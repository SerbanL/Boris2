#include "stdafx.h"
#include "MeshBaseCUDA.h"
#include "MeshBase.h"
#include "SuperMesh.h"
#include "BorisLib.h"

#if COMPILECUDA == 1

MeshBaseCUDA::MeshBaseCUDA(MeshBase* pMeshBase_) :
	meshRect(pMeshBase_->meshRect),
	n(pMeshBase_->n), h(pMeshBase_->h),
	n_e(pMeshBase_->n_e), h_e(pMeshBase_->h_e),
	n_t(pMeshBase_->n_t), h_t(pMeshBase_->h_t),
	n_m(pMeshBase_->n_m), h_m(pMeshBase_->h_m)
{
	pMeshBase = pMeshBase_;
}

//----------------------------------- MESH INFO GET/SET METHODS

int MeshBaseCUDA::GetMeshType(void)
{
	return (int)pMeshBase->GetMeshType();
}

//----------------------------------- VALUE GETTERS

//check if the ODECommon::available flag is true (ode step solved)
bool MeshBaseCUDA::CurrentTimeStepSolved(void)
{
	return pMeshBase->pSMesh->CurrentTimeStepSolved();
}

//check evaluation speedup flag in ODECommon
int MeshBaseCUDA::GetEvaluationSpeedup(void)
{
	return pMeshBase->pSMesh->GetEvaluationSpeedup();
}

//check in ODECommon the type of field update we need to do depending on the ODE evaluation step
int MeshBaseCUDA::Check_Step_Update(void)
{
	return pMeshBase->pSMesh->Check_Step_Update();
}

cuBReal MeshBaseCUDA::GetStageTime(void)
{
	return pMeshBase->pSMesh->GetStageTime();
}

int MeshBaseCUDA::GetStageStep(void)
{
	return pMeshBase->pSMesh->stage_step.minor;
}

#endif