#include "stdafx.h"
#include "MeshBase.h"
#include "SuperMesh.h"

//set cone_angle value depending on required traget acceptance rate (acceptance_rate) and current value of mc_acceptance_rate
void MeshBase::MonteCarlo_AdaptiveAngle(double& cone_angle, double acceptance_rate)
{
	if (acceptance_rate < 0 || acceptance_rate > 1) acceptance_rate = MONTECARLO_TARGETACCEPTANCE;

	//adaptive cone angle - nice and simple, efficient for all temperatures; much better than MCM fixed cone angle variants, or cone angle set using a formula.
	if (mc_acceptance_rate < acceptance_rate) {

		//acceptance probability too low : decrease cone angle
		cone_angle -= MONTECARLO_CONEANGLEDEG_DELTA;
		if (cone_angle < pSMesh->Get_MonteCarlo_ConeAngleLimits().i) cone_angle = pSMesh->Get_MonteCarlo_ConeAngleLimits().i;
	}
	else if (mc_acceptance_rate > acceptance_rate) {

		//acceptance probability too high : increase cone angle
		cone_angle += MONTECARLO_CONEANGLEDEG_DELTA;
		if (cone_angle > pSMesh->Get_MonteCarlo_ConeAngleLimits().j) cone_angle = pSMesh->Get_MonteCarlo_ConeAngleLimits().j;
	}
}

void MeshBase::Set_MonteCarlo_Constrained(DBL3 cmc_n_)
{
	if (cmc_n_.IsNull()) mc_constrain = false;
	else { mc_constrain = true; cmc_n = cmc_n_; }

#if COMPILECUDA == 1
	if (pMeshBaseCUDA) pMeshBaseCUDA->Set_MonteCarlo_Constrained(cmc_n);
#endif
}