#include "stdafx.h"
#include "DiffEqCUDA.h"
#include "DiffEq.h"

#include "Mesh.h"
#include "MeshCUDA.h"

#if COMPILECUDA == 1

DifferentialEquationCUDA::DifferentialEquationCUDA(DifferentialEquation *pmeshODE_) :
	ODECommonCUDA()
{
	pmeshODE = pmeshODE_;
	
	pMesh = pmeshODE->pMesh;
	pMeshCUDA = pmeshODE->pMesh->pMeshCUDA;
}

//called when using stochastic equations
void DifferentialEquationCUDA::GenerateThermalField(void)
{
	//if not in linked dTstoch mode, then only generate stochastic field at a minimum of dTstoch spacing
	if (!pmeshODE->link_dTstoch && pmeshODE->GetTime() < pmeshODE->time_stoch + pmeshODE->dTstoch) return;

	if (!pmeshODE->link_dTstoch) {

		double deltaT = pmeshODE->GetTime() - pmeshODE->time_stoch;
		pmeshODE->time_stoch = pmeshODE->GetTime();

		deltaTstoch.from_cpu(deltaT);

		GenerateThermalField_CUDA(deltaTstoch);
	}
	else GenerateThermalField_CUDA(*pdT);
}

void DifferentialEquationCUDA::GenerateThermalField_and_Torque(void)
{
	//if not in linked dTstoch mode, then only generate stochastic field at a minimum of dTstoch spacing
	if (!pmeshODE->link_dTstoch && pmeshODE->GetTime() < pmeshODE->time_stoch + pmeshODE->dTstoch) return;

	if (!pmeshODE->link_dTstoch) {

		double deltaT = pmeshODE->GetTime() - pmeshODE->time_stoch;
		pmeshODE->time_stoch = pmeshODE->GetTime();

		deltaTstoch.from_cpu(deltaT);

		GenerateThermalField_and_Torque_CUDA(deltaTstoch);
	}
	else GenerateThermalField_and_Torque_CUDA(*pdT);
}

#endif