#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//------------------------------------------------------------------------------------------------------ THERMAL VECs GENERATIONS

//NOT APPLICABLE FOR DIAMAGNETIC MESHES

void DifferentialEquationDM::GenerateThermalField(void)
{
}

void DifferentialEquationDM::GenerateThermalField_and_Torque(void)
{
}


//------------------------------------------------------------------------------------------------------ STOCHASTIC EQUATIONS

DBL3 DifferentialEquationDM::SLLG(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::SLLGSTT(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::SLLB(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::SLLBSTT(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

#endif