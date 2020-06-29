#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::LLG(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

DBL3 DifferentialEquationDM::LLGStatic(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::LLGSTT(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::LLB(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

//------------------------------------------------------------------------------------------------------

DBL3 DifferentialEquationDM::LLBSTT(int idx)
{
	double susrel = pMesh->susrel;
	pMesh->update_parameters_mcoarse(idx, pMesh->susrel, susrel);

	return susrel * pMesh->Heff[idx];
}

#endif