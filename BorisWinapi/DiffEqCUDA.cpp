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

#endif