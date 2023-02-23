#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh.h"

///////////////////////////////////////////////////////////////////////////////

DifferentialEquation::DifferentialEquation(Mesh *pMesh) :
	ODECommon(true),
	prng(GetSystemTickCount())
{
	//store this pointer in the pODE vector and get a unique id so it can be erased later in the destructor
	odeId.minor = pODE.push_back(this, odeId.major);

	this->pMesh = pMesh;
}

DifferentialEquation::~DifferentialEquation() 
{
	//erase the pODE entry
	if(pODE.is_id_set(odeId))
		pODE.erase(odeId);

	//if this mesh (ferromagnetic mesh) had the moving mesh trigger set, then it is now gone so must set the moving_mesh trigger to false in the ODE common too.
	if (pMesh->GetMoveMeshTrigger()) moving_mesh = false;

#if COMPILECUDA == 1
	if (pmeshODECUDA) {

		called_from_destructor = true;
		delete pmeshODECUDA;
	}
	pmeshODECUDA = nullptr;
#endif
}