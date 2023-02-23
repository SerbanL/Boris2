#include "stdafx.h"
#include "Atom_DiffEq.h"
#include "Atom_Mesh.h"

///////////////////////////////////////////////////////////////////////////////

Atom_DifferentialEquation::Atom_DifferentialEquation(Atom_Mesh *paMesh) :
	Atom_ODECommon(true),
	prng(GetSystemTickCount())
{
	//store this pointer in the pODE vector and get a unique id so it can be erased later in the destructor
	odeId.minor = pODE.push_back(this, odeId.major);

	this->paMesh = paMesh;
}

Atom_DifferentialEquation::~Atom_DifferentialEquation()
{
	//erase the pODE entry
	if(pODE.is_id_set(odeId))
		pODE.erase(odeId);

	//if this mesh had the moving mesh trigger set, then it is now gone so must set the moving_mesh trigger to false in the ODE common too.
	if (paMesh->GetMoveMeshTrigger()) moving_mesh = false;

#if COMPILECUDA == 1
	if (pameshODECUDA) {

		called_from_destructor = true;
		delete pameshODECUDA;
	}
	pameshODECUDA = nullptr;
#endif
}