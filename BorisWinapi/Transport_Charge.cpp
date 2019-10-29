#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

//-------------------Calculation Methods

DBL2 Transport::IterateChargeSolver_SOR(double damping)
{
	return pMesh->V.IteratePoisson_SOR<Transport>(&Transport::Evaluate_ChargeSolver_delsqV_RHS, *this, damping);
}

double Transport::Evaluate_ChargeSolver_delsqV_RHS(int idx) const
{
	//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma
	//The Poisson solver calls this method to evaluate the RHS of this equation
	return -(pMesh->V.grad_diri(idx) * pMesh->elC.grad_sided(idx)) / pMesh->elC[idx];
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of a potential and flux

//Charge transport only : V

//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
double Transport::afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

double Transport::afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

double Transport::bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return -(1.5 * pMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * pMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

double Transport::bfunc_V_pri(int cell1_idx, int cell2_idx) const
{
	return -(1.5 * pMesh->elC[cell1_idx] - 0.5 * pMesh->elC[cell2_idx]);
}

//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
double Transport::diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = pMesh->V.position_to_cellidx(relpos_m1);

	return Evaluate_ChargeSolver_delsqV_RHS(cellm1_idx);
}

double Transport::diff2_V_pri(int cell1_idx, DBL3 shift) const
{
	return Evaluate_ChargeSolver_delsqV_RHS(cell1_idx);
}

#endif