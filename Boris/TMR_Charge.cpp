#include "stdafx.h"
#include "TMR.h"

#ifdef MODULE_COMPILATION_TMR

#include "Mesh.h"
#include "SuperMesh.h"

#if COMPILECUDA == 1
#include "TMRCUDA.h"
#endif

//-------------------Calculation Methods

DBL2 TMR::IterateChargeSolver_SOR(double damping)
{
	return pMesh->V.IteratePoisson_SOR<TMR>(&TMR::Evaluate_ChargeSolver_delsqV_RHS, *this, damping);
}

double TMR::Evaluate_ChargeSolver_delsqV_RHS(int idx) const
{
	//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma
	//The Poisson solver calls this method to evaluate the RHS of this equation
	return -(pMesh->V.grad_diri(idx) * pMesh->elC.grad_sided(idx)) / pMesh->elC[idx];
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of a potential and flux

//Charge transport only : V

//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface.
double TMR::afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

double TMR::afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

double TMR::bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return -(1.5 * pMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * pMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

//Here sigma is constant along thickness, so no need to extrapolate.
double TMR::bfunc_V_pri(int cell1_idx, int cell2_idx) const
{
	return -(1.5 * pMesh->elC[cell1_idx] - 0.5 * pMesh->elC[cell2_idx]);
}

//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
double TMR::diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = pMesh->V.position_to_cellidx(relpos_m1);
	return diff2_V_pri(cellm1_idx, shift);
}

double TMR::diff2_V_pri(int cell1_idx, DBL3 shift) const
{
	//normalized, positive shift: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
	DBL3 nshift = mod(normalize(shift));

	return -((pMesh->V.grad_diri(cell1_idx) * nshift) * (pMesh->elC.grad_sided(cell1_idx) * nshift)) / pMesh->elC[cell1_idx];
}

#endif