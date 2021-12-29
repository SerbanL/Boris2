#include "stdafx.h"
#include "Atom_Transport.h"

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"

#if COMPILECUDA == 1
#include "Atom_TransportCUDA.h"
#endif

//-------------------Calculation Methods

DBL2 Atom_Transport::IterateChargeSolver_SOR(double damping)
{
	return paMesh->V.IteratePoisson_SOR<Atom_Transport>(&Atom_Transport::Evaluate_ChargeSolver_delsqV_RHS, *this, damping);
}

double Atom_Transport::Evaluate_ChargeSolver_delsqV_RHS(int idx) const
{
	//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma
	//The Poisson solver calls this method to evaluate the RHS of this equation
	return -(paMesh->V.grad_diri(idx) * paMesh->elC.grad_sided(idx)) / paMesh->elC[idx];
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of a potential and flux

//Charge transport only : V

//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
double Atom_Transport::afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

double Atom_Transport::afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

double Atom_Transport::bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return -(1.5 * paMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * paMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

double Atom_Transport::bfunc_V_pri(int cell1_idx, int cell2_idx) const
{
	return -(1.5 * paMesh->elC[cell1_idx] - 0.5 * paMesh->elC[cell2_idx]);
}

//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
double Atom_Transport::diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = paMesh->V.position_to_cellidx(relpos_m1);

	return Evaluate_ChargeSolver_delsqV_RHS(cellm1_idx);
}

double Atom_Transport::diff2_V_pri(int cell1_idx, DBL3 shift) const
{
	return Evaluate_ChargeSolver_delsqV_RHS(cell1_idx);
}

#endif