#include "stdafx.h"
#include "TMR.h"

#ifdef MODULE_COMPILATION_TMR

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TMRCUDA.h"
#endif

//-------------------Calculation Methods used by Spin Current Solver only

//---------------Charge part

DBL2 TMR::IterateSpinSolver_Charge_SOR(double damping)
{
	return pMesh->V.IteratePoisson_SOR<TMR>(&TMR::Evaluate_SpinSolver_delsqV_RHS, *this, damping);
}

//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes. Not needed for TMR.
void TMR::PrimeSpinSolver_Charge(void)
{
}

//call-back method for Poisson equation for V. For TMR this is also a Laplace equation, so return 0 here.
//VERIFIED - CORRECT
double TMR::Evaluate_SpinSolver_delsqV_RHS(int idx) const
{
	return -(pMesh->V.grad_diri(idx) * pMesh->elC.grad_sided(idx)) / pMesh->elC[idx];
}

//---------------Spin part

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
DBL2 TMR::IterateSpinSolver_Spin_SOR(double damping)
{
	return pMesh->S.IteratePoisson_SOR<TMR>(&TMR::Evaluate_SpinSolver_delsqS_RHS, *this, damping);
}

//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes. Not needed for TMR.
void TMR::PrimeSpinSolver_Spin(void)
{
}

//call-back method for Poisson equation for S
//VERIFIED - CORRECT
DBL3 TMR::Evaluate_SpinSolver_delsqS_RHS(int idx) const
{
	DBL3 delsq_S_RHS;

	double l_sf = pMesh->l_sf;
	double elecCond = pMesh->elecCond;
	pMesh->update_parameters_ecoarse(idx, pMesh->l_sf, l_sf, pMesh->elecCond, elecCond);

	//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes
	//longitudinal S decay term
	//Use elecCond to mark out metallic pinholes.
	//tunelling : l_sf tends to infinity
	if (elecCond > 0.0) delsq_S_RHS = (pMesh->S[idx] / (l_sf * l_sf));
	
	return delsq_S_RHS;
}

//-------------------

//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V. Not needed for TMR.
//VERIFIED - CORRECT
DBL3 TMR::NHNeumann_Vdiff(int idx) const
{
	return DBL3();
}

//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S. Not needed for TMR.
//VERIFIED - CORRECT
DBL33 TMR::NHNeumann_Sdiff(int idx) const
{
	return DBL33();
}

//-------------------CMBND computation methods

//Spin transport : V and S
//CMBND for V
//flux = a + b V' at the interface, b = -sigma, a = 0 for TMR
//VERIFIED - CORRECT
double TMR::afunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

//flux = a + b V' at the interface, b = -sigma, a = 0 for TMR
//VERIFIED - CORRECT
double TMR::afunc_st_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

//VERIFIED - CORRECT
double TMR::bfunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return -(1.5 * pMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * pMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

//VERIFIED - CORRECT
double TMR::bfunc_st_V_pri(int cell1_idx, int cell2_idx) const
{
	return -(1.5 * pMesh->elC[cell1_idx] - 0.5 * pMesh->elC[cell2_idx]);
}

//second order differential of V along the shift axis
//this is simply Evaluate_SpinSolver_delsqV_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
double TMR::diff2_st_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = pMesh->V.position_to_cellidx(relpos_m1);
	return diff2_st_V_pri(cellm1_idx, shift);
}

//second order differential of V along the shift axis
//this is simply Evaluate_SpinSolver_delsqV_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
double TMR::diff2_st_V_pri(int cell1_idx, DBL3 shift) const
{
	//normalized, positive shift: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
	DBL3 nshift = mod(normalize(shift));

	return -((pMesh->V.grad_diri(cell1_idx) * nshift) * (pMesh->elC.grad_sided(cell1_idx) * nshift)) / pMesh->elC[cell1_idx];
}

//CMBND for S
//flux = a + b S' at the interface, b = -De, a = 0 for TMR
//VERIFIED - CORRECT
DBL3 TMR::afunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

//flux = a + b S' at the interface, b = -De, a = 0 for TMR
//VERIFIED - CORRECT
DBL3 TMR::afunc_st_S_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

//flux = a + b S' at the interface, b = -De
//VERIFIED - CORRECT
double TMR::bfunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	double De = pMesh->De;
	double elecCond = pMesh->elecCond;
	pMesh->update_parameters_atposition(relpos_m1, pMesh->De, De, pMesh->elecCond, elecCond);

	//metallic pinholes : use De
	if (elecCond > 0.0) return -De;
	//tunelling : use 1
	else return -1;
}

//flux = a + b S' at the interface, b = -De
//VERIFIED - CORRECT
double TMR::bfunc_st_S_pri(int cell1_idx, int cell2_idx) const
{
	double De = pMesh->De;
	double elecCond = pMesh->elecCond;
	pMesh->update_parameters_ecoarse(cell1_idx, pMesh->De, De, pMesh->elecCond, elecCond);

	//metallic pinholes : use De
	if (elecCond > 0.0) return -De;
	//tunelling : use 1
	else return -1;
}

//second order differential of S along the shift axis
//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
DBL3 TMR::diff2_st_S_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = pMesh->S.position_to_cellidx(relpos_m1);
	return diff2_st_S_pri(cellm1_idx, shift);
}

//second order differential of S along the shift axis
//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
DBL3 TMR::diff2_st_S_pri(int cell1_idx, DBL3 shift) const
{
	DBL3 delsq_S_RHS;

	double l_sf = pMesh->l_sf;
	double elecCond = pMesh->elecCond;
	pMesh->update_parameters_ecoarse(cell1_idx, pMesh->l_sf, l_sf, pMesh->elecCond, elecCond);

	//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes
	//longitudinal S decay term
	//Use elecCond to mark out metallic pinholes.
	//tunelling : l_sf tends to infinity
	if (elecCond > 0.0) delsq_S_RHS = (pMesh->S[cell1_idx] / (l_sf * l_sf));

	return delsq_S_RHS;
}

//-------------------

//VERIFIED - CORRECT
double TMR::cfunc_sec(DBL3 relpos, DBL3 stencil) const
{
	double De = pMesh->De;
	double elecCond = pMesh->elecCond;
	pMesh->update_parameters_atposition(relpos, pMesh->De, De, pMesh->elecCond, elecCond);

	//metallic pinholes : use De
	if (elecCond > 0.0) return De / (pMesh->elC.weighted_average(relpos, stencil) * MUB_E);
	//tunelling : use 1
	else return 1.0 / (pMesh->elC.weighted_average(relpos, stencil) * MUB_E);
}

//VERIFIED - CORRECT
double TMR::cfunc_pri(int cell_idx) const
{
	double De = pMesh->De;
	double elecCond = pMesh->elecCond;
	pMesh->update_parameters_ecoarse(cell_idx, pMesh->De, De, pMesh->elecCond, elecCond);

	//metallic pinholes : use De
	if (elecCond > 0.0) return De / (pMesh->elC[cell_idx] * MUB_E);
	//tunelling : use 1
	else return 1.0 / (pMesh->elC[cell_idx] * MUB_E);
}

//-------------------

//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations. Not used for TMR module.
//VERIFIED - CORRECT
void TMR::CalculateSAField(void)
{
}

//Calculate interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set). Not used for TMR module.
//VERIFIED - CORRECT
void TMR::CalculateSAInterfaceField(TransportBase* ptrans_sec, CMBNDInfo& contact)
{
}

#endif