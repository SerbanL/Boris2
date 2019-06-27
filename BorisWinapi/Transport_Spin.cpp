#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

//-------------------Calculation Methods used by Spin Current Solver only

//---------------Charge part

DBL2 Transport::IterateSpinSolver_Charge_aSOR(bool start_iters, double conv_error)
{
	if (IsZ((double)pMesh->iSHA) || pMesh->M.linear_size()) {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		return pMesh->V.IteratePoisson_aSOR<Transport>(&Transport::Evaluate_SpinSolver_delsqV_RHS, *this, start_iters, conv_error);
	}
	else {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		return pMesh->V.IteratePoisson_aSOR<Transport>(&Transport::Evaluate_SpinSolver_delsqV_RHS, &Transport::NHNeumann_Vdiff, *this, start_iters, conv_error);
	}
}

DBL2 Transport::IterateSpinSolver_Charge_SOR(double damping)
{
	if (IsZ((double)pMesh->iSHA) || pMesh->M.linear_size()) {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		return pMesh->V.IteratePoisson_SOR<Transport>(&Transport::Evaluate_SpinSolver_delsqV_RHS, *this, damping);
	}
	else {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		return pMesh->V.IteratePoisson_SOR<Transport>(&Transport::Evaluate_SpinSolver_delsqV_RHS, &Transport::NHNeumann_Vdiff, *this, damping);
	}
}

double Transport::Evaluate_SpinSolver_delsqV_RHS(int idx) const
{
	//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma + betaD*De*e * del ((grad S) m) / sigma * muB
	
	double De = pMesh->De;
	double iSHA = pMesh->iSHA;
	double betaD = pMesh->betaD;
	pMesh->update_parameters_ecoarse(idx, pMesh->De, De, pMesh->iSHA, iSHA, pMesh->betaD, betaD);

	//The Poisson solver calls this method to evaluate the RHS of this equation
	double value = 0.0;
	
	if (IsZ((double)iSHA) || pMesh->M.linear_size()) {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		value = -(pMesh->V.grad_diri(idx) * pMesh->elC.grad_sided(idx)) / pMesh->elC[idx];
	}
	else {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		value = -(pMesh->V.grad_diri_nneu(idx, (iSHA * De / (MUB_E * pMesh->elC[idx])) * pMesh->S.curl_neu(idx)) * pMesh->elC.grad_sided(idx)) / pMesh->elC[idx];
	}
	
	//CPP-GMR contribution in magnetic meshes
	if (pMesh->M.linear_size() && IsNZ((double)betaD)) {

		double Ms = pMesh->Ms;
		pMesh->update_parameters_ecoarse(idx, pMesh->Ms, Ms);

		int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(idx));

		DBL33 grad_M = pMesh->M.grad_neu(idx_M);
		DBL33 grad_S = pMesh->S.grad_neu(idx);
		DBL3 delsq_S = pMesh->S.delsq_neu(idx);
		double div_grad_S_M = (grad_S.i * grad_M.i) + (grad_S.j * grad_M.j) + (grad_S.k * grad_M.k) + (pMesh->M[idx_M] * delsq_S);

		value += div_grad_S_M * betaD * De / (MUB_E * pMesh->elC[idx] * Ms);
	}
	
	return value;
}

//---------------Spin part

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
DBL2 Transport::IterateSpinSolver_Spin_aSOR(bool start_iters, double conv_error)
{
	if (IsZ((double)pMesh->SHA) || pMesh->M.linear_size()) {

		//no SHE contribution. Note, SHE is not included in magnetic meshes.
		return pMesh->S.IteratePoisson_aSOR<Transport>(&Transport::Evaluate_SpinSolver_delsqS_RHS, *this, start_iters, conv_error);
	}
	else {

		//SHE enabled, must use non-homogeneous Neumann boundary condition for grad S
		return pMesh->S.IteratePoisson_aSOR<Transport>(&Transport::Evaluate_SpinSolver_delsqS_RHS, &Transport::NHNeumann_Sdiff, *this, start_iters, conv_error);
	}
}

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
DBL2 Transport::IterateSpinSolver_Spin_SOR(double damping)
{
	if (IsZ((double)pMesh->SHA) || pMesh->M.linear_size()) {

		//no SHE contribution. Note, SHE is not included in magnetic meshes.
		return pMesh->S.IteratePoisson_SOR<Transport>(&Transport::Evaluate_SpinSolver_delsqS_RHS, *this, damping);
	}
	else {

		//SHE enabled, must use non-homogeneous Neumann boundary condition for grad S
		return pMesh->S.IteratePoisson_SOR<Transport>(&Transport::Evaluate_SpinSolver_delsqS_RHS, &Transport::NHNeumann_Sdiff, *this, damping);
	}
}

//call-back method for Poisson equation for S
DBL3 Transport::Evaluate_SpinSolver_delsqS_RHS(int idx) const
{
	DBL3 delsq_S_RHS;

	double De = pMesh->De;
	double SHA = pMesh->SHA;
	double l_sf = pMesh->l_sf;
	pMesh->update_parameters_ecoarse(idx, pMesh->De, De, pMesh->SHA, SHA, pMesh->l_sf, l_sf);

	//Terms occuring only in ferromagnetic meshes (where SHA = 0)
	if (pMesh->M.linear_size()) {

		double Ms = pMesh->Ms;
		double P = pMesh->P;
		double betaD = pMesh->betaD;
		double l_ex = pMesh->l_ex;
		double l_ph = pMesh->l_ph;
		pMesh->update_parameters_ecoarse(idx, pMesh->Ms, Ms, pMesh->P, P, pMesh->betaD, betaD, pMesh->l_ex, l_ex, pMesh->l_ph, l_ph);

		//1. term due to non-uniformity of M

		//find grad M and M at the M cell in which the current S cell center is
		int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(idx));

		DBL3 Mval = pMesh->M[idx_M];
		DBL33 grad_M = pMesh->M.grad_neu(idx_M);
		DBL3 Jc_dot_del_M = grad_M | pMesh->Jc[idx];

		delsq_S_RHS -= P * MUB_E * Jc_dot_del_M / (Ms * De);

		//2. transverse S decay terms

		delsq_S_RHS += ((pMesh->S[idx] ^ Mval) / (Ms * l_ex * l_ex) + (Mval ^ (pMesh->S[idx] ^ Mval)) / (Ms * Ms * l_ph * l_ph));

		//3. CPP-GMR term
		if (IsNZ((double)betaD)) {
				
			DBL33 grad_S = pMesh->S.grad_neu(idx);
			DBL3 delsq_S = pMesh->S.delsq_neu(idx);

			DBL3 grad_S_M_dot_del_M = grad_M | (grad_S * Mval);

			delsq_S_RHS += (grad_S_M_dot_del_M + Mval * ((grad_M.i * grad_S.i) + (grad_M.j * grad_S.j) + (grad_M.k * grad_S.k) + (Mval * delsq_S))) * betaD * P / (Ms * Ms);
		}
	}
	//terms occuring only in non-ferromagnetic meshes (i.e. those with SHA not zero)
	else {

		//1. SHA term. 

		//Check boundary conditions for this term : should be Dirichlet with 0 Jc value normal to the boundary except for electrodes.
		delsq_S_RHS += (SHA * MUB_E / De) * pMesh->Jc.diveps3_sided(idx);
	}

	//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes

	//1. longitudinal S decay term

	delsq_S_RHS += (pMesh->S[idx] / (l_sf * l_sf));

	return delsq_S_RHS;
}

//-------------------

//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V
DBL3 Transport::NHNeumann_Vdiff(int idx) const
{
	if (pMesh->M.linear_size()) return DBL3();

	double iSHA = pMesh->iSHA;
	double De = pMesh->De;
	pMesh->update_parameters_ecoarse(idx, pMesh->iSHA, iSHA, pMesh->De, De);

	return (iSHA * De / (MUB_E * pMesh->elC[idx])) * pMesh->S.curl_neu(idx);
}

//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S
DBL33 Transport::NHNeumann_Sdiff(int idx) const
{
	if (pMesh->M.linear_size()) return DBL33();

	double SHA = pMesh->SHA;
	double De = pMesh->De;
	pMesh->update_parameters_ecoarse(idx, pMesh->SHA, SHA, pMesh->De, De);

	return epsilon3(pMesh->Jc[idx]) * (SHA * MUB_E / De);
}

//-------------------CMBND computation methods

//Spin transport : V and S
//CMBND for V
double Transport::afunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	double a = 0.0;

	DBL3 u = shift.normalized() * -1;

	int idx_S1 = pMesh->S.position_to_cellidx(relpos_m1);
	int idx_S2 = pMesh->S.position_to_cellidx(relpos_m1 + shift);

	double betaD = pMesh->betaD;
	double iSHA = pMesh->iSHA;
	pMesh->update_parameters_atposition(relpos_m1, pMesh->betaD, betaD, pMesh->iSHA, iSHA);

	//CPP-GMR term
	if (pMesh->M.linear_size() && IsNZ((double)betaD)) {

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->Ms, Ms, pMesh->De, De);

		//need to find value at boundary so use interpolation
		
		//value a1
		DBL3 M1 = pMesh->M.weighted_average(relpos_m1, stencil);
		DBL33 grad_S1 = pMesh->S.grad_neu(idx_S1);

		double a1 = ((grad_S1 * M1) * betaD * De / (MUB_E * Ms)) * u;

		//value a2
		DBL3 M2 = pMesh->M.weighted_average(relpos_m1 + shift, stencil);
		DBL33 grad_S2 = pMesh->S.grad_neu(idx_S2);

		double a2 = ((grad_S2 * M2) * betaD * De / (MUB_E * Ms)) * u;

		//final interpolated a value
		a = (1.5 * a1 - 0.5 * a2);
	}

	//iSHE contribution
	if (IsNZ((double)iSHA) && !pMesh->M.linear_size()) {

		double SHA = pMesh->SHA;
		double De = pMesh->De;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->SHA, SHA, pMesh->De, De);

		//need to find value at boundary so use interpolation

		//value a1
		DBL3 Jc1 = pMesh->Jc.weighted_average(relpos_m1, stencil);

		double a1 = (iSHA * De / MUB_E) * pMesh->S.curl_nneu(idx_S1, epsilon3(Jc1) * (SHA * MUB_E / De)) * u;

		//value a2
		DBL3 Jc2 = pMesh->Jc.weighted_average(relpos_m1 + shift, stencil);

		double a2 = (iSHA * De / MUB_E) * pMesh->S.curl_nneu(idx_S2, epsilon3(Jc2) * (SHA * MUB_E / De)) * u;

		//final interpolated a value
		a += (1.5 * a1 - 0.5 * a2);
	}

	return a;
}

double Transport::afunc_st_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	double a = 0.0;

	DBL3 u = shift.normalized() * -1;

	double betaD = pMesh->betaD;
	double iSHA = pMesh->iSHA;
	pMesh->update_parameters_ecoarse(cell1_idx, pMesh->betaD, betaD, pMesh->iSHA, iSHA);

	if (pMesh->M.linear_size() && IsNZ((double)betaD)) {

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->Ms, Ms, pMesh->De, De);

		//need to find value at boundary so use interpolation

		//value a1
		int idx_M1 = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(cell1_idx));
		DBL3 M1 = pMesh->M[idx_M1];
		DBL33 grad_S1 = pMesh->S.grad_neu(cell1_idx);

		double a1 = ((grad_S1 * M1) * betaD * De / (MUB_E * Ms)) * u;

		//value a2
		int idx_M2 = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(cell2_idx));
		DBL3 M2 = pMesh->M[idx_M2];
		DBL33 grad_S2 = pMesh->S.grad_neu(cell2_idx);

		double a2 = ((grad_S2 * M2) * betaD * De / (MUB_E * Ms)) * u;

		//final interpolated a value
		a = (1.5 * a1 - 0.5 * a2);
	}
	
	//iSHE contribution
	if (IsNZ((double)iSHA) && !pMesh->M.linear_size()) {

		double SHA = pMesh->SHA;
		double De = pMesh->De;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->SHA, SHA, pMesh->De, De);

		//need to find value at boundary so use interpolation

		//value a1
		DBL3 Jc1 = pMesh->Jc[cell1_idx];

		double a1 = (iSHA * De / MUB_E) * pMesh->S.curl_nneu(cell1_idx, epsilon3(Jc1) * (SHA * MUB_E / De)) * u;

		//value a2
		DBL3 Jc2 = pMesh->Jc[cell2_idx];

		double a2 = (iSHA * De / MUB_E) * pMesh->S.curl_nneu(cell2_idx, epsilon3(Jc2) * (SHA * MUB_E / De)) * u;

		//final interpolated a value
		a += (1.5 * a1 - 0.5 * a2);
	}

	return a;
}

double Transport::bfunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return -(1.5 * pMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * pMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

double Transport::bfunc_st_V_pri(int cell1_idx, int cell2_idx) const
{
	return -(1.5 * pMesh->elC[cell1_idx] - 0.5 * pMesh->elC[cell2_idx]);
}

double Transport::diff2_st_V_sec(DBL3 relpos_m1, DBL3 stencil) const
{
	double iSHA = pMesh->iSHA;
	double De = pMesh->De;
	double betaD = pMesh->betaD;
	pMesh->update_parameters_atposition(relpos_m1, pMesh->iSHA, iSHA, pMesh->De, De, pMesh->betaD, betaD);

	int cellm1_idx = pMesh->V.position_to_cellidx(relpos_m1);

	double value = 0.0;
	
	if (IsZ((double)iSHA) || pMesh->M.linear_size()) {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		value = -(pMesh->V.grad_diri(cellm1_idx) * pMesh->elC.grad_sided(cellm1_idx)) / pMesh->elC[cellm1_idx];
	}
	else {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		value = -(pMesh->V.grad_diri_nneu(cellm1_idx, (iSHA * De / (MUB_E * pMesh->elC[cellm1_idx])) * pMesh->S.curl_neu(cellm1_idx)) * pMesh->elC.grad_sided(cellm1_idx)) / pMesh->elC[cellm1_idx];
	}

	if (pMesh->M.linear_size() && IsNZ((double)betaD)) {

		double Ms = pMesh->Ms;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->Ms, Ms);

		int idx_M = pMesh->M.position_to_cellidx(relpos_m1);

		DBL33 grad_S = pMesh->S.grad_neu(cellm1_idx);
		DBL33 grad_M = pMesh->M.grad_neu(idx_M);
		DBL3 delsq_S = pMesh->S.delsq_neu(cellm1_idx);

		double div_grad_S_M = (grad_S.i * grad_M.i) + (grad_S.j * grad_M.j) + (grad_S.k * grad_M.k) + (pMesh->M[idx_M] * delsq_S);

		value += div_grad_S_M * betaD * De / (MUB_E * Ms * pMesh->elC[cellm1_idx]);
	}

	return value;
}

double Transport::diff2_st_V_pri(int cell1_idx) const
{
	double iSHA = pMesh->iSHA;
	double De = pMesh->De;
	double betaD = pMesh->betaD;
	pMesh->update_parameters_ecoarse(cell1_idx, pMesh->iSHA, iSHA, pMesh->De, De, pMesh->betaD, betaD);

	double value = 0.0;
	
	if (IsZ((double)iSHA) || pMesh->M.linear_size()) {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		value = -(pMesh->V.grad_diri(cell1_idx) * pMesh->elC.grad_sided(cell1_idx)) / pMesh->elC[cell1_idx];
	}
	else {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		value = -(pMesh->V.grad_diri_nneu(cell1_idx, (iSHA * De / (MUB_E * pMesh->elC[cell1_idx])) * pMesh->S.curl_neu(cell1_idx)) * pMesh->elC.grad_sided(cell1_idx)) / pMesh->elC[cell1_idx];
	}

	if (pMesh->M.linear_size() && IsNZ((double)betaD)) {

		double Ms = pMesh->Ms;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->Ms, Ms);

		int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(cell1_idx));

		DBL33 grad_S = pMesh->S.grad_neu(cell1_idx);
		DBL33 grad_M = pMesh->M.grad_neu(idx_M);
		DBL3 delsq_S = pMesh->S.delsq_neu(cell1_idx);

		double div_grad_S_M = (grad_S.i * grad_M.i) + (grad_S.j * grad_M.j) + (grad_S.k * grad_M.k) + (pMesh->M[idx_M] * delsq_S);

		value += div_grad_S_M * betaD * De / (MUB_E * Ms * pMesh->elC[cell1_idx]);
	}

	return value;
}

//CMBND for S
DBL3 Transport::afunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	//need to find value at boundary so use interpolation

	//unit vector perpendicular to interface (pointing from secondary to primary mesh)
	DBL3 u = shift.normalized() * -1;

	//values on secondary side
	if (pMesh->M.linear_size()) {

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		double P = pMesh->P;
		double betaD = pMesh->betaD;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->Ms, Ms, pMesh->De, De, pMesh->P, P, pMesh->betaD, betaD);

		//a1 value
		DBL3 M1 = pMesh->M.weighted_average(relpos_m1, stencil);
		DBL3 Jc1 = pMesh->Jc.weighted_average(relpos_m1, stencil);
		int idx_S1 = pMesh->S.position_to_cellidx(relpos_m1);
		DBL33 grad_S1 = pMesh->S.grad_neu(idx_S1);

		DBL3 a1 = 
			-MUB_E * ((Jc1 | M1) | u) * (P / Ms)
			+ ((((grad_S1 * M1) | M1) * betaD * P * De / (Ms * Ms)) | u);

		//a2 value
		DBL3 M2 = pMesh->M.weighted_average(relpos_m1 + shift, stencil);
		DBL3 Jc2 = pMesh->Jc.weighted_average(relpos_m1 + shift, stencil);
		int idx_S2 = pMesh->S.position_to_cellidx(relpos_m1 + shift);
		DBL33 grad_S2 = pMesh->S.grad_neu(idx_S2);

		DBL3 a2 = 
			-MUB_E * ((Jc2 | M2) | u) * (P / Ms)
			+ ((((grad_S2 * M2) | M2) * betaD * P * De / (Ms * Ms)) | u);

		return (1.5 * a1 - 0.5 * a2);
	}
	else {

		double SHA = pMesh->SHA;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->SHA, SHA);

		//non-magnetic mesh
		DBL3 a1 = (epsilon3(pMesh->Jc.weighted_average(relpos_m1, stencil)) | u) * SHA * MUB_E;
		DBL3 a2 = (epsilon3(pMesh->Jc.weighted_average(relpos_m1 + shift, stencil)) | u) * SHA * MUB_E;

		return (1.5 * a1 - 0.5 * a2);
	}
}

DBL3 Transport::afunc_st_S_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	//need to find value at boundary so use interpolation

	DBL3 u = shift.normalized() * -1;

	//values on primary side
	if (pMesh->M.linear_size()) {

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		double P = pMesh->P;
		double betaD = pMesh->betaD;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->Ms, Ms, pMesh->De, De, pMesh->P, P, pMesh->betaD, betaD);

		//magnetic mesh
		int idx_M1 = pMesh->M.position_to_cellidx(pMesh->Jc.cellidx_to_position(cell1_idx));
		DBL3 M1 = pMesh->M[idx_M1];
		DBL3 Jc1 = pMesh->Jc[cell1_idx];

		DBL3 a1 = ((Jc1 | M1) | u) * (P / Ms) * (-MUB_E);

		DBL33 grad_S1 = pMesh->S.grad_neu(cell1_idx);

		a1 += (((grad_S1 * M1) | M1) * betaD * P * De / (Ms * Ms)) | u;

		//magnetic mesh
		int idx_M2 = pMesh->M.position_to_cellidx(pMesh->Jc.cellidx_to_position(cell2_idx));
		DBL3 M2 = pMesh->M[idx_M2];
		DBL3 Jc2 = pMesh->Jc[cell2_idx];

		DBL3 a2 = ((Jc2 | M2) | u) * (P / Ms) * (-MUB_E);

		DBL33 grad_S2 = pMesh->S.grad_neu(cell2_idx);

		a2 += (((grad_S2 * M2) | M2) * betaD * P * De / (Ms * Ms)) | u;

		return (1.5 * a1 - 0.5 * a2);
	}
	else {

		double SHA = pMesh->SHA;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->SHA, SHA);

		//non-magnetic mesh
		DBL3 a1 = (epsilon3(pMesh->Jc[cell1_idx]) | u) * SHA * MUB_E;
		DBL3 a2 = (epsilon3(pMesh->Jc[cell2_idx]) | u) * SHA * MUB_E;

		return (1.5 * a1 - 0.5 * a2);
	}
}

double Transport::bfunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	double De = pMesh->De;
	pMesh->update_parameters_atposition(relpos_m1, pMesh->De, De);

	return -1.0 * De;
}

double Transport::bfunc_st_S_pri(int cell1_idx, int cell2_idx) const
{
	double De = pMesh->De;
	pMesh->update_parameters_ecoarse(cell1_idx, pMesh->De, De);

	return -1.0 * De;
}

DBL3 Transport::diff2_st_S_sec(DBL3 relpos_m1, DBL3 stencil) const
{
	if (pMesh->M.linear_size()) {

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		double P = pMesh->P;
		double betaD = pMesh->betaD;
		double l_ex = pMesh->l_ex;
		double l_ph = pMesh->l_ph;
		double l_sf = pMesh->l_sf;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->Ms, Ms, pMesh->De, De, pMesh->P, P, pMesh->betaD, betaD, pMesh->l_ex, l_ex, pMesh->l_ph, l_ph, pMesh->l_sf, l_sf);

		DBL3 value = pMesh->S.weighted_average(relpos_m1, stencil) / (l_sf * l_sf);

		//add to value for magnetic mesh

		DBL3 Sval = pMesh->S.weighted_average(relpos_m1, stencil);
		DBL3 Mval = pMesh->M.weighted_average(relpos_m1, stencil);

		value += (Sval ^ Mval) / (Ms * l_ex * l_ex) + (Mval ^ (Sval ^ Mval)) / (Ms * Ms * l_ph * l_ph);

		//find grad M and M at the M cell in which the current S cell center is
		int idx_M = pMesh->M.position_to_cellidx(relpos_m1);
		DBL33 grad_M = pMesh->M.grad_neu(idx_M);
		DBL3 Jc_dot_del_M = grad_M | pMesh->Jc.weighted_average(relpos_m1, stencil);

		value -= P * MUB_E * Jc_dot_del_M / (Ms * De);

		if (IsNZ((double)betaD)) {

			int idx_S = pMesh->S.position_to_cellidx(relpos_m1);
			DBL33 grad_S = pMesh->S.grad_neu(idx_S);
			DBL3 delsq_S = pMesh->S.delsq_neu(idx_S);

			DBL3 grad_S_M_dot_del_M = grad_M | (grad_S * Mval);

			value += (grad_S_M_dot_del_M + Mval * ((grad_M.i * grad_S.i) + (grad_M.j * grad_S.j) + (grad_M.k * grad_S.k) + (Mval * delsq_S))) * betaD * P / (Ms * Ms);
		}

		return value;
	}
	else {

		double SHA = pMesh->SHA;
		double De = pMesh->De;
		double l_sf = pMesh->l_sf;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->SHA, SHA, pMesh->De, De, pMesh->l_sf, l_sf);

		DBL3 value = pMesh->S.weighted_average(relpos_m1, stencil) / (l_sf * l_sf);

		//add to value for non-magnetic mesh
		int idx_Jc = pMesh->Jc.position_to_cellidx(relpos_m1);
		value += (SHA * MUB_E / De) * pMesh->Jc.diveps3_sided(idx_Jc);

		return value;
	}
}

DBL3 Transport::diff2_st_S_pri(int cell1_idx) const
{
	if (pMesh->M.linear_size()) {

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		double P = pMesh->P;
		double betaD = pMesh->betaD;
		double l_ex = pMesh->l_ex;
		double l_ph = pMesh->l_ph;
		double l_sf = pMesh->l_sf;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->Ms, Ms, pMesh->De, De, pMesh->P, P, pMesh->betaD, betaD, pMesh->l_ex, l_ex, pMesh->l_ph, l_ph, pMesh->l_sf, l_sf);

		DBL3 value = pMesh->S[cell1_idx] / (l_sf * l_sf);

		//add to value for magnetic mesh

		int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(cell1_idx));

		DBL3 Sval = pMesh->S[cell1_idx];
		DBL3 Mval = pMesh->M[idx_M];

		value += (Sval ^ Mval) / (Ms * l_ex * l_ex) + (Mval ^ (Sval ^ Mval)) / (Ms * Ms * l_ph * l_ph);

		//find grad M and M at the M cell in which the current S cell center is
		DBL33 grad_M = pMesh->M.grad_neu(idx_M);
		DBL3 Jc_dot_del_M = grad_M | pMesh->Jc[cell1_idx];

		value -= P * MUB_E * Jc_dot_del_M / (Ms * De);

		if (IsNZ((double)betaD)) {

			DBL33 grad_S = pMesh->S.grad_neu(cell1_idx);
			DBL3 delsq_S = pMesh->S.delsq_neu(cell1_idx);

			DBL3 grad_S_M_dot_del_M = grad_M | (grad_S * Mval);

			value += (grad_S_M_dot_del_M + Mval * ((grad_M.i * grad_S.i) + (grad_M.j * grad_S.j) + (grad_M.k * grad_S.k) + (Mval * delsq_S))) * betaD * P / (Ms * Ms);
		}

		return value;
	}
	else {

		double SHA = pMesh->SHA;
		double De = pMesh->De;
		double l_sf = pMesh->l_sf;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->SHA, SHA, pMesh->De, De, pMesh->l_sf, l_sf);

		DBL3 value = pMesh->S[cell1_idx] / (l_sf * l_sf);

		//add to value for non-magnetic mesh
		value += (SHA * MUB_E / De) * pMesh->Jc.diveps3_sided(cell1_idx);

		return value;
	}
}

//-------------------

double Transport::cfunc_sec(DBL3 relpos, DBL3 stencil) const
{
	double De = pMesh->De;
	pMesh->update_parameters_atposition(relpos, pMesh->De, De);

	return De / (pMesh->elC.weighted_average(relpos, stencil) * MUB_E);
}

double Transport::cfunc_pri(int cell_idx) const
{
	double De = pMesh->De;
	pMesh->update_parameters_ecoarse(cell_idx, pMesh->De, De);

	return De / (pMesh->elC[cell_idx] * MUB_E);
}

//-------------------

//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
void Transport::CalculateSAField(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->M.linear_size(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			double grel = pMesh->grel;
			double De = pMesh->De;
			double ts_eff = pMesh->ts_eff;
			double l_ex = pMesh->l_ex;
			double l_ph = pMesh->l_ph;
			pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->grel, grel, pMesh->De, De, pMesh->ts_eff, ts_eff, pMesh->l_ex, l_ex, pMesh->l_ph, l_ph);

			if (IsNZ((double)grel)) {

				DBL3 Sa = pMesh->S.weighted_average(pMesh->M.cellidx_to_position(idx), pMesh->h);

				pMesh->Heff[idx] += (De * ts_eff / (GAMMA * grel * Ms)) * (Sa / (l_ex * l_ex) + (pMesh->M[idx] ^ Sa) / (l_ph * l_ph * Ms));
			}
		}
	}
}

//Calculate interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
void Transport::CalculateSAInterfaceField(Transport* ptrans_sec, CMBNDInfo& contact)
{
	//the top contacting mesh sets G values
	bool GInterface_Enabled = ((contact.IsPrimaryTop() && pMesh->GInterface_Enabled()) || (!contact.IsPrimaryTop() && ptrans_sec->pMesh->GInterface_Enabled()));

	if (pMesh->MComputation_Enabled() && !ptrans_sec->pMesh->Magnetisation_Enabled() && GInterface_Enabled) {

		//interface conductance method with F being the primary mesh : calculate and set spin torque

		//convert the cells box from S mesh to M mesh
		INT3 mbox_start = pMesh->M.cellidx_from_position(pMesh->S.cellidx_to_position(contact.cells_box.s) + pMesh->meshRect.s);
		INT3 mbox_end = pMesh->M.cellidx_from_position(pMesh->S.cellidx_to_position(contact.cells_box.e - INT3(1)) + pMesh->meshRect.s);

		if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
		if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
		if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

		INT3 box_sizes = mbox_end - mbox_start;

		//the cellsize perpendicular to the contact (in the M mesh)
		double dh = (DBL3(contact.cell_shift) & pMesh->h).norm();

		Mesh* pMesh_sec = ptrans_sec->pMesh;

		//primary cells in this contact
#pragma omp parallel for
		for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

			int i = (box_idx % box_sizes.x) + mbox_start.i;
			int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
			int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;
			
			//index of magnetic cell 1
			int mcell1_idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

			if (pMesh->M.is_empty(mcell1_idx)) continue;
			
			double Ms = pMesh->Ms;
			double grel = pMesh->grel;
			double tsi_eff = pMesh->tsi_eff;
			pMesh->update_parameters_mcoarse(mcell1_idx, pMesh->grel, grel, pMesh->Ms, Ms, pMesh->tsi_eff, tsi_eff);

			if (IsNZ((double)grel)) {

				//position at interface relative to primary mesh
				DBL3 mhshift_primary = contact.hshift_primary.normalized() & pMesh->h;
				DBL3 relpos_interf = ((DBL3(i, j, k) + DBL3(0.5)) & pMesh->h) + mhshift_primary / 2;

				DBL3 relpos_1 = relpos_interf - contact.hshift_primary / 2;

				DBL3 relpos_m1 = pMesh->meshRect.s - pMesh_sec->meshRect.s + relpos_interf + contact.hshift_secondary / 2;

				DBL3 stencil = pMesh->h - mod(mhshift_primary) + mod(contact.hshift_secondary);

				//S values
				DBL3 S_1 = pMesh->S.weighted_average(relpos_1, stencil);
				DBL3 S_2 = pMesh->S.weighted_average(relpos_1 - contact.hshift_primary, stencil);
				DBL3 S_m1 = pMesh_sec->S.weighted_average(relpos_m1, stencil);
				DBL3 S_m2 = pMesh_sec->S.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

				//c values
				double c_m1 = ptrans_sec->cfunc_sec(relpos_m1, stencil);
				double c_m2 = ptrans_sec->cfunc_sec(relpos_m1 + contact.hshift_secondary, stencil);
				double c_1 = cfunc_sec(relpos_1, stencil);
				double c_2 = cfunc_sec(relpos_1 - contact.hshift_primary, stencil);

				//Calculate S drop at the interface
				DBL3 Vs_F = 1.5 * c_1 * S_1 - 0.5 * c_2 * S_2;
				DBL3 Vs_N = 1.5 * c_m1 * S_m1 - 0.5 * c_m2 * S_m2;
				DBL3 dVs = Vs_F - Vs_N;

				//Get G values from top contacting mesh
				DBL2 Gmix;
				if (contact.IsPrimaryTop()) {

					Gmix = pMesh->Gmix;
					pMesh->update_parameters_mcoarse(mcell1_idx, pMesh->Gmix, Gmix);
				}
				else {

					Gmix = pMesh_sec->Gmix;
					pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gmix, Gmix);
				}

				double gI = (2.0 * GMUB_2E / dh) * Gmix.j / (-GAMMA * grel * Ms);
				double gR = (2.0 * GMUB_2E / dh) * Gmix.i / (-GAMMA * grel * Ms);

				pMesh->Heff[mcell1_idx] += tsi_eff * (gI * dVs + gR * (pMesh->M[mcell1_idx] ^ dVs) / Ms);
			}
		}
	}
}

//-------------------Display Calculation Methods

//prepare displayVEC ready for calculation of display quantity
bool Transport::PrepareDisplayVEC(DBL3 cellsize)
{
	if (pSMesh->SolveSpinCurrent() && pMesh->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC.assign(cellsize, pMesh->meshRect, DBL3(0.0));

		return true;
	}
	else displayVEC.clear();

	return false;
}

//return x, y, or z component of spin current (component = 0, 1, or 2)
VEC<DBL3>& Transport::GetSpinCurrent(int component)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC(pMesh->h_e)) return displayVEC;

	//compute spin current and store result in displayVEC depending on required component

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->S.linear_size(); idx++) {

		DBL33 Js = DBL33();

		if (pMesh->S.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			double P = pMesh->P;
			double De = pMesh->De;
			double betaD = pMesh->betaD;
			double SHA = pMesh->SHA;
			pMesh->update_parameters_ecoarse(idx, pMesh->Ms, Ms, pMesh->P, P, pMesh->De, De, pMesh->betaD, betaD, pMesh->SHA, SHA);

			if (pMesh->M.linear_size()) {

				//magnetic mesh terms

				//1. drift
				int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(idx));
						
				DBL3 Mval = pMesh->M[idx_M];
				DBL33 grad_S = pMesh->S.grad_neu(idx);

				Js = (pMesh->Jc[idx] | Mval) * (P / Ms) * (-MUB_E);

				//2. diffusion with homogeneous Neumann boundary condition
				Js -= grad_S * De;

				//3. CPP-GMR term
				if (IsNZ((double)betaD)) {

					DBL3 delsq_S = pMesh->S.delsq_neu(idx);

					Js += ((grad_S * Mval) | Mval) * P * betaD * De / (Ms * Ms);
				}
			}
			else {

				//non-magnetic mesh terms

				//1. SHE contribution
				Js = epsilon3(pMesh->Jc[idx]) * SHA * MUB_E;

				//2. diffusion with non-homogeneous Neumann boundary condition
				Js -= pMesh->S.grad_nneu(idx, epsilon3(pMesh->Jc[idx]) * (SHA * MUB_E / De)) * De;
			}
		}

		switch (component) {

		case 0:
			displayVEC[idx] = Js.x;
			break;
		case 1:
			displayVEC[idx] = Js.y;
			break;
		case 2:
			displayVEC[idx] = Js.z;
			break;
		}
	}

	return displayVEC;
}

DBL3 Transport::GetAverageSpinCurrent(int component, Rect rectangle)
{
#if COMPILECUDA == 1
	//update displayVEC with spin current in TransportCUDA
	GetSpinCurrentCUDA(component);

	//average spin current in displayVEC in TransportCUDA
	if (pSMesh->SolveSpinCurrent()) return cuReal3(reinterpret_cast<TransportCUDA*>(pModuleCUDA)->displayVEC()->average_nonempty(pMesh->n_e.dim(), rectangle));
	else return cuReal3(0.0);
#endif

	//update displayVEC with spin current
	GetSpinCurrent(component);

	//average spin current in displayVEC
	return displayVEC.average_nonempty_omp(rectangle);
}

//return spin torque computed from spin accumulation
VEC<DBL3>& Transport::GetSpinTorque(void)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC(pMesh->h)) return displayVEC;

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->M.linear_size(); idx++) {

		if (pMesh->M.is_empty(idx)) {

			displayVEC[idx] = DBL3();
			continue;
		}

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		double ts_eff = pMesh->ts_eff;
		double l_ex = pMesh->l_ex;
		double l_ph = pMesh->l_ph;
		pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->De, De, pMesh->ts_eff, ts_eff, pMesh->l_ex, l_ex, pMesh->l_ph, l_ph);

		//average S in the magnetization cell with index idx
		DBL3 Sav = pMesh->S.weighted_average(pMesh->M.cellidx_to_position(idx), pMesh->h);
		DBL3 M = pMesh->M[idx];

		displayVEC[idx] = ts_eff * ((Sav ^ M) * De / (Ms * l_ex * l_ex) + (M ^ (Sav ^ M)) * De / (Ms * Ms * l_ph * l_ph));
	}

	return displayVEC;
}

#if COMPILECUDA == 1
cu_obj<cuVEC<cuReal3>>& Transport::GetSpinCurrentCUDA(int component)
{ 
	return reinterpret_cast<TransportCUDA*>(pModuleCUDA)->GetSpinCurrent(component); 
}

cu_obj<cuVEC<cuReal3>>& Transport::GetSpinTorqueCUDA(void) 
{ 
	return reinterpret_cast<TransportCUDA*>(pModuleCUDA)->GetSpinTorque(); 
}
#endif

//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
void Transport::CalculateDisplaySAInterfaceTorque(Transport* ptrans_sec, CMBNDInfo& contact)
{
	//the top contacting mesh sets G values
	bool GInterface_Enabled = ((contact.IsPrimaryTop() && pMesh->GInterface_Enabled()) || (!contact.IsPrimaryTop() && ptrans_sec->pMesh->GInterface_Enabled()));

	if (pMesh->MComputation_Enabled() && !ptrans_sec->pMesh->Magnetisation_Enabled() && GInterface_Enabled) {

		//interface conductance method with F being the primary mesh : calculate and set spin torque

		//convert the cells box from S mesh to M mesh
		INT3 mbox_start = pMesh->M.cellidx_from_position(pMesh->S.cellidx_to_position(contact.cells_box.s) + pMesh->meshRect.s);
		INT3 mbox_end = pMesh->M.cellidx_from_position(pMesh->S.cellidx_to_position(contact.cells_box.e - INT3(1)) + pMesh->meshRect.s);

		if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
		if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
		if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

		INT3 box_sizes = mbox_end - mbox_start;

		//the cellsize perpendicular to the contact (in the M mesh)
		double dh = (DBL3(contact.cell_shift) & pMesh->h).norm();

		Mesh* pMesh_sec = ptrans_sec->pMesh;

		//primary cells in this contact
#pragma omp parallel for
		for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

			int i = (box_idx % box_sizes.x) + mbox_start.i;
			int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
			int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

			//index of magnetic cell 1
			int mcell1_idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

			if (pMesh->M.is_empty(mcell1_idx)) continue;

			double Ms = pMesh->Ms;
			double tsi_eff = pMesh->tsi_eff;
			pMesh->update_parameters_mcoarse(mcell1_idx, pMesh->Ms, Ms, pMesh->tsi_eff, tsi_eff);

			//position at interface relative to primary mesh
			DBL3 mhshift_primary = contact.hshift_primary.normalized() & pMesh->h;
			DBL3 relpos_interf = ((DBL3(i, j, k) + DBL3(0.5)) & pMesh->h) + mhshift_primary / 2;

			DBL3 relpos_1 = relpos_interf - contact.hshift_primary / 2;

			DBL3 relpos_m1 = pMesh->meshRect.s - pMesh_sec->meshRect.s + relpos_interf + contact.hshift_secondary / 2;

			DBL3 stencil = pMesh->h - mod(mhshift_primary) + mod(contact.hshift_secondary);

			//S values
			DBL3 S_1 = pMesh->S.weighted_average(relpos_1, stencil);
			DBL3 S_2 = pMesh->S.weighted_average(relpos_1 - contact.hshift_primary, stencil);
			DBL3 S_m1 = pMesh_sec->S.weighted_average(relpos_m1, stencil);
			DBL3 S_m2 = pMesh_sec->S.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

			//c values
			double c_m1 = ptrans_sec->cfunc_sec(relpos_m1, stencil);
			double c_m2 = ptrans_sec->cfunc_sec(relpos_m1 + contact.hshift_secondary, stencil);
			double c_1 = cfunc_sec(relpos_1, stencil);
			double c_2 = cfunc_sec(relpos_1 - contact.hshift_primary, stencil);

			//Calculate S drop at the interface
			DBL3 Vs_F = 1.5 * c_1 * S_1 - 0.5 * c_2 * S_2;
			DBL3 Vs_N = 1.5 * c_m1 * S_m1 - 0.5 * c_m2 * S_m2;
			DBL3 dVs = Vs_F - Vs_N;

			//Get G values from top contacting mesh
			DBL2 Gmix;
			if (contact.IsPrimaryTop()) {

				Gmix = pMesh->Gmix;
				pMesh->update_parameters_mcoarse(mcell1_idx, pMesh->Gmix, Gmix);
			}
			else {

				Gmix = pMesh_sec->Gmix;
				pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gmix, Gmix);
			}

			double gI = (2.0 * GMUB_2E / dh) * Gmix.j / Ms;
			double gR = (2.0 * GMUB_2E / dh) * Gmix.i / Ms;

			displayVEC[mcell1_idx] += tsi_eff * (gI * (pMesh->M[mcell1_idx] ^ dVs) + gR * (pMesh->M[mcell1_idx] ^ (pMesh->M[mcell1_idx] ^ dVs)) / Ms);
		}
	}
}

#endif