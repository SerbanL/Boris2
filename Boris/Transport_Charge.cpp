#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

//-------------------Calculation Methods

DBL2 Transport::IterateChargeSolver_SOR(double damping)
{
	if (!is_thermoelectric_mesh) {

		//no thermoelectric effect
		return pMesh->V.IteratePoisson_SOR<Transport>(&Transport::Evaluate_ChargeSolver_delsqV_RHS, *this, damping);
	}
	else {

		//include thermoelectric effect
		return pMesh->V.IteratePoisson_SOR<Transport>(&Transport::Evaluate_ChargeSolver_delsqV_Thermoelectric_RHS, &Transport::NHNeumann_Vdiff_Thermoelectric, *this, damping);
	}
}

//call-back method for Poisson equation to evaluate RHS
double Transport::Evaluate_ChargeSolver_delsqV_RHS(int idx) const
{
	//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma
	//The Poisson solver calls this method to evaluate the RHS of this equation

	return -(pMesh->V.grad_diri(idx) * pMesh->elC.grad_sided(idx)) / pMesh->elC[idx];
}

//call-back method for Poisson equation to evaluate RHS (with thermoelectric effect included)
double Transport::Evaluate_ChargeSolver_delsqV_Thermoelectric_RHS(int idx) const
{
	//include thermoelectric effect with Seebeck coefficient : delsq V = -S delsq T, obtained from div J = 0, where J = -sigma(grad V + S * grad T).
	//here we ignore gradients in sigma and S

	double Sc = pMesh->Sc;
	double thermCond = pMesh->thermCond;
	pMesh->update_parameters_ecoarse(idx, pMesh->Sc, Sc, pMesh->thermCond, thermCond);

	//corresponding index in Temp
	int idx_temp = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(idx));

	return -Sc * pMesh->Temp.delsq_robin(idx_temp, thermCond);
}

//Non-homogeneous Neumann boundary condition for V' when solving with thermoelectric effect included - call-back method for Poisson equation for V, charge solver only
DBL3 Transport::NHNeumann_Vdiff_Thermoelectric(int idx) const
{
	//Gradient of V normal to boundary is -S * grad T normal to boundary.
	//grad T normal to boundary obtained using Robin boundary condition as:

	//-K * grad T . n = heat flux normal to boundary = alpha(Tb - Ta), where alpha is Robin value, Ta ambient temperature and Tb temperature at boundary, K thermal conductivity
	//Thus grad V . n = S*alpha*(Tb-Ta) / K

	DBL3 bdiff = DBL3();
	DBL3 shift = pMesh->V.get_shift_to_emptycell(idx);

	if (!shift.IsNull()) {

		//corresponding index in Temp
		int idx_temp = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(idx));

		if (pMesh->Temp.is_cmbnd(idx_temp)) {

			//at composite media boundaries (for Temp) cannot use Robin boundary for heat flux
			//instead we can approximate it using a sided differential in the cell just next to the boundary
			double Sc = pMesh->Sc;
			pMesh->update_parameters_ecoarse(idx, pMesh->Sc, Sc);

			bdiff = Sc * pMesh->Temp.grad_sided(idx_temp);
		}
		else {

			//boundary, not a cmbnd. Use grad V . n = S*alpha*(Tb-Ta) / K
			double Sc = pMesh->Sc;
			double thermCond = pMesh->thermCond;
			pMesh->update_parameters_ecoarse(idx, pMesh->Sc, Sc, pMesh->thermCond, thermCond);

			bdiff = pMesh->Temp.get_robin_value(pMesh->V.cellidx_to_position(idx), shift) * Sc / thermCond;
			//at negative boundary inverse sign since heat flux normal also has opposite sign
			if (shift.x < 0) bdiff.x *= -1;
			if (shift.y < 0) bdiff.y *= -1;
			if (shift.z < 0) bdiff.z *= -1;
		}
	}

	return bdiff;
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of a potential and flux

//Charge transport only : V

//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
//With thermoelectric effect this becomes:
//1. No net current (e.g. no electrodes attached):
//Jc = -sigma * grad V - sigma * Sc * grad T = a + b * grad V -> a = -sigma * Sc * grad T and b = -sigma taken at the interface
//grad T normal to interface is given by Robin boundary condition : grad T . n = -alpha(Tb-Ta)/K
//2. Net current generated (open potential condition):
//Jc = sigma * Sc * grad T, so a = sigma * Sc * grad T, b = 0
double Transport::afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	if (!is_thermoelectric_mesh) return 0.0;
	else {

		//include thermoelectric effect

		int cellm1_idx = pMesh->V.position_to_cellidx(relpos_m1);

		//corresponding index in Temp
		int idx_temp1 = pMesh->Temp.position_to_cellidx(relpos_m1);
		int idx_temp2 = pMesh->Temp.position_to_cellidx(relpos_m1 + shift);

		double Sc = pMesh->Sc;
		pMesh->update_parameters_ecoarse(cellm1_idx, pMesh->Sc, Sc);

		//normalized shift in normal direction to boundary: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
		//do not use mod here as we need the temperature gradient to point in normal direction to boundary
		DBL3 nshift = normalize(shift);

		double T_grad1 = pMesh->Temp.grad_sided(idx_temp1) * nshift;
		double T_grad2 = pMesh->Temp.grad_sided(idx_temp2) * nshift;
		double T_grad = 1.5 * T_grad1 - 0.5 * T_grad2;

		//shift is from m1 to m2 cell, so use minus sign here if open potential mode
		if (pSTrans->IsOpenPotential()) return -Sc * pMesh->elC[cellm1_idx] * T_grad;
		else return Sc * pMesh->elC[cellm1_idx] * T_grad;
	}
}

double Transport::afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	if (!is_thermoelectric_mesh) return 0.0;
	else {

		//include thermoelectric effect

		double Sc = pMesh->Sc;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->Sc, Sc);

		//normalized shift in normal direction to boundary: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
		//do not use mod here as we need the temperature gradient to point in normal direction to boundary
		DBL3 nshift = normalize(shift);

		//corresponding index in Temp
		int idx_temp1 = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(cell1_idx));
		int idx_temp2 = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(cell2_idx));

		double T_grad1 = pMesh->Temp.grad_sided(idx_temp1) * nshift;
		double T_grad2 = pMesh->Temp.grad_sided(idx_temp2) * nshift;
		double T_grad = 1.5 * T_grad1 - 0.5 * T_grad2;

		//shift is from cell1 to cell2 so no need for minus sign adjustment
		return Sc * pMesh->elC[cell1_idx] * T_grad;
	}	
}

double Transport::bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	if (is_thermoelectric_mesh && pSTrans->IsOpenPotential()) return 0.0;
	else return -(1.5 * pMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * pMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

double Transport::bfunc_V_pri(int cell1_idx, int cell2_idx) const
{
	if (is_thermoelectric_mesh && pSTrans->IsOpenPotential()) return 0.0;
	else return -(1.5 * pMesh->elC[cell1_idx] - 0.5 * pMesh->elC[cell2_idx]);
}

//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
//With thermoelectric effect this becomes (Thomson effect not included):
//delsq V = -grad V * grad elC / elC - Sc * grad T * grad elC / elC - Sc * delsq Temp
double Transport::diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = pMesh->V.position_to_cellidx(relpos_m1);
	return diff2_V_pri(cellm1_idx, shift);
}

double Transport::diff2_V_pri(int cell1_idx, DBL3 shift) const
{
	if (!is_thermoelectric_mesh) {

		//normalized, positive shift: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
		DBL3 nshift = mod(normalize(shift));

		//no thermoelectric effect
		return -((pMesh->V.grad_diri(cell1_idx) * nshift) * (pMesh->elC.grad_sided(cell1_idx) * nshift)) / pMesh->elC[cell1_idx];
	}
	else {

		//include thermoelectric effect with Seebeck coefficient

		double Sc = pMesh->Sc;
		pMesh->update_parameters_ecoarse(cell1_idx, pMesh->Sc, Sc);

		//corresponding index in Temp
		int idx_temp = pMesh->Temp.position_to_cellidx(pMesh->V.cellidx_to_position(cell1_idx));

		return -Sc * pMesh->Temp.delsq_neu(idx_temp);
	}
}

#endif