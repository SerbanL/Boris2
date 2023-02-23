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
	if (!is_thermoelectric_mesh) {

		//no thermoelectric effect
		return paMesh->V.IteratePoisson_SOR<Atom_Transport>(&Atom_Transport::Evaluate_ChargeSolver_delsqV_RHS, *this, damping);
	}
	else {

		//include thermoelectric effect
		return paMesh->V.IteratePoisson_SOR<Atom_Transport>(&Atom_Transport::Evaluate_ChargeSolver_delsqV_Thermoelectric_RHS, &Atom_Transport::NHNeumann_Vdiff_Thermoelectric, *this, damping);
	}
}

double Atom_Transport::Evaluate_ChargeSolver_delsqV_RHS(int idx) const
{
	//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma
	//The Poisson solver calls this method to evaluate the RHS of this equation
	return -(paMesh->V.grad_diri(idx) * paMesh->elC.grad_sided(idx)) / paMesh->elC[idx];
}

//call-back method for Poisson equation to evaluate RHS (with thermoelectric effect included)
double Atom_Transport::Evaluate_ChargeSolver_delsqV_Thermoelectric_RHS(int idx) const
{
	//include thermoelectric effect with Seebeck coefficient : delsq V = -S delsq T, obtained from div J = 0, where J = -sigma(grad V + S * grad T).
	//here we ignore gradients in sigma and S

	double Sc = paMesh->Sc;
	double thermCond = paMesh->thermCond;
	paMesh->update_parameters_ecoarse(idx, paMesh->Sc, Sc, paMesh->thermCond, thermCond);

	//corresponding index in Temp
	int idx_temp = paMesh->Temp.position_to_cellidx(paMesh->V.cellidx_to_position(idx));

	return -Sc * paMesh->Temp.delsq_robin(idx_temp, thermCond);
}

//Non-homogeneous Neumann boundary condition for V' when solving with thermoelectric effect included - call-back method for Poisson equation for V, charge solver only
DBL3 Atom_Transport::NHNeumann_Vdiff_Thermoelectric(int idx) const
{
	//Gradient of V normal to boundary is -S * grad T normal to boundary.
	//grad T normal to boundary obtained using Robin boundary condition as:

	//-K * grad T . n = heat flux normal to boundary = alpha(Tb - Ta), where alpha is Robin value, Ta ambient temperature and Tb temperature at boundary, K thermal conductivity
	//Thus grad V . n = S*alpha*(Tb-Ta) / K

	DBL3 bdiff = DBL3();
	DBL3 shift = paMesh->V.get_shift_to_emptycell(idx);

	if (!shift.IsNull()) {

		//corresponding index in Temp
		int idx_temp = paMesh->Temp.position_to_cellidx(paMesh->V.cellidx_to_position(idx));

		if (paMesh->Temp.is_cmbnd(idx_temp)) {

			//at composite media boundaries (for Temp) cannot use Robin boundary for heat flux
			//instead we can approximate it using a sided differential in the cell just next to the boundary
			double Sc = paMesh->Sc;
			paMesh->update_parameters_ecoarse(idx, paMesh->Sc, Sc);

			bdiff = Sc * paMesh->Temp.grad_sided(idx_temp);
		}
		else {

			//boundary, not a cmbnd. Use grad V . n = S*alpha*(Tb-Ta) / K
			double Sc = paMesh->Sc;
			double thermCond = paMesh->thermCond;
			paMesh->update_parameters_ecoarse(idx, paMesh->Sc, Sc, paMesh->thermCond, thermCond);

			bdiff = paMesh->Temp.get_robin_value(paMesh->V.cellidx_to_position(idx), shift) * Sc / thermCond;
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
double Atom_Transport::afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	if (!is_thermoelectric_mesh) return 0.0;
	else {

		//include thermoelectric effect

		int cellm1_idx = paMesh->V.position_to_cellidx(relpos_m1);

		//corresponding index in Temp
		int idx_temp1 = paMesh->Temp.position_to_cellidx(relpos_m1);
		int idx_temp2 = paMesh->Temp.position_to_cellidx(relpos_m1 + shift);

		double Sc = paMesh->Sc;
		paMesh->update_parameters_ecoarse(cellm1_idx, paMesh->Sc, Sc);

		//normalized shift in normal direction to boundary: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
		//do not use mod here as we need the temperature gradient to point in normal direction to boundary
		DBL3 nshift = normalize(shift);

		double T_grad1 = paMesh->Temp.grad_sided(idx_temp1) * nshift;
		double T_grad2 = paMesh->Temp.grad_sided(idx_temp2) * nshift;
		double T_grad = 1.5 * T_grad1 - 0.5 * T_grad2;

		//shift is from m1 to m2 cell, so use minus sign here if open potential mode
		if (pSTrans->IsOpenPotential()) return -Sc * paMesh->elC[cellm1_idx] * T_grad;
		else return Sc * paMesh->elC[cellm1_idx] * T_grad;
	}
}

double Atom_Transport::afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	if (!is_thermoelectric_mesh) return 0.0;
	else {

		//include thermoelectric effect

		double Sc = paMesh->Sc;
		paMesh->update_parameters_ecoarse(cell1_idx, paMesh->Sc, Sc);

		//normalized shift in normal direction to boundary: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
		//do not use mod here as we need the temperature gradient to point in normal direction to boundary
		DBL3 nshift = normalize(shift);

		//corresponding index in Temp
		int idx_temp1 = paMesh->Temp.position_to_cellidx(paMesh->V.cellidx_to_position(cell1_idx));
		int idx_temp2 = paMesh->Temp.position_to_cellidx(paMesh->V.cellidx_to_position(cell2_idx));

		double T_grad1 = paMesh->Temp.grad_sided(idx_temp1) * nshift;
		double T_grad2 = paMesh->Temp.grad_sided(idx_temp2) * nshift;
		double T_grad = 1.5 * T_grad1 - 0.5 * T_grad2;

		//shift is from cell1 to cell2 so no need for minus sign adjustment
		return Sc * paMesh->elC[cell1_idx] * T_grad;
	}
}

double Atom_Transport::bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	if (is_thermoelectric_mesh && pSTrans->IsOpenPotential()) return 0.0;
	else return -(1.5 * paMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * paMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

double Atom_Transport::bfunc_V_pri(int cell1_idx, int cell2_idx) const
{
	if (is_thermoelectric_mesh && pSTrans->IsOpenPotential()) return 0.0;
	else return -(1.5 * paMesh->elC[cell1_idx] - 0.5 * paMesh->elC[cell2_idx]);
}

//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
double Atom_Transport::diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = paMesh->V.position_to_cellidx(relpos_m1);
	return diff2_V_pri(cellm1_idx, shift);
}

double Atom_Transport::diff2_V_pri(int cell1_idx, DBL3 shift) const
{
	if (!is_thermoelectric_mesh) {

		//normalized, positive shift: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
		DBL3 nshift = mod(normalize(shift));

		//no thermoelectric effect
		return -((paMesh->V.grad_diri(cell1_idx) * nshift) * (paMesh->elC.grad_sided(cell1_idx) * nshift)) / paMesh->elC[cell1_idx];
	}
	else {

		//include thermoelectric effect with Seebeck coefficient

		double Sc = paMesh->Sc;
		paMesh->update_parameters_ecoarse(cell1_idx, paMesh->Sc, Sc);

		//corresponding index in Temp
		int idx_temp = paMesh->Temp.position_to_cellidx(paMesh->V.cellidx_to_position(cell1_idx));

		return -Sc * paMesh->Temp.delsq_neu(idx_temp);
	}
}

#endif