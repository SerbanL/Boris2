#include "Atom_HeatCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_HEAT) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

//-------------------Calculation Methods

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 1-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void IterateHeatEquation_1TM_Kernel(ManagedAtom_MeshCUDA& cuaMesh, cuBReal* heatEq_RHS)
{
	cuVEC_VC<cuBReal>& Temp = *cuaMesh.pTemp;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		cuBReal density = *cuaMesh.pdensity;
		cuBReal shc = *cuaMesh.pshc;
		cuBReal thermCond = *cuaMesh.pthermCond;
		cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pdensity, density, *cuaMesh.pshc, shc, *cuaMesh.pthermCond, thermCond);

		cuBReal cro = density * shc;
		cuBReal K = thermCond;

		//heat equation with Robin boundaries (based on Newton's law of cooling)
		heatEq_RHS[idx] = Temp.delsq_robin(idx, K) * K / cro;

		//add Joule heating if set
		if (E.linear_size()) {

			cuReal3 position = Temp.cellidx_to_position(idx);

			cuReal3 E_value = E.weighted_average(position, Temp.h);
			cuBReal elC_value = elC.weighted_average(position, Temp.h);

			//add Joule heating source term
			heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro;
		}

		//add heat source contribution if set
		if (cuIsNZ(cuaMesh.pQ->get0())) {

			cuBReal Q = *cuaMesh.pQ;
			cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pQ, Q);

			heatEq_RHS[idx] += Q / cro;
		}
	}
}

__global__ void IterateHeatEquation_1TM_Equation_Kernel(
	ManagedAtom_MeshCUDA& cuaMesh, cuBReal* heatEq_RHS,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Q_equation,
	cuBReal time)
{
	cuVEC_VC<cuBReal>& Temp = *cuaMesh.pTemp;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		cuBReal density = *cuaMesh.pdensity;
		cuBReal shc = *cuaMesh.pshc;
		cuBReal thermCond = *cuaMesh.pthermCond;
		cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pdensity, density, *cuaMesh.pshc, shc, *cuaMesh.pthermCond, thermCond);

		cuBReal cro = density * shc;
		cuBReal K = thermCond;

		//heat equation with Robin boundaries (based on Newton's law of cooling)
		heatEq_RHS[idx] = Temp.delsq_robin(idx, K) * K / cro;

		//add Joule heating if set
		if (E.linear_size()) {

			cuReal3 position = Temp.cellidx_to_position(idx);

			cuReal3 E_value = E.weighted_average(position, Temp.h);
			cuBReal elC_value = elC.weighted_average(position, Temp.h);

			//add Joule heating source term
			heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro;
		}

		//add heat source contribution
		cuReal3 relpos = Temp.cellidx_to_position(idx);
		cuBReal Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

		heatEq_RHS[idx] += Q / cro;
	}
}

__global__ void TemperatureFTCS_Atom_Kernel(cuVEC_VC<cuBReal>& Temp, cuBReal* heatEq_RHS, cuBReal dT)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		Temp[idx] += dT * heatEq_RHS[idx];
	}
}

void Atom_HeatCUDA::IterateHeatEquation_1TM(cuBReal dT)
{

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_1TM_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, heatEq_RHS);
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_1TM_Equation_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
			paMeshCUDA->cuaMesh, heatEq_RHS,
			Q_equation.get_x(), paMeshCUDA->GetStageTime());
	}

	//2. Now use forward time to advance by dT
	TemperatureFTCS_Atom_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->Temp, heatEq_RHS, dT);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 2-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void IterateHeatEquation_2TM_Kernel(ManagedAtom_MeshCUDA& cuaMesh, cuBReal* heatEq_RHS, cuBReal dT)
{
	cuVEC_VC<cuBReal>& Temp = *cuaMesh.pTemp;
	cuVEC_VC<cuBReal>& Temp_l = *cuaMesh.pTemp_l;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx)) return;

		cuBReal density = *cuaMesh.pdensity;
		cuBReal shc = *cuaMesh.pshc;
		cuBReal shc_e = *cuaMesh.pshc_e;
		cuBReal G_el = *cuaMesh.pG_e;
		cuBReal thermCond = *cuaMesh.pthermCond;
		cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pdensity, density, *cuaMesh.pshc, shc, *cuaMesh.pshc_e, shc_e, *cuaMesh.pG_e, G_el, *cuaMesh.pthermCond, thermCond);

		cuBReal cro_e = density * shc_e;
		cuBReal K = thermCond;

		//1. Itinerant Electrons Temperature

		if (Temp.is_not_cmbnd(idx)) {

			//heat equation with Robin boundaries (based on Newton's law of cooling) and coupling to lattice
			heatEq_RHS[idx] = (Temp.delsq_robin(idx, K) * K - G_el * (Temp[idx] - Temp_l[idx])) / cro_e;

			//add Joule heating if set
			if (E.linear_size()) {

				cuReal3 position = Temp.cellidx_to_position(idx);

				cuBReal elC_value = elC.weighted_average(position, Temp.h);
				cuReal3 E_value = E.weighted_average(position, Temp.h);

				//add Joule heating source term
				heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro_e;
			}

			//add heat source contribution if set
			if (cuIsNZ(cuaMesh.pQ->get0())) {

				cuBReal Q = *cuaMesh.pQ;
				cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pQ, Q);

				heatEq_RHS[idx] += Q / cro_e;
			}
		}

		//2. Lattice Temperature

		//lattice specific heat capacity + electron specific heat capacity gives the total specific heat capacity
		cuBReal cro_l = density * (shc - shc_e);

		Temp_l[idx] += dT * G_el * (Temp[idx] - Temp_l[idx]) / cro_l;
	}
}

__global__ void IterateHeatEquation_2TM_Equation_Kernel(
	ManagedAtom_MeshCUDA& cuaMesh, cuBReal* heatEq_RHS,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Q_equation,
	cuBReal time, cuBReal dT)
{
	cuVEC_VC<cuBReal>& Temp = *cuaMesh.pTemp;
	cuVEC_VC<cuBReal>& Temp_l = *cuaMesh.pTemp_l;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx)) return;

		cuBReal density = *cuaMesh.pdensity;
		cuBReal shc = *cuaMesh.pshc;
		cuBReal shc_e = *cuaMesh.pshc_e;
		cuBReal G_el = *cuaMesh.pG_e;
		cuBReal thermCond = *cuaMesh.pthermCond;
		cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pdensity, density, *cuaMesh.pshc, shc, *cuaMesh.pshc_e, shc_e, *cuaMesh.pG_e, G_el, *cuaMesh.pthermCond, thermCond);

		cuBReal cro_e = density * shc_e;
		cuBReal K = thermCond;

		//1. Itinerant Electrons Temperature

		if (Temp.is_not_cmbnd(idx)) {

			//heat equation with Robin boundaries (based on Newton's law of cooling) and coupling to lattice
			heatEq_RHS[idx] = (Temp.delsq_robin(idx, K) * K - G_el * (Temp[idx] - Temp_l[idx])) / cro_e;

			//add Joule heating if set
			if (E.linear_size()) {

				cuReal3 position = Temp.cellidx_to_position(idx);

				cuBReal elC_value = elC.weighted_average(position, Temp.h);
				cuReal3 E_value = E.weighted_average(position, Temp.h);

				//add Joule heating source term
				heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro_e;
			}

			//add heat source contribution
			cuReal3 relpos = Temp.cellidx_to_position(idx);
			cuBReal Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

			heatEq_RHS[idx] += Q / cro_e;
		}

		//2. Lattice Temperature

		//lattice specific heat capacity + electron specific heat capacity gives the total specific heat capacity
		cuBReal cro_l = density * (shc - shc_e);

		Temp_l[idx] += dT * G_el * (Temp[idx] - Temp_l[idx]) / cro_l;
	}
}

void Atom_HeatCUDA::IterateHeatEquation_2TM(cuBReal dT)
{

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_2TM_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, heatEq_RHS, dT);
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_2TM_Equation_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
			paMeshCUDA->cuaMesh, heatEq_RHS,
			Q_equation.get_x(), paMeshCUDA->GetStageTime(), dT);
	}

	//2. Now use forward time to advance by dT
	TemperatureFTCS_Atom_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->Temp, heatEq_RHS, dT);
}

//-------------------Setters

//non-uniform temperature setting
__global__ void SetBaseTemperature_Nonuniform_Kernel(ManagedAtom_MeshCUDA& cuaMesh, cuBReal Temperature)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC_VC<cuBReal>& Temp = *cuaMesh.pTemp;
	cuVEC_VC<cuBReal>& Temp_l = *cuaMesh.pTemp_l;

	if (idx < Temp.linear_size()) {

		if (Temp.is_not_empty(idx)) {

			cuBReal cT = *cuaMesh.pcT;
			cuaMesh.update_parameters_tcoarse(idx, *cuaMesh.pcT, cT);

			Temp[idx] = cT * Temperature;

			if (Temp_l.linear_size()) Temp_l[idx] = cT * Temperature;
		}
	}
}

//set Temp non-uniformly as specified through the cT mesh parameter
void Atom_HeatCUDA::SetBaseTemperature_Nonuniform(cuBReal Temperature)
{
	SetBaseTemperature_Nonuniform_Kernel <<< (paMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, Temperature);
}

//set Temp uniformly to base temperature
void Atom_HeatCUDA::SetBaseTemperature(cuBReal Temperature)
{
	paMeshCUDA->Temp()->setnonempty(Temperature);
	paMeshCUDA->Temp_l()->setnonempty(Temperature);
}

#endif

#endif