#include "HeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_HEAT

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

//-------------------Calculation Methods

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 1-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void IterateHeatEquation_1TM_Kernel(ManagedMeshCUDA& cuMesh, cuBReal* heatEq_RHS)
{
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp; 
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		cuBReal density = *cuMesh.pdensity;
		cuBReal shc = *cuMesh.pshc;
		cuBReal thermCond = *cuMesh.pthermCond;
		cuMesh.update_parameters_tcoarse(idx, *cuMesh.pdensity, density, *cuMesh.pshc, shc, *cuMesh.pthermCond, thermCond);

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
		if (cuIsNZ(cuMesh.pQ->get0())) {

			cuBReal Q = *cuMesh.pQ;
			cuMesh.update_parameters_tcoarse(idx, *cuMesh.pQ, Q);

			heatEq_RHS[idx] += Q / cro;
		}
	}
}

__global__ void IterateHeatEquation_1TM_Equation_Kernel(
	ManagedMeshCUDA& cuMesh, cuBReal* heatEq_RHS,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Q_equation,
	cuBReal time)
{
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		cuBReal density = *cuMesh.pdensity;
		cuBReal shc = *cuMesh.pshc;
		cuBReal thermCond = *cuMesh.pthermCond;
		cuMesh.update_parameters_tcoarse(idx, *cuMesh.pdensity, density, *cuMesh.pshc, shc, *cuMesh.pthermCond, thermCond);

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

__global__ void TemperatureFTCS_Kernel(cuVEC_VC<cuBReal>& Temp, cuBReal* heatEq_RHS, cuBReal dT)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		Temp[idx] += dT * heatEq_RHS[idx];
	}
}

void HeatCUDA::IterateHeatEquation_1TM(cuBReal dT)
{

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_1TM_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, heatEq_RHS);
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {
	
		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_1TM_Equation_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
			pMeshCUDA->cuMesh, heatEq_RHS,
			Q_equation.get_x(), pMeshCUDA->GetStageTime());
	}

	//2. Now use forward time to advance by dT
	TemperatureFTCS_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->Temp, heatEq_RHS, dT);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 2-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void IterateHeatEquation_2TM_Kernel(ManagedMeshCUDA& cuMesh, cuBReal* heatEq_RHS, cuBReal dT)
{
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;
	cuVEC_VC<cuBReal>& Temp_l = *cuMesh.pTemp_l;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx)) return;

		cuBReal density = *cuMesh.pdensity;
		cuBReal shc = *cuMesh.pshc;
		cuBReal shc_e = *cuMesh.pshc_e;
		cuBReal G_el = *cuMesh.pG_e;
		cuBReal thermCond = *cuMesh.pthermCond;
		cuMesh.update_parameters_tcoarse(idx, *cuMesh.pdensity, density, *cuMesh.pshc, shc, *cuMesh.pshc_e, shc_e, *cuMesh.pG_e, G_el, *cuMesh.pthermCond, thermCond);

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
			if (cuIsNZ(cuMesh.pQ->get0())) {

				cuBReal Q = *cuMesh.pQ;
				cuMesh.update_parameters_tcoarse(idx, *cuMesh.pQ, Q);

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
	ManagedMeshCUDA& cuMesh, cuBReal* heatEq_RHS,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Q_equation,
	cuBReal time, cuBReal dT)
{
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;
	cuVEC_VC<cuBReal>& Temp_l = *cuMesh.pTemp_l;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx)) return;

		cuBReal density = *cuMesh.pdensity;
		cuBReal shc = *cuMesh.pshc;
		cuBReal shc_e = *cuMesh.pshc_e;
		cuBReal G_el = *cuMesh.pG_e;
		cuBReal thermCond = *cuMesh.pthermCond;
		cuMesh.update_parameters_tcoarse(idx, *cuMesh.pdensity, density, *cuMesh.pshc, shc, *cuMesh.pshc_e, shc_e, *cuMesh.pG_e, G_el, *cuMesh.pthermCond, thermCond);

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

void HeatCUDA::IterateHeatEquation_2TM(cuBReal dT)
{

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_2TM_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, heatEq_RHS, dT);
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		IterateHeatEquation_2TM_Equation_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
			pMeshCUDA->cuMesh, heatEq_RHS,
			Q_equation.get_x(), pMeshCUDA->GetStageTime(), dT);
	}

	//2. Now use forward time to advance by dT
	TemperatureFTCS_Kernel << < (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->Temp, heatEq_RHS, dT);
}

//-------------------Setters

//non-uniform temperature setting
__global__ void SetBaseTemperature_Nonuniform_Kernel(ManagedMeshCUDA& cuMesh, cuBReal Temperature)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;
	cuVEC_VC<cuBReal>& Temp_l = *cuMesh.pTemp_l;

	if (idx < Temp.linear_size()) {

		if (Temp.is_not_empty(idx)) {

			cuBReal cT = *cuMesh.pcT;
			cuMesh.update_parameters_tcoarse(idx, *cuMesh.pcT, cT);

			Temp[idx] = cT * Temperature;

			if (Temp_l.linear_size()) Temp_l[idx] = cT * Temperature;
		}
	}
}

//set Temp non-uniformly as specified through the cT mesh parameter
void HeatCUDA::SetBaseTemperature_Nonuniform(cuBReal Temperature)
{
	SetBaseTemperature_Nonuniform_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Temperature);
}

//set Temp uniformly to base temperature
void HeatCUDA::SetBaseTemperature(cuBReal Temperature)
{
	pMeshCUDA->Temp()->setnonempty(Temperature);
	pMeshCUDA->Temp_l()->setnonempty(Temperature);
}

#endif

#endif