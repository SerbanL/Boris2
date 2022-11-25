#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "SuperMeshCUDA.h"

#include "BorisCUDALib.cuh"

#include "Atom_MeshParamsControlCUDA.h"

//--------------------------------------------------------------- Electrical Conductivity with AMR

__global__ void Atom_CalculateElectricalConductivity_AMR_Kernel(ManagedAtom_MeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal amrPercentage = *cuMesh.pamrPercentage;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.pamrPercentage, amrPercentage);

			//get current density value at this conductivity cell
			cuReal3 jc_value = cu_normalize(elC[idx] * E[idx]);

			//get M value (M is on n, h mesh so could be different)
			cuReal3 m_value = cu_normalize(M1[elC.cellidx_to_position(idx)]);

			cuBReal dotproduct = jc_value * m_value;

			elC[idx] = elecCond / (1 + amrPercentage * dotproduct*dotproduct / 100);
		}
	}
}

//calculate electrical conductivity with AMR present
void Atom_TransportCUDA::CalculateElectricalConductivity_AMR(void)
{
	Atom_CalculateElectricalConductivity_AMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
}

//--------------------------------------------------------------- Electrical Conductivity with TAMR

__global__ void Atom_CalculateElectricalConductivity_TAMR_Kernel(ManagedAtom_MeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal tamrPercentage = *cuMesh.ptamrPercentage;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.ptamrPercentage, tamrPercentage);

			//get M value (M is on n, h mesh so could be different)
			cuReal3 m_value = cu_normalize(M1[elC.cellidx_to_position(idx)]);

			cuBReal dotprod1 = m_value * mcanis_ea1;
			cuBReal dotprod2 = m_value * mcanis_ea2;
			cuBReal dotprod3 = m_value * mcanis_ea3;

			//default formula : ro = ro0 * (1 + TAMR * sin^2(theta)), TAMR is the ratio.
			elC[idx] = elecCond / (1 + (tamrPercentage / 100) * (1 - dotprod1 * dotprod1));
		}
	}
}

__global__ void Atom_CalculateElectricalConductivity_TAMR_Equation_Kernel(ManagedAtom_MeshCUDA& cuMesh, ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& TAMR_conductivity_equation)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal tamrPercentage = *cuMesh.ptamrPercentage;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.ptamrPercentage, tamrPercentage);

			//get M value (M is on n, h mesh so could be different)
			cuReal3 m_value = cu_normalize(M1[elC.cellidx_to_position(idx)]);

			cuBReal dotprod1 = m_value * mcanis_ea1;
			cuBReal dotprod2 = m_value * mcanis_ea2;
			cuBReal dotprod3 = m_value * mcanis_ea3;

			elC[idx] = elecCond * TAMR_conductivity_equation.evaluate(tamrPercentage / 100, dotprod1, dotprod2, dotprod3);
		}
	}
}

//calculate electrical conductivity with AMR present
void Atom_TransportCUDA::CalculateElectricalConductivity_TAMR(void)
{
	if (TAMR_conductivity_equation.is_set()) {

		Atom_CalculateElectricalConductivity_TAMR_Equation_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, TAMR_conductivity_equation.get_x());
	}
	else {

		//default formula set
		Atom_CalculateElectricalConductivity_TAMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
	}
}

//--------------------------------------------------------------- Electrical Conductivity with TAMR and AMR

__global__ void Atom_CalculateElectricalConductivity_TAMR_and_AMR_Kernel(ManagedAtom_MeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal amrPercentage = *cuMesh.pamrPercentage;
			cuBReal tamrPercentage = *cuMesh.ptamrPercentage;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.pamrPercentage, amrPercentage, *cuMesh.ptamrPercentage, tamrPercentage, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

			//get current density value at this conductivity cell
			cuReal3 jc_value = cu_normalize(elC[idx] * E[idx]);

			//get M value (M is on n, h mesh so could be different)
			cuReal3 m_value = cu_normalize(M1[elC.cellidx_to_position(idx)]);

			cuBReal dotproduct = jc_value * m_value;
			cuBReal dotprod1 = m_value * mcanis_ea1;
			cuBReal dotprod2 = m_value * mcanis_ea2;
			cuBReal dotprod3 = m_value * mcanis_ea3;

			cuBReal resistivity_ratio = (1 + (tamrPercentage / 100) * (1 - dotprod1 * dotprod1)) + (1 + amrPercentage * dotproduct*dotproduct / 100);
			elC[idx] = elecCond / resistivity_ratio;
		}
	}
}

__global__ void Atom_CalculateElectricalConductivity_TAMR_Equation_and_AMR_Kernel(ManagedAtom_MeshCUDA& cuMesh, ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& TAMR_conductivity_equation)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal amrPercentage = *cuMesh.pamrPercentage;
			cuBReal tamrPercentage = *cuMesh.ptamrPercentage;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.pamrPercentage, amrPercentage, *cuMesh.ptamrPercentage, tamrPercentage, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

			//get current density value at this conductivity cell
			cuReal3 jc_value = cu_normalize(elC[idx] * E[idx]);

			//get M value (M is on n, h mesh so could be different)
			cuReal3 m_value = cu_normalize(M1[elC.cellidx_to_position(idx)]);

			cuBReal dotproduct = jc_value * m_value;
			cuBReal dotprod1 = m_value * mcanis_ea1;
			cuBReal dotprod2 = m_value * mcanis_ea2;
			cuBReal dotprod3 = m_value * mcanis_ea3;

			cuBReal resistivity_ratio = (1 / TAMR_conductivity_equation.evaluate(tamrPercentage / 100, dotprod1, dotprod2, dotprod3) + (1 + amrPercentage * dotproduct*dotproduct / 100));
			elC[idx] = elecCond / resistivity_ratio;
		}
	}
}

//calculate electrical conductivity with AMR present
void Atom_TransportCUDA::CalculateElectricalConductivity_TAMR_and_AMR(void)
{
	if (TAMR_conductivity_equation.is_set()) {

		Atom_CalculateElectricalConductivity_TAMR_Equation_and_AMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, TAMR_conductivity_equation.get_x());
	}
	else {

		//default formula set
		Atom_CalculateElectricalConductivity_TAMR_and_AMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
	}
}

//--------------------------------------------------------------- Electrical Conductivity without AMR

__global__ void Atom_CalculateElectricalConductivity_NoAMR_Kernel(ManagedAtom_MeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond);

			elC[idx] = elecCond;
		}
	}
}

//calculate electrical conductivity without AMR
void Atom_TransportCUDA::CalculateElectricalConductivity_NoAMR(void)
{
	Atom_CalculateElectricalConductivity_NoAMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
}

//--------------------------------------------------------------- Electric Field

//No thermoelectric effect

//Current density when only charge solver is used
__global__ void Atom_CalculateElectricField_Charge_Kernel(cuVEC<cuReal3>& E, cuVEC_VC<cuBReal>& V)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < V.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			E[idx] = -1.0 * V.grad_diri(idx);
		}
		else E[idx] = cuReal3(0.0);
	}
}

//Thermoelectric effect but no net current generated

__global__ void Atom_CalculateElectricField_Charge_ThermoElectric_Kernel(ManagedAtom_MeshCUDA& cuMesh)
{
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& V = *cuMesh.pV;
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < E.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			cuBReal Sc = *cuMesh.pSc;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pSc, Sc);

			//corresponding index in Temp
			int idx_temp = Temp.position_to_cellidx(V.cellidx_to_position(idx));

			//include thermoelectric effect, but no net current generated

			cuReal3 shift = V.get_shift_to_emptycell(idx);
			if (!shift.IsNull()) {

				//use sided differentials if both neighbors not available
				E[idx] = -1 * (V.grad_sided(idx) + Sc * Temp.grad_sided(idx_temp));
			}
			else E[idx] = -1 * (V.grad_neu(idx) + Sc * Temp.grad_neu(idx_temp));
		}
		else E[idx] = cuReal3(0);
	}
}

//Thermoelectric effect with a net current generated

__global__ void Atom_CalculateElectricField_Charge_ThermoElectric_OpenPotential_Kernel(ManagedAtom_MeshCUDA& cuMesh, cuBReal& mesh_thermoelectric_net_current)
{
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& V = *cuMesh.pV;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal mesh_thermoelectric_net_current_ = 0.0;

	if (idx < E.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			cuBReal Sc = *cuMesh.pSc;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pSc, Sc);

			//corresponding index in Temp
			int idx_temp = Temp.position_to_cellidx(V.cellidx_to_position(idx));

			E[idx] = Sc * Temp.grad_sided(idx_temp);

			if (V.is_cmbnd(idx) || V.is_dirichlet(idx)) {

				cuReal3 therm_current_density = elC[idx] * E[idx];

				//below divide by 2 since we want net current, not twice the current come into or out of the mesh

				//x
				if (V.is_cmbnd_x(idx) || V.is_dirichlet_x(idx)) {

					cuBReal cmbnd_current_density = therm_current_density.x;
					mesh_thermoelectric_net_current_ -= cmbnd_current_density * V.h.y * V.h.z / 2;
				}

				//y
				else if (V.is_cmbnd_y(idx) || V.is_dirichlet_y(idx)) {

					cuBReal cmbnd_current_density = therm_current_density.y;
					mesh_thermoelectric_net_current_ -= cmbnd_current_density * V.h.x * V.h.z / 2;
				}

				//z
				if (V.is_cmbnd_z(idx) || V.is_dirichlet_z(idx)) {

					cuBReal cmbnd_current_density = therm_current_density.z;
					mesh_thermoelectric_net_current_ -= cmbnd_current_density * V.h.x * V.h.y / 2;
				}
			}
		}
		else E[idx] = cuReal3(0);
	}

	reduction_sum(0, 1, &mesh_thermoelectric_net_current_, mesh_thermoelectric_net_current);
}

//-------------------Calculation Methods : Electric Field

//calculate electric field as the negative gradient of V
void Atom_TransportCUDA::CalculateElectricField(bool open_potential)
{
	if (!is_thermoelectric_mesh) {

		//no thermoelectric effect in this mesh
		Atom_CalculateElectricField_Charge_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->E, paMeshCUDA->V);
	}
	else {

		//include thermoeletric effect

		if (!open_potential) {

			//no net thermoelectric current generated
			Atom_CalculateElectricField_Charge_ThermoElectric_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
		}
		else {

			//thermoelectric effect but with open potential, so a net current is generated - must count total current coming out of this mesh
			mesh_thermoelectric_net_current.from_cpu(0.0);

			Atom_CalculateElectricField_Charge_ThermoElectric_OpenPotential_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, mesh_thermoelectric_net_current);
		}
	}
}

#endif

#endif