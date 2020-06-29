#include "Atom_ZeemanCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

__global__ void Atom_ZeemanCUDA_UpdateField_Cubic(ManagedAtom_MeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		Heff1[idx] += (cHA * Ha);

		if (do_reduction) {

			//energy density
			int non_empty_cells = M1.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MUB * M1[idx] * (cuBReal)MU0 * (cHA * Ha) / (non_empty_cells * M1.h.dim());
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void Atom_ZeemanCUDA_UpdateField_Equation_Cubic(
	ManagedAtom_MeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		cuReal3 relpos = M1.cellidx_to_position(idx);
		cuReal3 H = cuReal3(
			H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

		Heff1[idx] += (cHA * H);

		if (do_reduction) {

			//energy density
			int non_empty_cells = M1.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MUB * M1[idx] * (cuBReal)MU0 * (cHA * H) / (non_empty_cells * M1.h.dim());
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Atom_ZeemanCUDA::UpdateField(void)
{
	/////////////////////////////////////////
	// Fixed set field
	/////////////////////////////////////////

	if (!H_equation.is_set()) {

		if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

			if (paMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				Atom_ZeemanCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, Ha, energy, true);
			}
			else Atom_ZeemanCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, Ha, energy, false);
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

			if (paMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				Atom_ZeemanCUDA_UpdateField_Equation_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
					paMeshCUDA->cuaMesh,
					H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
					paMeshCUDA->GetStageTime(),
					energy, true);
			}
			else Atom_ZeemanCUDA_UpdateField_Equation_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				paMeshCUDA->cuaMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				paMeshCUDA->GetStageTime(),
				energy, false);
		}
	}
}

//-------------------Energy density methods

__global__ void Atom_ZeemanCUDA_GetEnergy_Cubic(ManagedAtom_MeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && avRect.contains(M1.cellidx_to_position(idx))) {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			//energy, not energy density : will convert to energy density later
			energy_ = -(cuBReal)MUB * M1[idx] * (cuBReal)MU0 * (cHA * Ha);
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void Atom_ZeemanCUDA_GetEnergy_Equation_Cubic(
	ManagedAtom_MeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && avRect.contains(M1.cellidx_to_position(idx))) {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			cuReal3 relpos = M1.cellidx_to_position(idx);
			cuReal3 H = cuReal3(
				H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

			//energy, not energy density : will convert to energy density later
			energy_ = -(cuBReal)MUB * M1[idx] * (cuBReal)MU0 * (cHA * H);
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

cuBReal Atom_ZeemanCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	/////////////////////////////////////////
	// Fixed set field
	/////////////////////////////////////////

	if (!H_equation.is_set()) {

		if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

			Atom_ZeemanCUDA_GetEnergy_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, Ha, energy, points_count, avRect);
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

			Atom_ZeemanCUDA_GetEnergy_Equation_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				paMeshCUDA->cuaMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				paMeshCUDA->GetStageTime(),
				energy, points_count, avRect);
		}
	}

	size_t points_count_cpu = points_count.to_cpu();

	if (points_count_cpu) return energy.to_cpu() / (points_count_cpu * paMeshCUDA->h.dim());
	else return 0.0;
}

//-------------------Others

BError Atom_ZeemanCUDA::SetFieldEquation(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	if (!H_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif

#endif