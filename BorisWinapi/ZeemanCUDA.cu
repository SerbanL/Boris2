#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ZEEMAN

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void ZeemanCUDA_UpdateField_FMDM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		Heff[idx] += (cHA * Ha);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * (cHA * Ha) / non_empty_cells;
		}
	}

	if(do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void ZeemanCUDA_UpdateField_Equation_FMDM(
	ManagedMeshCUDA& cuMesh, 
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		cuReal3 relpos = M.cellidx_to_position(idx);
		cuReal3 H = cuReal3(
			H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

		Heff[idx] += (cHA * H);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * (cHA * H) / non_empty_cells;
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void ZeemanCUDA_UpdateField_AFM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		Heff[idx] += (cHA * Ha);
		Heff2[idx] += (cHA * Ha);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * (cHA * Ha) / (2 * non_empty_cells);
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void ZeemanCUDA_UpdateField_Equation_AFM(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		cuReal3 relpos = M.cellidx_to_position(idx);
		cuReal3 H = cuReal3(
			H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

		Heff[idx] += (cHA * H);
		Heff2[idx] += (cHA * H);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * (cHA * H) / (2 * non_empty_cells);
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void ZeemanCUDA::UpdateField(void)
{
	/////////////////////////////////////////
	// Fixed set field
	/////////////////////////////////////////

	if (!H_equation.is_set()) {

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_AFM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, true);
			}
			else ZeemanCUDA_UpdateField_AFM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, false);
		}

		else {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_FMDM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, true);
			}
			else ZeemanCUDA_UpdateField_FMDM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, false);
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_Equation_AFM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
					pMeshCUDA->cuMesh,
					H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
					pMeshCUDA->GetStageTime(),
					energy, true);
			}
			else ZeemanCUDA_UpdateField_Equation_AFM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				energy, false);
		}

		else {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_Equation_FMDM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
					pMeshCUDA->cuMesh,
					H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
					pMeshCUDA->GetStageTime(),
					energy, true);
			}
			else ZeemanCUDA_UpdateField_Equation_FMDM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				energy, false);
		}
	}
}

//-------------------Energy density methods

__global__ void ZeemanCUDA_GetEnergy_FMDM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			energy_ = -(cuBReal)MU0 * M[idx] * (cHA * Ha);
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void ZeemanCUDA_GetEnergy_Equation_FMDM(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			cuReal3 relpos = M.cellidx_to_position(idx);
			cuReal3 H = cuReal3(
				H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

			energy_ = -(cuBReal)MU0 * M[idx] * (cHA * H);
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void ZeemanCUDA_GetEnergy_AFM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * (cHA * Ha) / 2;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void ZeemanCUDA_GetEnergy_Equation_AFM(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			cuReal3 relpos = M.cellidx_to_position(idx);
			cuReal3 H = cuReal3(
				H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

			energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * (cHA * H) / 2;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

cuBReal ZeemanCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	/////////////////////////////////////////
	// Fixed set field
	/////////////////////////////////////////

	if (!H_equation.is_set()) {

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			ZeemanCUDA_GetEnergy_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Ha, energy, points_count, avRect);
		}

		else {

			ZeemanCUDA_GetEnergy_FMDM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Ha, energy, points_count, avRect);
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			ZeemanCUDA_GetEnergy_Equation_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				energy, points_count, avRect);
		}

		else {

			ZeemanCUDA_GetEnergy_Equation_FMDM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				energy, points_count, avRect);
		}
	}

	size_t points_count_cpu = points_count.to_cpu();

	if (points_count_cpu) return energy.to_cpu() / points_count_cpu;
	else return 0.0;
}

//-------------------Others

BError ZeemanCUDA::SetFieldEquation(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
{
	BError error(CLASS_STR(ZeemanCUDA));

	if (!H_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif

#endif