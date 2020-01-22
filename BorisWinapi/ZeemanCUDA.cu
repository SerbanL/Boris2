#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ZEEMAN

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Mesh_FerromagneticCUDA.h"
#include "Mesh_AntiFerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void ZeemanCUDA_UpdateField_FM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, bool do_reduction)
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

__global__ void ZeemanCUDA_UpdateField_Equation_FM(
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

		if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_FM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, true);
			}
			else ZeemanCUDA_UpdateField_FM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, false);
		}

		else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_AFM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, true);
			}
			else ZeemanCUDA_UpdateField_AFM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, false);
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_Equation_FM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
					pMeshCUDA->cuMesh, 
					H_equation.get_x(), H_equation.get_y(), H_equation.get_z(), 
					pMeshCUDA->GetStageTime(),
					energy, true);
			}
			else ZeemanCUDA_UpdateField_Equation_FM << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				energy, false);
		}

		else if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

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
	}
}

BError ZeemanCUDA::SetFieldEquation(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
{
	BError error(CLASS_STR(ZeemanCUDA));

	if (!H_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif

#endif