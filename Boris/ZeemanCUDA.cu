#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ZEEMAN

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

//----------------------- Initialization

__global__ void set_ZeemanCUDA_pointers_kernel(
	ManagedMeshCUDA& cuMesh, cuVEC<cuReal3>& Havec, cuVEC<cuReal3>& globalField)
{
	if (threadIdx.x == 0) cuMesh.pHavec = &Havec;
	if (threadIdx.x == 1) cuMesh.pglobalField = &globalField;
}

void ZeemanCUDA::set_ZeemanCUDA_pointers(void)
{
	set_ZeemanCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(pMeshCUDA->cuMesh, Havec, globalField);
}

//----------------------- Computation

__global__ void ZeemanCUDA_UpdateField_FM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuVEC<cuReal3>& Havec, cuVEC<cuReal3>& globalField, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hext = cuReal3();

		if (Havec.linear_size()) Hext = Havec[idx];
		else {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			Hext = cHA * Ha;
		}

		if (globalField.linear_size()) Hext += globalField[idx];

		Heff[idx] = Hext;

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * Hext / non_empty_cells;
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hext;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * M[idx] * Hext;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void ZeemanCUDA_UpdateField_Equation_FM(
	ManagedMeshCUDA& cuMesh, 
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuVEC<cuReal3>& globalField, 
	ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hext = cuReal3();

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		cuReal3 relpos = M.cellidx_to_position(idx);
		cuReal3 H = cuReal3(
			H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

		Hext = cHA * H;
		if (globalField.linear_size()) Hext += globalField[idx];

		Heff[idx] = Hext;

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * Hext / non_empty_cells;
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hext;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * M[idx] * Hext;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void ZeemanCUDA_UpdateField_AFM(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuVEC<cuReal3>& Havec, cuVEC<cuReal3>& globalField, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hext = cuReal3();

		if (Havec.linear_size()) Hext = Havec[idx];
		else {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			Hext = cHA * Ha;
		}

		if (globalField.linear_size()) Hext += globalField[idx];

		Heff[idx] = Hext;
		Heff2[idx] = Hext;

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * Hext / (2 * non_empty_cells);
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hext;
		if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[idx] = Hext;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -MU0 * M[idx] * Hext;
		if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[idx] = -MU0 * M2[idx] * Hext;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void ZeemanCUDA_UpdateField_Equation_AFM(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuVEC<cuReal3>& globalField,
	ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hext = cuReal3();

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		cuReal3 relpos = M.cellidx_to_position(idx);
		cuReal3 H = cuReal3(
			H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

		Hext = cHA * H;
		if (globalField.linear_size()) Hext += globalField[idx];

		Heff[idx] = Hext;
		Heff2[idx] = Hext;

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * Hext / (2 * non_empty_cells);
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hext;
		if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[idx] = Hext;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -MU0 * M[idx] * Hext;
		if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[idx] = -MU0 * M2[idx] * Hext;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
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

				ZeemanCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Ha, Havec, globalField, cuModule, true);
			}
			else ZeemanCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Ha, Havec, globalField, cuModule, false);
		}

		else {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Ha, Havec, globalField, cuModule, true);
			}
			else ZeemanCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Ha, Havec, globalField, cuModule, false);
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_Equation_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
					pMeshCUDA->cuMesh,
					H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
					pMeshCUDA->GetStageTime(),
					globalField,
					cuModule, true);
			}
			else ZeemanCUDA_UpdateField_Equation_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				globalField,
				cuModule, false);
		}

		else {

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				ZeroEnergy();

				ZeemanCUDA_UpdateField_Equation_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
					pMeshCUDA->cuMesh,
					H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
					pMeshCUDA->GetStageTime(),
					globalField,
					cuModule, true);
			}
			else ZeemanCUDA_UpdateField_Equation_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				pMeshCUDA->cuMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				pMeshCUDA->GetStageTime(),
				globalField,
				cuModule, false);
		}
	}
}

//-------------------Others

BError ZeemanCUDA::SetFieldEquation(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
{
	BError error(CLASS_STR(ZeemanCUDA));

	if (!H_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (Havec()->size_cpu().dim()) Havec()->clear();

	return error;
}

#endif

#endif