#include "Atom_ZeemanCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

//----------------------- Initialization

__global__ void set_Atom_ZeemanCUDA_pointers_kernel(
	ManagedAtom_MeshCUDA& cuaMesh, cuVEC<cuReal3>& Havec, cuVEC<cuReal3>& globalField)
{
	if (threadIdx.x == 0) cuaMesh.pHavec = &Havec;
	if (threadIdx.x == 1) cuaMesh.pglobalField = &globalField;
}

void Atom_ZeemanCUDA::set_Atom_ZeemanCUDA_pointers(void)
{
	set_Atom_ZeemanCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(paMeshCUDA->cuaMesh, Havec, globalField);
}

//----------------------- Computation

__global__ void Atom_ZeemanCUDA_UpdateField_Cubic(ManagedAtom_MeshCUDA& cuMesh, cuReal3& Ha, cuVEC<cuReal3>& Havec, cuVEC<cuReal3>& globalField, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuReal3 Hext = cuReal3();

		if (Havec.linear_size()) Hext = Havec[idx];
		else {

			cuBReal cHA = *cuMesh.pcHA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

			Hext = cHA * Ha;
		}

		if (globalField.linear_size()) Hext += globalField[idx];

		Heff1[idx] = Hext;

		if (do_reduction) {

			//energy density
			int non_empty_cells = M1.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MUB_MU0 * M1[idx] * Hext / (non_empty_cells * M1.h.dim());
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hext;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MUB_MU0 * M1[idx] * Hext / M1.h.dim();
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void Atom_ZeemanCUDA_UpdateField_Equation_Cubic(
	ManagedAtom_MeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& H_equation_z,
	cuBReal time,
	cuVEC<cuReal3>& globalField,
	ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuReal3 Hext = cuReal3();

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		cuReal3 relpos = M1.cellidx_to_position(idx);
		cuReal3 H = cuReal3(
			H_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
			H_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

		Hext = cHA * H;
		if (globalField.linear_size()) Hext += globalField[idx];

		Heff1[idx] = Hext;

		if (do_reduction) {

			//energy density
			int non_empty_cells = M1.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MUB_MU0 * M1[idx] * Hext / (non_empty_cells * M1.h.dim());
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hext;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MUB_MU0 * M1[idx] * Hext / M1.h.dim();
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
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

				Atom_ZeemanCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, Ha, Havec, globalField, cuModule, true);
			}
			else Atom_ZeemanCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, Ha, Havec, globalField, cuModule, false);
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
					globalField,
					cuModule, true);
			}
			else Atom_ZeemanCUDA_UpdateField_Equation_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
				paMeshCUDA->cuaMesh,
				H_equation.get_x(), H_equation.get_y(), H_equation.get_z(),
				paMeshCUDA->GetStageTime(),
				globalField,
				cuModule, false);
		}
	}
}

//-------------------Others

BError Atom_ZeemanCUDA::SetFieldEquation(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	if (!H_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (Havec()->size_cpu().dim()) Havec()->clear();

	return error;
}

#endif

#endif