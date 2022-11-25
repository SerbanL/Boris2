#include "MElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"

__global__ void Set_Strain_From_Formula_Sd_Sod_Kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sd_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sd_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sd_equation_z,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sod_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sod_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sod_equation_z,
	cuBReal time)
{
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < strain_diag.linear_size()) {

		if (strain_diag.is_not_empty(idx)) {

			cuReal3 relpos = strain_diag.cellidx_to_position(idx);
			strain_diag[idx] = cuReal3(
				Sd_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sd_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sd_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

			strain_odiag[idx] = cuReal3(
				Sod_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sod_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sod_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));
		}
	}
}

__global__ void Set_Strain_From_Formula_Sd_Kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sd_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sd_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sd_equation_z,
	cuBReal time)
{
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < strain_diag.linear_size()) {

		if (strain_diag.is_not_empty(idx)) {

			cuReal3 relpos = strain_diag.cellidx_to_position(idx);
			strain_diag[idx] = cuReal3(
				Sd_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sd_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sd_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));

			strain_odiag[idx] = cuReal3();
		}
	}
}

__global__ void Set_Strain_From_Formula_Sod_Kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sod_equation_x,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sod_equation_y,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sod_equation_z,
	cuBReal time)
{
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < strain_diag.linear_size()) {

		if (strain_diag.is_not_empty(idx)) {

			cuReal3 relpos = strain_diag.cellidx_to_position(idx);
			strain_diag[idx] = cuReal3();

			strain_odiag[idx] = cuReal3(
				Sod_equation_x.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sod_equation_y.evaluate(relpos.x, relpos.y, relpos.z, time),
				Sod_equation_z.evaluate(relpos.x, relpos.y, relpos.z, time));
		}
	}
}

//----------------------------------------------- Auxiliary

//Run-time auxiliary to set strain directly from user supplied text formulas
void MElasticCUDA::Set_Strain_From_Formula(void)
{
	if (Sd_equation.is_set() && Sod_equation.is_set()) {

		Set_Strain_From_Formula_Sd_Sod_Kernel <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh,
			Sd_equation.get_x(), Sd_equation.get_y(), Sd_equation.get_z(),
			Sod_equation.get_x(), Sod_equation.get_y(), Sod_equation.get_z(),
			pMeshCUDA->GetStageTime());
	}
	else if (Sd_equation.is_set()) {

		Set_Strain_From_Formula_Sd_Kernel <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh,
			Sd_equation.get_x(), Sd_equation.get_y(), Sd_equation.get_z(),
			pMeshCUDA->GetStageTime());
	}
	else if (Sod_equation.is_set()) {

		Set_Strain_From_Formula_Sod_Kernel <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh,
			Sod_equation.get_x(), Sod_equation.get_y(), Sod_equation.get_z(),
			pMeshCUDA->GetStageTime());
	}
}

//----------------------- UpdateField LAUNCHER

void MElasticCUDA::UpdateField(void)
{
	if (Sd_equation.is_set_vector() || Sod_equation.is_set_vector()) {

		//strain specified using a formula
		Set_Strain_From_Formula();
	}
}

#endif

#endif