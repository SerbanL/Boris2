#include "MElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void MElasticCUDA_UpdateField_FM(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		if (M.is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal2 MEc = *cuMesh.pMEc;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pMEc, MEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			cuReal3 position = M.cellidx_to_position(idx);
			//xx, yy, zz
			cuReal3 Sd = strain_diag[position];
			//yz, xz, xy
			cuReal3 Sod = strain_odiag[position];

			//normalised magnetisation
			//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

			cuReal3 m = cuReal3(M[idx] * mcanis_ea1, M[idx] * mcanis_ea2, M[idx] * mcanis_ea3) / Ms;
			Sd = cuReal3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
			Sod = cuReal3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

			cuReal3 Hmel_1 = (-2.0 * MEc.i / (MU0 * Ms)) * cuReal3(
				m.x*Sd.x*mcanis_ea1.x + m.y*Sd.y*mcanis_ea2.x + m.z*Sd.z*mcanis_ea3.x,
				m.x*Sd.x*mcanis_ea1.y + m.y*Sd.y*mcanis_ea2.y + m.z*Sd.z*mcanis_ea3.y,
				m.x*Sd.x*mcanis_ea1.z + m.y*Sd.y*mcanis_ea2.z + m.z*Sd.z*mcanis_ea3.z);

			cuReal3 Hmel_2 = (-2.0 * MEc.j / (MU0 * Ms)) * cuReal3(
				Sod.z * (mcanis_ea1.x*m.y + mcanis_ea2.x*m.x) + Sod.y * (mcanis_ea1.x*m.z + mcanis_ea3.x*m.x) + Sod.x * (mcanis_ea2.x*m.z + mcanis_ea3.x*m.y),
				Sod.z * (mcanis_ea1.y*m.y + mcanis_ea2.y*m.x) + Sod.y * (mcanis_ea1.y*m.z + mcanis_ea3.y*m.x) + Sod.x * (mcanis_ea2.y*m.z + mcanis_ea3.y*m.y),
				Sod.z * (mcanis_ea1.z*m.y + mcanis_ea2.z*m.x) + Sod.y * (mcanis_ea1.z*m.z + mcanis_ea3.z*m.x) + Sod.x * (mcanis_ea2.z*m.z + mcanis_ea3.z*m.y));

			Heff[idx] += Hmel_1 + Hmel_2;

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * (Hmel_1 + Hmel_2) / (2 * non_empty_cells);
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void MElasticCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		MElasticCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, true);
	}
	else MElasticCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, false);
}

#endif

#endif