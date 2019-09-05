#include "ExchangeBaseCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void CalculateExchangeCoupling_kernel(
	ManagedMeshCUDA& mesh_sec, ManagedMeshCUDA& mesh_pri,
	CMBNDInfoCUDA& contact,
	cuReal& energy, bool do_reduction)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal energy_ = 0.0;

	cuINT3 box_sizes = contact.cells_box.size();

	cuVEC_VC<cuReal3>& M_pri = *mesh_pri.pM;
	cuVEC<cuReal3>& Heff_pri = *mesh_pri.pHeff;
	cuVEC_VC<cuReal3>& M_sec = *mesh_sec.pM;

	if (box_idx < box_sizes.dim()) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		cuReal hRsq = contact.hshift_primary.norm();
		hRsq *= hRsq;

		int cell1_idx = i + j * M_pri.n.x + k * M_pri.n.x*M_pri.n.y;

		if (M_pri.is_not_empty(cell1_idx) && M_pri.is_cmbnd(cell1_idx)) {

			//calculate second primary cell index
			int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * M_pri.n.x + (k + contact.cell_shift.k) * M_pri.n.x*M_pri.n.y;

			//relative position of cell -1 in secondary mesh
			cuReal3 relpos_m1 = M_pri.rect.s - M_sec.rect.s + ((cuReal3(i, j, k) + cuReal3(0.5)) & M_pri.h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

			//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
			cuReal3 stencil = M_pri.h - cu_mod(contact.hshift_primary) + cu_mod(contact.hshift_secondary);

			cuReal Ms = *mesh_pri.pMs;
			cuReal A = *mesh_pri.pA;
			mesh_pri.update_parameters_mcoarse(cell1_idx, *mesh_pri.pA, A, *mesh_pri.pMs, Ms);

			cuReal3 Hexch;

			//values at cells -1, 1
			cuReal3 M_1 = M_pri[cell1_idx];
			cuReal3 M_m1 = M_sec.weighted_average(relpos_m1, stencil);

			if (cell2_idx < M_pri.n.dim() && M_pri.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				cuReal3 M_2 = M_pri[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_2 + M_m1 - 2 * M_1) / hRsq;
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_m1 - M_1) / hRsq;
			}

			Heff_pri[cell1_idx] += Hexch;

			if (do_reduction) {

				int non_empty_cells = M_pri.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuReal)MU0 * M_1 * Hexch / (2 * non_empty_cells);
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- CalculateExchangeCoupling LAUNCHER

//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
void ExchangeBaseCUDA::CalculateExchangeCoupling(cu_obj<cuReal>& energy)
{
	for (int contact_idx = 0; contact_idx < CMBNDcontacts.size(); contact_idx++) {

		size_t size = CMBNDcontacts[contact_idx].cells_box.size().dim();

		//the contacting meshes indexes : secondary mesh index is the one in contact with this one (the primary)
		int idx_sec = CMBNDcontacts[contact_idx].mesh_idx.i;
		int idx_pri = CMBNDcontacts[contact_idx].mesh_idx.j;

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			CalculateExchangeCoupling_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*pContactingManagedMeshes[idx_sec], *pContactingManagedMeshes[idx_pri], CMBNDcontactsCUDA[contact_idx], energy, true);
		}
		else {

			CalculateExchangeCoupling_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*pContactingManagedMeshes[idx_sec], *pContactingManagedMeshes[idx_pri], CMBNDcontactsCUDA[contact_idx], energy, false);
		}
	}
}

#endif