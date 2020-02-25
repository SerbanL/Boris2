#include "iDMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_IDMEXCHANGE

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

// both contacting meshes are ferromagnetic
__global__ void CalculateiDMExchangeCoupling_FM_kernel(
	ManagedMeshCUDA& mesh_sec, ManagedMeshCUDA& mesh_pri,
	CMBNDInfoCUDA& contact,
	cuBReal& energy, bool do_reduction)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	cuINT3 box_sizes = contact.cells_box.size();

	cuVEC_VC<cuReal3>& M_pri = *mesh_pri.pM;
	cuVEC<cuReal3>& Heff_pri = *mesh_pri.pHeff;
	cuVEC_VC<cuReal3>& M_sec = *mesh_sec.pM;

	if (box_idx < box_sizes.dim()) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		cuReal3 hshift_primary = contact.hshift_primary;
		cuBReal hRsq = hshift_primary.norm();
		hRsq *= hRsq;

		int cell1_idx = i + j * M_pri.n.x + k * M_pri.n.x*M_pri.n.y;

		if (M_pri.is_not_empty(cell1_idx) && M_pri.is_cmbnd(cell1_idx)) {

			//calculate second primary cell index
			int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * M_pri.n.x + (k + contact.cell_shift.k) * M_pri.n.x*M_pri.n.y;

			//relative position of cell -1 in secondary mesh
			cuReal3 relpos_m1 = M_pri.rect.s - M_sec.rect.s + ((cuReal3(i, j, k) + cuReal3(0.5)) & M_pri.h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

			//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
			cuReal3 stencil = M_pri.h - cu_mod(contact.hshift_primary) + cu_mod(contact.hshift_secondary);

			cuBReal Ms = *mesh_pri.pMs;
			cuBReal A = *mesh_pri.pA;
			cuBReal D = *mesh_pri.pD;
			mesh_pri.update_parameters_mcoarse(cell1_idx, *mesh_pri.pA, A, *mesh_pri.pD, D, *mesh_pri.pMs, Ms);

			cuReal3 Hexch;

			//values at cells -1, 1
			cuReal3 M_1 = M_pri[cell1_idx];
			cuReal3 M_m1 = M_sec.weighted_average(relpos_m1, stencil);

			if (cell2_idx < M_pri.n.dim() && M_pri.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				cuReal3 M_2 = M_pri[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_2 + M_m1 - 2 * M_1) / hRsq;

				//add iDMI contributions at CMBND cells, correcting for the sided differentials already applied here
				//the contributions are different depending on the CMBND coupling direction
				if (cuIsNZ(hshift_primary.x)) {

					//along x
					Hexch += (-2 * D / ((cuBReal)MU0*Ms*Ms)) * cuReal3(M_2.z + M_m1.z - 2 * M_1.z, 0, -M_2.x - M_m1.x + 2 * M_1.x) / (2 * hshift_primary.x);
				}
				else if (cuIsNZ(hshift_primary.y)) {

					//along y
					Hexch += (-2 * D / ((cuBReal)MU0*Ms*Ms)) * cuReal3(0, M_2.z + M_m1.z - 2 * M_1.z, -M_2.y - M_m1.y + 2 * M_1.y) / (2 * hshift_primary.y);
				}
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_m1 - M_1) / hRsq;
			}

			Heff_pri[cell1_idx] += Hexch;

			if (do_reduction) {

				int non_empty_cells = M_pri.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * M_1 * Hexch / (2 * non_empty_cells);
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

// both contacting meshes are antiferromagnetic
__global__ void CalculateiDMExchangeCoupling_AFM_kernel(
	ManagedMeshCUDA& mesh_sec, ManagedMeshCUDA& mesh_pri,
	CMBNDInfoCUDA& contact,
	cuBReal& energy, bool do_reduction)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	cuINT3 box_sizes = contact.cells_box.size();

	cuVEC_VC<cuReal3>& M_pri = *mesh_pri.pM;
	cuVEC_VC<cuReal3>& M2_pri = *mesh_pri.pM2;

	cuVEC<cuReal3>& Heff_pri = *mesh_pri.pHeff;
	cuVEC<cuReal3>& Heff2_pri = *mesh_pri.pHeff2;

	cuVEC_VC<cuReal3>& M_sec = *mesh_sec.pM;
	cuVEC_VC<cuReal3>& M2_sec = *mesh_sec.pM2;

	if (box_idx < box_sizes.dim()) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		cuReal3 hshift_primary = contact.hshift_primary;
		cuBReal hRsq = hshift_primary.norm();
		hRsq *= hRsq;

		int cell1_idx = i + j * M_pri.n.x + k * M_pri.n.x*M_pri.n.y;

		if (M_pri.is_not_empty(cell1_idx) && M_pri.is_cmbnd(cell1_idx)) {

			//calculate second primary cell index
			int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * M_pri.n.x + (k + contact.cell_shift.k) * M_pri.n.x*M_pri.n.y;

			//relative position of cell -1 in secondary mesh
			cuReal3 relpos_m1 = M_pri.rect.s - M_sec.rect.s + ((cuReal3(i, j, k) + cuReal3(0.5)) & M_pri.h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

			//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
			cuReal3 stencil = M_pri.h - cu_mod(contact.hshift_primary) + cu_mod(contact.hshift_secondary);

			cuReal2 Ms_AFM = *mesh_pri.pMs_AFM;
			cuReal2 A_AFM = *mesh_pri.pA_AFM;
			cuReal2 D_AFM = *mesh_pri.pD_AFM;
			mesh_pri.update_parameters_mcoarse(cell1_idx, *mesh_pri.pA_AFM, A_AFM, *mesh_pri.pD_AFM, D_AFM, *mesh_pri.pMs_AFM, Ms_AFM);

			cuReal3 Hexch, Hexch_B;

			//values at cells -1, 1
			cuReal3 M_1 = M_pri[cell1_idx];
			cuReal3 M_m1 = M_sec.weighted_average(relpos_m1, stencil);

			cuReal3 M_1_B = M2_pri[cell1_idx];
			cuReal3 M_m1_B = M2_sec.weighted_average(relpos_m1, stencil);

			if (cell2_idx < M_pri.n.dim() && M_pri.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				cuReal3 M_2 = M_pri[cell2_idx];
				cuReal3 M_2_B = M2_pri[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * (M_2 + M_m1 - 2 * M_1) / hRsq;
				Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * (M_2_B + M_m1_B - 2 * M_1_B) / hRsq;

				//add iDMI contributions at CMBND cells, correcting for the sided differentials already applied here
				//the contributions are different depending on the CMBND coupling direction
				if (cuIsNZ(hshift_primary.x)) {

					//along x
					Hexch += (-2 * D_AFM.i / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.i)) * cuReal3(M_2.z + M_m1.z - 2 * M_1.z, 0, -M_2.x - M_m1.x + 2 * M_1.x) / (2 * hshift_primary.x);
					Hexch_B += (-2 * D_AFM.j / ((cuBReal)MU0*Ms_AFM.j*Ms_AFM.j)) * cuReal3(M_2_B.z + M_m1_B.z - 2 * M_1_B.z, 0, -M_2_B.x - M_m1_B.x + 2 * M_1_B.x) / (2 * hshift_primary.x);
				}
				else if (cuIsNZ(hshift_primary.y)) {

					//along y
					Hexch += (-2 * D_AFM.i / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.i)) * cuReal3(0, M_2.z + M_m1.z - 2 * M_1.z, -M_2.y - M_m1.y + 2 * M_1.y) / (2 * hshift_primary.y);
					Hexch_B += (-2 * D_AFM.j / ((cuBReal)MU0*Ms_AFM.j*Ms_AFM.j)) * cuReal3(0, M_2_B.z + M_m1_B.z - 2 * M_1_B.z, -M_2_B.y - M_m1_B.y + 2 * M_1_B.y) / (2 * hshift_primary.y);
				}
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * (M_m1 - M_1) / hRsq;
				Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * (M_m1_B - M_1_B) / hRsq;
			}

			Heff_pri[cell1_idx] += Hexch;
			Heff2_pri[cell1_idx] += Hexch_B;

			if (do_reduction) {

				int non_empty_cells = M_pri.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M_1 * Hexch + M_1_B * Hexch_B) / (4 * non_empty_cells);
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- CalculateExchangeCoupling LAUNCHER

//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
void iDMExchangeCUDA::CalculateExchangeCoupling(cu_obj<cuBReal>& energy)
{
	for (int contact_idx = 0; contact_idx < CMBNDcontacts.size(); contact_idx++) {

		size_t size = CMBNDcontacts[contact_idx].cells_box.size().dim();

		//the contacting meshes indexes : secondary mesh index is the one in contact with this one (the primary)
		int idx_sec = CMBNDcontacts[contact_idx].mesh_idx.i;
		int idx_pri = CMBNDcontacts[contact_idx].mesh_idx.j;

		if (pContactingMeshes[idx_pri]->GetMeshType() == MESH_ANTIFERROMAGNETIC && pContactingMeshes[idx_sec]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			// both contacting meshes are antiferromagnetic

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				CalculateiDMExchangeCoupling_AFM_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*pContactingManagedMeshes[idx_sec], *pContactingManagedMeshes[idx_pri], CMBNDcontactsCUDA[contact_idx], energy, true);
			}
			else {

				CalculateiDMExchangeCoupling_AFM_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*pContactingManagedMeshes[idx_sec], *pContactingManagedMeshes[idx_pri], CMBNDcontactsCUDA[contact_idx], energy, false);
			}
		}
		else if (pContactingMeshes[idx_pri]->GetMeshType() == MESH_FERROMAGNETIC && pContactingMeshes[idx_sec]->GetMeshType() == MESH_FERROMAGNETIC) {

			//both meshes are ferromagnetic

			if (pMeshCUDA->CurrentTimeStepSolved()) {

				CalculateiDMExchangeCoupling_FM_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*pContactingManagedMeshes[idx_sec], *pContactingManagedMeshes[idx_pri], CMBNDcontactsCUDA[contact_idx], energy, true);
			}
			else {

				CalculateiDMExchangeCoupling_FM_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*pContactingManagedMeshes[idx_sec], *pContactingManagedMeshes[idx_pri], CMBNDcontactsCUDA[contact_idx], energy, false);
			}
		}
	}
}

#endif

#endif