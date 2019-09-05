#include "stdafx.h"
#include "ExchangeBase.h"
#include "Supermesh.h"
#include "Mesh_Ferromagnetic.h"

ExchangeBase::ExchangeBase(Mesh *pMesh_)
{
	pMesh = dynamic_cast<FMesh*>(pMesh_);
}

BError ExchangeBase::Initialize(void)
{
	BError error(CLASS_STR(ExchangeBase));

	//clear everything then rebuild
	pM.clear();
	pFMeshes.clear();
	CMBNDcontacts.clear();

	if (pMesh->GetMeshExchangeCoupling()) {

		//ExchangeBase given friend access to Mesh
		SuperMesh* pSMesh = pMesh->pSMesh;

		//now build pM
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			//use MComputation_Enabled check, not Magnetisation_Enabled check as we only want to couple to other ferromagnetic meshes, not dipole meshes.
			if ((*pSMesh)[idx]->MComputation_Enabled()) {

				pFMeshes.push_back(dynamic_cast<FMesh*>((*pSMesh)[idx]));
				pM.push_back(&(*pSMesh)[idx]->M);
			}
		}

		//set cmbnd flags
		for (int idx = 0; idx < pFMeshes.size(); idx++) {

			if (pFMeshes[idx] == pMesh) {

				//set CMBND flags, even for 1 cell thickness in cmbnd direction
				CMBNDcontacts = pM[idx]->set_cmbnd_flags(idx, pM, false);
				break;
			}
		}
	}

	return error;
}

//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
void ExchangeBase::CalculateExchangeCoupling(double& energy)
{
	for (int contact_idx = 0; contact_idx < CMBNDcontacts.size(); contact_idx++) {

		//the contacting meshes indexes : secondary mesh index is the one in contact with this one (the primary)
		int idx_sec = CMBNDcontacts[contact_idx].mesh_idx.i;
		int idx_pri = CMBNDcontacts[contact_idx].mesh_idx.j;

		//secondary and primary ferromagnetic meshes
		FMesh& FMesh_pri = *pFMeshes[idx_pri];
		FMesh& FMesh_sec = *pFMeshes[idx_sec];

		SZ3 n = FMesh_pri.M.n;
		DBL3 h = FMesh_pri.M.h;

		double hRsq = CMBNDcontacts[contact_idx].hshift_primary.norm();
		hRsq *= hRsq;

		INT3 box_sizes = CMBNDcontacts[contact_idx].cells_box.size();

		double energy_ = 0.0;

		//primary cells in this contact
#pragma omp parallel for reduction (+:energy_)
		for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

			int i = (box_idx % box_sizes.x) + CMBNDcontacts[contact_idx].cells_box.s.i;
			int j = ((box_idx / box_sizes.x) % box_sizes.y) + CMBNDcontacts[contact_idx].cells_box.s.j;
			int k = (box_idx / (box_sizes.x * box_sizes.y)) + CMBNDcontacts[contact_idx].cells_box.s.k;

			int cell1_idx = i + j * n.x + k * n.x*n.y;

			if (FMesh_pri.M.is_empty(cell1_idx) || FMesh_pri.M.is_not_cmbnd(cell1_idx)) continue;

			//calculate second primary cell index
			int cell2_idx = (i + CMBNDcontacts[contact_idx].cell_shift.i) + (j + CMBNDcontacts[contact_idx].cell_shift.j) * n.x + (k + CMBNDcontacts[contact_idx].cell_shift.k) * n.x*n.y;

			//relative position of cell -1 in secondary mesh
			DBL3 relpos_m1 = FMesh_pri.M.rect.s - FMesh_sec.M.rect.s + ((DBL3(i, j, k) + DBL3(0.5)) & h) + (CMBNDcontacts[contact_idx].hshift_primary + CMBNDcontacts[contact_idx].hshift_secondary) / 2;

			//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
			DBL3 stencil = h - mod(CMBNDcontacts[contact_idx].hshift_primary) + mod(CMBNDcontacts[contact_idx].hshift_secondary);

			double Ms = FMesh_pri.Ms;
			double A = FMesh_pri.A;
			FMesh_pri.update_parameters_mcoarse(cell1_idx, FMesh_pri.A, A, FMesh_pri.Ms, Ms);

			DBL3 Hexch;

			//values at cells -1, 1
			DBL3 M_1 = FMesh_pri.M[cell1_idx];
			DBL3 M_m1 = FMesh_sec.M.weighted_average(relpos_m1, stencil);

			if (cell2_idx < n.dim() && FMesh_pri.M.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				DBL3 M_2 = FMesh_pri.M[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_2 + M_m1 - 2 * M_1) / hRsq;
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_m1 - M_1) / hRsq;
			}
			
			FMesh_pri.Heff[cell1_idx] += Hexch;

			energy_ += M_1 * Hexch;
		}

		energy += energy_;
	}
}
