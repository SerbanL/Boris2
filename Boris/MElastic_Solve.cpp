#include "stdafx.h"
#include "MElastic.h"

#ifdef MODULE_COMPILATION_MELASTIC

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#include "SuperMesh.h"

#include "MElastic_Boundaries.h"

#if COMPILECUDA == 1
#include "MElasticCUDA.h"
#endif

double MElastic::Calculate_MElastic_Field(void)
{
	double energy = 0;

	//Energy density formula:

	//Let e1, e2, e3 be the orthogonal cubic axes, e.g. x, y, z in the simplest case

	//diagonal terms in stress tensor :
	//let Sd = (exx, eyy, ezz) be the diagonal terms.

	//off-diagonal terms in stress tensor (remember it is symmetric so we only have 3 of these):
	//let Sod = (eyz, exz, exy)

	//Then:

	//energy density from diagonal terms:
	//emel_d = B1 * [ (m.e1)^2*(Sd.e1) + (m.e2)^2*(Sd.e2) + (m.e3)^2*(Sd.e3) ]

	//energy density from off-diagonal terms (but remember this is zero for uniform strains):
	//emel_od = 2 * B2 * [ (m.e1)*(m.e2)*(Sod.e3) + (m.e1)*(m.e3)*(Sod.e2) + (m.e2)*(m.e3)*(Sod.e1) ]

	//Obtain Hmel as usual : Hmel = -1/mu0Ms * de/dm, where exchange stiffness-related terms are ignored here.

	if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
				DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
				DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
				DBL2 MEc = pMesh->MEc;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM, pMesh->MEc, MEc, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

				DBL3 position = pMesh->M.cellidx_to_position(idx);
				//xx, yy, zz
				DBL3 Sd = pMesh->strain_diag[position];
				//yz, xz, xy
				DBL3 Sod = pMesh->strain_odiag[position];

				//normalised magnetization
				//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

				DBL3 mA = DBL3(pMesh->M[idx] * mcanis_ea1, pMesh->M[idx] * mcanis_ea2, pMesh->M[idx] * mcanis_ea3) / Ms_AFM.i;
				DBL3 mB = DBL3(pMesh->M2[idx] * mcanis_ea1, pMesh->M2[idx] * mcanis_ea2, pMesh->M2[idx] * mcanis_ea3) / Ms_AFM.j;
				
				Sd = DBL3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
				Sod = DBL3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

				DBL3 Hmel_1_A = (-2.0 * MEc.i / (MU0 * Ms_AFM.i)) * DBL3(
					mA.x*Sd.x*mcanis_ea1.x + mA.y*Sd.y*mcanis_ea2.x + mA.z*Sd.z*mcanis_ea3.x,
					mA.x*Sd.x*mcanis_ea1.y + mA.y*Sd.y*mcanis_ea2.y + mA.z*Sd.z*mcanis_ea3.y,
					mA.x*Sd.x*mcanis_ea1.z + mA.y*Sd.y*mcanis_ea2.z + mA.z*Sd.z*mcanis_ea3.z);

				DBL3 Hmel_2_A = (-2.0 * MEc.j / (MU0 * Ms_AFM.i)) * DBL3(
					Sod.z * (mcanis_ea1.x*mA.y + mcanis_ea2.x*mA.x) + Sod.y * (mcanis_ea1.x*mA.z + mcanis_ea3.x*mA.x) + Sod.x * (mcanis_ea2.x*mA.z + mcanis_ea3.x*mA.y),
					Sod.z * (mcanis_ea1.y*mA.y + mcanis_ea2.y*mA.x) + Sod.y * (mcanis_ea1.y*mA.z + mcanis_ea3.y*mA.x) + Sod.x * (mcanis_ea2.y*mA.z + mcanis_ea3.y*mA.y),
					Sod.z * (mcanis_ea1.z*mA.y + mcanis_ea2.z*mA.x) + Sod.y * (mcanis_ea1.z*mA.z + mcanis_ea3.z*mA.x) + Sod.x * (mcanis_ea2.z*mA.z + mcanis_ea3.z*mA.y));

				DBL3 Hmel_1_B = (-2.0 * MEc.i / (MU0 * Ms_AFM.j)) * DBL3(
					mB.x*Sd.x*mcanis_ea1.x + mB.y*Sd.y*mcanis_ea2.x + mB.z*Sd.z*mcanis_ea3.x,
					mB.x*Sd.x*mcanis_ea1.y + mB.y*Sd.y*mcanis_ea2.y + mB.z*Sd.z*mcanis_ea3.y,
					mB.x*Sd.x*mcanis_ea1.z + mB.y*Sd.y*mcanis_ea2.z + mB.z*Sd.z*mcanis_ea3.z);

				DBL3 Hmel_2_B = (-2.0 * MEc.j / (MU0 * Ms_AFM.j)) * DBL3(
					Sod.z * (mcanis_ea1.x*mB.y + mcanis_ea2.x*mB.x) + Sod.y * (mcanis_ea1.x*mB.z + mcanis_ea3.x*mB.x) + Sod.x * (mcanis_ea2.x*mB.z + mcanis_ea3.x*mB.y),
					Sod.z * (mcanis_ea1.y*mB.y + mcanis_ea2.y*mB.x) + Sod.y * (mcanis_ea1.y*mB.z + mcanis_ea3.y*mB.x) + Sod.x * (mcanis_ea2.y*mB.z + mcanis_ea3.y*mB.y),
					Sod.z * (mcanis_ea1.z*mB.y + mcanis_ea2.z*mB.x) + Sod.y * (mcanis_ea1.z*mB.z + mcanis_ea3.z*mB.x) + Sod.x * (mcanis_ea2.z*mB.z + mcanis_ea3.z*mB.y));

				pMesh->Heff[idx] += Hmel_1_A + Hmel_2_A;
				pMesh->Heff2[idx] += Hmel_1_B + Hmel_2_B;

				energy += pMesh->M[idx] * (Hmel_1_A + Hmel_2_A) / 2 + pMesh->M2[idx] * (Hmel_1_B + Hmel_2_B) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Hmel_1_A + Hmel_2_A;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Hmel_1_B + Hmel_2_B;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * (Hmel_1_A + Hmel_2_A) / 2;
				if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * pMesh->M2[idx] * (Hmel_1_B + Hmel_2_B) / 2;
			}
		}
	}

	else {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
				DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
				DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
				DBL2 MEc = pMesh->MEc;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->MEc, MEc, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

				DBL3 position = pMesh->M.cellidx_to_position(idx);
				//xx, yy, zz
				DBL3 Sd = pMesh->strain_diag[position];
				//yz, xz, xy
				DBL3 Sod = pMesh->strain_odiag[position];

				//normalised magnetization
				//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

				DBL3 m = DBL3(pMesh->M[idx] * mcanis_ea1, pMesh->M[idx] * mcanis_ea2, pMesh->M[idx] * mcanis_ea3) / Ms;
				Sd = DBL3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
				Sod = DBL3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

				DBL3 Hmel_1 = (-2.0 * MEc.i / (MU0 * Ms)) * DBL3(
					m.x*Sd.x*mcanis_ea1.x + m.y*Sd.y*mcanis_ea2.x + m.z*Sd.z*mcanis_ea3.x,
					m.x*Sd.x*mcanis_ea1.y + m.y*Sd.y*mcanis_ea2.y + m.z*Sd.z*mcanis_ea3.y,
					m.x*Sd.x*mcanis_ea1.z + m.y*Sd.y*mcanis_ea2.z + m.z*Sd.z*mcanis_ea3.z);

				DBL3 Hmel_2 = (-2.0 * MEc.j / (MU0 * Ms)) * DBL3(
					Sod.z * (mcanis_ea1.x*m.y + mcanis_ea2.x*m.x) + Sod.y * (mcanis_ea1.x*m.z + mcanis_ea3.x*m.x) + Sod.x * (mcanis_ea2.x*m.z + mcanis_ea3.x*m.y),
					Sod.z * (mcanis_ea1.y*m.y + mcanis_ea2.y*m.x) + Sod.y * (mcanis_ea1.y*m.z + mcanis_ea3.y*m.x) + Sod.x * (mcanis_ea2.y*m.z + mcanis_ea3.y*m.y),
					Sod.z * (mcanis_ea1.z*m.y + mcanis_ea2.z*m.x) + Sod.y * (mcanis_ea1.z*m.z + mcanis_ea3.z*m.x) + Sod.x * (mcanis_ea2.z*m.z + mcanis_ea3.z*m.y));

				pMesh->Heff[idx] += Hmel_1 + Hmel_2;

				energy += pMesh->M[idx] * (Hmel_1 + Hmel_2);

				if (Module_Heff.linear_size()) Module_Heff[idx] = Hmel_1 + Hmel_2;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * (Hmel_1 + Hmel_2) / 2;
			}
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//----------------------------------------------- Computational Helpers

void MElastic::Iterate_Elastic_Solver_Velocity(double dT)
{
	DBL3 h_m = pMesh->u_disp.h;
	SZ3 n_m = pMesh->u_disp.n;

	double time = pSMesh->GetStageTime();

	//1a. Update velocity

	//loop over vertices
	for (int k = 0; k < n_m.k + 1; k++) {
#pragma omp parallel for
		for (int j = 0; j < n_m.j + 1; j++) {
			for (int i = 0; i < n_m.i + 1; i++) {

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

				double density = pMesh->density;
				double mdamping = pMesh->mdamping;
				pMesh->update_parameters_scoarse(idx_u, pMesh->density, density, pMesh->mdamping, mdamping);

				INT3 ijk = INT3(i, j, k);

				//external forces on different faces (keep track separately in case an edge cell is excited simultaneously by 2 or more external forces
				DBL3 Fext_xface = DBL3(), Fext_yface = DBL3(), Fext_zface = DBL3();

				//is there an external force? If so, get it, otherwise it will be zero
				if (
					((i == 0 || i == n_m.i) && pMesh->strain_diag.is_dirichlet_x(idx_u)) || 
					((j == 0 || j == n_m.j) && pMesh->strain_diag.is_dirichlet_y(idx_u)) || 
					((k == 0 || k == n_m.k) && pMesh->strain_diag.is_dirichlet_z(idx_u))) {

					//search through all available surfaces to get external force
					for (int sidx = 0; sidx < external_stress_surfaces.size(); sidx++) {

						int orientation = external_stress_surfaces[sidx].contains(ijk_u);
						if (orientation) {

							switch (abs(orientation)) {

								//x face
							case 1:
								Fext_xface = external_stress_surfaces[sidx].get_ext_force_edges(ijk, time);
								break;

								//y face
							case 2:
								Fext_yface = external_stress_surfaces[sidx].get_ext_force_edges(ijk, time);
								break;

								//z face
							case 3:
								Fext_zface = external_stress_surfaces[sidx].get_ext_force_edges(ijk, time);
								break;
							};
						}
					}
				}

				//update vx
				if (i < n_m.i) {

					//set zero at fixed faces (for vx only y and z faces are applicable)
					if (((j == 0 || j == n_m.j) && pMesh->u_disp.is_dirichlet_y(idx_u)) || ((k == 0 || k == n_m.k) && pMesh->u_disp.is_dirichlet_z(idx_u))) {

						vx[ijk] = 0.0;
					}
					else {

						int njend = (j < n_m.j);
						int nkend = (k < n_m.k);

						//check for required axis normal faces being present
						bool zface_u = 
							j < n_m.j &&
							(pMesh->u_disp.is_not_empty(idx_u) || 
							(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));
							
						bool zface_l = 
							j > 0 && 
							(pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x) || 
							(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));
							
						bool yface_u =
							k < n_m.k &&
							(pMesh->u_disp.is_not_empty(idx_u) ||
							(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x)));
								
						bool yface_l = 
							k > 0 &&
							(pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
							(j > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

						//at least one face is required, otherwise velocity must be zero
						if (zface_u || zface_l || yface_u || yface_l) {

							double dsxx_dx = 0.0, dsxy_dy = 0.0, dsxz_dz = 0.0;

							//always interior
							dsxx_dx = (sdd[INT3(i + 1, j, k)].x - sdd[ijk].x) / h_m.x;

							//interior
							if (zface_u && zface_l) dsxy_dy = (sxy[ijk] - sxy[INT3(i, j - 1, k)]) / h_m.y;
							else if (zface_l) dsxy_dy = (Fext_yface.x - sxy[INT3(i, j - 1, k)]) / (h_m.y / 2);
							else if (zface_u) dsxy_dy = (sxy[ijk] - Fext_yface.x) / (h_m.y / 2);

							//interior
							if (yface_u && yface_l) dsxz_dz = (sxz[ijk] - sxz[INT3(i, j, k - 1)]) / h_m.z;
							else if (yface_l) dsxz_dz = (Fext_zface.x - sxz[INT3(i, j, k - 1)]) / (h_m.z / 2);
							else if (yface_u) dsxz_dz = (sxz[ijk] - Fext_zface.x) / (h_m.z / 2);

							vx[ijk] += dT * (dsxx_dx + dsxy_dy + dsxz_dz - mdamping * vx[ijk]) / density;
						}
						else vx[ijk] = 0.0;
					}
				}

				//update vy
				if (j < n_m.j) {

					//set zero at fixed faces (for vy only x and z faces are applicable)
					if (((i == 0 || i == n_m.i) && pMesh->u_disp.is_dirichlet_x(idx_u)) || ((k == 0 || k == n_m.k) && pMesh->u_disp.is_dirichlet_z(idx_u))) {

						vy[ijk] = 0.0;
					}
					else {

						int niend = (i < n_m.i);
						int nkend = (k < n_m.k);

						//check for required axis normal faces being present
						bool zface_u =
							i < n_m.i &&
							(pMesh->u_disp.is_not_empty(idx_u) ||
							(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

						bool zface_l =
							i > 0 &&
							(pMesh->u_disp.is_not_empty(idx_u - niend) ||
							(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

						bool xface_u =
							k < n_m.k &&
							(pMesh->u_disp.is_not_empty(idx_u) ||
							(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend)));

						bool xface_l =
							k > 0 &&
							(pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
							(i > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

						//at least one face is required, otherwise velocity must be zero
						if (zface_u || zface_l || xface_u || xface_l) {

							double dsxy_dx = 0.0, dsyy_dy = 0.0, dsyz_dz = 0.0;

							//always interior
							dsyy_dy = (sdd[INT3(i, j + 1, k)].y - sdd[ijk].y) / h_m.y;

							//interior
							if (zface_u && zface_l) dsxy_dx = (sxy[ijk] - sxy[INT3(i - 1, j, k)]) / h_m.x;
							else if (zface_l) dsxy_dx = (Fext_xface.y - sxy[INT3(i - 1, j, k)]) / (h_m.x / 2);
							else if (zface_u) dsxy_dx = (sxy[ijk] - Fext_xface.y) / (h_m.x / 2);

							//interior
							if (xface_u && xface_l) dsyz_dz = (syz[ijk] - syz[INT3(i, j, k - 1)]) / h_m.z;
							else if (xface_l) dsyz_dz = (Fext_zface.y - syz[INT3(i, j, k - 1)]) / (h_m.z / 2);
							else if (xface_u) dsyz_dz = (syz[ijk] - Fext_zface.y) / (h_m.z / 2);

							vy[ijk] += dT * (dsxy_dx + dsyy_dy + dsyz_dz - mdamping * vy[ijk]) / density;
						}
						else vy[ijk] = 0.0;
					}
				}

				//update vz
				if (k < n_m.k) {

					//set zero at fixed faces (for vz only x and y faces are applicable)
					if (((i == 0 || i == n_m.i) && pMesh->u_disp.is_dirichlet_x(idx_u)) || ((j == 0 || j == n_m.j) && pMesh->u_disp.is_dirichlet_y(idx_u))) {

						vz[ijk] = 0.0;
					}
					else {

						int niend = (i < n_m.i);
						int njend = (j < n_m.j);

						//check for required axis normal faces being present
						bool yface_u =
							i < n_m.i &&
							(pMesh->u_disp.is_not_empty(idx_u) ||
							(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x)));

						bool yface_l =
							i > 0 &&
							(pMesh->u_disp.is_not_empty(idx_u - niend) ||
							(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

						bool xface_u =
							j < n_m.j &&
							(pMesh->u_disp.is_not_empty(idx_u) ||
							(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend)));

						bool xface_l =
							j > 0 &&
							(pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x) ||
							(i > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

						//at least one face is required, otherwise velocity must be zero
						if (yface_u || yface_l || xface_u || xface_l) {

							double dsxz_dx = 0.0, dsyz_dy = 0.0, dszz_dz = 0.0;

							//always interior
							dszz_dz = (sdd[INT3(i, j, k + 1)].z - sdd[ijk].z) / h_m.z;

							//interior
							if (yface_u && yface_l) dsxz_dx = (sxz[ijk] - sxz[INT3(i - 1, j, k)]) / h_m.x;
							else if (yface_l) dsxz_dx = (Fext_xface.z - sxz[INT3(i - 1, j, k)]) / (h_m.x / 2);
							else if (yface_u) dsxz_dx = (sxz[ijk] - Fext_xface.z) / (h_m.x / 2);

							//interior
							if (xface_u && xface_l) dsyz_dy = (syz[ijk] - syz[INT3(i, j - 1, k)]) / h_m.y;
							else if (xface_l) dsyz_dy = (Fext_yface.z - syz[INT3(i, j - 1, k)]) / (h_m.y / 2);
							else if (xface_u) dsyz_dy = (syz[ijk] - Fext_yface.z) / (h_m.y / 2);

							vz[ijk] += dT * (dsxz_dx + dsyz_dy + dszz_dz - mdamping * vz[ijk]) / density;
						}
						else vz[ijk] = 0.0;
					}
				}

				//show u_disp approximation for visualization (not used for further computations so it's fine)
				if (i < n_m.i && j < n_m.j && k < n_m.k) {

					pMesh->u_disp[idx_u] += dT * DBL3(vx[ijk], vy[ijk], vz[ijk]);
				}
			}
		}
	}
}

void MElastic::Iterate_Elastic_Solver_Stress(double dT)
{
	DBL3 h_m = pMesh->u_disp.h;
	SZ3 n_m = pMesh->u_disp.n;

	double time = pSMesh->GetStageTime();

	//1b. Update stress

	for (int k = 0; k < n_m.k + 1; k++) {
#pragma omp parallel for
		for (int j = 0; j < n_m.j + 1; j++) {
			for (int i = 0; i < n_m.i + 1; i++) {

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

				DBL3 cC = pMesh->cC;
				pMesh->update_parameters_scoarse(idx_u, pMesh->cC, cC);

				INT3 ijk = INT3(i, j, k);

				//external forces on different faces (keep track separately in case an edge cell is excited simultaneously by 2 or more external forces
				//also only need the component normal to the surface, at the vertex
				double Fext_xface = 0.0, Fext_yface = 0.0, Fext_zface = 0.0;

				//is there an external force? If so, get it, otherwise it will be zero
				if (
					((i == 0 || i == n_m.i) && pMesh->strain_diag.is_dirichlet_x(idx_u)) ||
					((j == 0 || j == n_m.j) && pMesh->strain_diag.is_dirichlet_y(idx_u)) ||
					((k == 0 || k == n_m.k) && pMesh->strain_diag.is_dirichlet_z(idx_u))) {

					//search through all available surfaces to get external force
					for (int sidx = 0; sidx < external_stress_surfaces.size(); sidx++) {

						int orientation = external_stress_surfaces[sidx].contains(ijk_u);
						if (orientation) {

							switch (abs(orientation)) {

								//x face
							case 1:
								Fext_xface = external_stress_surfaces[sidx].get_ext_force_vertices(ijk, time);
								break;

								//y face
							case 2:
								Fext_yface = external_stress_surfaces[sidx].get_ext_force_vertices(ijk, time);
								break;

								//z face
							case 3:
								Fext_zface = external_stress_surfaces[sidx].get_ext_force_vertices(ijk, time);
								break;
							};
						}
					}
				}

				//update sxx, syy, szz
				int niend = (i < n_m.i);
				int njend = (j < n_m.j);
				int nkend = (k < n_m.k);

				//check if required edges are present
				bool xedge_u = 
					i < n_m.i &&
					(pMesh->u_disp.is_not_empty(idx_u) ||
					(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
					(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x)) || 
					(j > 0 && k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

				bool xedge_l =
					i > 0 &&
					(pMesh->u_disp.is_not_empty(idx_u - niend) ||
					(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)) ||
					(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x - niend)) ||
					(j > 0 && k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x - niend)));

				bool yedge_u =
					j < n_m.j &&
					(pMesh->u_disp.is_not_empty(idx_u) ||
					(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
					(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend)) ||
					(i > 0 && k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

				bool yedge_l =
					j > 0 &&
					(pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x) ||
					(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)) ||
					(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend - njend * n_m.x)) ||
					(i > 0 && k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend - njend * n_m.x)));

				bool zedge_u =
					k < n_m.k &&
					(pMesh->u_disp.is_not_empty(idx_u) ||
					(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
					(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend)) ||
					(i > 0 && j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

				bool zedge_l =
					k > 0 &&
					(pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
					(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x - nkend * n_m.x*n_m.y)) ||
					(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend - nkend * n_m.x*n_m.y)) ||
					(i > 0 && j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x - niend - nkend * n_m.x*n_m.y)));

				//check for fixed faces at ends
				bool xfixed_l = (i == 0 && pMesh->u_disp.is_dirichlet_px(idx_u));
				bool xfixed_u = (i == n_m.i && pMesh->u_disp.is_dirichlet_nx(idx_u));

				bool yfixed_l = (j == 0 && pMesh->u_disp.is_dirichlet_py(idx_u));
				bool yfixed_u = (j == n_m.j && pMesh->u_disp.is_dirichlet_ny(idx_u));

				bool zfixed_l = (k == 0 && pMesh->u_disp.is_dirichlet_pz(idx_u));
				bool zfixed_u = (k == n_m.k && pMesh->u_disp.is_dirichlet_nz(idx_u));

				double dvx_dx = 0.0;

				//interior
				if (xedge_u && xedge_l) dvx_dx = (vx[ijk] - vx[INT3(i - 1, j, k)]) / h_m.x;
				else if (xedge_l) {

					//is it a fixed face or free?
					if (xfixed_u) dvx_dx = -vx[INT3(i - 1, j, k)] / (h_m.x / 2);
					else {

						//approximate zero
						dvx_dx = 0.0;
					}
				}
				else if (xedge_u) {

					//is it a fixed face or free?
					if (xfixed_l) dvx_dx = vx[ijk] / (h_m.x / 2);
					else {

						//approximate zero
						dvx_dx = 0.0;
					}
				}

				double dvy_dy = 0.0;

				//interior
				if (yedge_u && yedge_l) dvy_dy = (vy[ijk] - vy[INT3(i, j - 1, k)]) / h_m.y;
				//at +y face
				else if (yedge_l) {

					//is it a fixed face or free?
					if (yfixed_u) dvy_dy = -vy[INT3(i, j - 1, k)] / (h_m.y / 2);
					else {

						//approximate zero
						dvy_dy = 0.0;
					}
				}
				//at -y face
				else if (yedge_u) {

					//is it a fixed face or free?
					if (yfixed_l) dvy_dy = vy[ijk] / (h_m.y / 2);
					else {

						//approximate zero
						dvy_dy = 0.0;
					}
				}

				double dvz_dz = 0.0;

				//interior
				if (zedge_u && zedge_l) dvz_dz = (vz[ijk] - vz[INT3(i, j, k - 1)]) / h_m.z;
				//at +z face
				else if (zedge_l) {

					//is it a fixed face or free?
					if (zfixed_u) dvz_dz = -vz[INT3(i, j, k - 1)] / (h_m.z / 2);
					else {

						//approximate zero
						dvz_dz = 0.0;
					}
				}
				//at -z face
				else if (zedge_u) {

					//is it a fixed face or free?
					if (zfixed_l) dvz_dz = vz[ijk] / (h_m.z / 2);
					else {

						//approximate zero
						dvz_dz = 0.0;
					}
				}

				if ((!xedge_u && !xfixed_u) || (!xedge_l && !xfixed_l)) sdd[ijk].x = Fext_xface;
				else sdd[ijk].x += dT * (cC.i * dvx_dx + cC.j * (dvy_dy + dvz_dz));

				if ((!yedge_u && !yfixed_u) || (!yedge_l && !yfixed_l)) sdd[ijk].y = Fext_yface;
				else sdd[ijk].y += dT * (cC.i * dvy_dy + cC.j * (dvx_dx + dvz_dz));

				if ((!zedge_u && !zfixed_u) || (!zedge_l && !zfixed_l)) sdd[ijk].z = Fext_zface;
				else sdd[ijk].z += dT * (cC.i * dvz_dz + cC.j * (dvx_dx + dvy_dy));

				//update sxy
				if (i < n_m.i && j < n_m.j) {

					bool zface =
						(pMesh->u_disp.is_not_empty(idx_u) ||
						(k > 0 && pMesh->u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

					if (zface) {

						double dvx_dy = (vx[INT3(i, j + 1, k)] - vx[ijk]) / h_m.y;
						double dvy_dx = (vy[INT3(i + 1, j, k)] - vy[ijk]) / h_m.x;

						sxy[ijk] += dT * cC.k * (dvx_dy + dvy_dx) / 2;
					}
					else sxy[ijk] = 0.0;
				}

				//update sxz
				if (i < n_m.i && k < n_m.k) {

					bool yface =
						(pMesh->u_disp.is_not_empty(idx_u) ||
						(j > 0 && pMesh->u_disp.is_not_empty(idx_u - njend * n_m.x)));

					if (yface) {

						double dvx_dz = (vx[INT3(i, j, k + 1)] - vx[ijk]) / h_m.z;
						double dvz_dx = (vz[INT3(i + 1, j, k)] - vz[ijk]) / h_m.x;

						sxz[ijk] += dT * cC.k * (dvx_dz + dvz_dx) / 2;
					}
					else sxz[ijk] = 0.0;
				}

				//update syz
				if (j < n_m.j && k < n_m.k) {

					bool xface =
						(pMesh->u_disp.is_not_empty(idx_u) ||
						(i > 0 && pMesh->u_disp.is_not_empty(idx_u - niend)));

					if (xface) {

						double dvy_dz = (vy[INT3(i, j, k + 1)] - vy[ijk]) / h_m.z;
						double dvz_dy = (vz[INT3(i, j + 1, k)] - vz[ijk]) / h_m.y;

						syz[ijk] += dT * cC.k * (dvy_dz + dvz_dy) / 2;
					}
					else syz[ijk] = 0.0;
				}
			}
		}
	}
}

//updatre strain from stress
void MElastic::Calculate_Strain_From_Stress(void)
{
	SZ3 n_m = pMesh->n_m;

	//1c. Get strain from stress (cannot lump it in with step b since need to average over several vertices and faces to get cell-centred stress values)
	for (int k = 0; k < n_m.k; k++) {
#pragma omp parallel for
		for (int j = 0; j < n_m.j; j++) {
			for (int i = 0; i < n_m.i; i++) {

				int idx_u = i + j * n_m.x + k * n_m.x * n_m.y;

				//no need to update stress values for empty cells (so keep them zero)
				if (pMesh->u_disp.is_empty(idx_u)) continue;

				DBL3 cC = pMesh->cC;
				pMesh->update_parameters_scoarse(idx_u, pMesh->cC, cC);

				INT3 ijk = INT3(i, j, k);
				//index for diagonal stress (sxx, syy, szz)
				int idx_sd = i + j * (n_m.i + 1) + k * (n_m.i + 1) * (n_m.j + 1);

				//invert the elastic coefficients matrix for diagonal stress-strain terms
				double det = cC.i*cC.i*cC.i + 2 * cC.j*cC.j*cC.j - 3 * cC.i*cC.j*cC.j;
				double diag = (cC.i*cC.i - cC.j*cC.j) / det;
				double odiag = -(cC.i*cC.j - cC.j*cC.j) / det;

				//strain - diagonal. Use average of the 8 vertex diagonal stress values to get cell-centre value.
				double sig_d_x = (
					sdd[idx_sd].x + sdd[idx_sd + 1].x + sdd[idx_sd + n_m.x + 1].x + sdd[idx_sd + n_m.x + 2].x +
					sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1)].x + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + 1].x + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 1].x + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 2].x) / 8;

				double sig_d_y = (
					sdd[idx_sd].y + sdd[idx_sd + 1].y + sdd[idx_sd + n_m.x + 1].y + sdd[idx_sd + n_m.x + 2].y +
					sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1)].y + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + 1].y + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 1].y + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 2].y) / 8;

				double sig_d_z = (
					sdd[idx_sd].z + sdd[idx_sd + 1].z + sdd[idx_sd + n_m.x + 1].z + sdd[idx_sd + n_m.x + 2].z +
					sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1)].z + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + 1].z + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 1].z + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 2].z) / 8;

				pMesh->strain_diag[idx_u].x = diag * sig_d_x + odiag * (sig_d_y + sig_d_z);
				pMesh->strain_diag[idx_u].y = diag * sig_d_y + odiag * (sig_d_x + sig_d_z);
				pMesh->strain_diag[idx_u].z = diag * sig_d_z + odiag * (sig_d_x + sig_d_y);

				//strain - off-diagonal (yz, xz, xy). Use average of 2 face off-diagonal stress values to get cell-centre value.
				pMesh->strain_odiag[idx_u].x = (syz[ijk] + syz[ijk + INT3(1, 0, 0)]) / (2 * cC.k);
				pMesh->strain_odiag[idx_u].y = (sxz[ijk] + sxz[ijk + INT3(0, 1, 0)]) / (2 * cC.k);
				pMesh->strain_odiag[idx_u].z = (sxy[ijk] + sxy[ijk + INT3(0, 0, 1)]) / (2 * cC.k);
			}
		}
	}
}

//---------------------------------------------- CMBND

void MElastic::make_velocity_continuous(
	CMBNDInfo& contact,
	VEC<double>& vx_sec, VEC<double>& vy_sec, VEC<double>& vz_sec, VEC_VC<DBL3>& u_disp_sec)
{
	const Box& cb = contact.cells_box;
	int axis;
	if (contact.cell_shift.x) axis = 1;
	else if (contact.cell_shift.y) axis = 2;
	else axis = 3;

	//+/-x normal face
	if (axis == 1) {

		int i = (contact.IsPrimaryTop() ? cb.s.i : cb.e.i);

		double spacing = pMesh->u_disp.h.x + u_disp_sec.h.x;
		double wpri = 1.0 - pMesh->u_disp.h.x / spacing;
		double wsec = 1.0 - u_disp_sec.h.x / spacing;

		for (int k = cb.s.k; k < cb.e.k + 1; k++) {
#pragma omp parallel for
			for (int j = cb.s.j; j < cb.e.j + 1; j++) {

				INT3 ijk = INT3(i, j, k);

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < pMesh->u_disp.n.i ? i : pMesh->u_disp.n.i - 1, j < pMesh->u_disp.n.j ? j : pMesh->u_disp.n.j - 1, k < pMesh->u_disp.n.k ? k : pMesh->u_disp.n.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * pMesh->u_disp.n.x + ijk_u.k * pMesh->u_disp.n.x * pMesh->u_disp.n.y;

				if (pMesh->u_disp.is_empty(idx_u) || pMesh->u_disp.is_not_cmbnd(idx_u)) continue;

				//absolute position of interface vertex
				DBL3 abs_pos = (pMesh->u_disp.h & ijk) + pMesh->u_disp.rect.s;

				if (contact.IsPrimaryTop()) {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3(-u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//vy
					if (j < cb.e.j) {

						//value at interface : interpolate between primary and secondary
						vy[INT3(0, j, k)] = vy[INT3(1, j, k)] * wpri + vy_sec[ijk_sec + INT3(0, 0, k == cb.e.k)] * wsec;
					}

					//vz
					if (k < cb.e.k) {

						//value at interface : interpolate between primary and secondary
						vz[INT3(0, j, k)] = vz[INT3(1, j, k)] * wpri + vz_sec[ijk_sec + INT3(0, j == cb.e.j, 0)] * wsec;
					}
				}
				else {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3(u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//vy
					if (j < cb.e.j) {

						vy[INT3(vy.n.x - 1, j, k)] = vy[INT3(vy.n.x - 2, j, k)] * wpri + vy_sec[ijk_sec + INT3(1, 0, k == cb.e.k)] * wsec;
					}

					//vz
					if (k < cb.e.k) {

						vz[INT3(vz.n.x - 1, j, k)] = vz[INT3(vz.n.x - 2, j, k)] * wpri + vz_sec[ijk_sec + INT3(1, j == cb.e.j, 0)] * wsec;
					}
				}
			}
		}
	}

	//+/-y normal face
	else if (axis == 2) {

		int j = (contact.IsPrimaryTop() ? cb.s.j : cb.e.j);

		double spacing = pMesh->u_disp.h.y + u_disp_sec.h.y;
		double wpri = 1.0 - pMesh->u_disp.h.y / spacing;
		double wsec = 1.0 - u_disp_sec.h.y / spacing;

		for (int k = cb.s.k; k < cb.e.k + 1; k++) {
#pragma omp parallel for
			for (int i = cb.s.i; i < cb.e.i + 1; i++) {

				INT3 ijk = INT3(i, j, k);

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < pMesh->u_disp.n.i ? i : pMesh->u_disp.n.i - 1, j < pMesh->u_disp.n.j ? j : pMesh->u_disp.n.j - 1, k < pMesh->u_disp.n.k ? k : pMesh->u_disp.n.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * pMesh->u_disp.n.x + ijk_u.k * pMesh->u_disp.n.x * pMesh->u_disp.n.y;

				if (pMesh->u_disp.is_empty(idx_u) || pMesh->u_disp.is_not_cmbnd(idx_u)) continue;

				//absolute position of interface vertex
				DBL3 abs_pos = (pMesh->u_disp.h & ijk) + pMesh->u_disp.rect.s;

				if (contact.IsPrimaryTop()) {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, -u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//vx
					if (i < cb.e.i) {

						vx[INT3(i, 0, k)] = vx[INT3(i, 1, k)] * wpri + vx_sec[ijk_sec + INT3(0, 0, k == cb.e.k)] * wsec;
					}

					//vz
					if (k < cb.e.k) {

						vz[INT3(i, 0, k)] = vz[INT3(i, 1, k)] * wpri + vz_sec[ijk_sec + INT3(i == cb.e.i, 0, 0)] * wsec;
					}
				}
				else {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//vx
					if (i < cb.e.i) {

						vx[INT3(i, vx.n.y - 1, k)] = vx[INT3(i, vx.n.y - 2, k)] * wpri + vx_sec[ijk_sec + INT3(0, 1, k == cb.e.k)] * wsec;
					}

					//vz
					if (k < cb.e.k) {

						vz[INT3(i, vz.n.y - 1, k)] = vz[INT3(i, vz.n.y - 1, k)] * wpri + vz_sec[ijk_sec + INT3(i == cb.e.i, 1, 0)] * wsec;
					}
				}
			}
		}
	}

	//+/-z normal face
	else if (axis == 3) {

		int k = (contact.IsPrimaryTop() ? cb.s.k : cb.e.k);

		double spacing = pMesh->u_disp.h.z + u_disp_sec.h.z;
		double wpri = 1.0 - pMesh->u_disp.h.z / spacing;
		double wsec = 1.0 - u_disp_sec.h.z / spacing;

	#pragma omp parallel for
		for (int j = cb.s.j; j < cb.e.j + 1; j++) {
			for (int i = cb.s.i; i < cb.e.i + 1; i++) {

				INT3 ijk = INT3(i, j, k);

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < pMesh->u_disp.n.i ? i : pMesh->u_disp.n.i - 1, j < pMesh->u_disp.n.j ? j : pMesh->u_disp.n.j - 1, k < pMesh->u_disp.n.k ? k : pMesh->u_disp.n.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * pMesh->u_disp.n.x + ijk_u.k * pMesh->u_disp.n.x * pMesh->u_disp.n.y;

				if (pMesh->u_disp.is_empty(idx_u) || pMesh->u_disp.is_not_cmbnd(idx_u)) continue;

				//absolute position of interface vertex
				DBL3 abs_pos = (pMesh->u_disp.h & ijk) + pMesh->u_disp.rect.s;

				if (contact.IsPrimaryTop()) {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, -u_disp_sec.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//vx
					if (i < cb.e.i) {

						vx[INT3(i, j, 0)] = vx[INT3(i, j, 1)] * wpri + vx_sec[ijk_sec + INT3(0, j == cb.e.j, 0)] * wsec;
					}

					//vy
					if (j < cb.e.j) {

						vy[INT3(i, j, 0)] = vy[INT3(i, j, 1)] * wpri + vy_sec[ijk_sec + INT3(i == cb.e.i, 0, 0)] * wsec;
					}
				}
				else {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, u_disp_sec.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//vx
					if (i < cb.e.i) {

						vx[INT3(i, j, vx.n.z - 1)] = vx[INT3(i, j, vx.n.z - 2)] * wpri + vx_sec[ijk_sec + INT3(0, j == cb.e.j, 1)] * wsec;
					}

					//vy
					if (j < cb.e.j) {

						vy[INT3(i, j, vy.n.z - 1)] = vy[INT3(i, j, vy.n.z - 2)] * wpri + vy_sec[ijk_sec + INT3(i == cb.e.i, 0, 1)] * wsec;
					}
				}
			}
		}
	}
}

void MElastic::make_stress_continuous(
	CMBNDInfo& contact,
	VEC<DBL3>& sdd_sec, VEC<double>& sxy_sec, VEC<double>& sxz_sec, VEC<double>& syz_sec,
	VEC_VC<DBL3>& u_disp_sec)
{
	const Box& cb = contact.cells_box;
	int axis;
	if (contact.cell_shift.x) axis = 1;
	else if (contact.cell_shift.y) axis = 2;
	else axis = 3;

	//+/-x normal face
	if (axis == 1) {

		int i = (contact.IsPrimaryTop() ? cb.s.i : cb.e.i);

		double spacing = pMesh->u_disp.h.x + u_disp_sec.h.x;
		double wpri = 1.0 - pMesh->u_disp.h.x / spacing;
		double wsec = 1.0 - u_disp_sec.h.x / spacing;

		for (int k = cb.s.k; k < cb.e.k + 1; k++) {
#pragma omp parallel for
			for (int j = cb.s.j; j < cb.e.j + 1; j++) {

				INT3 ijk = INT3(i, j, k);

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < pMesh->u_disp.n.i ? i : pMesh->u_disp.n.i - 1, j < pMesh->u_disp.n.j ? j : pMesh->u_disp.n.j - 1, k < pMesh->u_disp.n.k ? k : pMesh->u_disp.n.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * pMesh->u_disp.n.x + ijk_u.k * pMesh->u_disp.n.x * pMesh->u_disp.n.y;

				if (pMesh->u_disp.is_empty(idx_u) || pMesh->u_disp.is_not_cmbnd(idx_u)) continue;

				//absolute position of interface vertex
				DBL3 abs_pos = (pMesh->u_disp.h & ijk) + pMesh->u_disp.rect.s;

				if (contact.IsPrimaryTop()) {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3(-u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					sdd[INT3(0, j, k)] = sdd[INT3(1, j, k)] * wpri + sdd_sec[ijk_sec + INT3(0, j == cb.e.j, k == cb.e.k)] * wsec;

					//syz
					if (j < cb.e.j && k < cb.e.k) {

						syz[INT3(0, j, k)] = syz[INT3(1, j, k)] * wpri + syz_sec[ijk_sec] * wsec;
					}
				}
				else {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3(u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					sdd[INT3(sdd.n.x - 1, j, k)] = sdd[INT3(sdd.n.x - 2, j, k)] * wpri + sdd_sec[ijk_sec + INT3(1, j == cb.e.j, k == cb.e.k)] * wsec;

					//syz
					if (j < cb.e.j && k < cb.e.k) {

						syz[INT3(syz.n.x - 1, j, k)] = syz[INT3(syz.n.x - 2, j, k)] * wpri + syz_sec[ijk_sec + INT3(1, 0, 0)] * wsec;
					}
				}
			}
		}
	}

	//+/-y normal face
	if (axis == 2) {

		int j = (contact.IsPrimaryTop() ? cb.s.j : cb.e.j);

		double spacing = pMesh->u_disp.h.y + u_disp_sec.h.y;
		double wpri = 1.0 - pMesh->u_disp.h.y / spacing;
		double wsec = 1.0 - u_disp_sec.h.y / spacing;

		for (int k = cb.s.k; k < cb.e.k + 1; k++) {
#pragma omp parallel for
			for (int i = cb.s.i; i < cb.e.i + 1; i++) {

				INT3 ijk = INT3(i, j, k);

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < pMesh->u_disp.n.i ? i : pMesh->u_disp.n.i - 1, j < pMesh->u_disp.n.j ? j : pMesh->u_disp.n.j - 1, k < pMesh->u_disp.n.k ? k : pMesh->u_disp.n.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * pMesh->u_disp.n.x + ijk_u.k * pMesh->u_disp.n.x * pMesh->u_disp.n.y;

				if (pMesh->u_disp.is_empty(idx_u) || pMesh->u_disp.is_not_cmbnd(idx_u)) continue;

				//absolute position of interface vertex
				DBL3 abs_pos = (pMesh->u_disp.h & ijk) + pMesh->u_disp.rect.s;

				if (contact.IsPrimaryTop()) {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, -u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//sxx, syy, szz
					sdd[INT3(i, 0, k)] = sdd[INT3(i, 1, k)] * wpri + sdd_sec[ijk_sec + INT3(i == cb.e.i, 0, k == cb.e.k)] * wsec;

					//sxz
					if (i < cb.e.i && k < cb.e.k) {

						sxz[INT3(i, 0, k)] = sxz[INT3(i, 1, k)] * wpri + sxz_sec[ijk_sec] * wsec;
					}
				}
				else {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * pMesh->u_disp.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//sxx, syy, szz
					sdd[INT3(i, sdd.n.y - 1, k)] = sdd[INT3(i, sdd.n.y - 2, k)] * wpri + sdd_sec[ijk_sec + INT3(i == cb.e.i, 1, k == cb.e.k)] * wsec;

					//sxz
					if (i < cb.e.i && k < cb.e.k) {

						sxz[INT3(i, sxz.n.y - 1, k)] = sxz[INT3(i, sxz.n.y - 2, k)] * wpri + sxz_sec[ijk_sec + INT3(0, 1, 0)] * wsec;
					}
				}
			}
		}
	}

	//+/-z normal face
	if (axis == 3) {

		int k = (contact.IsPrimaryTop() ? cb.s.k : cb.e.k);

		double spacing = pMesh->u_disp.h.z + u_disp_sec.h.z;
		double wpri = 1.0 - pMesh->u_disp.h.z / spacing;
		double wsec = 1.0 - u_disp_sec.h.z / spacing;

#pragma omp parallel for
		for (int j = cb.s.j; j < cb.e.j + 1; j++) {
			for (int i = cb.s.i; i < cb.e.i + 1; i++) {

				INT3 ijk = INT3(i, j, k);

				//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
				INT3 ijk_u = INT3(i < pMesh->u_disp.n.i ? i : pMesh->u_disp.n.i - 1, j < pMesh->u_disp.n.j ? j : pMesh->u_disp.n.j - 1, k < pMesh->u_disp.n.k ? k : pMesh->u_disp.n.k - 1);
				int idx_u = ijk_u.i + ijk_u.j * pMesh->u_disp.n.x + ijk_u.k * pMesh->u_disp.n.x * pMesh->u_disp.n.y;

				if (pMesh->u_disp.is_empty(idx_u) || pMesh->u_disp.is_not_cmbnd(idx_u)) continue;

				//absolute position of interface vertex
				DBL3 abs_pos = (pMesh->u_disp.h & ijk) + pMesh->u_disp.rect.s;

				if (contact.IsPrimaryTop()) {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, -u_disp_sec.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//sxx, syy, szz
					sdd[INT3(i, j, 0)] = sdd[INT3(i, j, 1)] * wpri + sdd_sec[ijk_sec + INT3(i == cb.e.i, j == cb.e.j, 0)] * wsec;

					//sxy
					if (i < cb.e.i && j < cb.e.j) {

						sxy[INT3(i, j, 0)] = sxy[INT3(i, j, 1)] * wpri + sxy_sec[ijk_sec] * wsec;
					}
				}
				else {

					//absolute position in secondary, half a cellsize into it, and middle of primary cell face
					DBL3 abs_pos_sec = abs_pos + DBL3((i == cb.e.i ? -1 : +1) * pMesh->u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * pMesh->u_disp.h.y / 2, u_disp_sec.h.z / 2);
					//index of secondary cell just next to boundary
					INT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

					//sxx, syy, szz
					sdd[INT3(i, j, sdd.n.z - 1)] = sdd[INT3(i, j, sdd.n.z - 2)] * wpri + sdd_sec[ijk_sec + INT3(i == cb.e.i, j == cb.e.j, 1)] * wsec;

					//sxy
					if (i < cb.e.i && j < cb.e.j) {

						sxy[INT3(i, j, sxy.n.z - 1)] = sxy[INT3(i, j, sxy.n.z - 2)] * wpri + sxy_sec[ijk_sec + INT3(0, 0, 1)] * wsec;
					}
				}
			}
		}
	}
}

#endif
