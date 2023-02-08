#include "MElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "ManagedDiffEqFMCUDA.h"
#include "ManagedDiffEqAFMCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "MElastic_BoundariesCUDA.h"

//----------------------- Calculate_MElastic_Field KERNELS

__global__ void MElasticCUDA_UpdateField_FM(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
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
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuReal2 MEc = *cuMesh.pMEc;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pMEc, MEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

			cuReal3 position = M.cellidx_to_position(idx);
			//xx, yy, zz
			cuReal3 Sd = strain_diag[position];
			//yz, xz, xy
			cuReal3 Sod = strain_odiag[position];

			//normalised magnetization
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

			if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hmel_1 + Hmel_2;
			if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * M[idx] * (Hmel_1 + Hmel_2) / 2;
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void MElasticCUDA_UpdateField_AFM(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		if (M.is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuReal2 MEc = *cuMesh.pMEc;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pMEc, MEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

			cuReal3 position = M.cellidx_to_position(idx);
			//xx, yy, zz
			cuReal3 Sd = strain_diag[position];
			//yz, xz, xy
			cuReal3 Sod = strain_odiag[position];

			//normalised magnetization
			//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

			cuReal3 mA = cuReal3(M[idx] * mcanis_ea1, M[idx] * mcanis_ea2, M[idx] * mcanis_ea3) / Ms_AFM.i;
			cuReal3 mB = cuReal3(M2[idx] * mcanis_ea1, M2[idx] * mcanis_ea2, M2[idx] * mcanis_ea3) / Ms_AFM.j;

			Sd = cuReal3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
			Sod = cuReal3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

			cuReal3 Hmel_1_A = (-2.0 * MEc.i / (MU0 * Ms_AFM.i)) * cuReal3(
				mA.x*Sd.x*mcanis_ea1.x + mA.y*Sd.y*mcanis_ea2.x + mA.z*Sd.z*mcanis_ea3.x,
				mA.x*Sd.x*mcanis_ea1.y + mA.y*Sd.y*mcanis_ea2.y + mA.z*Sd.z*mcanis_ea3.y,
				mA.x*Sd.x*mcanis_ea1.z + mA.y*Sd.y*mcanis_ea2.z + mA.z*Sd.z*mcanis_ea3.z);

			cuReal3 Hmel_2_A = (-2.0 * MEc.j / (MU0 * Ms_AFM.i)) * cuReal3(
				Sod.z * (mcanis_ea1.x*mA.y + mcanis_ea2.x*mA.x) + Sod.y * (mcanis_ea1.x*mA.z + mcanis_ea3.x*mA.x) + Sod.x * (mcanis_ea2.x*mA.z + mcanis_ea3.x*mA.y),
				Sod.z * (mcanis_ea1.y*mA.y + mcanis_ea2.y*mA.x) + Sod.y * (mcanis_ea1.y*mA.z + mcanis_ea3.y*mA.x) + Sod.x * (mcanis_ea2.y*mA.z + mcanis_ea3.y*mA.y),
				Sod.z * (mcanis_ea1.z*mA.y + mcanis_ea2.z*mA.x) + Sod.y * (mcanis_ea1.z*mA.z + mcanis_ea3.z*mA.x) + Sod.x * (mcanis_ea2.z*mA.z + mcanis_ea3.z*mA.y));

			cuReal3 Hmel_1_B = (-2.0 * MEc.i / (MU0 * Ms_AFM.j)) * cuReal3(
				mB.x*Sd.x*mcanis_ea1.x + mB.y*Sd.y*mcanis_ea2.x + mB.z*Sd.z*mcanis_ea3.x,
				mB.x*Sd.x*mcanis_ea1.y + mB.y*Sd.y*mcanis_ea2.y + mB.z*Sd.z*mcanis_ea3.y,
				mB.x*Sd.x*mcanis_ea1.z + mB.y*Sd.y*mcanis_ea2.z + mB.z*Sd.z*mcanis_ea3.z);

			cuReal3 Hmel_2_B = (-2.0 * MEc.j / (MU0 * Ms_AFM.j)) * cuReal3(
				Sod.z * (mcanis_ea1.x*mB.y + mcanis_ea2.x*mB.x) + Sod.y * (mcanis_ea1.x*mB.z + mcanis_ea3.x*mB.x) + Sod.x * (mcanis_ea2.x*mB.z + mcanis_ea3.x*mB.y),
				Sod.z * (mcanis_ea1.y*mB.y + mcanis_ea2.y*mB.x) + Sod.y * (mcanis_ea1.y*mB.z + mcanis_ea3.y*mB.x) + Sod.x * (mcanis_ea2.y*mB.z + mcanis_ea3.y*mB.y),
				Sod.z * (mcanis_ea1.z*mB.y + mcanis_ea2.z*mB.x) + Sod.y * (mcanis_ea1.z*mB.z + mcanis_ea3.z*mB.x) + Sod.x * (mcanis_ea2.z*mB.z + mcanis_ea3.z*mB.y));

			Heff[idx] += Hmel_1_A + Hmel_2_A;
			Heff2[idx] += Hmel_1_B + Hmel_2_B;

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] * (Hmel_1_A + Hmel_2_A) + M2[idx] * (Hmel_1_B + Hmel_2_B)) / (2 * non_empty_cells);
			}

			if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hmel_1_A + Hmel_2_A;
			if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[idx] = Hmel_1_B + Hmel_2_B;
			if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * M[idx] * (Hmel_1_A + Hmel_2_A) / 2;
			if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[idx] = -(cuBReal)MU0 * M2[idx] * (Hmel_1_B + Hmel_2_B) / 2;
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- Calculate_MElastic_Field LAUNCHER

//compute magnetoelastic effective field to use in magnetization equation.
void MElasticCUDA::Calculate_MElastic_Field(void)
{
	//disabled by setting magnetoelastic coefficient to zero (also disabled in non-magnetic meshes)
	if (melastic_field_disabled) return;

	ZeroEnergy();

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			MElasticCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else {

			MElasticCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
		}
	}
	else if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			MElasticCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else {

			MElasticCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
		}
	}
}

//----------------------- Iterate_Elastic_Solver KERNELS

__global__ void Iterate_Elastic_Solver_Velocity_Kernel(
	ManagedMeshCUDA& cuMesh,
	MElastic_BoundaryCUDA* external_stress_surfaces, size_t num_surfaces,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuReal3>& sdd,
	cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuBReal time, cuBReal dT)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	//kernel launch with size (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1) 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int i = idx % (n_m.i + 1);
	int j = (idx / (n_m.i + 1)) % (n_m.j + 1);
	int k = idx / ((n_m.i + 1)*(n_m.j + 1));

	if (idx < (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1)) {

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

		cuBReal density = *cuMesh.pdensity;
		cuBReal mdamping = *cuMesh.pmdamping;
		cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pdensity, density, *cuMesh.pmdamping, mdamping);

		cuINT3 ijk = cuINT3(i, j, k);

		//external forces on different faces (keep track separately in case an edge cell is excited simultaneously by 2 or more external forces
		cuReal3 Fext_xface = cuReal3(), Fext_yface = cuReal3(), Fext_zface = cuReal3();
		
		//is there an external force? If so, get it, otherwise it will be zero
		if (
			((i == 0 || i == n_m.i) && strain_diag.is_dirichlet_x(idx_u)) ||
			((j == 0 || j == n_m.j) && strain_diag.is_dirichlet_y(idx_u)) ||
			((k == 0 || k == n_m.k) && strain_diag.is_dirichlet_z(idx_u))) {

			//search through all available surfaces to get external force
			for (int sidx = 0; sidx < num_surfaces; sidx++) {

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
			if (((j == 0 || j == n_m.j) && u_disp.is_dirichlet_y(idx_u)) || ((k == 0 || k == n_m.k) && u_disp.is_dirichlet_z(idx_u))) {

				vx[ijk] = 0.0;
			}
			else {

				int njend = (j < n_m.j);
				int nkend = (k < n_m.k);

				//check for required axis normal faces being present
				bool zface_u =
					j < n_m.j &&
					(u_disp.is_not_empty(idx_u) ||
					(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

				bool zface_l =
					j > 0 &&
					(u_disp.is_not_empty(idx_u - njend * n_m.x) ||
					(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

				bool yface_u =
					k < n_m.k &&
					(u_disp.is_not_empty(idx_u) ||
					(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)));

				bool yface_l =
					k > 0 &&
					(u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
					(j > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

				//at least one face is required, otherwise velocity must be zero
				if (zface_u || zface_l || yface_u || yface_l) {

					cuBReal dsxx_dx = 0.0, dsxy_dy = 0.0, dsxz_dz = 0.0;

					//always interior
					dsxx_dx = (sdd[cuINT3(i + 1, j, k)].x - sdd[ijk].x) / h_m.x;

					//interior
					if (zface_u && zface_l) dsxy_dy = (sxy[ijk] - sxy[cuINT3(i, j - 1, k)]) / h_m.y;
					else if (zface_l) dsxy_dy = (Fext_yface.x - sxy[cuINT3(i, j - 1, k)]) / (h_m.y / 2);
					else if (zface_u) dsxy_dy = (sxy[ijk] - Fext_yface.x) / (h_m.y / 2);

					//interior
					if (yface_u && yface_l) dsxz_dz = (sxz[ijk] - sxz[cuINT3(i, j, k - 1)]) / h_m.z;
					else if (yface_l) dsxz_dz = (Fext_zface.x - sxz[cuINT3(i, j, k - 1)]) / (h_m.z / 2);
					else if (yface_u) dsxz_dz = (sxz[ijk] - Fext_zface.x) / (h_m.z / 2);

					vx[ijk] += dT * (dsxx_dx + dsxy_dy + dsxz_dz - mdamping * vx[ijk]) / density;
				}
				else vx[ijk] = 0.0;
			}
		}

		//update vy
		if (j < n_m.j) {

			//set zero at fixed faces (for vy only x and z faces are applicable)
			if (((i == 0 || i == n_m.i) && u_disp.is_dirichlet_x(idx_u)) || ((k == 0 || k == n_m.k) && u_disp.is_dirichlet_z(idx_u))) {

				vy[ijk] = 0.0;
			}
			else {

				int niend = (i < n_m.i);
				int nkend = (k < n_m.k);

				//check for required axis normal faces being present
				bool zface_u =
					i < n_m.i &&
					(u_disp.is_not_empty(idx_u) ||
					(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

				bool zface_l =
					i > 0 &&
					(u_disp.is_not_empty(idx_u - niend) ||
					(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

				bool xface_u =
					k < n_m.k &&
					(u_disp.is_not_empty(idx_u) ||
					(i > 0 && u_disp.is_not_empty(idx_u - niend)));

				bool xface_l =
					k > 0 &&
					(u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
					(i > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

				//at least one face is required, otherwise velocity must be zero
				if (zface_u || zface_l || xface_u || xface_l) {

					cuBReal dsxy_dx = 0.0, dsyy_dy = 0.0, dsyz_dz = 0.0;

					//always interior
					dsyy_dy = (sdd[cuINT3(i, j + 1, k)].y - sdd[ijk].y) / h_m.y;

					//interior
					if (zface_u && zface_l) dsxy_dx = (sxy[ijk] - sxy[cuINT3(i - 1, j, k)]) / h_m.x;
					else if (zface_l) dsxy_dx = (Fext_xface.y - sxy[cuINT3(i - 1, j, k)]) / (h_m.x / 2);
					else if (zface_u) dsxy_dx = (sxy[ijk] - Fext_xface.y) / (h_m.x / 2);

					//interior
					if (xface_u && xface_l) dsyz_dz = (syz[ijk] - syz[cuINT3(i, j, k - 1)]) / h_m.z;
					else if (xface_l) dsyz_dz = (Fext_zface.y - syz[cuINT3(i, j, k - 1)]) / (h_m.z / 2);
					else if (xface_u) dsyz_dz = (syz[ijk] - Fext_zface.y) / (h_m.z / 2);

					vy[ijk] += dT * (dsxy_dx + dsyy_dy + dsyz_dz - mdamping * vy[ijk]) / density;
				}
				else vy[ijk] = 0.0;
			}
		}

		//update vz
		if (k < n_m.k) {

			//set zero at fixed faces (for vz only x and y faces are applicable)
			if (((i == 0 || i == n_m.i) && u_disp.is_dirichlet_x(idx_u)) || ((j == 0 || j == n_m.j) && u_disp.is_dirichlet_y(idx_u))) {

				vz[ijk] = 0.0;
			}
			else {

				int niend = (i < n_m.i);
				int njend = (j < n_m.j);

				//check for required axis normal faces being present
				bool yface_u =
					i < n_m.i &&
					(u_disp.is_not_empty(idx_u) ||
					(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)));

				bool yface_l =
					i > 0 &&
					(u_disp.is_not_empty(idx_u - niend) ||
					(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

				bool xface_u =
					j < n_m.j &&
					(u_disp.is_not_empty(idx_u) ||
					(i > 0 && u_disp.is_not_empty(idx_u - niend)));

				bool xface_l =
					j > 0 &&
					(u_disp.is_not_empty(idx_u - njend * n_m.x) ||
					(i > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

				//at least one face is required, otherwise velocity must be zero
				if (yface_u || yface_l || xface_u || xface_l) {

					cuBReal dsxz_dx = 0.0, dsyz_dy = 0.0, dszz_dz = 0.0;

					//always interior
					dszz_dz = (sdd[cuINT3(i, j, k + 1)].z - sdd[ijk].z) / h_m.z;

					//interior
					if (yface_u && yface_l) dsxz_dx = (sxz[ijk] - sxz[cuINT3(i - 1, j, k)]) / h_m.x;
					else if (yface_l) dsxz_dx = (Fext_xface.z - sxz[cuINT3(i - 1, j, k)]) / (h_m.x / 2);
					else if (yface_u) dsxz_dx = (sxz[ijk] - Fext_xface.z) / (h_m.x / 2);

					//interior
					if (xface_u && xface_l) dsyz_dy = (syz[ijk] - syz[cuINT3(i, j - 1, k)]) / h_m.y;
					else if (xface_l) dsyz_dy = (Fext_yface.z - syz[cuINT3(i, j - 1, k)]) / (h_m.y / 2);
					else if (xface_u) dsyz_dy = (syz[ijk] - Fext_yface.z) / (h_m.y / 2);

					vz[ijk] += dT * (dsxz_dx + dsyz_dy + dszz_dz - mdamping * vz[ijk]) / density;
				}
				else vz[ijk] = 0.0;
			}
		}

		//IMPORTANT: don't update mechanical displacement now. Do it in the stress update routine instead.
		//The reason is velocity components on CMBND will be incorrect here, but these will be set from continuity condition correctly after all meshes have updated, and before stress is calculated.
		//Thus mechnical displacement computed here will be incorrect, but when computed in stress update routine it will be correct
	}
}

__device__ void Iterate_Elastic_Solver_Stress_CUDA(
	cuINT3 ijk, cuINT3 ijk_u, int idx_u,
	ManagedMeshCUDA& cuMesh,
	MElastic_BoundaryCUDA* external_stress_surfaces, size_t num_surfaces,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuBReal time, cuBReal dT,
	bool thermoelasticity_enabled,
	cuBReal* Temp_previous, cuBReal magnetic_dT,
	cuReal3 dsdd_dt_ms, cuBReal dsxy_dt_ms, cuBReal dsxz_dt_ms, cuBReal dsyz_dt_ms)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	cuReal3 cC = *cuMesh.pcC;
	cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC);
	cuBReal cr = cC.j / cC.i;

	//needed for thermoelasticity (includes time derivative of temperature)
	cuBReal dsdd_dt_te = 0.0;
	if (thermoelasticity_enabled) {

		cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;
		cuVEC_VC<cuBReal>& Temp_l = *cuMesh.pTemp_l;

		int idx_T = Temp.position_to_cellidx(u_disp.cellidx_to_position(idx_u));

		if (Temp.is_not_empty(idx_T)) {

			cuBReal thalpha = *cuMesh.pthalpha;
			cuReal3 cC = *cuMesh.pcC;
			cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC, *cuMesh.pthalpha, thalpha);

			cuBReal Temperature = 0.0;
			//for 2TM we need to use the lattice temperature
			if (Temp_l.linear_size()) Temperature = Temp_l[idx_T];
			else Temperature = Temp[idx_T];

			dsdd_dt_te = (cC.i + 2 * cC.j) * thalpha * (Temperature - Temp_previous[idx_T]) / magnetic_dT;
		}
	}

	cuBReal adT = dsdd_dt_te / cC.i;
	cuReal3 dms = dsdd_dt_ms / cC.i;

	//external forces on different faces (keep track separately in case an edge cell is excited simultaneously by 2 or more external forces
	cuBReal Fext_xface = 0.0, Fext_yface = 0.0, Fext_zface = 0.0;
	//time derivatives of forces on the different faces, divided by c11
	cuBReal dFx = 0.0, dFy = 0.0, dFz = 0.0;
			
	//is there an external force? If so, get it, otherwise it will be zero
	if (
		((ijk.i == 0 || ijk.i == n_m.i) && strain_diag.is_dirichlet_x(idx_u)) ||
		((ijk.j == 0 || ijk.j == n_m.j) && strain_diag.is_dirichlet_y(idx_u)) ||
		((ijk.k == 0 || ijk.k == n_m.k) && strain_diag.is_dirichlet_z(idx_u))) {

		//search through all available surfaces to get external force
		for (int sidx = 0; sidx < num_surfaces; sidx++) {

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
	int niend = (ijk.i < n_m.i);
	int njend = (ijk.j < n_m.j);
	int nkend = (ijk.k < n_m.k);

	//check if required edges are present
	bool xedge_u =
		ijk.i < n_m.i &&
		(u_disp.is_not_empty(idx_u) ||
		(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
		(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
		(ijk.j > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

	bool xedge_l =
		ijk.i > 0 &&
		(u_disp.is_not_empty(idx_u - niend) ||
		(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)) ||
		(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)) ||
		(ijk.j > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x - niend)));

	bool yedge_u =
		ijk.j < n_m.j &&
		(u_disp.is_not_empty(idx_u) ||
		(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
		(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend)) ||
		(ijk.i > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

	bool yedge_l =
		ijk.j > 0 &&
		(u_disp.is_not_empty(idx_u - njend * n_m.x) ||
		(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)) ||
		(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend - njend * n_m.x)) ||
		(ijk.i > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend - njend * n_m.x)));

	bool zedge_u =
		ijk.k < n_m.k &&
		(u_disp.is_not_empty(idx_u) ||
		(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
		(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend)) ||
		(ijk.i > 0 && ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

	bool zedge_l =
		ijk.k > 0 &&
		(u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
		(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - nkend * n_m.x*n_m.y)) ||
		(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend - nkend * n_m.x*n_m.y)) ||
		(ijk.i > 0 && ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend - nkend * n_m.x*n_m.y)));

	//check for fixed faces at ends
	bool xfixed_l = (ijk.i == 0 && u_disp.is_dirichlet_px(idx_u));
	bool xfixed_u = (ijk.i == n_m.i && u_disp.is_dirichlet_nx(idx_u));

	bool yfixed_l = (ijk.j == 0 && u_disp.is_dirichlet_py(idx_u));
	bool yfixed_u = (ijk.j == n_m.j && u_disp.is_dirichlet_ny(idx_u));

	bool zfixed_l = (ijk.k == 0 && u_disp.is_dirichlet_pz(idx_u));
	bool zfixed_u = (ijk.k == n_m.k && u_disp.is_dirichlet_nz(idx_u));

	cuBReal dvx_dx = 0.0;

	//interior
	if (xedge_u && xedge_l) dvx_dx = (vx[ijk] - vx[cuINT3(ijk.i - 1, ijk.j, ijk.k)]) / h_m.x;
	//fixed face : Dirichlet value of zero for velocity derivative
	else if (xedge_l && xfixed_u) {

		dvx_dx = -vx[cuINT3(ijk.i - 1, ijk.j, ijk.k)] / (h_m.x / 2);
	}
	else if (xedge_u && xfixed_l) {

		dvx_dx = vx[ijk] / (h_m.x / 2);
	}
	//free face
	else {

		//both side derivatives
		if (yedge_l && yedge_u && zedge_l && zedge_u) {

			dvx_dx = -cr * ((vy[ijk] - vy[cuINT3(ijk.i, ijk.j - 1, ijk.k)]) / h_m.y + (vz[ijk] - vz[cuINT3(ijk.i, ijk.j, ijk.k - 1)]) / h_m.z) + adT + dms.x + dFx;
		}
		//only z derivative
		else if (zedge_l && zedge_u) {

			//dvx = (dFx - cr*dFy + dmsx - cr*dmsy) / (1 - cr^2) + (adT - cr*dvz) / (1 + cr)
			dvx_dx = (adT - cr * (vz[ijk] - vz[cuINT3(ijk.i, ijk.j, ijk.k - 1)]) / h_m.z) / (1 + cr) + (dms.x - cr * dms.y + dFx - cr * dFy) / (1 - cr * cr);
		}
		//only y derivative
		else if (yedge_l && yedge_u) {

			//dvx = (dFx - cr*dFz + dmsx - cr*dmsz) / (1 - cr^2) + (adT - cr*dvy) / (1 + cr)
			dvx_dx = (adT - cr * (vy[ijk] - vy[cuINT3(ijk.i, ijk.j - 1, ijk.k)]) / h_m.y) / (1 + cr) + (dms.x - cr * dms.z + dFx - cr * dFz) / (1 - cr * cr);
		}
		//no side derivatives : corner point. In this case all diagonal stress components set from external conditions, so derivatives not needed (set zero)
		else dvx_dx = 0.0;
	}

	cuBReal dvy_dy = 0.0;

	//interior
	if (yedge_u && yedge_l) dvy_dy = (vy[ijk] - vy[cuINT3(ijk.i, ijk.j - 1, ijk.k)]) / h_m.y;
	//fixed face : Dirichlet value of zero for velocity derivative
	else if (yedge_l && yfixed_u) {

		dvy_dy = -vy[cuINT3(ijk.i, ijk.j - 1, ijk.k)] / (h_m.y / 2);
	}
	else if (yedge_u && yfixed_l) {

		dvy_dy = vy[ijk] / (h_m.y / 2);
	}
	//free face
	else {

		//both side derivatives
		if (xedge_l && xedge_u && zedge_l && zedge_u) {

			dvy_dy = -cr * ((vx[ijk] - vx[cuINT3(ijk.i - 1, ijk.j, ijk.k)]) / h_m.x + (vz[ijk] - vz[cuINT3(ijk.i, ijk.j, ijk.k - 1)]) / h_m.z) + adT + dms.y + dFy;
		}
		//only z derivative
		else if (zedge_l && zedge_u) {

			//dvy = (dFy - cr*dFx + dmsy - cr*dmsx) / (1 - cr^2) + (adT - cr*dvz) / (1 + cr)
			dvy_dy = (adT - cr * (vz[ijk] - vz[cuINT3(ijk.i, ijk.j, ijk.k - 1)]) / h_m.z) / (1 + cr) + (dms.y - cr * dms.x + dFy - cr * dFx) / (1 - cr * cr);
		}
		//only x derivative
		else if (xedge_l && xedge_u) {

			//dvy = (dFy - cr*dFz + dmsy - cr*dmsz) / (1 - cr^2) + (adT - cr*dvx) / (1 + cr)
			dvy_dy = (adT - cr * (vx[ijk] - vx[cuINT3(ijk.i - 1, ijk.j, ijk.k)]) / h_m.x) / (1 + cr) + (dms.y - cr * dms.z + dFy - cr * dFz) / (1 - cr * cr);
		}
		//no side derivatives : corner point. In this case all diagonal stress components set from external conditions, so derivatives not needed (set zero)
		else dvy_dy = 0.0;
	}

	cuBReal dvz_dz = 0.0;

	//interior
	if (zedge_u && zedge_l) dvz_dz = (vz[ijk] - vz[cuINT3(ijk.i, ijk.j, ijk.k - 1)]) / h_m.z;
	//fixed face : Dirichlet value of zero for velocity derivative
	else if (zedge_l && zfixed_u) {

		dvz_dz = -vz[cuINT3(ijk.i, ijk.j, ijk.k - 1)] / (h_m.z / 2);
	}
	//fixed face : Dirichlet value of zero for velocity derivative
	else if (zedge_u && zfixed_l) {

		dvz_dz = vz[ijk] / (h_m.z / 2);
	}
	//free face
	else {

		//both side derivatives
		if (xedge_l && xedge_u && yedge_l && yedge_u) {

			dvz_dz = -cr * ((vx[ijk] - vx[cuINT3(ijk.i - 1, ijk.j, ijk.k)]) / h_m.x + (vy[ijk] - vy[cuINT3(ijk.i, ijk.j - 1, ijk.k)]) / h_m.y) + adT + dms.z + dFz;
		}
		//only y derivative
		else if (yedge_l && yedge_u) {

			//dvz = (dFz - cr*dFx + dmsz - cr*dmsx) / (1 - cr^2) + (adT - cr*dvy) / (1 + cr)
			dvz_dz = (adT - cr * (vy[ijk] - vy[cuINT3(ijk.i, ijk.j - 1, ijk.k)]) / h_m.y) / (1 + cr) + (dms.z - cr * dms.x + dFz - cr * dFx) / (1 - cr * cr);
		}
		//only x derivative
		else if (xedge_l && xedge_u) {

			//dvz = (dFz - cr*dFy + dmsz - cr*dmsy) / (1 - cr^2) + (adT - cr*dvx) / (1 + cr)
			dvz_dz = (adT - cr * (vx[ijk] - vx[cuINT3(ijk.i - 1, ijk.j, ijk.k)]) / h_m.x) / (1 + cr) + (dms.z - cr * dms.y + dFz - cr * dFy) / (1 - cr * cr);
		}
		//no side derivatives : corner point. In this case all diagonal stress components set from external conditions, so derivatives not needed (set zero)
		else dvz_dz = 0.0;
	}

	//update sdd if not empty
	if ((xedge_u || xedge_l) && (yedge_u || yedge_l) && (zedge_u || zedge_l)) {

		if ((!xedge_u && !xfixed_u) || (!xedge_l && !xfixed_l)) sdd[ijk].x = Fext_xface;
		else sdd[ijk].x += dT * (cC.i * dvx_dx + cC.j * (dvy_dy + dvz_dz) - dsdd_dt_ms.x - dsdd_dt_te);

		if ((!yedge_u && !yfixed_u) || (!yedge_l && !yfixed_l)) sdd[ijk].y = Fext_yface;
		else sdd[ijk].y += dT * (cC.i * dvy_dy + cC.j * (dvx_dx + dvz_dz) - dsdd_dt_ms.y - dsdd_dt_te);

		if ((!zedge_u && !zfixed_u) || (!zedge_l && !zfixed_l)) sdd[ijk].z = Fext_zface;
		else sdd[ijk].z += dT * (cC.i * dvz_dz + cC.j * (dvx_dx + dvy_dy) - dsdd_dt_ms.z - dsdd_dt_te);
	}
	else sdd[ijk] = cuReal3();

	//update sxy
	if (ijk.i < n_m.i && ijk.j < n_m.j) {

		bool zface =
			(u_disp.is_not_empty(idx_u) ||
			(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

		if (zface) {

			cuBReal dvx_dy = (vx[cuINT3(ijk.i, ijk.j + 1, ijk.k)] - vx[ijk]) / h_m.y;
			cuBReal dvy_dx = (vy[cuINT3(ijk.i + 1, ijk.j, ijk.k)] - vy[ijk]) / h_m.x;

			sxy[ijk] += dT * (cC.k * (dvx_dy + dvy_dx) / 2 - dsxy_dt_ms);
		}
		else sxy[ijk] = 0.0;
	}

	//update sxz
	if (ijk.i < n_m.i && ijk.k < n_m.k) {

		bool yface =
			(u_disp.is_not_empty(idx_u) ||
			(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)));

		if (yface) {

			cuBReal dvx_dz = (vx[cuINT3(ijk.i, ijk.j, ijk.k + 1)] - vx[ijk]) / h_m.z;
			cuBReal dvz_dx = (vz[cuINT3(ijk.i + 1, ijk.j, ijk.k)] - vz[ijk]) / h_m.x;

			sxz[ijk] += dT * (cC.k * (dvx_dz + dvz_dx) / 2 - dsxz_dt_ms);
		}
		else sxz[ijk] = 0.0;
	}

	//update syz
	if (ijk.j < n_m.j && ijk.k < n_m.k) {

		bool xface =
			(u_disp.is_not_empty(idx_u) ||
			(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend)));

		if (xface) {

			cuBReal dvy_dz = (vy[cuINT3(ijk.i, ijk.j, ijk.k + 1)] - vy[ijk]) / h_m.z;
			cuBReal dvz_dy = (vz[cuINT3(ijk.i, ijk.j + 1, ijk.k)] - vz[ijk]) / h_m.y;

			syz[ijk] += dT * (cC.k * (dvy_dz + dvz_dy) / 2 - dsyz_dt_ms);
		}
		else syz[ijk] = 0.0;
	}

	//update mechanical displacement using velocity (remember u is cell-centred)
	if (ijk.i < n_m.i && ijk.j < n_m.j && ijk.k < n_m.k) {

		if (u_disp.is_not_empty(idx_u)) {

			//find velocity values cell-centred
			cuBReal vx_cc = (vx[ijk] + vx[ijk + cuINT3(0, 1, 0)] + vx[ijk + cuINT3(0, 0, 1)] + vx[ijk + cuINT3(0, 1, 1)]) / 4;
			cuBReal vy_cc = (vy[ijk] + vy[ijk + cuINT3(1, 0, 0)] + vy[ijk + cuINT3(0, 0, 1)] + vy[ijk + cuINT3(1, 0, 1)]) / 4;
			cuBReal vz_cc = (vz[ijk] + vz[ijk + cuINT3(1, 0, 0)] + vz[ijk + cuINT3(0, 1, 0)] + vz[ijk + cuINT3(1, 1, 0)]) / 4;

			u_disp[idx_u] += dT * cuReal3(vx_cc, vy_cc, vz_cc);
		}
		else u_disp[idx_u] = cuReal3();
	}
}

__global__ void Iterate_Elastic_Solver_Stress_FM_Kernel(
	ManagedMeshCUDA& cuMesh,
	MElastic_BoundaryCUDA* external_stress_surfaces, size_t num_surfaces,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuBReal time, cuBReal dT,
	bool magnetostriction_enabled, bool thermoelasticity_enabled,
	cuBReal* Temp_previous, cuBReal magnetic_dT,
	ManagedDiffEqFMCUDA& cuDiffEq_FM)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;

	cuSZ3& n_m = u_disp.n;

	//kernel launch with size (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1) 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int i = idx % (n_m.i + 1);
	int j = (idx / (n_m.i + 1)) % (n_m.j + 1);
	int k = idx / ((n_m.i + 1)*(n_m.j + 1));

	if (idx < (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1)) {

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

		//needed for magnetostriction (time derivatives of stress due to magnetostriction)
		cuReal3 dsdd_dt_ms = cuReal3();
		cuBReal dsxy_dt_ms = 0.0, dsxz_dt_ms = 0.0, dsyz_dt_ms = 0.0;
		if (magnetostriction_enabled) {

			cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
			cuVEC_VC<cuReal3>& M = *cuMesh.pM;

			int idx_M = M.position_to_cellidx(u_disp.cellidx_to_position(idx_u));

			if (M.is_not_empty(idx_M)) {

				cuBReal Ms = *cuMesh.pMs;
				cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
				cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
				cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
				cuReal2 mMEc = *cuMesh.pmMEc;
				cuMesh.update_parameters_mcoarse(idx_M, *cuMesh.pMs, Ms, *cuMesh.pmMEc, mMEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

				cuReal3 m = cuReal3(M[idx_M] * mcanis_ea1, M[idx_M] * mcanis_ea2, M[idx_M] * mcanis_ea3) / Ms;
				cuReal3 dM_dt = (M[idx_M] - (*cuDiffEq_FM.psM1)[idx_M]) / magnetic_dT;
				cuReal3 dm_dt = cuReal3(dM_dt * mcanis_ea1, dM_dt * mcanis_ea2, dM_dt * mcanis_ea3) / Ms;

				dsdd_dt_ms = 2 * mMEc.i * cuReal3(
					m.x*dm_dt.x*mcanis_ea1.x + m.y*dm_dt.y*mcanis_ea2.x + m.z*dm_dt.z*mcanis_ea3.x,
					m.x*dm_dt.x*mcanis_ea1.y + m.y*dm_dt.y*mcanis_ea2.y + m.z*dm_dt.z*mcanis_ea3.y,
					m.x*dm_dt.x*mcanis_ea1.z + m.y*dm_dt.y*mcanis_ea2.z + m.z*dm_dt.z*mcanis_ea3.z);

				dsxy_dt_ms = 2 * mMEc.j * ((m.x*dm_dt.y + m.y*dm_dt.x)*mcanis_ea3.z + (m.x*dm_dt.z + m.z*dm_dt.x)*mcanis_ea2.z + (m.y*dm_dt.z + m.z*dm_dt.y)*mcanis_ea1.z);
				dsxz_dt_ms = 2 * mMEc.j * ((m.x*dm_dt.y + m.y*dm_dt.x)*mcanis_ea3.y + (m.x*dm_dt.z + m.z*dm_dt.x)*mcanis_ea2.y + (m.y*dm_dt.z + m.z*dm_dt.y)*mcanis_ea1.y);
				dsyz_dt_ms = 2 * mMEc.j * ((m.x*dm_dt.y + m.y*dm_dt.x)*mcanis_ea3.x + (m.x*dm_dt.z + m.z*dm_dt.x)*mcanis_ea2.x + (m.y*dm_dt.z + m.z*dm_dt.y)*mcanis_ea1.x);
			}
		}

		//now solve the main part, with the possible addition of magnetostriction contribution
		Iterate_Elastic_Solver_Stress_CUDA(
			ijk, ijk_u, idx_u,
			cuMesh,
			external_stress_surfaces, num_surfaces,
			vx, vy, vz,
			sdd, sxy, sxz, syz,
			time, dT,
			thermoelasticity_enabled,
			Temp_previous, magnetic_dT,
			dsdd_dt_ms, dsxy_dt_ms, dsxz_dt_ms, dsyz_dt_ms);
	}
}

__global__ void Iterate_Elastic_Solver_Stress_AFM_Kernel(
	ManagedMeshCUDA& cuMesh,
	MElastic_BoundaryCUDA* external_stress_surfaces, size_t num_surfaces,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuBReal time, cuBReal dT,
	bool magnetostriction_enabled, bool thermoelasticity_enabled,
	cuBReal* Temp_previous, cuBReal magnetic_dT,
	ManagedDiffEqAFMCUDA& cuDiffEq_AFM)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;

	cuSZ3& n_m = u_disp.n;

	//kernel launch with size (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1) 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int i = idx % (n_m.i + 1);
	int j = (idx / (n_m.i + 1)) % (n_m.j + 1);
	int k = idx / ((n_m.i + 1)*(n_m.j + 1));

	if (idx < (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1)) {

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

		//needed for magnetostriction (time derivatives of stress due to magnetostriction)
		cuReal3 dsdd_dt_ms = cuReal3();
		cuBReal dsxy_dt_ms = 0.0, dsxz_dt_ms = 0.0, dsyz_dt_ms = 0.0;
		if (magnetostriction_enabled) {

			cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
			cuVEC_VC<cuReal3>& M = *cuMesh.pM;
			cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

			int idx_M = M.position_to_cellidx(u_disp.cellidx_to_position(idx_u));

			if (M.is_not_empty(idx_M)) {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
				cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
				cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
				cuReal2 mMEc = *cuMesh.pmMEc;
				cuMesh.update_parameters_mcoarse(idx_M, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pmMEc, mMEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

				cuReal3 mA = cuReal3(M[idx_M] * mcanis_ea1, M[idx_M] * mcanis_ea2, M[idx_M] * mcanis_ea3) / Ms_AFM.i;
				cuReal3 mB = cuReal3(M2[idx_M] * mcanis_ea1, M2[idx_M] * mcanis_ea2, M2[idx_M] * mcanis_ea3) / Ms_AFM.j;
				cuReal3 dM_dtA = (M[idx_M] - (*cuDiffEq_AFM.psM1)[idx_M]) / magnetic_dT;
				cuReal3 dm_dtA = cuReal3(dM_dtA * mcanis_ea1, dM_dtA * mcanis_ea2, dM_dtA * mcanis_ea3) / Ms_AFM.i;

				cuReal3 dM_dtB = (M2[idx_M] - (*cuDiffEq_AFM.psM1_2)[idx_M]) / magnetic_dT;
				cuReal3 dm_dtB = cuReal3(dM_dtB * mcanis_ea1, dM_dtB * mcanis_ea2, dM_dtB * mcanis_ea3) / Ms_AFM.j;

				dsdd_dt_ms = mMEc.i * cuReal3(
					(mA.x*dm_dtA.x + mB.x*dm_dtB.x)*mcanis_ea1.x + (mA.y*dm_dtA.y + mB.y*dm_dtB.y)*mcanis_ea2.x + (mA.z*dm_dtA.z + mB.z*dm_dtB.z)*mcanis_ea3.x,
					(mA.x*dm_dtA.x + mB.x*dm_dtB.x)*mcanis_ea1.y + (mA.y*dm_dtA.y + mB.y*dm_dtB.y)*mcanis_ea2.y + (mA.z*dm_dtA.z + mB.z*dm_dtB.z)*mcanis_ea3.y,
					(mA.x*dm_dtA.x + mB.x*dm_dtB.x)*mcanis_ea1.z + (mA.y*dm_dtA.y + mB.y*dm_dtB.y)*mcanis_ea2.z + (mA.z*dm_dtA.z + mB.z*dm_dtB.z)*mcanis_ea3.z);

				dsxy_dt_ms = mMEc.j * ((mA.x*dm_dtA.y + mA.y*dm_dtA.x + mB.x*dm_dtB.y + mB.y*dm_dtB.x)*mcanis_ea3.z + (mA.x*dm_dtA.z + mA.z*dm_dtA.x + mB.x*dm_dtB.z + mB.z*dm_dtB.x)*mcanis_ea2.z + (mA.y*dm_dtA.z + mA.z*dm_dtA.y + mB.y*dm_dtB.z + mB.z*dm_dtB.y)*mcanis_ea1.z);
				dsxz_dt_ms = mMEc.j * ((mA.x*dm_dtA.y + mA.y*dm_dtA.x + mB.x*dm_dtB.y + mB.y*dm_dtB.x)*mcanis_ea3.y + (mA.x*dm_dtA.z + mA.z*dm_dtA.x + mB.x*dm_dtB.z + mB.z*dm_dtB.x)*mcanis_ea2.y + (mA.y*dm_dtA.z + mA.z*dm_dtA.y + mB.y*dm_dtB.z + mB.z*dm_dtB.y)*mcanis_ea1.y);
				dsyz_dt_ms = mMEc.j * ((mA.x*dm_dtA.y + mA.y*dm_dtA.x + mB.x*dm_dtB.y + mB.y*dm_dtB.x)*mcanis_ea3.x + (mA.x*dm_dtA.z + mA.z*dm_dtA.x + mB.x*dm_dtB.z + mB.z*dm_dtB.x)*mcanis_ea2.x + (mA.y*dm_dtA.z + mA.z*dm_dtA.y + mB.y*dm_dtB.z + mB.z*dm_dtB.y)*mcanis_ea1.x);
			}
		}

		//now solve the main part, with the possible addition of magnetostriction contribution
		Iterate_Elastic_Solver_Stress_CUDA(
			ijk, ijk_u, idx_u,
			cuMesh,
			external_stress_surfaces, num_surfaces,
			vx, vy, vz,
			sdd, sxy, sxz, syz,
			time, dT,
			thermoelasticity_enabled,
			Temp_previous, magnetic_dT,
			dsdd_dt_ms, dsxy_dt_ms, dsxz_dt_ms, dsyz_dt_ms);
	}
}

__global__ void Iterate_Elastic_Solver_Stress_NoMS_Kernel(
	ManagedMeshCUDA& cuMesh,
	MElastic_BoundaryCUDA* external_stress_surfaces, size_t num_surfaces,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuBReal time, cuBReal dT,
	bool thermoelasticity_enabled,
	cuBReal* Temp_previous, cuBReal magnetic_dT)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;

	cuSZ3& n_m = u_disp.n;

	//kernel launch with size (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1) 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int i = idx % (n_m.i + 1);
	int j = (idx / (n_m.i + 1)) % (n_m.j + 1);
	int k = idx / ((n_m.i + 1)*(n_m.j + 1));

	if (idx < (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1)) {

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

		//now solve the main part without magnetostriction
		Iterate_Elastic_Solver_Stress_CUDA(
			ijk, ijk_u, idx_u,
			cuMesh,
			external_stress_surfaces, num_surfaces,
			vx, vy, vz,
			sdd, sxy, sxz, syz,
			time, dT,
			thermoelasticity_enabled,
			Temp_previous, magnetic_dT,
			cuReal3(), 0.0, 0.0, 0.0);
	}
}

__global__ void Calculate_Strain_Kernel(
	ManagedMeshCUDA& cuMesh,
	cuVEC<cuReal3>& sdd,
	cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	//kernel launch with size n_m.i * n_m.j * n_m.k 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < u_disp.linear_size()) {
		
		if (u_disp.is_not_empty(idx)) {

				//get all 9 first-order differentials of u
				cuReal33 grad_u = u_disp.grad_sided(idx);

				//diagonal components
				strain_diag[idx] = cuReal3(grad_u.x.x, grad_u.y.y, grad_u.z.z);

				//off-diagonal components (yz, xz, xy)
				strain_odiag[idx] = 0.5 * cuReal3(grad_u.y.z + grad_u.z.y, grad_u.x.z + grad_u.z.x, grad_u.x.y + grad_u.y.x);
		}
		else {

			strain_diag[idx] = cuReal3();
			strain_odiag[idx] = cuReal3();
		}
	}
}

//----------------------- Iterate_Elastic_Solver LAUNCHERS

//update velocity for dT time increment (also updating displacement)
void MElasticCUDA::Iterate_Elastic_Solver_Velocity(double dT)
{
	size_t size = (pMeshCUDA->n_m.i + 1) * (pMeshCUDA->n_m.j + 1) * (pMeshCUDA->n_m.k + 1);

	//1a. Update velocity
	Iterate_Elastic_Solver_Velocity_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(pMeshCUDA->cuMesh, external_stress_surfaces_arr, external_stress_surfaces.size(),
		vx, vy, vz, sdd, sxy, sxz, syz,
		pMeshCUDA->GetStageTime(), dT);
}

//update stress for dT time increment
void MElasticCUDA::Iterate_Elastic_Solver_Stress(double dT, double magnetic_dT)
{
	size_t size = (pMeshCUDA->n_m.i + 1) * (pMeshCUDA->n_m.j + 1) * (pMeshCUDA->n_m.k + 1);

	//1b. Update stress
	if (magnetostriction_enabled) {

		if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			Iterate_Elastic_Solver_Stress_AFM_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, external_stress_surfaces_arr, external_stress_surfaces.size(),
					vx, vy, vz, sdd, sxy, sxz, syz,
					pMeshCUDA->GetStageTime(), dT,
					magnetostriction_enabled, thermoelasticity_enabled,
					Temp_previous, magnetic_dT,
					reinterpret_cast<AFMeshCUDA*>(pMeshCUDA)->Get_ManagedDiffEqCUDA());
		}
		else if (pMeshCUDA->GetMeshType() == MESH_FERROMAGNETIC) {

			Iterate_Elastic_Solver_Stress_FM_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, external_stress_surfaces_arr, external_stress_surfaces.size(),
					vx, vy, vz, sdd, sxy, sxz, syz,
					pMeshCUDA->GetStageTime(), dT,
					magnetostriction_enabled, thermoelasticity_enabled,
					Temp_previous, magnetic_dT,
					reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->Get_ManagedDiffEqCUDA());
		}
	}
	else {

		Iterate_Elastic_Solver_Stress_NoMS_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, external_stress_surfaces_arr, external_stress_surfaces.size(),
				vx, vy, vz, sdd, sxy, sxz, syz,
				pMeshCUDA->GetStageTime(), dT,
				thermoelasticity_enabled,
				Temp_previous, magnetic_dT);
	}
}

void MElasticCUDA::Calculate_Strain(void)
{
	Calculate_Strain_Kernel <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(pMeshCUDA->cuMesh, sdd, sxy, sxz, syz);
}

//---------------------------------------------- Initial Conditions Launchers and Kernels

__global__ void Set_Initial_Stress_Kernel(
	ManagedMeshCUDA& cuMesh,
	cuVEC<cuReal3>& sdd,
	cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	bool magnetostriction_enabled, bool thermoelasticity_enabled, cuBReal& T_ambient)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& strain_diag = *cuMesh.pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *cuMesh.pstrain_odiag;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	//kernel launch with size (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1) 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int i = idx % (n_m.i + 1);
	int j = (idx / (n_m.i + 1)) % (n_m.j + 1);
	int k = idx / ((n_m.i + 1)*(n_m.j + 1));

	if (idx < (n_m.i + 1) * (n_m.j + 1) * (n_m.k + 1)) {

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < n_m.i ? i : n_m.i - 1, j < n_m.j ? j : n_m.j - 1, k < n_m.k ? k : n_m.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * n_m.x + ijk_u.k * n_m.x * n_m.y;

		cuINT3 ijk = cuINT3(i, j, k);

		//update sxx, syy, szz
		int niend = (ijk.i < n_m.i);
		int njend = (ijk.j < n_m.j);
		int nkend = (ijk.k < n_m.k);

		//check if required edges are present
		bool xedge_u =
			ijk.i < n_m.i &&
			(u_disp.is_not_empty(idx_u) ||
			(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
				(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
				(ijk.j > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

		bool xedge_l =
			ijk.i > 0 &&
			(u_disp.is_not_empty(idx_u - niend) ||
			(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)) ||
				(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)) ||
				(ijk.j > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x - niend)));

		bool yedge_u =
			ijk.j < n_m.j &&
			(u_disp.is_not_empty(idx_u) ||
			(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
				(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend)) ||
				(ijk.i > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

		bool yedge_l =
			ijk.j > 0 &&
			(u_disp.is_not_empty(idx_u - njend * n_m.x) ||
			(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)) ||
				(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend - njend * n_m.x)) ||
				(ijk.i > 0 && ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend - njend * n_m.x)));

		bool zedge_u =
			ijk.k < n_m.k &&
			(u_disp.is_not_empty(idx_u) ||
			(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
				(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend)) ||
				(ijk.i > 0 && ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

		bool zedge_l =
			ijk.k > 0 &&
			(u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
			(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - nkend * n_m.x*n_m.y)) ||
				(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend - nkend * n_m.x*n_m.y)) ||
				(ijk.i > 0 && ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend - nkend * n_m.x*n_m.y)));

		cuReal3 Stress_MS_dd = cuReal3();
		cuBReal Stress_MS_xy = 0.0, Stress_MS_xz = 0.0, Stress_MS_yz = 0.0;
		if (magnetostriction_enabled) {

			cuVEC_VC<cuReal3>& M = *cuMesh.pM;
			cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

			int idx_M = M.position_to_cellidx(u_disp.cellidx_to_position(idx_u));

			if (M.is_not_empty(idx_M)) {

				//MESH_ANTIFERROMAGNETIC
				if (M2.linear_size()) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
					cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
					cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
					cuReal2 mMEc = *cuMesh.pmMEc;
					cuMesh.update_parameters_mcoarse(idx_M, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pmMEc, mMEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

					cuReal3 mA = cuReal3(M[idx_M] * mcanis_ea1, M[idx_M] * mcanis_ea2, M[idx_M] * mcanis_ea3) / Ms_AFM.i;
					cuReal3 mB = cuReal3(M2[idx_M] * mcanis_ea1, M2[idx_M] * mcanis_ea2, M2[idx_M] * mcanis_ea3) / Ms_AFM.j;

					Stress_MS_dd = (mMEc.i / 2) * cuReal3(
						(mA.x*mA.x + mB.x*mB.x)*mcanis_ea1.x + (mA.y*mA.y + mB.y*mB.y)*mcanis_ea2.x + (mA.z*mA.z + mB.z*mB.z)*mcanis_ea3.x,
						(mA.x*mA.x + mB.x*mB.x)*mcanis_ea1.y + (mA.y*mA.y + mB.y*mB.y)*mcanis_ea2.y + (mA.z*mA.z + mB.z*mB.z)*mcanis_ea3.y,
						(mA.x*mA.x + mB.x*mB.x)*mcanis_ea1.z + (mA.y*mA.y + mB.y*mB.y)*mcanis_ea2.z + (mA.z*mA.z + mB.z*mB.z)*mcanis_ea3.z);

					Stress_MS_xy = mMEc.j * ((mA.x*mA.y + mB.x*mB.y)*mcanis_ea3.z + (mA.x*mA.z + mB.x*mB.z)*mcanis_ea2.z + (mA.y*mA.z + mB.y*mB.z)*mcanis_ea1.z);
					Stress_MS_xz = mMEc.j * ((mA.x*mA.y + mB.x*mB.y)*mcanis_ea3.y + (mA.x*mA.z + mB.x*mB.z)*mcanis_ea2.y + (mA.y*mA.z + mB.y*mB.z)*mcanis_ea1.y);
					Stress_MS_yz = mMEc.j * ((mA.x*mA.y + mB.x*mB.y)*mcanis_ea3.x + (mA.x*mA.z + mB.x*mB.z)*mcanis_ea2.x + (mA.y*mA.z + mB.y*mB.z)*mcanis_ea1.x);
				}
				//MESH_FERROMAGNETIC
				else {

					cuBReal Ms = *cuMesh.pMs;
					cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
					cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
					cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
					cuReal2 mMEc = *cuMesh.pmMEc;
					cuMesh.update_parameters_mcoarse(idx_M, *cuMesh.pMs, Ms, *cuMesh.pmMEc, mMEc, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

					cuReal3 m = cuReal3(M[idx_M] * mcanis_ea1, M[idx_M] * mcanis_ea2, M[idx_M] * mcanis_ea3) / Ms;

					Stress_MS_dd = mMEc.i * cuReal3(
						m.x*m.x*mcanis_ea1.x + m.y*m.y*mcanis_ea2.x + m.z*m.z*mcanis_ea3.x,
						m.x*m.x*mcanis_ea1.y + m.y*m.y*mcanis_ea2.y + m.z*m.z*mcanis_ea3.y,
						m.x*m.x*mcanis_ea1.z + m.y*m.y*mcanis_ea2.z + m.z*m.z*mcanis_ea3.z);

					Stress_MS_xy = 2 * mMEc.j * (m.x*m.y*mcanis_ea3.z + m.x*m.z*mcanis_ea2.z + m.y*m.z*mcanis_ea1.z);
					Stress_MS_xz = 2 * mMEc.j * (m.x*m.y*mcanis_ea3.y + m.x*m.z*mcanis_ea2.y + m.y*m.z*mcanis_ea1.y);
					Stress_MS_yz = 2 * mMEc.j * (m.x*m.y*mcanis_ea3.x + m.x*m.z*mcanis_ea2.x + m.y*m.z*mcanis_ea1.x);
				}
			}
		}

		cuBReal Stress_Temp = 0.0;
		if (thermoelasticity_enabled) {

			cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;
			cuVEC_VC<cuBReal>& Temp_l = *cuMesh.pTemp_l;

			int idx_T = Temp.position_to_cellidx(u_disp.cellidx_to_position(idx_u));

			if (Temp.is_not_empty(idx_T)) {

				cuBReal thalpha = *cuMesh.pthalpha;
				cuReal3 cC = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC, *cuMesh.pthalpha, thalpha);

				cuBReal Temperature = 0.0;
				//for 2TM we need to use the lattice temperature
				if (Temp_l.linear_size()) Temperature = Temp_l[idx_T];
				else Temperature = Temp[idx_T];

				Stress_Temp = (cC.i + 2 * cC.j) * thalpha * (Temperature - T_ambient);
			}
		}

		//update sdd if not empty
		if ((xedge_u || xedge_l) && (yedge_u || yedge_l) && (zedge_u || zedge_l)) {

			sdd[ijk].x = -(Stress_MS_dd.x + Stress_Temp);
			sdd[ijk].y = -(Stress_MS_dd.y + Stress_Temp);
			sdd[ijk].z = -(Stress_MS_dd.z + Stress_Temp);
		}
		else sdd[ijk] = cuReal3();

		//update sxy
		if (ijk.i < n_m.i && ijk.j < n_m.j) {

			bool zface =
				(u_disp.is_not_empty(idx_u) ||
				(ijk.k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

			if (zface) {

				sxy[ijk] = -Stress_MS_xy;
			}
			else sxy[ijk] = 0.0;
		}

		//update sxz
		if (ijk.i < n_m.i && ijk.k < n_m.k) {

			bool yface =
				(u_disp.is_not_empty(idx_u) ||
				(ijk.j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)));

			if (yface) {

				sxz[ijk] = -Stress_MS_xz;
			}
			else sxz[ijk] = 0.0;
		}

		//update syz
		if (ijk.j < n_m.j && ijk.k < n_m.k) {

			bool xface =
				(u_disp.is_not_empty(idx_u) ||
				(ijk.i > 0 && u_disp.is_not_empty(idx_u - niend)));

			if (xface) {

				syz[ijk] = -Stress_MS_yz;
			}
			else syz[ijk] = 0.0;
		}
	}
}

//if thermoelasticity or magnetostriction is enabled, then initial stress must be set correctly
void MElasticCUDA::Set_Initial_Stress(void)
{
	if (!magnetostriction_enabled && !thermoelasticity_enabled) {

		sdd()->set(cuReal3());
		sxy()->set(0.0); sxz()->set(0.0); syz()->set(0.0);
	}
	else {

		//reset for dT / dt computation
		if (thermoelasticity_enabled) {

			if (Temp_previous.resize(pMeshCUDA->n_t.dim())) Save_Current_Temperature();
		}

		//reset for dm / dt computation
		if (magnetostriction_enabled) pMeshCUDA->SaveMagnetization();

		size_t size = (pMeshCUDA->n_m.i + 1) * (pMeshCUDA->n_m.j + 1) * (pMeshCUDA->n_m.k + 1);

		Set_Initial_Stress_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, sdd, sxy, sxz, syz, magnetostriction_enabled, thermoelasticity_enabled, T_ambient);
	}
}

__global__ void Save_Current_Temperature_Kernel(cuBReal* Temp_previous, cuVEC_VC<cuBReal>& Temp, cuVEC_VC<cuBReal>& Temp_l)
{
	//kernel launch with size n_m.i * n_m.j * n_m.k 
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		//2TM
		if (Temp_l.linear_size()) Temp_previous[idx] = Temp_l[idx];
		//1TM
		else Temp_previous[idx] = Temp[idx];
	}
}

//if thermoelasticity is enabled then save current temperature values in Temp_previous (called after elastic solver fully incremented by magnetic_dT)
void MElasticCUDA::Save_Current_Temperature(void)
{
	if (thermoelasticity_enabled) {

		Save_Current_Temperature_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (Temp_previous, pMeshCUDA->Temp, pMeshCUDA->Temp_l);
	}
}

//---------------------------------------------- CMBND routines (verlocity) KERNELS

__global__ void make_velocity_continuous_x_Kernel(
	ManagedMeshCUDA& cuMesh, ManagedMeshCUDA& cuMesh_sec,
	CMBNDInfoCUDA& contact,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuBReal>& vx_sec, cuVEC<cuBReal>& vy_sec, cuVEC<cuBReal>& vz_sec)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& u_disp_sec = *cuMesh_sec.pu_disp;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	const cuBox& cb = contact.cells_box;

	//kernel launch with size (cb.e.k + 1 - cb.s.k) * (cb.e.j + 1 - cb.s.j)
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int ydim = cb.e.j + 1 - cb.s.j;
	int zdim = cb.e.k + 1 - cb.s.k;

	if (idx < ydim * zdim) {

		int i = (contact.IsPrimaryTop() ? cb.s.i : cb.e.i);
		int j = (idx % ydim) + cb.s.j;
		int k = (idx / ydim) + cb.s.k;

		cuBReal spacing = (h_m.x + u_disp_sec.h.x);
		cuBReal wpri = 1.0 - h_m.x / spacing;
		cuBReal wsec = 1.0 - u_disp_sec.h.x / spacing;

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < u_disp.n.i ? i : u_disp.n.i - 1, j < u_disp.n.j ? j : u_disp.n.j - 1, k < u_disp.n.k ? k : u_disp.n.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * u_disp.n.x + ijk_u.k * u_disp.n.x * u_disp.n.y;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		int niend = (i < n_m.i);
		int njend = (j < n_m.j);
		int nkend = (k < n_m.k);

		//check for required faces being present (used for vx and vy components)
		bool xface_u = u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u);
		bool xface_l_vz = j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x) && u_disp.is_cmbnd(idx_u - njend * n_m.x);
		bool xface_l_vy = k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y);

		//check for required edge being present (used for vy component)
		bool xedge_u = xface_u || xface_l_vy || xface_l_vz || (j > 0 && k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y - njend * n_m.x) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y - njend * n_m.x));

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(-u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);
			int idx_u_sec = ijk_sec.i + ijk_sec.j * u_disp_sec.n.x + ijk_sec.k * u_disp_sec.n.x * u_disp_sec.n.y;

			//vy
			if (j < cb.e.j) {

				//value at interface : interpolate between primary and secondary
				if (xface_u || xface_l_vy) vy[cuINT3(0, j, k)] = vy[cuINT3(1, j, k)] * wpri + vy_sec[ijk_sec + cuINT3(0, 0, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				//value at interface : interpolate between primary and secondary
				if (xface_u || xface_l_vz) vz[cuINT3(0, j, k)] = vz[cuINT3(1, j, k)] * wpri + vz_sec[ijk_sec + cuINT3(0, j == cb.e.j, 0)] * wsec;
			}

			//set vx continuity obtained from sxx continuity
			if (xedge_u) {

				cuReal3 cC_p = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC_p);
				cuReal3 cC_s = *cuMesh_sec.pcC;
				cuMesh_sec.update_parameters_scoarse(idx_u_sec, *cuMesh_sec.pcC, cC_s);
				cuBReal rdu = (cC_p.i / cC_s.i) * (u_disp_sec.h.x / h_m.x);

				//simplest case just interpolate values either side (see derivation in notes)
				//in theory there are further contributions here due to 1) side derivatives, 2) thermoelastic contribution, 3) magnetostriction contribution
				//these contributions will be zero if there's no change in cC, alphaT, B1 coefficients across the interface since they are proportional to respective coefficient differences (see notes)
				//however even if there's a material mismatch at the interface, these contributions are still virtually zero since scaled by (cellsize / c11). complicated to include them and no real extra accuracy.
				vx[cuINT3(0, j, k)] = (vx[cuINT3(1, j, k)] * (1 + 3 * rdu) + 2 * vx_sec[ijk_sec + cuINT3(-1, 0, 0)]) / (3 * (1 + rdu));
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);
			int idx_u_sec = ijk_sec.i + ijk_sec.j * u_disp_sec.n.x + ijk_sec.k * u_disp_sec.n.x * u_disp_sec.n.y;

			//vy
			if (j < cb.e.j) {

				if (xface_u || xface_l_vy) vy[cuINT3(vy.n.x - 1, j, k)] = vy[cuINT3(vy.n.x - 2, j, k)] * wpri + vy_sec[ijk_sec + cuINT3(1, 0, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				if (xface_u || xface_l_vz) vz[cuINT3(vz.n.x - 1, j, k)] = vz[cuINT3(vz.n.x - 2, j, k)] * wpri + vz_sec[ijk_sec + cuINT3(1, j == cb.e.j, 0)] * wsec;
			}

			//set vx continuity obtained from sxx continuity
			if (xedge_u) {

				cuReal3 cC_p = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC_p);
				cuReal3 cC_s = *cuMesh_sec.pcC;
				cuMesh_sec.update_parameters_scoarse(idx_u_sec, *cuMesh_sec.pcC, cC_s);
				cuBReal rdu = (cC_p.i / cC_s.i) * (u_disp_sec.h.x / h_m.x);

				//simplest case just interpolate values either side (see derivation in notes)
				//in theory there are further contributions here due to 1) side derivatives, 2) thermoelastic contribution, 3) magnetostriction contribution
				//these contributions will be zero if there's no change in cC, alphaT, B1 coefficients across the interface since they are proportional to respective coefficient differences (see notes)
				//however even if there's a material mismatch at the interface, these contributions are still virtually zero since scaled by (cellsize / c11). complicated to include them and no real extra accuracy.
				vx[cuINT3(vx.n.x - 1, j, k)] = (vx[cuINT3(vx.n.x - 2, j, k)] * (1 + 3 * rdu) + 2 * vx_sec[ijk_sec + cuINT3(1, 0, 0)]) / (3 * (1 + rdu));
			}
		}
	}
}

__global__ void make_velocity_continuous_y_Kernel(
	ManagedMeshCUDA& cuMesh, ManagedMeshCUDA& cuMesh_sec,
	CMBNDInfoCUDA& contact,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuBReal>& vx_sec, cuVEC<cuBReal>& vy_sec, cuVEC<cuBReal>& vz_sec)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& u_disp_sec = *cuMesh_sec.pu_disp;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	const cuBox& cb = contact.cells_box;

	//kernel launch with size (cb.e.k + 1 - cb.s.k) * (cb.e.i + 1 - cb.s.i)
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int xdim = cb.e.i + 1 - cb.s.i;
	int zdim = cb.e.k + 1 - cb.s.k;

	if (idx < xdim * zdim) {

		int i = (idx % xdim) + cb.s.i;
		int j = (contact.IsPrimaryTop() ? cb.s.j : cb.e.j);
		int k = (idx / xdim) + cb.s.k;

		cuBReal spacing = (h_m.y + u_disp_sec.h.y);
		cuBReal wpri = 1.0 - h_m.y / spacing;
		cuBReal wsec = 1.0 - u_disp_sec.h.y / spacing;

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < u_disp.n.i ? i : u_disp.n.i - 1, j < u_disp.n.j ? j : u_disp.n.j - 1, k < u_disp.n.k ? k : u_disp.n.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * u_disp.n.x + ijk_u.k * u_disp.n.x * u_disp.n.y;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		int niend = (i < n_m.i);
		int njend = (j < n_m.j);
		int nkend = (k < n_m.k);

		//check for required faces being present (used for vx and vy components)
		bool yface_u = u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u);
		bool yface_l_vz = i > 0 && u_disp.is_not_empty(idx_u - niend) && u_disp.is_cmbnd(idx_u - niend);
		bool yface_l_vx = k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y);

		//check for required edge being present (used for vy component)
		bool yedge_u = yface_u || yface_l_vx || yface_l_vz || (i > 0 && k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y - niend) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y - niend));

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, -u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);
			int idx_u_sec = ijk_sec.i + ijk_sec.j * u_disp_sec.n.x + ijk_sec.k * u_disp_sec.n.x * u_disp_sec.n.y;

			//vx
			if (i < cb.e.i) {

				if (yface_u || yface_l_vx) vx[cuINT3(i, 0, k)] = vx[cuINT3(i, 1, k)] * wpri + vx_sec[ijk_sec + cuINT3(0, 0, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				if (yface_u || yface_l_vz) vz[cuINT3(i, 0, k)] = vz[cuINT3(i, 1, k)] * wpri + vz_sec[ijk_sec + cuINT3(i == cb.e.i, 0, 0)] * wsec;
			}

			//set vy continuity obtained from syy continuity
			if (yedge_u) {

				cuReal3 cC_p = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC_p);
				cuReal3 cC_s = *cuMesh_sec.pcC;
				cuMesh_sec.update_parameters_scoarse(idx_u_sec, *cuMesh_sec.pcC, cC_s);
				cuBReal rdu = (cC_p.i / cC_s.i) * (u_disp_sec.h.y / h_m.y);

				//simplest case just interpolate values either side (see derivation in notes)
				//in theory there are further contributions here due to 1) side derivatives, 2) thermoelastic contribution, 3) magnetostriction contribution
				//these contributions will be zero if there's no change in cC, alphaT, B1 coefficients across the interface since they are proportional to respective coefficient differences (see notes)
				//however even if there's a material mismatch at the interface, these contributions are still virtually zero since scaled by (cellsize / c11). complicated to include them and no real extra accuracy.
				vy[cuINT3(i, 0, k)] = (vy[cuINT3(i, 1, k)] * (1 + 3 * rdu) + 2 * vy_sec[ijk_sec + cuINT3(0, -1, 0)]) / (3 * (1 + rdu));
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);
			int idx_u_sec = ijk_sec.i + ijk_sec.j * u_disp_sec.n.x + ijk_sec.k * u_disp_sec.n.x * u_disp_sec.n.y;

			//vx
			if (i < cb.e.i) {

				vx[cuINT3(i, vx.n.y - 1, k)] = vx[cuINT3(i, vx.n.y - 2, k)] * wpri + vx_sec[ijk_sec + cuINT3(0, 1, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				vz[cuINT3(i, vz.n.y - 1, k)] = vz[cuINT3(i, vz.n.y - 1, k)] * wpri + vz_sec[ijk_sec + cuINT3(i == cb.e.i, 1, 0)] * wsec;
			}

			//set vy continuity obtained from syy continuity
			if (yedge_u) {

				cuReal3 cC_p = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC_p);
				cuReal3 cC_s = *cuMesh_sec.pcC;
				cuMesh_sec.update_parameters_scoarse(idx_u_sec, *cuMesh_sec.pcC, cC_s);
				cuBReal rdu = (cC_p.i / cC_s.i) * (u_disp_sec.h.y / h_m.y);

				//simplest case just interpolate values either side (see derivation in notes)
				//in theory there are further contributions here due to 1) side derivatives, 2) thermoelastic contribution, 3) magnetostriction contribution
				//these contributions will be zero if there's no change in cC, alphaT, B1 coefficients across the interface since they are proportional to respective coefficient differences (see notes)
				//however even if there's a material mismatch at the interface, these contributions are still virtually zero since scaled by (cellsize / c11). complicated to include them and no real extra accuracy.
				vy[cuINT3(i, vy.n.y - 1, k)] = (vy[cuINT3(i, vy.n.y - 2, k)] * (1 + 3 * rdu) + 2 * vy_sec[ijk_sec + cuINT3(0, 1, 0)]) / (3 * (1 + rdu));
			}
		}
	}
}

__global__ void make_velocity_continuous_z_Kernel(
	ManagedMeshCUDA& cuMesh, ManagedMeshCUDA& cuMesh_sec,
	CMBNDInfoCUDA& contact,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuBReal>& vx_sec, cuVEC<cuBReal>& vy_sec, cuVEC<cuBReal>& vz_sec)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	cuVEC_VC<cuReal3>& u_disp_sec = *cuMesh_sec.pu_disp;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	const cuBox& cb = contact.cells_box;

	//kernel launch with size (cb.e.i + 1 - cb.s.i) * (cb.e.j + 1 - cb.s.j)
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int xdim = cb.e.i + 1 - cb.s.i;
	int ydim = cb.e.j + 1 - cb.s.j;

	if (idx < xdim * ydim) {

		int i = (idx % xdim) + cb.s.i;
		int j = (idx / xdim) + cb.s.j;
		int k = (contact.IsPrimaryTop() ? cb.s.k : cb.e.k);

		cuBReal spacing = (h_m.z + u_disp_sec.h.z);
		cuBReal wpri = 1.0 - h_m.z / spacing;
		cuBReal wsec = 1.0 - u_disp_sec.h.z / spacing;

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < u_disp.n.i ? i : u_disp.n.i - 1, j < u_disp.n.j ? j : u_disp.n.j - 1, k < u_disp.n.k ? k : u_disp.n.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * u_disp.n.x + ijk_u.k * u_disp.n.x * u_disp.n.y;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		int niend = (i < n_m.i);
		int njend = (j < n_m.j);
		int nkend = (k < n_m.k);

		//check for required faces being present (used for vx and vy components)
		bool zface_u = u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u);
		bool zface_l_vy = i > 0 && u_disp.is_not_empty(idx_u - niend) && u_disp.is_cmbnd(idx_u - niend);
		bool zface_l_vx = j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x) && u_disp.is_cmbnd(idx_u - njend * n_m.x);

		//check for required edge being present (used for vz component)
		bool zedge_u = zface_u || zface_l_vx || zface_l_vy || (i > 0 && j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend) && u_disp.is_cmbnd(idx_u - njend * n_m.x - niend));

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, -u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);
			int idx_u_sec = ijk_sec.i + ijk_sec.j * u_disp_sec.n.x + ijk_sec.k * u_disp_sec.n.x * u_disp_sec.n.y;

			//vx
			if (i < cb.e.i) {

				if (zface_u || zface_l_vx) vx[cuINT3(i, j, 0)] = vx[cuINT3(i, j, 1)] * wpri + vx_sec[ijk_sec + cuINT3(0, j == cb.e.j, 0)] * wsec;
			}

			//vy
			if (j < cb.e.j) {

				if (zface_u || zface_l_vy) vy[cuINT3(i, j, 0)] = vy[cuINT3(i, j, 1)] * wpri + vy_sec[ijk_sec + cuINT3(i == cb.e.i, 0, 0)] * wsec;
			}

			//set vz continuity obtained from szz continuity
			if (zedge_u) {

				cuReal3 cC_p = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC_p);
				cuReal3 cC_s = *cuMesh_sec.pcC;
				cuMesh_sec.update_parameters_scoarse(idx_u_sec, *cuMesh_sec.pcC, cC_s);
				cuBReal rdu = (cC_p.i / cC_s.i) * (u_disp_sec.h.z / h_m.z);

				//simplest case just interpolate values either side (see derivation in notes)
				//in theory there are further contributions here due to 1) side derivatives, 2) thermoelastic contribution, 3) magnetostriction contribution
				//these contributions will be zero if there's no change in cC, alphaT, B1 coefficients across the interface since they are proportional to respective coefficient differences (see notes)
				//however even if there's a material mismatch at the interface, these contributions are still virtually zero since scaled by (cellsize / c11). complicated to include them and no real extra accuracy.
				vz[cuINT3(i, j, 0)] = (vz[cuINT3(i, j, 1)] * (1 + 3 * rdu) + 2 * vz_sec[ijk_sec + cuINT3(0, 0, -1)]) / (3 * (1 + rdu));
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);
			int idx_u_sec = ijk_sec.i + ijk_sec.j * u_disp_sec.n.x + ijk_sec.k * u_disp_sec.n.x * u_disp_sec.n.y;

			//vx
			if (i < cb.e.i) {

				if (zface_u || zface_l_vx) vx[cuINT3(i, j, vx.n.z - 1)] = vx[cuINT3(i, j, vx.n.z - 2)] * wpri + vx_sec[ijk_sec + cuINT3(0, j == cb.e.j, 1)] * wsec;
			}

			//vy
			if (j < cb.e.j) {

				if (zface_u || zface_l_vy) vy[cuINT3(i, j, vy.n.z - 1)] = vy[cuINT3(i, j, vy.n.z - 2)] * wpri + vy_sec[ijk_sec + cuINT3(i == cb.e.i, 0, 1)] * wsec;
			}

			//set vz continuity obtained from szz continuity
			if (zedge_u) {

				cuReal3 cC_p = *cuMesh.pcC;
				cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC_p);
				cuReal3 cC_s = *cuMesh_sec.pcC;
				cuMesh_sec.update_parameters_scoarse(idx_u_sec, *cuMesh_sec.pcC, cC_s);
				cuBReal rdu = (cC_p.i / cC_s.i) * (u_disp_sec.h.z / h_m.z);

				//simplest case just interpolate values either side (see derivation in notes)
				//in theory there are further contributions here due to 1) side derivatives, 2) thermoelastic contribution, 3) magnetostriction contribution
				//these contributions will be zero if there's no change in cC, alphaT, B1 coefficients across the interface since they are proportional to respective coefficient differences (see notes)
				//however even if there's a material mismatch at the interface, these contributions are still virtually zero since scaled by (cellsize / c11). complicated to include them and no real extra accuracy.
				vz[cuINT3(i, j, vz.n.z - 1)] = (vz[cuINT3(i, j, vz.n.z - 2)] * (1 + 3 * rdu) + 2 * vz_sec[ijk_sec + cuINT3(0, 0, 1)]) / (3 * (1 + rdu));
			}
		}
	}
}

//---------------------------------------------- CMBND routines (stress) LAUNCHER

void MElasticCUDA::make_velocity_continuous(
	cuSZ3 box_dims, int axis,
	cu_obj<CMBNDInfoCUDA>& contact,
	MElasticCUDA* pMElastic_sec)
{
	//+/-x normal face
	if (axis == 1) {

		make_velocity_continuous_x_Kernel <<< ((box_dims.j + 1)*(box_dims.k + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, pMElastic_sec->pMeshCUDA->cuMesh,
			contact,
			vx, vy, vz,
			pMElastic_sec->vx, pMElastic_sec->vy, pMElastic_sec->vz);
	}

	//+/-y normal face
	else if (axis == 2) {

		make_velocity_continuous_y_Kernel <<< ((box_dims.i + 1)*(box_dims.k + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, pMElastic_sec->pMeshCUDA->cuMesh, 
			contact,
			vx, vy, vz,
			pMElastic_sec->vx, pMElastic_sec->vy, pMElastic_sec->vz);
	}

	//+/-z normal face
	else if (axis == 3) {

		make_velocity_continuous_z_Kernel <<< ((box_dims.i + 1)*(box_dims.j + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, pMElastic_sec->pMeshCUDA->cuMesh, 
			contact,
			vx, vy, vz,
			pMElastic_sec->vx, pMElastic_sec->vy, pMElastic_sec->vz);
	}
}

//---------------------------------------------- CMBND routines (stress) KERNELS

__global__ void make_stress_continuous_x_Kernel(
	ManagedMeshCUDA& cuMesh,
	CMBNDInfoCUDA& contact,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuVEC<cuReal3>& sdd_sec, cuVEC<cuBReal>& sxy_sec, cuVEC<cuBReal>& sxz_sec, cuVEC<cuBReal>& syz_sec,
	cuVEC_VC<cuReal3>& u_disp_sec)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;
	
	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	const cuBox& cb = contact.cells_box;

	//kernel launch with size (cb.e.k + 1 - cb.s.k) * (cb.e.j + 1 - cb.s.j)
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int ydim = cb.e.j + 1 - cb.s.j;
	int zdim = cb.e.k + 1 - cb.s.k;

	if (idx < ydim * zdim) {

		int i = (contact.IsPrimaryTop() ? cb.s.i : cb.e.i);
		int j = (idx % ydim) + cb.s.j;
		int k = (idx / ydim) + cb.s.k;

		cuBReal spacing = (h_m.x + u_disp_sec.h.x);
		cuBReal wpri = 1.0 - h_m.x / spacing;
		cuBReal wsec = 1.0 - u_disp_sec.h.x / spacing;

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < u_disp.n.i ? i : u_disp.n.i - 1, j < u_disp.n.j ? j : u_disp.n.j - 1, k < u_disp.n.k ? k : u_disp.n.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * u_disp.n.x + ijk_u.k * u_disp.n.x * u_disp.n.y;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		int niend = (ijk.i < n_m.i);
		int njend = (ijk.j < n_m.j);
		int nkend = (ijk.k < n_m.k);

		//check if required edges are present
		bool yedge_u =
			ijk.j < n_m.j &&
			((u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) ||
			(ijk.k > 0 && (u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y))));

		bool yedge_l =
			ijk.j > 0 &&
			((u_disp.is_not_empty(idx_u - njend * n_m.x) && u_disp.is_cmbnd(idx_u - njend * n_m.x)) ||
			(ijk.k > 0 && (u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y - njend * n_m.x) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y - njend * n_m.x))));

		bool zedge_u =
			ijk.k < n_m.k &&
			((u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) ||
			(ijk.j > 0 && (u_disp.is_not_empty(idx_u - njend * n_m.x) && u_disp.is_cmbnd(idx_u - njend * n_m.x))));

		bool zedge_l =
			ijk.k > 0 &&
			((u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y)) ||
			(ijk.j > 0 && (u_disp.is_not_empty(idx_u - njend * n_m.x - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - njend * n_m.x - nkend * n_m.x * n_m.y))));

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(-u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			if ((yedge_u || yedge_l) && (zedge_u || zedge_l))
				sdd[cuINT3(0, j, k)] = sdd[cuINT3(1, j, k)] * wpri + sdd_sec[ijk_sec + cuINT3(0, j == cb.e.j, k == cb.e.k)] * wsec;

			//syz
			if (j < cb.e.j && k < cb.e.k && u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u))
				syz[cuINT3(0, j, k)] = syz[cuINT3(1, j, k)] * wpri + syz_sec[ijk_sec] * wsec;
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			if ((yedge_u || yedge_l) && (zedge_u || zedge_l))
				sdd[cuINT3(sdd.n.x - 1, j, k)] = sdd[cuINT3(sdd.n.x - 2, j, k)] * wpri + sdd_sec[ijk_sec + cuINT3(1, j == cb.e.j, k == cb.e.k)] * wsec;

			//syz
			if (j < cb.e.j && k < cb.e.k && u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u))
				syz[cuINT3(syz.n.x - 1, j, k)] = syz[cuINT3(syz.n.x - 2, j, k)] * wpri + syz_sec[ijk_sec + cuINT3(1, 0, 0)] * wsec;
		}
	}
}

__global__ void make_stress_continuous_y_Kernel(
	ManagedMeshCUDA& cuMesh,
	CMBNDInfoCUDA& contact,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuVEC<cuReal3>& sdd_sec, cuVEC<cuBReal>& sxy_sec, cuVEC<cuBReal>& sxz_sec, cuVEC<cuBReal>& syz_sec,
	cuVEC_VC<cuReal3>& u_disp_sec)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	const cuBox& cb = contact.cells_box;

	//kernel launch with size (cb.e.k + 1 - cb.s.k) * (cb.e.i + 1 - cb.s.i)
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int xdim = cb.e.i + 1 - cb.s.i;
	int zdim = cb.e.k + 1 - cb.s.k;

	if (idx < xdim * zdim) {

		int i = (idx % xdim) + cb.s.i;
		int j = (contact.IsPrimaryTop() ? cb.s.j : cb.e.j);
		int k = (idx / xdim) + cb.s.k;

		cuBReal spacing = (h_m.y + u_disp_sec.h.y);
		cuBReal wpri = 1.0 - h_m.y / spacing;
		cuBReal wsec = 1.0 - u_disp_sec.h.y / spacing;

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < u_disp.n.i ? i : u_disp.n.i - 1, j < u_disp.n.j ? j : u_disp.n.j - 1, k < u_disp.n.k ? k : u_disp.n.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * u_disp.n.x + ijk_u.k * u_disp.n.x * u_disp.n.y;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		int niend = (ijk.i < n_m.i);
		int njend = (ijk.j < n_m.j);
		int nkend = (ijk.k < n_m.k);

		//check if required edges are present
		bool xedge_u =
			ijk.i < n_m.i &&
			((u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) ||
			(ijk.k > 0 && (u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y))));

		bool xedge_l =
			ijk.i > 0 &&
			((u_disp.is_not_empty(idx_u - niend) && u_disp.is_cmbnd(idx_u - niend)) ||
			(ijk.k > 0 && (u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y - niend) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y - niend))));

		bool zedge_u =
			ijk.k < n_m.k &&
			((u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) ||
			(ijk.i > 0 && (u_disp.is_not_empty(idx_u - niend) && u_disp.is_cmbnd(idx_u - niend))));

		bool zedge_l =
			ijk.k > 0 &&
			((u_disp.is_not_empty(idx_u - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - nkend * n_m.x * n_m.y)) ||
			(ijk.i > 0 && (u_disp.is_not_empty(idx_u - niend - nkend * n_m.x * n_m.y) && u_disp.is_cmbnd(idx_u - niend - nkend * n_m.x * n_m.y))));

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, -u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			if ((xedge_u || xedge_l) && (zedge_u || zedge_l))
				sdd[cuINT3(i, 0, k)] = sdd[cuINT3(i, 1, k)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, 0, k == cb.e.k)] * wsec;

			//sxz
			if (i < cb.e.i && k < cb.e.k && u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u))
				sxz[cuINT3(i, 0, k)] = sxz[cuINT3(i, 1, k)] * wpri + sxz_sec[ijk_sec] * wsec;
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			if ((xedge_u || xedge_l) && (zedge_u || zedge_l))
				sdd[cuINT3(i, sdd.n.y - 1, k)] = sdd[cuINT3(i, sdd.n.y - 2, k)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, 1, k == cb.e.k)] * wsec;

			//sxz
			if (i < cb.e.j && k < cb.e.k && u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u))
				sxz[cuINT3(i, sxz.n.y - 1, k)] = sxz[cuINT3(i, sxz.n.y - 2, k)] * wpri + sxz_sec[ijk_sec + cuINT3(0, 1, 0)] * wsec;
		}
	}
}

__global__ void make_stress_continuous_z_Kernel(
	ManagedMeshCUDA& cuMesh,
	CMBNDInfoCUDA& contact,
	cuVEC<cuReal3>& sdd, cuVEC<cuBReal>& sxy, cuVEC<cuBReal>& sxz, cuVEC<cuBReal>& syz,
	cuVEC<cuReal3>& sdd_sec, cuVEC<cuBReal>& sxy_sec, cuVEC<cuBReal>& sxz_sec, cuVEC<cuBReal>& syz_sec,
	cuVEC_VC<cuReal3>& u_disp_sec)
{
	cuVEC_VC<cuReal3>& u_disp = *cuMesh.pu_disp;

	cuReal3& h_m = u_disp.h;
	cuSZ3& n_m = u_disp.n;

	const cuBox& cb = contact.cells_box;

	//kernel launch with size (cb.e.i + 1 - cb.s.i) * (cb.e.j + 1 - cb.s.j)
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int xdim = cb.e.i + 1 - cb.s.i;
	int ydim = cb.e.j + 1 - cb.s.j;

	if (idx < xdim * ydim) {

		int i = (idx % xdim) + cb.s.i;
		int j = (idx / xdim) + cb.s.j;
		int k = (contact.IsPrimaryTop() ? cb.s.k : cb.e.k);

		cuBReal spacing = (h_m.z + u_disp_sec.h.z);
		cuBReal wpri = 1.0 - h_m.z / spacing;
		cuBReal wsec = 1.0 - u_disp_sec.h.z / spacing;

		cuINT3 ijk = cuINT3(i, j, k);

		//convert vertex index to cell-center index by capping maximum index size (use this to index u_disp)
		cuINT3 ijk_u = cuINT3(i < u_disp.n.i ? i : u_disp.n.i - 1, j < u_disp.n.j ? j : u_disp.n.j - 1, k < u_disp.n.k ? k : u_disp.n.k - 1);
		int idx_u = ijk_u.i + ijk_u.j * u_disp.n.x + ijk_u.k * u_disp.n.x * u_disp.n.y;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		int niend = (ijk.i < n_m.i);
		int njend = (ijk.j < n_m.j);
		int nkend = (ijk.k < n_m.k);

		//check if required edges are present
		bool xedge_u =
			ijk.i < n_m.i &&
			((u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) ||
			(ijk.j > 0 && (u_disp.is_not_empty(idx_u - njend * n_m.x) && u_disp.is_cmbnd(idx_u - njend * n_m.x))));

		bool xedge_l =
			ijk.i > 0 &&
			((u_disp.is_not_empty(idx_u - niend) && u_disp.is_cmbnd(idx_u - niend)) ||
			(ijk.j > 0 && (u_disp.is_not_empty(idx_u - njend * n_m.x - niend) && u_disp.is_cmbnd(idx_u - njend * n_m.x - niend))));

		bool yedge_u =
			ijk.j < n_m.j &&
			((u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) ||
			(ijk.i > 0 && (u_disp.is_not_empty(idx_u - niend) && u_disp.is_cmbnd(idx_u - niend))));

		bool yedge_l =
			ijk.j > 0 &&
			((u_disp.is_not_empty(idx_u - njend * n_m.x) && u_disp.is_cmbnd(idx_u - njend * n_m.x)) ||
			(ijk.i > 0 && (u_disp.is_not_empty(idx_u - niend - njend * n_m.x) && u_disp.is_cmbnd(idx_u - niend - njend * n_m.x))));

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, -u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			if ((xedge_u || xedge_l) && (yedge_u || yedge_l))
				sdd[cuINT3(i, j, 0)] = sdd[cuINT3(i, j, 1)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, j == cb.e.j, 0)] * wsec;
	
			//sxy
			if (i < cb.e.i && j < cb.e.j && u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) 
				sxy[cuINT3(i, j, 0)] = sxy[cuINT3(i, j, 1)] * wpri + sxy_sec[ijk_sec] * wsec;
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			if ((xedge_u || xedge_l) && (yedge_u || yedge_l))
				sdd[cuINT3(i, j, sdd.n.z - 1)] = sdd[cuINT3(i, j, sdd.n.z - 2)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, j == cb.e.j, 1)] * wsec;

			//sxy
			if (i < cb.e.i && j < cb.e.j && u_disp.is_not_empty(idx_u) && u_disp.is_cmbnd(idx_u)) 
				sxy[cuINT3(i, j, sxy.n.z - 1)] = sxy[cuINT3(i, j, sxy.n.z - 2)] * wpri + sxy_sec[ijk_sec + cuINT3(0, 0, 1)] * wsec;
				
		}
	}
}

//---------------------------------------------- CMBND routines (stress) LAUNCHERS

void MElasticCUDA::make_stress_continuous(
	cuSZ3 box_dims, int axis,
	cu_obj<CMBNDInfoCUDA>& contact,
	cu_obj<cuVEC<cuReal3>>& sdd_sec, cu_obj<cuVEC<cuBReal>>& sxy_sec, cu_obj<cuVEC<cuBReal>>& sxz_sec, cu_obj<cuVEC<cuBReal>>& syz_sec,
	cu_obj<cuVEC_VC<cuReal3>>& u_disp_sec)
{
	//+/-x normal face
	if (axis == 1) {

		make_stress_continuous_x_Kernel <<< ((box_dims.j + 1)*(box_dims.k + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact, 
			sdd, sxy, sxz, syz,
			sdd_sec, sxy_sec, sxz_sec, syz_sec,
			u_disp_sec);
	}

	//+/-y normal face
	else if (axis == 2) {

		make_stress_continuous_y_Kernel <<< ((box_dims.i + 1)*(box_dims.k + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			sdd, sxy, sxz, syz,
			sdd_sec, sxy_sec, sxz_sec, syz_sec,
			u_disp_sec);
	}

	//+/-z normal face
	else if (axis == 3) {

		make_stress_continuous_z_Kernel <<< ((box_dims.i + 1)*(box_dims.j + 1) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			sdd, sxy, sxz, syz,
			sdd_sec, sxy_sec, sxz_sec, syz_sec,
			u_disp_sec);
	}
}

#endif

#endif
