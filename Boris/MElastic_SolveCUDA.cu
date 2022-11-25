#include "MElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

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
	ZeroEnergy();

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			MElasticCUDA_UpdateField_AFM <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else {

			MElasticCUDA_UpdateField_AFM <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
		}
	}
	else {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			MElasticCUDA_UpdateField_FM <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else {

			MElasticCUDA_UpdateField_FM <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
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

		//show u_disp approximation for visualization (not used for further computations so it's fine)
		if (i < n_m.i && j < n_m.j && k < n_m.k) {

			u_disp[idx_u] += dT * cuReal3(vx[ijk], vy[ijk], vz[ijk]);
		}
	}
}

__global__ void Iterate_Elastic_Solver_Stress_Kernel(
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

		cuReal3 cC = *cuMesh.pcC;
		cuMesh.update_parameters_scoarse(idx_u, *cuMesh.pcC, cC);

		cuINT3 ijk = cuINT3(i, j, k);

		//external forces on different faces (keep track separately in case an edge cell is excited simultaneously by 2 or more external forces
		cuBReal Fext_xface = 0.0, Fext_yface = 0.0, Fext_zface = 0.0;

			
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
			(u_disp.is_not_empty(idx_u) ||
			(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
				(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
				(j > 0 && k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)));

		bool xedge_l =
			i > 0 &&
			(u_disp.is_not_empty(idx_u - niend) ||
			(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)) ||
				(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)) ||
				(j > 0 && k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x - niend)));

		bool yedge_u =
			j < n_m.j &&
			(u_disp.is_not_empty(idx_u) ||
			(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)) ||
				(i > 0 && u_disp.is_not_empty(idx_u - niend)) ||
				(i > 0 && k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend)));

		bool yedge_l =
			j > 0 &&
			(u_disp.is_not_empty(idx_u - njend * n_m.x) ||
			(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - njend * n_m.x)) ||
				(i > 0 && u_disp.is_not_empty(idx_u - niend - njend * n_m.x)) ||
				(i > 0 && k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y - niend - njend * n_m.x)));

		bool zedge_u =
			k < n_m.k &&
			(u_disp.is_not_empty(idx_u) ||
			(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)) ||
				(i > 0 && u_disp.is_not_empty(idx_u - niend)) ||
				(i > 0 && j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend)));

		bool zedge_l =
			k > 0 &&
			(u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y) ||
			(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - nkend * n_m.x*n_m.y)) ||
				(i > 0 && u_disp.is_not_empty(idx_u - niend - nkend * n_m.x*n_m.y)) ||
				(i > 0 && j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x - niend - nkend * n_m.x*n_m.y)));

		//check for fixed faces at ends
		bool xfixed_l = (i == 0 && u_disp.is_dirichlet_px(idx_u));
		bool xfixed_u = (i == n_m.i && u_disp.is_dirichlet_nx(idx_u));

		bool yfixed_l = (j == 0 && u_disp.is_dirichlet_py(idx_u));
		bool yfixed_u = (j == n_m.j && u_disp.is_dirichlet_ny(idx_u));

		bool zfixed_l = (k == 0 && u_disp.is_dirichlet_pz(idx_u));
		bool zfixed_u = (k == n_m.k && u_disp.is_dirichlet_nz(idx_u));

		cuBReal dvx_dx = 0.0;

		//interior
		if (xedge_u && xedge_l) dvx_dx = (vx[ijk] - vx[cuINT3(i - 1, j, k)]) / h_m.x;
		else if (xedge_l) {

			//is it a fixed face or free?
			if (xfixed_u) dvx_dx = -vx[cuINT3(i - 1, j, k)] / (h_m.x / 2);
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

		cuBReal dvy_dy = 0.0;

		//interior
		if (yedge_u && yedge_l) dvy_dy = (vy[ijk] - vy[cuINT3(i, j - 1, k)]) / h_m.y;
		//at +y face
		else if (yedge_l) {

			//is it a fixed face or free?
			if (yfixed_u) dvy_dy = -vy[cuINT3(i, j - 1, k)] / (h_m.y / 2);
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

		cuBReal dvz_dz = 0.0;

		//interior
		if (zedge_u && zedge_l) dvz_dz = (vz[ijk] - vz[cuINT3(i, j, k - 1)]) / h_m.z;
		//at +z face
		else if (zedge_l) {

			//is it a fixed face or free?
			if (zfixed_u) dvz_dz = -vz[cuINT3(i, j, k - 1)] / (h_m.z / 2);
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
				(u_disp.is_not_empty(idx_u) ||
				(k > 0 && u_disp.is_not_empty(idx_u - nkend * n_m.x*n_m.y)));

			if (zface) {

				cuBReal dvx_dy = (vx[cuINT3(i, j + 1, k)] - vx[ijk]) / h_m.y;
				cuBReal dvy_dx = (vy[cuINT3(i + 1, j, k)] - vy[ijk]) / h_m.x;

				sxy[ijk] += dT * cC.k * (dvx_dy + dvy_dx) / 2;
			}
			else sxy[ijk] = 0.0;
		}

		//update sxz
		if (i < n_m.i && k < n_m.k) {

			bool yface =
				(u_disp.is_not_empty(idx_u) ||
				(j > 0 && u_disp.is_not_empty(idx_u - njend * n_m.x)));

			if (yface) {

				cuBReal dvx_dz = (vx[cuINT3(i, j, k + 1)] - vx[ijk]) / h_m.z;
				cuBReal dvz_dx = (vz[cuINT3(i + 1, j, k)] - vz[ijk]) / h_m.x;

				sxz[ijk] += dT * cC.k * (dvx_dz + dvz_dx) / 2;
			}
			else sxz[ijk] = 0.0;
		}

		//update syz
		if (j < n_m.j && k < n_m.k) {

			bool xface =
				(u_disp.is_not_empty(idx_u) ||
				(i > 0 && u_disp.is_not_empty(idx_u - niend)));

			if (xface) {

				cuBReal dvy_dz = (vy[cuINT3(i, j, k + 1)] - vy[ijk]) / h_m.z;
				cuBReal dvz_dy = (vz[cuINT3(i, j + 1, k)] - vz[ijk]) / h_m.y;

				syz[ijk] += dT * cC.k * (dvy_dz + dvz_dy) / 2;
			}
			else syz[ijk] = 0.0;
		}
	}
}

//1c. Get strain from stress (cannot lump it in with step b since need to average over several vertices and faces to get cell-centred stress values)
__global__ void Calculate_Strain_from_Stress_Kernel(
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

	int i = idx % n_m.i;
	int j = (idx / n_m.i) % n_m.j;
	int k = idx / (n_m.i*n_m.j);

	if (idx < u_disp.linear_size()) {

		if (u_disp.is_not_empty(idx)) {

			cuReal3 cC = *cuMesh.pcC;
			cuMesh.update_parameters_scoarse(idx, *cuMesh.pcC, cC);

			cuINT3 ijk = cuINT3(i, j, k);

			//index for diagonal stress (sxx, syy, szz)
			int idx_sd = i + j * (n_m.i + 1) + k * (n_m.i + 1) * (n_m.j + 1);

			//invert the elastic coefficients matrix for diagonal stress-strain terms
			cuBReal det = cC.i*cC.i*cC.i + 2 * cC.j*cC.j*cC.j - 3 * cC.i*cC.j*cC.j;
			cuBReal diag = (cC.i*cC.i - cC.j*cC.j) / det;
			cuBReal odiag = -(cC.i*cC.j - cC.j*cC.j) / det;

			//strain - diagonal. Use average of the 8 vertex diagonal stress values to get cell-centre value.
			cuBReal sig_d_x = (
				sdd[idx_sd].x + sdd[idx_sd + 1].x + sdd[idx_sd + n_m.x + 1].x + sdd[idx_sd + n_m.x + 2].x +
				sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1)].x + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + 1].x + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 1].x + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 2].x) / 8;

			cuBReal sig_d_y = (
				sdd[idx_sd].y + sdd[idx_sd + 1].y + sdd[idx_sd + n_m.x + 1].y + sdd[idx_sd + n_m.x + 2].y +
				sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1)].y + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + 1].y + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 1].y + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 2].y) / 8;

			cuBReal sig_d_z = (
				sdd[idx_sd].z + sdd[idx_sd + 1].z + sdd[idx_sd + n_m.x + 1].z + sdd[idx_sd + n_m.x + 2].z +
				sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1)].z + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + 1].z + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 1].z + sdd[idx_sd + (n_m.x + 1)*(n_m.y + 1) + n_m.x + 2].z) / 8;

			strain_diag[idx].x = diag * sig_d_x + odiag * (sig_d_y + sig_d_z);
			strain_diag[idx].y = diag * sig_d_y + odiag * (sig_d_x + sig_d_z);
			strain_diag[idx].z = diag * sig_d_z + odiag * (sig_d_x + sig_d_y);

			//strain - off-diagonal (yz, xz, xy). Use average of 2 face off-diagonal stress values to get cell-centre value.
			strain_odiag[idx].x = (syz[ijk] + syz[ijk + cuINT3(1, 0, 0)]) / (2 * cC.k);
			strain_odiag[idx].y = (sxz[ijk] + sxz[ijk + cuINT3(0, 1, 0)]) / (2 * cC.k);
			strain_odiag[idx].z = (sxy[ijk] + sxy[ijk + cuINT3(0, 0, 1)]) / (2 * cC.k);
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
void MElasticCUDA::Iterate_Elastic_Solver_Stress(double dT)
{
	size_t size = (pMeshCUDA->n_m.i + 1) * (pMeshCUDA->n_m.j + 1) * (pMeshCUDA->n_m.k + 1);

	//1b. Update stress
	Iterate_Elastic_Solver_Stress_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(pMeshCUDA->cuMesh, external_stress_surfaces_arr, external_stress_surfaces.size(),
		vx, vy, vz, sdd, sxy, sxz, syz,
		pMeshCUDA->GetStageTime(), dT);
}

//update strain from stress
void MElasticCUDA::Calculate_Strain_From_Stress(void)
{
	size_t size = (pMeshCUDA->n_m.i + 1) * (pMeshCUDA->n_m.j + 1) * (pMeshCUDA->n_m.k + 1);

	//1c. Get strain from stress (cannot lump it in with step b since need to average over several vertices and faces to get cell-centred stress values)
	Calculate_Strain_from_Stress_Kernel <<< (pMeshCUDA->n_m.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(pMeshCUDA->cuMesh, sdd, sxy, sxz, syz);
}

//---------------------------------------------- CMBND routines (verlocity) KERNELS

__global__ void make_velocity_continuous_x_Kernel(
	ManagedMeshCUDA& cuMesh,
	CMBNDInfoCUDA& contact,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuBReal>& vx_sec, cuVEC<cuBReal>& vy_sec, cuVEC<cuBReal>& vz_sec,
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

		if (u_disp.is_empty(idx_u) || u_disp.is_not_cmbnd(idx_u)) return;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(-u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//vy
			if (j < cb.e.j) {

				//value at interface : interpolate between primary and secondary
				vy[cuINT3(0, j, k)] = vy[cuINT3(1, j, k)] * wpri + vy_sec[ijk_sec + cuINT3(0, 0, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				//value at interface : interpolate between primary and secondary
				vz[cuINT3(0, j, k)] = vz[cuINT3(1, j, k)] * wpri + vz_sec[ijk_sec + cuINT3(0, j == cb.e.j, 0)] * wsec;
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//vy
			if (j < cb.e.j) {

				vy[cuINT3(vy.n.x - 1, j, k)] = vy[cuINT3(vy.n.x - 2, j, k)] * wpri + vy_sec[ijk_sec + cuINT3(1, 0, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				vz[cuINT3(vz.n.x - 1, j, k)] = vz[cuINT3(vz.n.x - 2, j, k)] * wpri + vz_sec[ijk_sec + cuINT3(1, j == cb.e.j, 0)] * wsec;
			}
		}
	}
}

__global__ void make_velocity_continuous_y_Kernel(
	ManagedMeshCUDA& cuMesh,
	CMBNDInfoCUDA& contact,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuBReal>& vx_sec, cuVEC<cuBReal>& vy_sec, cuVEC<cuBReal>& vz_sec,
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

		if (u_disp.is_empty(idx_u) || u_disp.is_not_cmbnd(idx_u)) return;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, -u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//vx
			if (i < cb.e.i) {

				vx[cuINT3(i, 0, k)] = vx[cuINT3(i, 1, k)] * wpri + vx_sec[ijk_sec + cuINT3(0, 0, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				vz[cuINT3(i, 0, k)] = vz[cuINT3(i, 1, k)] * wpri + vz_sec[ijk_sec + cuINT3(i == cb.e.i, 0, 0)] * wsec;
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//vx
			if (i < cb.e.i) {

				vx[cuINT3(i, vx.n.y - 1, k)] = vx[cuINT3(i, vx.n.y - 2, k)] * wpri + vx_sec[ijk_sec + cuINT3(0, 1, k == cb.e.k)] * wsec;
			}

			//vz
			if (k < cb.e.k) {

				vz[cuINT3(i, vz.n.y - 1, k)] = vz[cuINT3(i, vz.n.y - 1, k)] * wpri + vz_sec[ijk_sec + cuINT3(i == cb.e.i, 1, 0)] * wsec;
			}
		}
	}
}

__global__ void make_velocity_continuous_z_Kernel(
	ManagedMeshCUDA& cuMesh,
	CMBNDInfoCUDA& contact,
	cuVEC<cuBReal>& vx, cuVEC<cuBReal>& vy, cuVEC<cuBReal>& vz,
	cuVEC<cuBReal>& vx_sec, cuVEC<cuBReal>& vy_sec, cuVEC<cuBReal>& vz_sec,
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

		if (u_disp.is_empty(idx_u) || u_disp.is_not_cmbnd(idx_u)) return;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, -u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//vx
			if (i < cb.e.i) {

				vx[cuINT3(i, j, 0)] = vx[cuINT3(i, j, 1)] * wpri + vx_sec[ijk_sec + cuINT3(0, j == cb.e.j, 0)] * wsec;
			}

			//vy
			if (j < cb.e.j) {

				vy[cuINT3(i, j, 0)] = vy[cuINT3(i, j, 1)] * wpri + vy_sec[ijk_sec + cuINT3(i == cb.e.i, 0, 0)] * wsec;
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//vx
			if (i < cb.e.i) {

				vx[cuINT3(i, j, vx.n.z - 1)] = vx[cuINT3(i, j, vx.n.z - 2)] * wpri + vx_sec[ijk_sec + cuINT3(0, j == cb.e.j, 1)] * wsec;
			}

			//vy
			if (j < cb.e.j) {

				vy[cuINT3(i, j, vy.n.z - 1)] = vy[cuINT3(i, j, vy.n.z - 2)] * wpri + vy_sec[ijk_sec + cuINT3(i == cb.e.i, 0, 1)] * wsec;
			}
		}
	}
}

//---------------------------------------------- CMBND routines (stress) LAUNCHER

void MElasticCUDA::make_velocity_continuous(
	size_t size, int axis,
	cu_obj<CMBNDInfoCUDA>& contact,
	cu_obj<cuVEC<cuBReal>>& vx_sec, cu_obj<cuVEC<cuBReal>>& vy_sec, cu_obj<cuVEC<cuBReal>>& vz_sec, cu_obj<cuVEC_VC<cuReal3>>& u_disp_sec)
{
	//+/-x normal face
	if (axis == 1) {

		make_velocity_continuous_x_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			vx, vy, vz,
			vx_sec, vy_sec, vz_sec,
			u_disp_sec);
	}

	//+/-y normal face
	else if (axis == 2) {

		make_velocity_continuous_y_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			vx, vy, vz,
			vx_sec, vy_sec, vz_sec,
			u_disp_sec);
	}

	//+/-z normal face
	else if (axis == 3) {

		make_velocity_continuous_z_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			vx, vy, vz,
			vx_sec, vy_sec, vz_sec,
			u_disp_sec);
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

		if (u_disp.is_empty(idx_u) || u_disp.is_not_cmbnd(idx_u)) return;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(-u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			sdd[cuINT3(0, j, k)] = sdd[cuINT3(1, j, k)] * wpri + sdd_sec[ijk_sec + cuINT3(0, j == cb.e.j, k == cb.e.k)] * wsec;

			//syz
			if (j < cb.e.j && k < cb.e.k) {

				syz[cuINT3(0, j, k)] = syz[cuINT3(1, j, k)] * wpri + syz_sec[ijk_sec] * wsec;
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3(u_disp_sec.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			sdd[cuINT3(sdd.n.x - 1, j, k)] = sdd[cuINT3(sdd.n.x - 2, j, k)] * wpri + sdd_sec[ijk_sec + cuINT3(1, j == cb.e.j, k == cb.e.k)] * wsec;

			//syz
			if (j < cb.e.j && k < cb.e.k) {

				syz[cuINT3(syz.n.x - 1, j, k)] = syz[cuINT3(syz.n.x - 2, j, k)] * wpri + syz_sec[ijk_sec + cuINT3(1, 0, 0)] * wsec;
			}
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

		if (u_disp.is_empty(idx_u) || u_disp.is_not_cmbnd(idx_u)) return;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, -u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			sdd[cuINT3(i, 0, k)] = sdd[cuINT3(i, 1, k)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, 0, k == cb.e.k)] * wsec;

			//sxz
			if (i < cb.e.i && k < cb.e.k) {

				sxz[cuINT3(i, 0, k)] = sxz[cuINT3(i, 1, k)] * wpri + sxz_sec[ijk_sec] * wsec;
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, u_disp_sec.h.y / 2, (k == cb.e.k ? -1 : +1) * u_disp.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			sdd[cuINT3(i, sdd.n.y - 1, k)] = sdd[cuINT3(i, sdd.n.y - 2, k)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, 1, k == cb.e.k)] * wsec;

			//sxz
			if (i < cb.e.i && k < cb.e.k) {

				sxz[cuINT3(i, sxz.n.y - 1, k)] = sxz[cuINT3(i, sxz.n.y - 2, k)] * wpri + sxz_sec[ijk_sec + cuINT3(0, 1, 0)] * wsec;
			}
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

		if (u_disp.is_empty(idx_u) || u_disp.is_not_cmbnd(idx_u)) return;

		//absolute position of interface vertex
		cuReal3 abs_pos = (u_disp.h & ijk) + u_disp.rect.s;

		if (contact.IsPrimaryTop()) {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, -u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			sdd[cuINT3(i, j, 0)] = sdd[cuINT3(i, j, 1)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, j == cb.e.j, 0)] * wsec;

			//sxy
			if (i < cb.e.i && j < cb.e.j) {

				sxy[cuINT3(i, j, 0)] = sxy[cuINT3(i, j, 1)] * wpri + sxy_sec[ijk_sec] * wsec;
			}
		}
		else {

			//absolute position in secondary, half a cellsize into it, and middle of primary cell face
			cuReal3 abs_pos_sec = abs_pos + cuReal3((i == cb.e.i ? -1 : +1) * u_disp.h.x / 2, (j == cb.e.j ? -1 : +1) * u_disp.h.y / 2, u_disp_sec.h.z / 2);
			//index of secondary cell just next to boundary
			cuINT3 ijk_sec = u_disp_sec.cellidx_from_position(abs_pos_sec);

			//sxx, syy, szz
			sdd[cuINT3(i, j, sdd.n.z - 1)] = sdd[cuINT3(i, j, sdd.n.z - 2)] * wpri + sdd_sec[ijk_sec + cuINT3(i == cb.e.i, j == cb.e.j, 1)] * wsec;

			//sxy
			if (i < cb.e.i && j < cb.e.j) {

				sxy[cuINT3(i, j, sxy.n.z - 1)] = sxy[cuINT3(i, j, sxy.n.z - 2)] * wpri + sxy_sec[ijk_sec + cuINT3(0, 0, 1)] * wsec;
			}
		}
	}
}

//---------------------------------------------- CMBND routines (stress) LAUNCHERS

void MElasticCUDA::make_stress_continuous(
	size_t size, int axis,
	cu_obj<CMBNDInfoCUDA>& contact,
	cu_obj<cuVEC<cuReal3>>& sdd_sec, cu_obj<cuVEC<cuBReal>>& sxy_sec, cu_obj<cuVEC<cuBReal>>& sxz_sec, cu_obj<cuVEC<cuBReal>>& syz_sec,
	cu_obj<cuVEC_VC<cuReal3>>& u_disp_sec)
{
	//+/-x normal face
	if (axis == 1) {

		make_stress_continuous_x_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact, 
			sdd, sxy, sxz, syz,
			sdd_sec, sxy_sec, sxz_sec, syz_sec,
			u_disp_sec);
	}

	//+/-y normal face
	else if (axis == 2) {

		make_stress_continuous_y_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			sdd, sxy, sxz, syz,
			sdd_sec, sxy_sec, sxz_sec, syz_sec,
			u_disp_sec);
	}

	//+/-z normal face
	else if (axis == 3) {

		make_stress_continuous_z_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, contact,
			sdd, sxy, sxz, syz,
			sdd_sec, sxy_sec, sxz_sec, syz_sec,
			u_disp_sec);
	}
}

#endif

#endif