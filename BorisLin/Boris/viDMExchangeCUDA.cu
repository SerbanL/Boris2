#include "viDMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_VIDMEXCHANGE

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"
#include "ModulesDefs.h"

//////////////////////////////////////////////////////////////////////// UPDATE FIELD

__global__ void viDMExchangeCUDA_FM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction, bool showtotalenergy)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hexch_A, Hexch_D;

		if (M.is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal A = *cuMesh.pA;
			cuBReal D = *cuMesh.pD;
			cuReal3 D_dir = *cuMesh.pD_dir;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A, *cuMesh.pD, D, *cuMesh.pD_dir, D_dir);

			cuBReal Aconst = 2 * A / ((cuBReal)MU0 * Ms * Ms);
			cuBReal Dconst = -2 * D / ((cuBReal)MU0 * Ms * Ms);

			if (M.is_plane_interior(idx)) {

				//interior point : can use cheaper neu versions

				//direct exchange contribution
				if (*cuMesh.pbase_temperature > 0.0 && *cuMesh.pT_Curie > 0.0) {

					//for finite temperature simulations the magnetization length may have a spatial variation
					//this will not affect the transverse torque (mxH), but will affect the longitudinal term in the sLLB equation (m.H) and cannot be neglected when close to Tc.

					cuReal33 Mg = M.grad_neu(idx);
					cuReal3 dMdx = Mg.x, dMdy = Mg.y, dMdz = Mg.z;

					cuBReal delsq_Msq = 2 * M[idx] * (M.dxx_neu(idx) + M.dyy_neu(idx) + M.dzz_neu(idx)) + 2 * (dMdx * dMdx + dMdy * dMdy + dMdz * dMdz);
					cuBReal Mnorm = M[idx].norm();
					Hexch_A = Aconst * (M.delsq_neu(idx) - M[idx] * delsq_Msq / (2 * Mnorm*Mnorm));
				}
				else {

					//zero temperature simulations : magnetization length could still vary but will only affect mxH term, so not needed for 0K simulations.
					Hexch_A = Aconst * M.delsq_neu(idx);
				}

				//Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				cuReal33 Mdiff = M.grad_neu(idx);

				cuReal3 hexch_D_x = cuReal3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
				cuReal3 hexch_D_y = cuReal3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
				cuReal3 hexch_D_z = cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

				Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
			}
			else {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				cuReal3 bnd_dm_dy_x = (D / (2 * A)) * cuReal3(-M[idx].y, M[idx].x, 0);
				cuReal3 bnd_dm_dz_x = (D / (2 * A)) * cuReal3(-M[idx].z, 0, M[idx].x);
				cuReal33 bnd_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

				cuReal3 bnd_dm_dx_y = (D / (2 * A)) * cuReal3(M[idx].y, -M[idx].x, 0);
				cuReal3 bnd_dm_dz_y = (D / (2 * A)) * cuReal3(0, -M[idx].z, M[idx].y);
				cuReal33 bnd_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

				cuReal3 bnd_dm_dx_z = (D / (2 * A)) * cuReal3(M[idx].z, 0, -M[idx].x);
				cuReal3 bnd_dm_dy_z = (D / (2 * A)) * cuReal3(0, M[idx].z, -M[idx].y);
				cuReal33 bnd_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

				cuReal33 bnd_nneu = D_dir.x * bnd_nneu_x + D_dir.y * bnd_nneu_y + D_dir.z * bnd_nneu_z;

				//direct exchange contribution
				if (*cuMesh.pbase_temperature > 0.0 && *cuMesh.pT_Curie > 0.0) {

					//for finite temperature simulations the magnetization length may have a spatial variation
					//this will not affect the transverse torque (mxH), but will affect the longitudinal term in the sLLB equation (m.H) and cannot be neglected when close to Tc.

					cuReal33 Mg = M.grad_nneu(idx, bnd_nneu);
					cuReal3 dMdx = Mg.x, dMdy = Mg.y, dMdz = Mg.z;

					cuBReal delsq_Msq = 2 * M[idx] * (M.dxx_nneu(idx, bnd_nneu) + M.dyy_nneu(idx, bnd_nneu) + M.dzz_nneu(idx, bnd_nneu)) + 2 * (dMdx * dMdx + dMdy * dMdy + dMdz * dMdz);
					cuBReal Mnorm = M[idx].norm();
					Hexch_A = Aconst * (M.delsq_nneu(idx, bnd_nneu) - M[idx] * delsq_Msq / (2 * Mnorm*Mnorm));
				}
				else {

					//zero temperature simulations : magnetization length could still vary but will only affect mxH term, so not needed for 0K simulations.
					Hexch_A = Aconst * M.delsq_nneu(idx, bnd_nneu);
				}

				//Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				cuReal33 Mdiff = M.grad_nneu(idx, bnd_nneu);

				cuReal3 hexch_D_x = cuReal3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
				cuReal3 hexch_D_y = cuReal3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
				cuReal3 hexch_D_z = cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

				Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
			}

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * (Hexch_A + Hexch_D) / (2 * non_empty_cells);
			}

			//spatial dependence display of effective field and energy density
			if (do_reduction && cuModule.pModule_Heff->linear_size() && cuModule.pModule_energy->linear_size()) {

				if (showtotalenergy) {

					//total : direct and DMI
					(*cuModule.pModule_Heff)[idx] = Hexch_A + Hexch_D;
					(*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * (M[idx] * (Hexch_A + Hexch_D)) / 2;
				}
				else {

					//just DMI
					(*cuModule.pModule_Heff)[idx] = Hexch_D;
					(*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * (M[idx] * Hexch_D) / 2;
				}
			}
		}

		Heff[idx] += (Hexch_A + Hexch_D);
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void viDMExchangeCUDA_AFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction, bool showtotalenergy)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hexch_A, Hexch_A2, Hexch_D, Hexch_D2;

		if (M.is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 A_AFM = *cuMesh.pA_AFM;
			cuReal2 Ah = *cuMesh.pAh;
			cuReal2 Anh = *cuMesh.pAnh;
			cuReal2 D_AFM = *cuMesh.pD_AFM;
			cuReal3 D_dir = *cuMesh.pD_dir;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pA_AFM, A_AFM, *cuMesh.pAh, Ah, *cuMesh.pAnh, Anh, *cuMesh.pD_AFM, D_AFM, *cuMesh.pD_dir, D_dir);

			if (M.is_plane_interior(idx)) {

				//interior point : can use cheaper neu versions

				//1. direct exchange contribution + AFM contribution
				cuReal3 delsq_M_A = M.delsq_neu(idx);
				cuReal3 delsq_M_B = M2.delsq_neu(idx);

				cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());

				Hexch_A = 2 * A_AFM.i * delsq_M_A / ((cuBReal)MU0 * Ms_AFM.i * Ms_AFM.i) + (-4 * Ah.i * (M[idx] ^ (M[idx] ^ M2[idx])) / (Mmag.i*Mmag.i) + Anh.i * delsq_M_B) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
				Hexch_A2 = 2 * A_AFM.j * delsq_M_B / ((cuBReal)MU0 * Ms_AFM.j * Ms_AFM.j) + (-4 * Ah.j * (M2[idx] ^ (M2[idx] ^ M[idx])) / (Mmag.j*Mmag.j) + Anh.j * delsq_M_A) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);

				//2. Dzyaloshinskii-Moriya interfacial exchange contribution

				cuReal33 Mdiff_A = M.grad_neu(idx);
				cuReal33 Mdiff_B = M2.grad_neu(idx);

				cuReal3 hexch_D_A_x = cuReal3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
				cuReal3 hexch_D_A_y = cuReal3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
				cuReal3 hexch_D_A_z = cuReal3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

				Hexch_D = -2 * D_AFM.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z) / ((cuBReal)MU0 * Ms_AFM.i * Ms_AFM.i);

				cuReal3 hexch_D_B_x = cuReal3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
				cuReal3 hexch_D_B_y = cuReal3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
				cuReal3 hexch_D_B_z = cuReal3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

				Hexch_D2 = -2 * D_AFM.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z) / ((cuBReal)MU0 * Ms_AFM.j * Ms_AFM.j);
			}
			else {

				cuReal33 bndA_nneu, bndB_nneu;

				cuReal2 nhconst = Anh / (2 * A_AFM);

				if (fabs(nhconst.i) != 1.0) {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					cuReal3 bnd_dm_dy_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(-M[idx].y + nhconst.i * M2[idx].y, M[idx].x - nhconst.i * M2[idx].x, 0);
					cuReal3 bnd_dm_dz_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(-M[idx].z + nhconst.i * M2[idx].z, 0, M[idx].x - nhconst.i * M2[idx].x);
					cuReal33 bndA_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

					cuReal3 bnd_dm_dx_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(M[idx].y - nhconst.i * M2[idx].y, -M[idx].x + nhconst.i * M2[idx].x, 0);
					cuReal3 bnd_dm_dz_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(0, -M[idx].z + nhconst.i * M2[idx].z, M[idx].y - nhconst.i * M2[idx].y);
					cuReal33 bndA_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

					cuReal3 bnd_dm_dx_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(M[idx].z - nhconst.i * M2[idx].z, 0, -M[idx].x + nhconst.i * M2[idx].x);
					cuReal3 bnd_dm_dy_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(0, M[idx].z - nhconst.i * M2[idx].z, -M[idx].y + nhconst.i * M2[idx].y);
					cuReal33 bndA_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

					bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
				}
				else {

					cuReal3 bnd_dm_dy_x = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(-M[idx].y, M[idx].x, 0);
					cuReal3 bnd_dm_dz_x = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(-M[idx].z, 0, M[idx].x);
					cuReal33 bndA_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

					cuReal3 bnd_dm_dx_y = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(M[idx].y, -M[idx].x, 0);
					cuReal3 bnd_dm_dz_y = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(0, -M[idx].z, M[idx].y);
					cuReal33 bndA_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

					cuReal3 bnd_dm_dx_z = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(M[idx].z, 0, -M[idx].x);
					cuReal3 bnd_dm_dy_z = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(0, M[idx].z, -M[idx].y);
					cuReal33 bndA_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

					bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
				}

				if (fabs(nhconst.j) != 1.0) {

					cuReal3 bnd_dm_dy_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(-M2[idx].y + nhconst.j * M[idx].y, M2[idx].x - nhconst.j * M[idx].x, 0);
					cuReal3 bnd_dm_dz_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(-M2[idx].z + nhconst.j * M[idx].z, 0, M2[idx].x - nhconst.j * M[idx].x);
					cuReal33 bndB_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

					cuReal3 bnd_dm_dx_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(M2[idx].y - nhconst.j * M[idx].y, -M2[idx].x + nhconst.j * M[idx].x, 0);
					cuReal3 bnd_dm_dz_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(0, -M2[idx].z + nhconst.j * M[idx].z, M2[idx].y - nhconst.j * M[idx].y);
					cuReal33 bndB_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

					cuReal3 bnd_dm_dx_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(M2[idx].z - nhconst.j * M[idx].z, 0, -M2[idx].x + nhconst.j * M[idx].x);
					cuReal3 bnd_dm_dy_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(0, M2[idx].z - nhconst.j * M[idx].z, -M2[idx].y + nhconst.j * M[idx].y);
					cuReal33 bndB_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

					bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
				}
				else {

					cuReal3 bnd_dm_dy_x = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(-M2[idx].y, M2[idx].x, 0);
					cuReal3 bnd_dm_dz_x = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(-M2[idx].z, 0, M2[idx].x);
					cuReal33 bndB_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

					cuReal3 bnd_dm_dx_y = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(M2[idx].y, -M2[idx].x, 0);
					cuReal3 bnd_dm_dz_y = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(0, -M2[idx].z, M2[idx].y);
					cuReal33 bndB_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

					cuReal3 bnd_dm_dx_z = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(M2[idx].z, 0, -M2[idx].x);
					cuReal3 bnd_dm_dy_z = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(0, M2[idx].z, -M2[idx].y);
					cuReal33 bndB_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

					bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
				}

				cuReal3 delsq_M_A = M.delsq_nneu(idx, bndA_nneu);
				cuReal3 delsq_M_B = M2.delsq_nneu(idx, bndB_nneu);

				cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());

				//1. direct exchange contribution + AFM contribution

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
				Hexch_A = 2 * A_AFM.i * delsq_M_A / ((cuBReal)MU0 * Ms_AFM.i * Ms_AFM.i) + (-4 * Ah.i * (M[idx] ^ (M[idx] ^ M2[idx])) / (Mmag.i*Mmag.i) + Anh.i * delsq_M_B) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
				Hexch_A2 = 2 * A_AFM.j * delsq_M_B / ((cuBReal)MU0 * Ms_AFM.j * Ms_AFM.j) + (-4 * Ah.j * (M2[idx] ^ (M2[idx] ^ M[idx])) / (Mmag.j*Mmag.j) + Anh.j * delsq_M_A) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);

				//2. Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
				cuReal33 Mdiff_A = M.grad_nneu(idx, bndA_nneu);
				cuReal33 Mdiff_B = M2.grad_nneu(idx, bndB_nneu);

				cuReal3 hexch_D_A_x = cuReal3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
				cuReal3 hexch_D_A_y = cuReal3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
				cuReal3 hexch_D_A_z = cuReal3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

				Hexch_D = -2 * D_AFM.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z) / ((cuBReal)MU0 * Ms_AFM.i * Ms_AFM.i);

				cuReal3 hexch_D_B_x = cuReal3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
				cuReal3 hexch_D_B_y = cuReal3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
				cuReal3 hexch_D_B_z = cuReal3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

				Hexch_D2 = -2 * D_AFM.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z) / ((cuBReal)MU0 * Ms_AFM.j * Ms_AFM.j);
			}

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] * (Hexch_A + Hexch_D) + M2[idx] * (Hexch_A2 + Hexch_D2)) / (4 * non_empty_cells);
			}

			//spatial dependence display of effective field and energy density
			if (do_reduction && cuModule.pModule_Heff->linear_size() && cuModule.pModule_energy->linear_size()) {

				if (showtotalenergy) {

					//total : direct and DMI
					(*cuModule.pModule_Heff)[idx] = Hexch_A + Hexch_D;
					(*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * (M[idx] * (Hexch_A + Hexch_D)) / 2;
				}
				else {

					//just DMI
					(*cuModule.pModule_Heff)[idx] = Hexch_D;
					(*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * (M[idx] * Hexch_D) / 2;
				}
			}

			if (do_reduction && cuModule.pModule_Heff2->linear_size() && cuModule.pModule_energy2->linear_size()) {

				if (showtotalenergy) {

					//total : direct and DMI
					(*cuModule.pModule_Heff2)[idx] = Hexch_A2 + Hexch_D2;
					(*cuModule.pModule_energy2)[idx] = -(cuBReal)MU0 * (M2[idx] * (Hexch_A2 + Hexch_D2)) / 2;
				}
				else {

					//just DMI
					(*cuModule.pModule_Heff2)[idx] = Hexch_D2;
					(*cuModule.pModule_energy2)[idx] = -(cuBReal)MU0 * (M2[idx] * Hexch_D2) / 2;
				}
			}
		}

		Heff[idx] += (Hexch_A + Hexch_D);
		Heff2[idx] += (Hexch_A2 + Hexch_D2);
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void viDMExchangeCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			viDMExchangeCUDA_AFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, cuModule, true, (MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_EXCHANGE);
		}
		else {

			viDMExchangeCUDA_AFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, cuModule, false, (MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_EXCHANGE);
		}
	}
	else {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			viDMExchangeCUDA_FM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, cuModule, true, (MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_EXCHANGE);
		}
		else {

			viDMExchangeCUDA_FM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, cuModule, false, (MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_EXCHANGE);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// COUPLING ACROSS MULTIPLE MESHES ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//Not implemented for this module
}

#endif

#endif