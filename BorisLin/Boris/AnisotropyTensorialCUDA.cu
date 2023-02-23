#include "AnisotropyTensorialCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANITENS

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

__global__ void Anisotropy_TensorialCUDA_FM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal4>& Kt = *cuMesh.pKt;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Heff_value = cuReal3();

		if (M.is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal K1 = *cuMesh.pK1;
			cuBReal K2 = *cuMesh.pK2;
			cuBReal K3 = *cuMesh.pK3;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pK1, K1, *cuMesh.pK2, K2, *cuMesh.pK3, K3, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);
			
			//calculate dot products
			cuBReal a = (M[idx] * mcanis_ea1) / Ms;
			cuBReal b = (M[idx] * mcanis_ea2) / Ms;
			cuBReal c = (M[idx] * mcanis_ea3) / Ms;

			for (int tidx = 0; tidx < Kt.linear_size(); tidx++) {

				//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
				//(-d / mu0Ms) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

				cuBReal ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
				cuBReal ap = 0.0, bp = 0.0, cp = 0.0;
				if (Kt[tidx].j > 0) { ap1 = pow(a, Kt[tidx].j - 1); ap = ap1 * a; }
				else ap = pow(a, Kt[tidx].j);
				if (Kt[tidx].k > 0) { bp1 = pow(b, Kt[tidx].k - 1); bp = bp1 * b; }
				else bp = pow(b, Kt[tidx].k);
				if (Kt[tidx].l > 0) { cp1 = pow(c, Kt[tidx].l - 1); cp = cp1 * c; }
				else cp = pow(c, Kt[tidx].l);

				cuBReal coeff;
				int order = Kt[tidx].j + Kt[tidx].k + Kt[tidx].l;
				if (order == 2) coeff = -K1 * Kt[tidx].i / ((cuBReal)MU0*Ms);
				else if (order == 4) coeff = -K2 * Kt[tidx].i / ((cuBReal)MU0*Ms);
				else if (order == 6) coeff = -K3 * Kt[tidx].i / ((cuBReal)MU0*Ms);
				else coeff = -Kt[tidx].i / ((cuBReal)MU0*Ms);

				Heff_value += coeff * (Kt[tidx].j * ap1*bp*cp * mcanis_ea1 + Kt[tidx].k * ap*bp1*cp * mcanis_ea2 + Kt[tidx].l * ap*bp*cp1 * mcanis_ea3);

				energy_ += -coeff * (cuBReal)MU0*Ms * ap*bp*cp;
			}

			if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Heff_value;
			if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = energy_;

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ /= non_empty_cells;
			}
		}

		Heff[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void Anisotropy_TensorialCUDA_AFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal4>& Kt = *cuMesh.pKt;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;
	cuVEC<cuReal4>& Kt2 = *cuMesh.pKt2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Heff_value = cuReal3();
		cuReal3 Heff2_value = cuReal3();

		if (M.is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 K1_AFM = *cuMesh.pK1_AFM;
			cuReal2 K2_AFM = *cuMesh.pK2_AFM;
			cuReal2 K3_AFM = *cuMesh.pK3_AFM;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;
			cuReal3 mcanis_ea3 = *cuMesh.pmcanis_ea3;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pK1_AFM, K1_AFM, *cuMesh.pK2_AFM, K2_AFM, *cuMesh.pK3_AFM, K3_AFM, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2, *cuMesh.pmcanis_ea3, mcanis_ea3);

			//calculate dot products
			cuBReal a = (M[idx] * mcanis_ea1) / Ms_AFM.i;
			cuBReal b = (M[idx] * mcanis_ea2) / Ms_AFM.i;
			cuBReal c = (M[idx] * mcanis_ea3) / Ms_AFM.i;

			cuBReal a2 = (M2[idx] * mcanis_ea1) / Ms_AFM.j;
			cuBReal b2 = (M2[idx] * mcanis_ea2) / Ms_AFM.j;
			cuBReal c2 = (M2[idx] * mcanis_ea3) / Ms_AFM.j;

			cuBReal energy1_ = 0.0, energy2_ = 0.0;

			for (int tidx = 0; tidx < Kt.linear_size(); tidx++) {

				//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
				//(-d / mu0Ms) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

				cuBReal ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
				cuBReal ap = 0.0, bp = 0.0, cp = 0.0;
				if (Kt[tidx].j > 0) { ap1 = pow(a, Kt[tidx].j - 1); ap = ap1 * a; }
				else ap = pow(a, Kt[tidx].j);
				if (Kt[tidx].k > 0) { bp1 = pow(b, Kt[tidx].k - 1); bp = bp1 * b; }
				else bp = pow(b, Kt[tidx].k);
				if (Kt[tidx].l > 0) { cp1 = pow(c, Kt[tidx].l - 1); cp = cp1 * c; }
				else cp = pow(c, Kt[tidx].l);

				cuBReal coeff;
				int order = Kt[tidx].j + Kt[tidx].k + Kt[tidx].l;
				if (order == 2) coeff = -K1_AFM.i * Kt[tidx].i / ((cuBReal)MU0*Ms_AFM.i);
				else if (order == 4) coeff = -K2_AFM.i * Kt[tidx].i / ((cuBReal)MU0*Ms_AFM.i);
				else if (order == 6) coeff = -K3_AFM.i * Kt[tidx].i / ((cuBReal)MU0*Ms_AFM.i);
				else coeff = -Kt[tidx].i / ((cuBReal)MU0*Ms_AFM.i);

				Heff_value += coeff * (Kt[tidx].j * ap1*bp*cp * mcanis_ea1 + Kt[tidx].k * ap*bp1*cp * mcanis_ea2 + Kt[tidx].l * ap*bp*cp1 * mcanis_ea3);

				energy1_ += -coeff * (cuBReal)MU0*Ms_AFM.i * ap*bp*cp;
			}

			for (int tidx = 0; tidx < Kt2.linear_size(); tidx++) {

				//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
				//(-d / mu0Ms) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

				cuBReal ap1 = 0.0, bp1 = 0.0, cp1 = 0.0;
				cuBReal ap = 0.0, bp = 0.0, cp = 0.0;
				if (Kt2[tidx].j > 0) { ap1 = pow(a2, Kt2[tidx].j - 1); ap = ap1 * a2; }
				else ap = pow(a2, Kt2[tidx].j);
				if (Kt2[tidx].k > 0) { bp1 = pow(b2, Kt2[tidx].k - 1); bp = bp1 * b2; }
				else bp = pow(b2, Kt2[tidx].k);
				if (Kt2[tidx].l > 0) { cp1 = pow(c2, Kt2[tidx].l - 1); cp = cp1 * c2; }
				else cp = pow(c2, Kt2[tidx].l);

				cuBReal coeff;
				int order = Kt2[tidx].j + Kt2[tidx].k + Kt2[tidx].l;
				if (order == 2) coeff = -K1_AFM.j * Kt2[tidx].i / ((cuBReal)MU0*Ms_AFM.j);
				else if (order == 4) coeff = -K2_AFM.j * Kt2[tidx].i / ((cuBReal)MU0*Ms_AFM.j);
				else if (order == 6) coeff = -K3_AFM.j * Kt2[tidx].i / ((cuBReal)MU0*Ms_AFM.j);
				else coeff = -Kt2[tidx].i / ((cuBReal)MU0*Ms_AFM.j);

				Heff2_value += coeff * (Kt2[tidx].j * ap1*bp*cp * mcanis_ea1 + Kt2[tidx].k * ap*bp1*cp * mcanis_ea2 + Kt2[tidx].l * ap*bp*cp1 * mcanis_ea3);

				energy2_ += -coeff * (cuBReal)MU0*Ms_AFM.j * ap*bp*cp;
			}

			if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Heff_value;
			if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[idx] = Heff2_value;
			if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = energy1_;
			if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[idx] = energy2_;

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = (energy1_ + energy2_) / (2 * non_empty_cells);
			}
		}

		Heff[idx] += Heff_value;
		Heff2[idx] += Heff2_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void Anisotropy_TensorialCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Anisotropy_TensorialCUDA_AFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else {

			Anisotropy_TensorialCUDA_AFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
		}
	}
	else {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Anisotropy_TensorialCUDA_FM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else {

			Anisotropy_TensorialCUDA_FM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
		}
	}
}

#endif

#endif
