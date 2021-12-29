#include "ManagedMeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"
#include "MeshParamsControlCUDA.h"

#include "ModulesDefs.h"

//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

//Energy Deltas

//Ferromagnetic

//switch function which adds all assigned energy contributions in this mesh to calculate energy change from current spin to Mnew spin : return energy change as new - old
//If Mnew is passed in as cuReal3(), then this function returns the current spin energy only - all functions in the switch statement below implement this eventuality.
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM(int spin_index, cuReal3 Mnew, int*& cuModules, int numModules, cuReal3& Ha)
{
	cuBReal energy = 0.0;

	for (int idx = 0; idx < numModules; idx++) {

		switch (cuModules[idx]) {

		case MOD_DEMAG_N:
			energy += Get_EnergyChange_FM_DemagNCUDA(spin_index, Mnew);
			break;

		case MOD_DEMAG:
			energy += Get_EnergyChange_FM_DemagCUDA(spin_index, Mnew);
			break;

		case MOD_SDEMAG_DEMAG:
			//same method as for MOD_DEMAG, but the effective field pointer now points to SDemag_DemagCUDA module effective field
			energy += Get_EnergyChange_FM_DemagCUDA(spin_index, Mnew);
			break;

		case MOD_EXCHANGE:
			energy += Get_EnergyChange_FM_ExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_DMEXCHANGE:
			energy += Get_EnergyChange_FM_DMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_EnergyChange_FM_iDMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_VIDMEXCHANGE:
			energy += Get_EnergyChange_FM_viDMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_SURFEXCHANGE:
			energy += Get_EnergyChange_FM_SurfExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_ZEEMAN:
			energy += Get_EnergyChange_FM_ZeemanCUDA(spin_index, Mnew, Ha);
			break;

		case MOD_MOPTICAL:
			energy += Get_EnergyChange_FM_MOpticalCUDA(spin_index, Mnew);
			break;

		case MOD_ANIUNI:
			energy += Get_EnergyChange_FM_AnisotropyCUDA(spin_index, Mnew);
			break;

		case MOD_ANICUBI:
			energy += Get_EnergyChange_FM_AnisotropyCubiCUDA(spin_index, Mnew);
			break;

		case MOD_ANIBI:
			energy += Get_EnergyChange_FM_AnisotropyBiaxialCUDA(spin_index, Mnew);
			break;

		case MOD_ANITENS:
			energy += Get_EnergyChange_FM_AnisotropyTensorialCUDA(spin_index, Mnew);
			break;

		case MOD_ROUGHNESS:
			energy += Get_EnergyChange_FM_RoughnessCUDA(spin_index, Mnew);
			break;

		case MOD_MELASTIC:
			energy += Get_EnergyChange_FM_MElasticCUDA(spin_index, Mnew);
			break;

		default:
			break;
		}
	}

	return energy;
}

//Demag_N
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_DemagNCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
	cuReal2 Nxy = *pNxy;
	cuBReal Nz = (1 - Nxy.x - Nxy.y);

	if (Mnew != cuReal3()) {

		return ((cuBReal)MU0 / 2) * M.h.dim() * (
			(Mnew * cuReal3(Nxy.x * Mnew.x, Nxy.y * Mnew.y, Nz * Mnew.z)) 
			- (M[spin_index] * cuReal3(Nxy.x * M[spin_index].x, Nxy.y * M[spin_index].y, Nz * M[spin_index].z)));
	}
	else return ((cuBReal)MU0 / 2) * M.h.dim() * (M[spin_index] * cuReal3(Nxy.x * M[spin_index].x, Nxy.y * M[spin_index].y, Nz * M[spin_index].z));
}

//Demag
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_DemagCUDA(int spin_index, cuReal3 Mnew)
{
	if (pDemag_Heff && pDemag_Heff->linear_size()) {

		cuVEC_VC<cuReal3>& M = *pM;

		if (Mnew != cuReal3()) return -(cuBReal)MU0 * M.h.dim() * (*pDemag_Heff)[M.cellidx_to_position(spin_index)] * (Mnew - M[spin_index]);
		else return -(cuBReal)MU0 * M.h.dim() * (*pDemag_Heff)[M.cellidx_to_position(spin_index)] * M[spin_index];
	}
	else return 0.0;
}

//Exch_6ngbr_Neu
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_ExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal Ms = *pMs;
	cuBReal A = *pA;
	update_parameters_mcoarse(spin_index, *pA, A, *pMs, Ms);

	cuReal3 Hexch = (2 * A / ((cuBReal)MU0*Ms*Ms)) * M.delsq_neu(spin_index);
	cuBReal energy_ = -(cuBReal)MU0 * M[spin_index] * Hexch / 2;

	if (Mnew != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew.renormalize(M[spin_index].norm());

		cuReal3 Mold = M[spin_index];
		M[spin_index] = Mnew;
		Hexch = (2 * A / ((cuBReal)MU0*Ms*Ms)) * M.delsq_neu(spin_index);
		cuBReal energynew_ = -(cuBReal)MU0 * M[spin_index] * Hexch / 2;
		M[spin_index] = Mold;

		//multiply by 2 as we are not double-counting here (same as for Demag)
		return M.h.dim() * (energynew_ - energy_) * 2;
	}
	else return M.h.dim() * energy_ * 2;
}

//DMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_DMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal Ms = *pMs;
	cuBReal A = *pA;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pA, A, *pD, D, *pMs, Ms);

	cuBReal Aconst = 2 * A / ((cuBReal)MU0 * Ms * Ms);
	cuBReal Dconst = -2 * D / ((cuBReal)MU0 * Ms * Ms);

	auto Get_Energy = [&](void) -> cuBReal
	{
		cuReal3 Hexch_A, Hexch_D;

		if (M.is_interior(spin_index)) {

			//interior point : can use cheaper neu versions

			//direct exchange contribution
			Hexch_A = Aconst * M.delsq_neu(spin_index);

			//Dzyaloshinskii-Moriya exchange contribution

			//Hdm, ex = -2D / (mu0*Ms) * curl m
			Hexch_D = Dconst * M.curl_neu(spin_index);
		}
		else {

			//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. equivalent to m x h -> 0 when relaxing.
			cuReal3 bnd_dm_dx = (D / (2 * A)) * cuReal3(0, -M[spin_index].z, M[spin_index].y);
			cuReal3 bnd_dm_dy = (D / (2 * A)) * cuReal3(M[spin_index].z, 0, -M[spin_index].x);
			cuReal3 bnd_dm_dz = (D / (2 * A)) * cuReal3(-M[spin_index].y, M[spin_index].x, 0);
			cuReal33 bnd_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, bnd_dm_dz);

			//direct exchange contribution
			//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
			Hexch_A = Aconst * M.delsq_nneu(spin_index, bnd_nneu);

			//Dzyaloshinskii-Moriya exchange contribution

			//Hdm, ex = -2D / (mu0*Ms) * curl m
			//For cmbnd cells curl_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
			Hexch_D = Dconst * M.curl_nneu(spin_index, bnd_nneu);
		}

		return M[spin_index] * (Hexch_A + Hexch_D);
	};

	cuBReal energy_ = Get_Energy();

	if (Mnew != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew.renormalize(M[spin_index].norm());

		//new spin energy
		cuReal3 Mold = M[spin_index];
		M[spin_index] = Mnew;
		cuBReal energynew_ = Get_Energy();
		M[spin_index] = Mold;

		//do not divide by 2 as we are not double-counting here
		return -MU0 * M.h.dim() * (energynew_ - energy_);
	}
	else return -MU0 * M.h.dim() * energy_;
}

//iDMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_iDMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal Ms = *pMs;
	cuBReal A = *pA;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pA, A, *pD, D, *pMs, Ms);

	cuBReal Aconst = 2 * A / ((cuBReal)MU0 * Ms * Ms);
	cuBReal Dconst = -2 * D / ((cuBReal)MU0 * Ms * Ms);

	auto Get_Energy = [&](void) -> cuBReal
	{
		cuReal3 Hexch_A, Hexch_D;

		if (M.is_interior(spin_index)) {

			//interior point : can use cheaper neu versions

			//direct exchange contribution
			Hexch_A = Aconst * M.delsq_neu(spin_index);

			//Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			cuReal33 Mdiff = M.grad_neu(spin_index);

			//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
			Hexch_D = Dconst * cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
		}
		else {

			//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
			cuReal3 bnd_dm_dx = (D / (2 * A)) * cuReal3(M[spin_index].z, 0, -M[spin_index].x);
			cuReal3 bnd_dm_dy = (D / (2 * A)) * cuReal3(0, M[spin_index].z, -M[spin_index].y);
			cuReal33 bnd_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, cuReal3());

			//direct exchange contribution
			//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
			Hexch_A = Aconst * M.delsq_nneu(spin_index, bnd_nneu);

			//Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
			cuReal33 Mdiff = M.grad_nneu(spin_index, bnd_nneu);

			//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
			Hexch_D = Dconst * cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
		}

		return M[spin_index] * (Hexch_A + Hexch_D);
	};

	cuBReal energy_ = Get_Energy();

	if (Mnew != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew.renormalize(M[spin_index].norm());

		//new spin energy
		cuReal3 Mold = M[spin_index];
		M[spin_index] = Mnew;
		cuBReal energynew_ = Get_Energy();
		M[spin_index] = Mold;

		//do not divide by 2 as we are not double-counting here
		return -MU0 * M.h.dim() * (energynew_ - energy_);
	}
	else return -MU0 * M.h.dim() * energy_;
}

//viDMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_viDMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal Ms = *pMs;
	cuBReal A = *pA;
	cuBReal D = *pD;
	cuReal3 D_dir = *pD_dir;
	update_parameters_mcoarse(spin_index, *pA, A, *pD, D, *pMs, Ms, *pD_dir, D_dir);

	cuBReal Aconst = 2 * A / ((cuBReal)MU0 * Ms * Ms);
	cuBReal Dconst = -2 * D / ((cuBReal)MU0 * Ms * Ms);

	auto Get_Energy = [&](void) -> cuBReal
	{
		cuReal3 Hexch_A, Hexch_D;

		if (M.is_plane_interior(spin_index)) {

			//interior point : can use cheaper neu versions

			//direct exchange contribution
			Hexch_A = Aconst * M.delsq_neu(spin_index);

			//Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			cuReal33 Mdiff = M.grad_neu(spin_index);

			cuReal3 hexch_D_x = cuReal3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
			cuReal3 hexch_D_y = cuReal3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
			cuReal3 hexch_D_z = cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

			Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
		}
		else {

			//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
			cuReal3 bnd_dm_dy_x = (D / (2 * A)) * cuReal3(-M[spin_index].y, M[spin_index].x, 0);
			cuReal3 bnd_dm_dz_x = (D / (2 * A)) * cuReal3(-M[spin_index].z, 0, M[spin_index].x);
			cuReal33 bnd_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

			cuReal3 bnd_dm_dx_y = (D / (2 * A)) * cuReal3(M[spin_index].y, -M[spin_index].x, 0);
			cuReal3 bnd_dm_dz_y = (D / (2 * A)) * cuReal3(0, -M[spin_index].z, M[spin_index].y);
			cuReal33 bnd_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

			cuReal3 bnd_dm_dx_z = (D / (2 * A)) * cuReal3(M[spin_index].z, 0, -M[spin_index].x);
			cuReal3 bnd_dm_dy_z = (D / (2 * A)) * cuReal3(0, M[spin_index].z, -M[spin_index].y);
			cuReal33 bnd_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

			cuReal33 bnd_nneu = D_dir.x * bnd_nneu_x + D_dir.y * bnd_nneu_y + D_dir.z * bnd_nneu_z;

			//direct exchange contribution
			//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
			Hexch_A = Aconst * M.delsq_nneu(spin_index, bnd_nneu);

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			cuReal33 Mdiff = M.grad_nneu(spin_index, bnd_nneu);

			cuReal3 hexch_D_x = cuReal3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
			cuReal3 hexch_D_y = cuReal3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
			cuReal3 hexch_D_z = cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

			Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
		}

		return M[spin_index] * (Hexch_A + Hexch_D);
	};

	cuBReal energy_ = Get_Energy();

	if (Mnew != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew.renormalize(M[spin_index].norm());

		//new spin energy
		cuReal3 Mold = M[spin_index];
		M[spin_index] = Mnew;
		cuBReal energynew_ = Get_Energy();
		M[spin_index] = Mold;

		//do not divide by 2 as we are not double-counting here
		return -MU0 * M.h.dim() * (energynew_ - energy_);
	}
	else return -MU0 * M.h.dim() * energy_;
}

//SurfExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_SurfExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuBReal energy_new = 0, energy_old = 0;

	cuVEC_VC<cuReal3>& M = *pM;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	//thickness of layer - SurfExchange applies for layers in the xy plane
	cuBReal thickness = M.rect.e.z - M.rect.s.z;

	//if spin is on top surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == n.z - 1 && (pMeshFM_Top_size + pMeshAFM_Top_size)) {

		if (!M.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			cuBReal Ms = *pMs;
			update_parameters_mcoarse(spin_index, *pMs, Ms);

			//check all meshes for coupling : FM meshes first
			for (int mesh_idx = 0; mesh_idx < (int)pMeshFM_Top_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Top = *(pMeshFM_Top[mesh_idx].pM);

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Top.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Top.rect.s.y,
					M_Top.h.z / 2);

				//can't couple to an empty cell
				if (!M_Top.rect.contains(cell_rel_pos + M_Top.rect.s) || M_Top.is_empty(cell_rel_pos)) continue;

				//Surface exchange field from a ferromagnetic mesh (RKKY)

				//Top mesh sets J1 and J2 values
				cuBReal J1 = *(pMeshFM_Top[mesh_idx].pJ1);
				cuBReal J2 = *(pMeshFM_Top[mesh_idx].pJ2);
				pMeshFM_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMeshFM_Top[mesh_idx].pJ1), J1, *(pMeshFM_Top[mesh_idx].pJ2), J2);

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Top[cell_rel_pos].normalized();
				cuReal3 m_i = M[spin_index] / Ms;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuBReal dot_prod = m_i * m_j;
				energy_old = (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / thickness;

				if (Mnew != cuReal3()) {

					cuReal3 mnew_i = Mnew / Ms;
					cuBReal dot_prod_new = mnew_i * m_j;
					energy_new = (-1 * J1 - 2 * J2 * dot_prod_new) * dot_prod_new / thickness;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}

			//next AFM meshes
			for (int mesh_idx = 0; mesh_idx < (int)pMeshAFM_Top_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Top = *(pMeshAFM_Top[mesh_idx].pM);
				cuVEC_VC<cuReal3>& M2_Top = *(pMeshAFM_Top[mesh_idx].pM2);

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Top.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Top.rect.s.y,
					M_Top.h.z / 2);

				//can't couple to an empty cell
				if (!M_Top.rect.contains(cell_rel_pos + M_Top.rect.s) || M_Top.is_empty(cell_rel_pos)) continue;

				//Surface exchange field from an antiferromagnetic mesh (exchange bias)

				//Top mesh sets J1 and J2 values
				cuBReal J1 = *(pMeshAFM_Top[mesh_idx].pJ1);
				cuBReal J2 = *(pMeshAFM_Top[mesh_idx].pJ2);
				pMeshAFM_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMeshAFM_Top[mesh_idx].pJ1), J1, *(pMeshAFM_Top[mesh_idx].pJ2), J2);

				//get magnetization values in top mesh cell to couple with
				cuReal3 m_j1 = M_Top[cell_rel_pos].normalized();
				cuReal3 m_j2 = M2_Top[cell_rel_pos].normalized();
				cuReal3 m_i = M[spin_index] / Ms;

				//total surface exchange field in coupling cells, including contributions from both sub-lattices
				energy_old += (-J1 * (m_i * m_j1) - J2 * (m_i * m_j2)) / thickness;

				if (Mnew != cuReal3()) {

					cuReal3 mnew_i = Mnew / Ms;
					energy_new += (-J1 * (mnew_i * m_j1) - J2 * (mnew_i * m_j2)) / thickness;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (spin_index / (n.x * n.y) == 0 && (pMeshFM_Bot_size + pMeshAFM_Bot_size)) {

		//surface exchange coupling at the bottom

		if (!M.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			cuBReal Ms = *pMs;
			cuBReal J1 = *pJ1;
			cuBReal J2 = *pJ2;
			update_parameters_mcoarse(spin_index, *pMs, Ms, *pJ1, J1, *pJ2, J2);

			//check all meshes for coupling : FM meshes first
			for (int mesh_idx = 0; mesh_idx < (int)pMeshFM_Bot_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMeshFM_Bot[mesh_idx].pM);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Bot.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Bot.rect.s.y,
					M_Bot.rect.e.z - M_Bot.rect.s.z - M_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M_Bot.rect.contains(cell_rel_pos + M_Bot.rect.s) || M_Bot.is_empty(cell_rel_pos)) continue;

				//Surface exchange field from a ferromagnetic mesh (RKKY)

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M[spin_index] / Ms;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuBReal dot_prod = m_i * m_j;
				energy_old += (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / thickness;

				if (Mnew != cuReal3()) {

					cuReal3 mnew_i = Mnew / Ms;
					cuBReal dot_prod_new = mnew_i * m_j;
					energy_new += (-1 * J1 - 2 * J2 * dot_prod_new) * dot_prod_new / thickness;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}

			//next AFM meshes
			for (int mesh_idx = 0; mesh_idx < (int)pMeshAFM_Bot_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMeshAFM_Bot[mesh_idx].pM);
				cuVEC_VC<cuReal3>& M2_Bot = *(pMeshAFM_Bot[mesh_idx].pM2);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Bot.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Bot.rect.s.y,
					M_Bot.rect.e.z - M_Bot.rect.s.z - M_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M_Bot.rect.contains(cell_rel_pos + M_Bot.rect.s) || M_Bot.is_empty(cell_rel_pos)) continue;

				//Surface exchange field from an antiferromagnetic mesh (exchange bias)

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j1 = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_j2 = M2_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M[spin_index] / Ms;

				//total surface exchange field in coupling cells, including contributions from both sub-lattices
				energy_old += (-J1 * (m_i * m_j1) - J2 * (m_i * m_j2)) / thickness;

				if (Mnew != cuReal3()) {

					cuReal3 mnew_i = Mnew / Ms;
					energy_new += (-J1 * (mnew_i * m_j1) - J2 * (mnew_i * m_j2)) / thickness;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	//multiply by n.z: the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
	if (Mnew != cuReal3()) return M.h.dim() * n.z * (energy_new - energy_old);
	else return M.h.dim() * n.z * energy_old;
}

//ZeemanCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_ZeemanCUDA(int spin_index, cuReal3 Mnew, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuReal3 Hext = cuReal3();

	if (pHavec && pHavec->linear_size()) {

		Hext = (*pHavec)[spin_index];
	}
	else {

		cuBReal cHA = *pcHA;
		update_parameters_mcoarse(spin_index, *pcHA, cHA);

		Hext = cHA * Ha;
	}

	if (Mnew != cuReal3()) return -M.h.dim() * (Mnew - M[spin_index]) * (cuBReal)MU0 * Hext;
	else return -M.h.dim() * M[spin_index] * (cuBReal)MU0 * Hext;
}

//MOpticalCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_MOpticalCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal cHmo = *pcHmo;
	update_parameters_mcoarse(spin_index, *pcHmo, cHmo);

	if (Mnew != cuReal3()) return -M.h.dim() * (Mnew - M[spin_index]) * (cuBReal)MU0 * cuReal3(0, 0, cHmo);
	else return -M.h.dim() * M[spin_index] * (cuBReal)MU0 * cuReal3(0, 0, cHmo);
}

//AnisotropyCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_AnisotropyCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal Ms = *pMs;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	update_parameters_mcoarse(spin_index, *pMs, Ms, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1);

	cuBReal dotprod = (M[spin_index] * mcanis_ea1) / Ms;
	cuBReal dpsq = dotprod * dotprod;

	if (Mnew != cuReal3()) {

		//calculate m.ea dot product
		cuBReal dotprod_new = Mnew * mcanis_ea1 / Ms;
		cuBReal dpsq_new = dotprod_new * dotprod_new;

		return M.h.dim() * ((K1 + K2 * (1 - dpsq_new)) * (1 - dpsq_new) - (K1 + K2 * (1 - dpsq)) * (1 - dpsq));
	}
	else return M.h.dim() * (K1 + K2 * (1 - dpsq)) * (1 - dpsq);
}

//AnisotropyCubiCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_AnisotropyCubiCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal Ms = *pMs;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pMs, Ms, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	cuReal3 S = M[spin_index] / Ms;
	cuReal3 S_new = Mnew / Ms;

	//calculate m.ea1, m.ea2 and m.ea3 dot products
	cuBReal d1 = S * mcanis_ea1;
	cuBReal d2 = S * mcanis_ea2;
	cuBReal d3 = S * mcanis_ea3;
	cuBReal d123 = d1 * d2 * d3;

	if (Mnew != cuReal3()) {

		cuBReal d1_new = S_new * mcanis_ea1;
		cuBReal d2_new = S_new * mcanis_ea2;
		cuBReal d3_new = S_new * mcanis_ea3;
		cuBReal d123_new = d1_new * d2_new * d3_new;

		return M.h.dim() * (
			(K1 * (d1_new*d1_new*d2_new*d2_new + d1_new*d1_new*d3_new*d3_new + d2_new*d2_new*d3_new*d3_new) + K2 * d123_new*d123_new) 
			- (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123));
	}
	else return M.h.dim() * (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123);
}

//AnisotropyBiaxialCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_AnisotropyBiaxialCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal Ms = *pMs;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pMs, Ms, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1 = M[spin_index] * mcanis_ea1 / Ms;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1 = M[spin_index] * mcanis_ea2 / Ms;
	cuBReal b2 = M[spin_index] * mcanis_ea3 / Ms;

	if (Mnew != cuReal3()) {

		//calculate m.ea1 dot product (uniaxial contribution)
		cuBReal u1new = Mnew * mcanis_ea1 / Ms;

		//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
		cuBReal b1new = Mnew * mcanis_ea2 / Ms;
		cuBReal b2new = Mnew * mcanis_ea3 / Ms;

		return M.h.dim() * ((K1 * (1 - u1new*u1new) + K2 * b1new*b1new*b2new*b2new) - (K1 * (1 - u1 * u1) + K2 * b1*b1*b2*b2));
	}
	else return M.h.dim() * (K1 * (1 - u1*u1) + K2 * b1*b1*b2*b2);
}

//AnisotropyTensorialCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_AnisotropyTensorialCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal K3 = *pK3;
	cuBReal Ms = *pMs;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pMs, Ms, *pK1, K1, *pK2, K2, *pK3, K3, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	auto Get_Energy = [&](cuBReal a, cuBReal b, cuBReal c) -> cuBReal {

		cuBReal energy_ = 0.0;

		for (int tidx = 0; tidx < pKt->linear_size(); tidx++) {

			cuBReal coeff;
			int order = (*pKt)[tidx].j + (*pKt)[tidx].k + (*pKt)[tidx].l;
			if (order == 2) coeff = K1 * (*pKt)[tidx].i;
			else if (order == 4) coeff = K2 * (*pKt)[tidx].i;
			else if (order == 6) coeff = K3 * (*pKt)[tidx].i;
			else coeff = (*pKt)[tidx].i;

			energy_ += coeff * pow(a, (*pKt)[tidx].j)*pow(b, (*pKt)[tidx].k)*pow(c, (*pKt)[tidx].l);
		}

		return energy_;
	};
	
	//calculate dot products
	cuBReal a = M[spin_index] * mcanis_ea1 / Ms;
	cuBReal b = M[spin_index] * mcanis_ea2 / Ms;
	cuBReal c = M[spin_index] * mcanis_ea3 / Ms;

	cuBReal energy_ = Get_Energy(a, b, c);

	if (Mnew != cuReal3()) {

		cuBReal anew = Mnew * mcanis_ea1 / Ms;
		cuBReal bnew = Mnew * mcanis_ea2 / Ms;
		cuBReal cnew = Mnew * mcanis_ea3 / Ms;

		cuBReal energynew_ = Get_Energy(anew, bnew, cnew);

		return M.h.dim() * (energynew_ - energy_);
	}
	else return M.h.dim() * energy_;
}

//RoughnessCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_RoughnessCUDA(int spin_index, cuReal3 Mnew)
{
	if (pFmul_rough && pFmul_rough->linear_size()) {

		cuVEC_VC<cuReal3>& M = *pM;

		cuReal33 Fmat = cuReal33(
			cuReal3((*pFmul_rough)[spin_index].x, (*pFomul_rough)[spin_index].x, (*pFomul_rough)[spin_index].y),
			cuReal3((*pFomul_rough)[spin_index].x, (*pFmul_rough)[spin_index].y, (*pFomul_rough)[spin_index].z),
			cuReal3((*pFomul_rough)[spin_index].y, (*pFomul_rough)[spin_index].z, (*pFmul_rough)[spin_index].z));

		cuReal3 Hrough = Fmat * M[spin_index];
		
		if (Mnew != cuReal3()) {

			cuReal3 Hrough_new = Fmat * Mnew;

			return -M.h.dim() * (cuBReal)MU0 * (Hrough_new * Mnew - Hrough * M[spin_index]);
		}
		else return -M.h.dim() * (cuBReal)MU0 * Hrough * M[spin_index];
	}
	else return 0.0;
}

//MElasticCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_MElasticCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& strain_diag = *pstrain_diag;
	cuVEC_VC<cuReal3>& strain_odiag = *pstrain_odiag;

	cuBReal Ms = *pMs;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	cuReal2 MEc = *pMEc;
	update_parameters_mcoarse(spin_index, *pMs, Ms, *pMEc, MEc, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	cuReal3 position = M.cellidx_to_position(spin_index);
	//xx, yy, zz
	cuReal3 Sd = strain_diag[position];
	//yz, xz, xy
	cuReal3 Sod = strain_odiag[position];

	//normalised magnetization
	//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

	Sd = cuReal3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
	Sod = cuReal3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

	auto Get_Energy = [&](cuReal3 M) -> cuBReal
	{
		cuReal3 m = cuReal3(M * mcanis_ea1, M * mcanis_ea2, M * mcanis_ea3) / Ms;

		cuReal3 Hmel_1 = (-2.0 * MEc.i / ((cuBReal)MU0 * Ms)) * cuReal3(
			m.x*Sd.x*mcanis_ea1.x + m.y*Sd.y*mcanis_ea2.x + m.z*Sd.z*mcanis_ea3.x,
			m.x*Sd.x*mcanis_ea1.y + m.y*Sd.y*mcanis_ea2.y + m.z*Sd.z*mcanis_ea3.y,
			m.x*Sd.x*mcanis_ea1.z + m.y*Sd.y*mcanis_ea2.z + m.z*Sd.z*mcanis_ea3.z);

		cuReal3 Hmel_2 = (-2.0 * MEc.j / ((cuBReal)MU0 * Ms)) * cuReal3(
			Sod.z * (mcanis_ea1.x*m.y + mcanis_ea2.x*m.x) + Sod.y * (mcanis_ea1.x*m.z + mcanis_ea3.x*m.x) + Sod.x * (mcanis_ea2.x*m.z + mcanis_ea3.x*m.y),
			Sod.z * (mcanis_ea1.y*m.y + mcanis_ea2.y*m.x) + Sod.y * (mcanis_ea1.y*m.z + mcanis_ea3.y*m.x) + Sod.x * (mcanis_ea2.y*m.z + mcanis_ea3.y*m.y),
			Sod.z * (mcanis_ea1.z*m.y + mcanis_ea2.z*m.x) + Sod.y * (mcanis_ea1.z*m.z + mcanis_ea3.z*m.x) + Sod.x * (mcanis_ea2.z*m.z + mcanis_ea3.z*m.y));

		return -(cuBReal)MU0 * M * (Hmel_1 + Hmel_2) / 2;
	};

	if (Mnew != cuReal3()) return M.h.dim() * (Get_Energy(Mnew) - Get_Energy(M[spin_index]));
	else return M.h.dim() * Get_Energy(M[spin_index]);
}

#endif