#include "ManagedMeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"
#include "MeshParamsControlCUDA.h"

#include "ModulesDefs.h"

//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

//Energy Deltas

//Antiferromagnetic

//switch function which adds all assigned energy contributions in this mesh to calculate energy change from current spin to Mnew spin : return energy change as new - old
//If Mnew is passed in as cuReal3(), then this function returns the current spin energy only - all functions in the switch statement below implement this eventuality.
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B, int*& cuModules, int numModules, cuReal3& Ha)
{
	cuReal2 energy = cuReal2();

	for (int idx = 0; idx < numModules; idx++) {

		switch (cuModules[idx]) {

		case MOD_DEMAG_N:
			energy += Get_EnergyChange_AFM_DemagNCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_DEMAG:
			energy += Get_EnergyChange_AFM_DemagCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_SDEMAG_DEMAG:
			//same method as for MOD_DEMAG, but the effective field pointer now points to SDemag_DemagCUDA module effective field
			energy += Get_EnergyChange_AFM_DemagCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_EXCHANGE:
			energy += Get_EnergyChange_AFM_ExchangeCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_DMEXCHANGE:
			energy += Get_EnergyChange_AFM_DMExchangeCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_EnergyChange_AFM_iDMExchangeCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_VIDMEXCHANGE:
			energy += Get_EnergyChange_AFM_viDMExchangeCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_SURFEXCHANGE:
			energy += Get_EnergyChange_AFM_SurfExchangeCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_ZEEMAN:
			energy += Get_EnergyChange_AFM_ZeemanCUDA(spin_index, Mnew_A, Mnew_B, Ha);
			break;

		case MOD_MOPTICAL:
			energy += Get_EnergyChange_AFM_MOpticalCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_ANIUNI:
			energy += Get_EnergyChange_AFM_AnisotropyCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_ANICUBI:
			energy += Get_EnergyChange_AFM_AnisotropyCubiCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_ANIBI:
			energy += Get_EnergyChange_AFM_AnisotropyBiaxialCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_ANITENS:
			energy += Get_EnergyChange_AFM_AnisotropyTensorialCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_ROUGHNESS:
			energy += Get_EnergyChange_AFM_RoughnessCUDA(spin_index, Mnew_A, Mnew_B);
			break;

		case MOD_MELASTIC:
			//Not available for AFM
			break;

		default:
			break;
		}
	}

	return energy;
}

//Demag_N
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_DemagNCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& Mvec = *pM;
	cuVEC_VC<cuReal3>& Mvec2 = *pM2;

	//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
	cuReal2 Nxy = *pNxy;
	cuBReal Nz = (1 - Nxy.x - Nxy.y);

	cuReal3 M = (Mvec[spin_index] + Mvec2[spin_index]) / 2;
	cuReal3 Mnew = (Mnew_A + Mnew_B) / 2;

	cuBReal energy_ = 0.0;

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		energy_ = ((cuBReal)MU0 / 2) * Mvec.h.dim() * (
			(Mnew * cuReal3(Nxy.x * Mnew.x, Nxy.y * Mnew.y, Nz * Mnew.z)) -
			(M * cuReal3(Nxy.x * M.x, Nxy.y * M.y, Nz * M.z)));
	}
	else energy_ = ((cuBReal)MU0 / 2) * Mvec.h.dim() * (M * cuReal3(Nxy.x * M.x, Nxy.y * M.y, Nz * M.z));

	return cuReal2(energy_, energy_);
}

//Demag
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_DemagCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	//Module_Heff needs to be calculated (done during a Monte Carlo simulation, where this method would be used)
	if (pDemag_Heff && pDemag_Heff->linear_size()) {

		cuVEC_VC<cuReal3>& Mvec = *pM;
		cuVEC_VC<cuReal3>& Mvec2 = *pM2;

		cuReal3 M = (Mvec[spin_index] + Mvec2[spin_index]) / 2;
		cuReal3 Mnew = (Mnew_A + Mnew_B) / 2;

		cuBReal energy_ = 0.0;

		//do not divide by 2 as we are not double-counting here
		if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

			energy_ = -(cuBReal)MU0 * Mvec.h.dim() * (*pDemag_Heff)[Mvec.cellidx_to_position(spin_index)] * (Mnew - M);
		}
		else {

			energy_ = -(cuBReal)MU0 * Mvec.h.dim() * (*pDemag_Heff)[Mvec.cellidx_to_position(spin_index)] * M;
		}

		return cuReal2(energy_, energy_);
	}
	else return cuReal2();
}

//Exch_6ngbr_Neu
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_ExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 A_AFM = *pA_AFM;
	cuReal2 Ah = *pAh;
	cuReal2 Anh = *pAnh;
	update_parameters_mcoarse(spin_index, *pA_AFM, A_AFM, *pMs_AFM, Ms_AFM, *pAh, Ah, *pAnh, Anh);

	auto Get_Energy = [&](void) -> cuReal2 {

		cuReal3 delsq_M_A = M.delsq_neu(spin_index);
		cuReal3 delsq_M_B = M2.delsq_neu(spin_index);

		cuReal3 Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
		cuReal3 Hexch2 = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

		return cuReal2(M[spin_index] * Hexch, M2[spin_index] * Hexch2);
	};

	cuReal2 energy_ = Get_Energy();

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew_A.renormalize(M[spin_index].norm());
		Mnew_B.renormalize(M2[spin_index].norm());

		cuReal3 Mold_A = M[spin_index];
		cuReal3 Mold_B = M2[spin_index];

		M[spin_index] = Mnew_A;
		M2[spin_index] = Mnew_B;

		cuReal2 energynew_ = Get_Energy();
		
		M[spin_index] = Mold_A;
		M2[spin_index] = Mold_B;

		//do not divide by 2 as we are not double-counting here
		return -(cuBReal)MU0 * M.h.dim() * (energynew_ - energy_);
	}
	//If Mnew is null then this method is used to obtain current energy only, not energy change
	else return -(cuBReal)MU0 * M.h.dim() * energy_;
}

//DMExchangeCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_DMExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 A_AFM = *pA_AFM;
	cuReal2 Ah = *pAh;
	cuReal2 Anh = *pAnh;
	cuReal2 D_AFM = *pD_AFM;
	cuBReal Dh = *pDh;
	cuReal3 dh_dir = *pdh_dir;
	update_parameters_mcoarse(spin_index, *pA_AFM, A_AFM, *pAh, Ah, *pAnh, Anh, *pD_AFM, D_AFM, *pMs_AFM, Ms_AFM, *pDh, Dh, *pdh_dir, dh_dir);

	cuReal2 Aconst = 2 * A_AFM / (MU0 * (Ms_AFM & Ms_AFM));
	cuReal2 Dconst = -2 * D_AFM / (MU0 * (Ms_AFM & Ms_AFM));
	cuBReal Dhconst = (Dh / (MU0*Ms_AFM.i*Ms_AFM.j));

	auto Get_Energy = [&](void) -> cuReal2
	{
		cuReal3 Hexch_A, Hexch_A2, Hexch_D, Hexch_D2;

		if (M.is_plane_interior(spin_index)) {

			//interior point : can use cheaper neu versions

			//1. direct exchange contribution + AFM contribution
			cuReal3 delsq_M_A = M.delsq_neu(spin_index);
			cuReal3 delsq_M_B = M2.delsq_neu(spin_index);

			Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			//2. Dzyaloshinskii-Moriya exchange contribution

			//Hdm, ex = -2D / (mu0*Ms) * curl m
			Hexch_D = Dconst.i * M.curl_neu(spin_index);
			Hexch_D2 = Dconst.j * M2.curl_neu(spin_index);

			//3. Homogeneous DMI contribution
			Hexch_D += Dhconst * (dh_dir ^ M2[spin_index]);
			Hexch_D2 += -Dhconst * (dh_dir ^ M[spin_index]);
		}
		else {

			//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
			cuReal3 bnd_dm_dx = (D_AFM.i / (2 * A_AFM.i)) * cuReal3(0, -M[spin_index].z, M[spin_index].y);
			cuReal3 bnd_dm_dy = (D_AFM.i / (2 * A_AFM.i)) * cuReal3(M[spin_index].z, 0, -M[spin_index].x);
			cuReal3 bnd_dm_dz = (D_AFM.i / (2 * A_AFM.i)) * cuReal3(-M[spin_index].y, M[spin_index].x, 0);
			cuReal33 bndA_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, bnd_dm_dz);

			bnd_dm_dx = (D_AFM.j / (2 * A_AFM.j)) * cuReal3(0, -M2[spin_index].z, M2[spin_index].y);
			bnd_dm_dy = (D_AFM.j / (2 * A_AFM.j)) * cuReal3(M2[spin_index].z, 0, -M2[spin_index].x);
			bnd_dm_dz = (D_AFM.j / (2 * A_AFM.j)) * cuReal3(-M2[spin_index].y, M2[spin_index].x, 0);
			cuReal33 bndB_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, bnd_dm_dz);

			cuReal3 delsq_M_A = M.delsq_nneu(spin_index, bndA_nneu);
			cuReal3 delsq_M_B = M2.delsq_nneu(spin_index, bndB_nneu);

			//1. direct exchange contribution + AFM contribution

			//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
			Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			//2. Dzyaloshinskii-Moriya exchange contribution

			//Hdm, ex = -2D / (mu0*Ms) * curl m
			//For cmbnd cells curl_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
			Hexch_D = Dconst.i * M.curl_nneu(spin_index, bndA_nneu);
			Hexch_D2 = Dconst.j * M2.curl_nneu(spin_index, bndB_nneu);

			//3. Homogeneous DMI contribution
			Hexch_D += Dhconst * (dh_dir ^ M2[spin_index]);
			Hexch_D2 += -Dhconst * (dh_dir ^ M[spin_index]);
		}

		return cuReal2(M[spin_index] * (Hexch_A + Hexch_D), M2[spin_index] * (Hexch_A2 + Hexch_D2));
	};

	cuReal2 energy_ = Get_Energy();

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew_A.renormalize(M[spin_index].norm());
		Mnew_B.renormalize(M2[spin_index].norm());

		cuReal3 Mold_A = M[spin_index];
		cuReal3 Mold_B = M2[spin_index];

		M[spin_index] = Mnew_A;
		M2[spin_index] = Mnew_B;

		cuReal2 energynew_ = Get_Energy();

		M[spin_index] = Mold_A;
		M2[spin_index] = Mold_B;

		//do not divide by 2 as we are not double-counting here
		return -MU0 * M.h.dim() * (energynew_ - energy_);
	}
	//If Mnew is null then this method is used to obtain current energy only, not energy change
	else return -MU0 * M.h.dim() * energy_;
}

//iDMExchangeCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_iDMExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 A_AFM = *pA_AFM;
	cuReal2 Ah = *pAh;
	cuReal2 Anh = *pAnh;
	cuReal2 D_AFM = *pD_AFM;
	update_parameters_mcoarse(spin_index, *pA_AFM, A_AFM, *pAh, Ah, *pAnh, Anh, *pD_AFM, D_AFM, *pMs_AFM, Ms_AFM);

	cuReal2 Aconst = 2 * A_AFM / (MU0 * (Ms_AFM & Ms_AFM));
	cuReal2 Dconst = -2 * D_AFM / (MU0 * (Ms_AFM & Ms_AFM));

	auto Get_Energy = [&](void) -> cuReal2
	{
		cuReal3 Hexch_A, Hexch_A2, Hexch_D, Hexch_D2;

		if (M.is_plane_interior(spin_index)) {

			//interior point : can use cheaper neu versions

			//1. direct exchange contribution + AFM contribution
			cuReal3 delsq_M_A = M.delsq_neu(spin_index);
			cuReal3 delsq_M_B = M2.delsq_neu(spin_index);

			Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			//2. Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			cuReal33 Mdiff_A = M.grad_neu(spin_index);
			cuReal33 Mdiff_B = M2.grad_neu(spin_index);

			//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
			Hexch_D = Dconst.i * cuReal3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);
			Hexch_D2 = Dconst.j * cuReal3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);
		}
		else {

			cuReal33 bndA_nneu, bndB_nneu;

			cuReal2 nhconst = Anh / (2 * A_AFM);

			if (fabs(nhconst.i) != 1.0) {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				cuReal3 bnd_dm_dx = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(M[spin_index].z - nhconst.i * M2[spin_index].z, 0, -M[spin_index].x + nhconst.i * M2[spin_index].x);
				cuReal3 bnd_dm_dy = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(0, M[spin_index].z - nhconst.i * M2[spin_index].z, -M[spin_index].y + nhconst.i * M2[spin_index].y);

				bndA_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, cuReal3());
			}
			else {

				cuReal3 bnd_dm_dx = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(M[spin_index].z, 0, -M[spin_index].x);
				cuReal3 bnd_dm_dy = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(0, M[spin_index].z, -M[spin_index].y);

				bndA_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, cuReal3());
			}

			if (fabs(nhconst.j) != 1.0) {

				cuReal3 bnd_dm_dx = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(M2[spin_index].z - nhconst.j * M[spin_index].z, 0, -M2[spin_index].x + nhconst.j * M[spin_index].x);
				cuReal3 bnd_dm_dy = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(0, M2[spin_index].z - nhconst.j * M[spin_index].z, -M2[spin_index].y + nhconst.j * M[spin_index].y);

				bndB_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, cuReal3());
			}
			else {

				cuReal3 bnd_dm_dx = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(M2[spin_index].z, 0, -M2[spin_index].x);
				cuReal3 bnd_dm_dy = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(0, M2[spin_index].z, -M2[spin_index].y);

				bndB_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, cuReal3());
			}

			cuReal3 delsq_M_A = M.delsq_nneu(spin_index, bndA_nneu);
			cuReal3 delsq_M_B = M2.delsq_nneu(spin_index, bndB_nneu);

			//1. direct exchange contribution + AFM contribution

			//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
			Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			//2. Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
			cuReal33 Mdiff_A = M.grad_nneu(spin_index, bndA_nneu);
			cuReal33 Mdiff_B = M2.grad_nneu(spin_index, bndB_nneu);

			//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
			Hexch_D = Dconst.i * cuReal3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);
			Hexch_D2 = Dconst.j * cuReal3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);
		}

		return cuReal2(M[spin_index] * (Hexch_A + Hexch_D), M2[spin_index] * (Hexch_A2 + Hexch_D2));
	};

	cuReal2 energy_ = Get_Energy();

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew_A.renormalize(M[spin_index].norm());
		Mnew_B.renormalize(M2[spin_index].norm());

		cuReal3 Mold_A = M[spin_index];
		cuReal3 Mold_B = M2[spin_index];

		M[spin_index] = Mnew_A;
		M2[spin_index] = Mnew_B;

		cuReal2 energynew_ = Get_Energy();

		M[spin_index] = Mold_A;
		M2[spin_index] = Mold_B;

		//do not divide by 2 as we are not double-counting here
		return -MU0 * M.h.dim() * (energynew_ - energy_);
	}
	//If Mnew is null then this method is used to obtain current energy only, not energy change
	else return -MU0 * M.h.dim() * energy_;
}

//viDMExchangeCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_viDMExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 A_AFM = *pA_AFM;
	cuReal2 Ah = *pAh;
	cuReal2 Anh = *pAnh;
	cuReal2 D_AFM = *pD_AFM;
	cuReal3 D_dir = *pD_dir;
	update_parameters_mcoarse(spin_index, *pA_AFM, A_AFM, *pAh, Ah, *pAnh, Anh, *pD_AFM, D_AFM, *pMs_AFM, Ms_AFM, *pD_dir, D_dir);

	cuReal2 Aconst = 2 * A_AFM / (MU0 * (Ms_AFM & Ms_AFM));
	cuReal2 Dconst = -2 * D_AFM / (MU0 * (Ms_AFM & Ms_AFM));

	auto Get_Energy = [&](void) -> cuReal2
	{
		cuReal3 Hexch_A, Hexch_A2, Hexch_D, Hexch_D2;

		if (M.is_plane_interior(spin_index)) {

			//interior point : can use cheaper neu versions

			//1. direct exchange contribution + AFM contribution
			cuReal3 delsq_M_A = M.delsq_neu(spin_index);
			cuReal3 delsq_M_B = M2.delsq_neu(spin_index);

			Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			//2. Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			cuReal33 Mdiff_A = M.grad_neu(spin_index);
			cuReal33 Mdiff_B = M2.grad_neu(spin_index);

			cuReal3 hexch_D_A_x = cuReal3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
			cuReal3 hexch_D_A_y = cuReal3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
			cuReal3 hexch_D_A_z = cuReal3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

			Hexch_D = Dconst.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z);

			cuReal3 hexch_D_B_x = cuReal3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
			cuReal3 hexch_D_B_y = cuReal3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
			cuReal3 hexch_D_B_z = cuReal3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

			Hexch_D2 = Dconst.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z);
		}
		else {

			cuReal33 bndA_nneu, bndB_nneu;

			cuReal2 nhconst = Anh / (2 * A_AFM);

			if (fabs(nhconst.i) != 1.0) {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				cuReal3 bnd_dm_dy_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(-M[spin_index].y + nhconst.i * M2[spin_index].y, M[spin_index].x - nhconst.i * M2[spin_index].x, 0);
				cuReal3 bnd_dm_dz_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(-M[spin_index].z + nhconst.i * M2[spin_index].z, 0, M[spin_index].x - nhconst.i * M2[spin_index].x);
				cuReal33 bndA_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

				cuReal3 bnd_dm_dx_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(M[spin_index].y - nhconst.i * M2[spin_index].y, -M[spin_index].x + nhconst.i * M2[spin_index].x, 0);
				cuReal3 bnd_dm_dz_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(0, -M[spin_index].z + nhconst.i * M2[spin_index].z, M[spin_index].y - nhconst.i * M2[spin_index].y);
				cuReal33 bndA_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

				cuReal3 bnd_dm_dx_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(M[spin_index].z - nhconst.i * M2[spin_index].z, 0, -M[spin_index].x + nhconst.i * M2[spin_index].x);
				cuReal3 bnd_dm_dy_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * cuReal3(0, M[spin_index].z - nhconst.i * M2[spin_index].z, -M[spin_index].y + nhconst.i * M2[spin_index].y);
				cuReal33 bndA_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

				bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
			}
			else {

				cuReal3 bnd_dm_dy_x = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(-M[spin_index].y, M[spin_index].x, 0);
				cuReal3 bnd_dm_dz_x = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(-M[spin_index].z, 0, M[spin_index].x);
				cuReal33 bndA_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

				cuReal3 bnd_dm_dx_y = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(M[spin_index].y, -M[spin_index].x, 0);
				cuReal3 bnd_dm_dz_y = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(0, -M[spin_index].z, M[spin_index].y);
				cuReal33 bndA_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

				cuReal3 bnd_dm_dx_z = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(M[spin_index].z, 0, -M[spin_index].x);
				cuReal3 bnd_dm_dy_z = (D_AFM.i / (4 * A_AFM.i)) * cuReal3(0, M[spin_index].z, -M[spin_index].y);
				cuReal33 bndA_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

				bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
			}

			if (fabs(nhconst.j) != 1.0) {

				cuReal3 bnd_dm_dy_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(-M2[spin_index].y + nhconst.j * M[spin_index].y, M2[spin_index].x - nhconst.j * M[spin_index].x, 0);
				cuReal3 bnd_dm_dz_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(-M2[spin_index].z + nhconst.j * M[spin_index].z, 0, M2[spin_index].x - nhconst.j * M[spin_index].x);
				cuReal33 bndB_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

				cuReal3 bnd_dm_dx_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(M2[spin_index].y - nhconst.j * M[spin_index].y, -M2[spin_index].x + nhconst.j * M[spin_index].x, 0);
				cuReal3 bnd_dm_dz_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(0, -M2[spin_index].z + nhconst.j * M[spin_index].z, M2[spin_index].y - nhconst.j * M[spin_index].y);
				cuReal33 bndB_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

				cuReal3 bnd_dm_dx_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(M2[spin_index].z - nhconst.j * M[spin_index].z, 0, -M2[spin_index].x + nhconst.j * M[spin_index].x);
				cuReal3 bnd_dm_dy_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * cuReal3(0, M2[spin_index].z - nhconst.j * M[spin_index].z, -M2[spin_index].y + nhconst.j * M[spin_index].y);
				cuReal33 bndB_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

				bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
			}
			else {

				cuReal3 bnd_dm_dy_x = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(-M2[spin_index].y, M2[spin_index].x, 0);
				cuReal3 bnd_dm_dz_x = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(-M2[spin_index].z, 0, M2[spin_index].x);
				cuReal33 bndB_nneu_x = cuReal33(cuReal3(), bnd_dm_dy_x, bnd_dm_dz_x);

				cuReal3 bnd_dm_dx_y = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(M2[spin_index].y, -M2[spin_index].x, 0);
				cuReal3 bnd_dm_dz_y = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(0, -M2[spin_index].z, M2[spin_index].y);
				cuReal33 bndB_nneu_y = cuReal33(bnd_dm_dx_y, cuReal3(), bnd_dm_dz_y);

				cuReal3 bnd_dm_dx_z = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(M2[spin_index].z, 0, -M2[spin_index].x);
				cuReal3 bnd_dm_dy_z = (D_AFM.j / (4 * A_AFM.j)) * cuReal3(0, M2[spin_index].z, -M2[spin_index].y);
				cuReal33 bndB_nneu_z = cuReal33(bnd_dm_dx_z, bnd_dm_dy_z, cuReal3());

				bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
			}

			cuReal3 delsq_M_A = M.delsq_nneu(spin_index, bndA_nneu);
			cuReal3 delsq_M_B = M2.delsq_nneu(spin_index, bndB_nneu);

			//1. direct exchange contribution + AFM contribution

			//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
			Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			//2. Dzyaloshinskii-Moriya interfacial exchange contribution

			//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
			//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
			cuReal33 Mdiff_A = M.grad_nneu(spin_index, bndA_nneu);
			cuReal33 Mdiff_B = M2.grad_nneu(spin_index, bndB_nneu);

			cuReal3 hexch_D_A_x = cuReal3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
			cuReal3 hexch_D_A_y = cuReal3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
			cuReal3 hexch_D_A_z = cuReal3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

			Hexch_D = Dconst.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z);

			cuReal3 hexch_D_B_x = cuReal3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
			cuReal3 hexch_D_B_y = cuReal3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
			cuReal3 hexch_D_B_z = cuReal3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

			Hexch_D2 = Dconst.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z);
		}

		return cuReal2(M[spin_index] * (Hexch_A + Hexch_D), M2[spin_index] * (Hexch_A2 + Hexch_D2));
	};

	cuReal2 energy_ = Get_Energy();

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
		Mnew_A.renormalize(M[spin_index].norm());
		Mnew_B.renormalize(M2[spin_index].norm());

		cuReal3 Mold_A = M[spin_index];
		cuReal3 Mold_B = M2[spin_index];

		M[spin_index] = Mnew_A;
		M2[spin_index] = Mnew_B;

		cuReal2 energynew_ = Get_Energy();

		M[spin_index] = Mold_A;
		M2[spin_index] = Mold_B;

		//do not divide by 2 as we are not double-counting here
		return -MU0 * M.h.dim() * (energynew_ - energy_);
	}
	//If Mnew is null then this method is used to obtain current energy only, not energy change
	else return -MU0 * M.h.dim() * energy_;
}

//SurfExchangeCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_SurfExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuReal2 energy_new = cuReal2(), energy_old = cuReal2();

	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	//thickness of layer - SurfExchange applies for layers in the xy plane
	cuBReal thickness = M.rect.e.z - M.rect.s.z;

	//if spin is on top surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == n.z - 1 && (pMeshFM_Top_size + pMeshAFM_Top_size)) {

		if (!M.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			cuReal2 Ms_AFM = *pMs_AFM;
			update_parameters_mcoarse(spin_index, *pMs_AFM, Ms_AFM);
			
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

				//Surface exchange field from a ferromagnetic mesh

				//Top mesh sets J1 and J2 values
				cuBReal J1 = *(pMeshFM_Top[mesh_idx].pJ1);
				cuBReal J2 = *(pMeshFM_Top[mesh_idx].pJ2);
				pMeshFM_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMeshFM_Top[mesh_idx].pJ1), J1, *(pMeshFM_Top[mesh_idx].pJ2), J2);
				
				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Top[cell_rel_pos].normalized();
				cuReal3 m_i1 = M[spin_index] / Ms_AFM.i;
				cuReal3 m_i2 = M2[spin_index] / Ms_AFM.j;

				energy_old.i = (-J1 * (m_i1 * m_j)) / thickness;
				energy_old.j = (-J2 * (m_i2 * m_j)) / thickness;

				if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

					cuReal3 mnew_i1 = Mnew_A / Ms_AFM.i;
					cuReal3 mnew_i2 = Mnew_B / Ms_AFM.j;

					energy_new.i = (-J1 * (mnew_i1 * m_j)) / thickness;
					energy_new.j = (-J2 * (mnew_i2 * m_j)) / thickness;
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
				cuReal3 m_i1 = M[spin_index] / Ms_AFM.i;
				cuReal3 m_i2 = M2[spin_index] / Ms_AFM.j;

				energy_old.i = (-J1 * (m_i1 * m_j1)) / thickness;
				energy_old.j = (-J2 * (m_i2 * m_j2)) / thickness;

				if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

					cuReal3 mnew_i1 = Mnew_A / Ms_AFM.i;
					cuReal3 mnew_i2 = Mnew_B / Ms_AFM.j;

					energy_new.i = (-J1 * (mnew_i1 * m_j1)) / thickness;
					energy_new.j = (-J2 * (mnew_i2 * m_j2)) / thickness;
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

			cuReal2 Ms_AFM = *pMs_AFM;
			cuBReal J1 = *pJ1;
			cuBReal J2 = *pJ2;
			update_parameters_mcoarse(spin_index, *pMs_AFM, Ms_AFM, *pJ1, J1, *pJ2, J2);

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

				//Surface exchange field from a ferromagnetic mesh

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_i1 = M[spin_index] / Ms_AFM.i;
				cuReal3 m_i2 = M2[spin_index] / Ms_AFM.j;

				energy_old.i += (-J1 * (m_i1 * m_j)) / thickness;
				energy_old.j += (-J2 * (m_i2 * m_j)) / thickness;

				if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

					cuReal3 mnew_i1 = Mnew_A / Ms_AFM.i;
					cuReal3 mnew_i2 = Mnew_B / Ms_AFM.j;

					energy_new.i += (-J1 * (mnew_i1 * m_j)) / thickness;
					energy_new.j += (-J2 * (mnew_i2 * m_j)) / thickness;
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
				cuReal3 m_i1 = M[spin_index] / Ms_AFM.i;
				cuReal3 m_i2 = M2[spin_index] / Ms_AFM.j;

				energy_old.i += (-J1 * (m_i1 * m_j1)) / thickness;
				energy_old.j += (-J2 * (m_i2 * m_j2)) / thickness;

				if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

					cuReal3 mnew_i1 = Mnew_A / Ms_AFM.i;
					cuReal3 mnew_i2 = Mnew_B / Ms_AFM.j;

					energy_new.i += (-J1 * (mnew_i1 * m_j1)) / thickness;
					energy_new.j += (-J2 * (mnew_i2 * m_j2)) / thickness;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	//multiply by n.z: the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) return M.h.dim() * n.z * (energy_new - energy_old);
	else return M.h.dim() * n.z * energy_old;
}

//ZeemanCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_ZeemanCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal3 Hext = cuReal3();

	if (pHavec && pHavec->linear_size()) {

		Hext = (*pHavec)[spin_index];
	}
	else {

		cuBReal cHA = *pcHA;
		update_parameters_mcoarse(spin_index, *pcHA, cHA);

		Hext = cHA * Ha;
	}

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) return -M.h.dim() * cuReal2((Mnew_A - M[spin_index]) * (cuBReal)MU0 * Hext, (Mnew_B - M2[spin_index]) * (cuBReal)MU0 * Hext);
	else return -M.h.dim() * cuReal2(M[spin_index] * (cuBReal)MU0 * Hext, M2[spin_index] * (cuBReal)MU0 * Hext);
}

//MOpticalCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_MOpticalCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuBReal cHmo = *pcHmo;
	update_parameters_mcoarse(spin_index, *pcHmo, cHmo);

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) return -MU0 * M.h.dim() * cuReal2((Mnew_A - M[spin_index]) * cuReal3(0, 0, cHmo), (Mnew_B - M2[spin_index]) * cuReal3(0, 0, cHmo));
	else return -MU0 * M.h.dim() * cuReal2(M[spin_index] * cuReal3(0, 0, cHmo), M2[spin_index] * cuReal3(0, 0, cHmo));
}

//AnisotropyCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_AnisotropyCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 K1_AFM = *pK1_AFM;
	cuReal2 K2_AFM = *pK2_AFM;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	update_parameters_mcoarse(spin_index, *pMs_AFM, Ms_AFM, *pK1_AFM, K1_AFM, *pK2_AFM, K2_AFM, *pmcanis_ea1, mcanis_ea1);

	//calculate m.ea dot product
	cuBReal dotprod = (M[spin_index] * mcanis_ea1) / Ms_AFM.i;
	cuBReal dotprod2 = (M2[spin_index] * mcanis_ea1) / Ms_AFM.j;

	cuBReal energyA = (K1_AFM.i + K2_AFM.i * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod);
	cuBReal energyB = (K1_AFM.j + K2_AFM.j * (1 - dotprod2 * dotprod2)) * (1 - dotprod2 * dotprod2);

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		dotprod = (Mnew_A * mcanis_ea1) / Ms_AFM.i;
		dotprod2 = (Mnew_B * mcanis_ea1) / Ms_AFM.j;

		cuBReal energynewA = (K1_AFM.i + K2_AFM.i * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod);
		cuBReal energynewB = (K1_AFM.j + K2_AFM.j * (1 - dotprod2 * dotprod2)) * (1 - dotprod2 * dotprod2);

		return M.h.dim() * cuReal2(energynewA - energyA, energynewB - energyB);
	}
	else return M.h.dim() * cuReal2(energyA, energyB);
}

//AnisotropyCubiCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_AnisotropyCubiCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 K1_AFM = *pK1_AFM;
	cuReal2 K2_AFM = *pK2_AFM;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pMs_AFM, Ms_AFM, *pK1_AFM, K1_AFM, *pK2_AFM, K2_AFM, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	//calculate m.ea1, m.ea2 and m.ea3 dot products
	cuBReal d1 = (M[spin_index] * mcanis_ea1) / Ms_AFM.i;
	cuBReal d2 = (M[spin_index] * mcanis_ea2) / Ms_AFM.i;
	cuBReal d3 = (M[spin_index] * mcanis_ea3) / Ms_AFM.i;
	cuBReal d123 = d1 * d2*d3;

	//same thing for sub-lattice B

	cuBReal d1B = (M2[spin_index] * mcanis_ea1) / Ms_AFM.j;
	cuBReal d2B = (M2[spin_index] * mcanis_ea2) / Ms_AFM.j;
	cuBReal d3B = (M2[spin_index] * mcanis_ea3) / Ms_AFM.j;
	cuBReal d123B = d1B * d2B*d3B;

	cuBReal energyA = K1_AFM.i * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3) + K2_AFM.i * d123*d123;
	cuBReal energyB = K1_AFM.j * (d1B*d1B*d2B*d2B + d1B * d1B*d3B*d3B + d2B * d2B*d3B*d3B) + K2_AFM.j * d123B*d123B;

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		//calculate m.ea1, m.ea2 and m.ea3 dot products
		cuBReal d1new = (M[spin_index] * mcanis_ea1) / Ms_AFM.i;
		cuBReal d2new = (M[spin_index] * mcanis_ea2) / Ms_AFM.i;
		cuBReal d3new = (M[spin_index] * mcanis_ea3) / Ms_AFM.i;
		cuBReal d123new = d1new * d2new*d3new;

		//same thing for sub-lattice B

		cuBReal d1Bnew = (M2[spin_index] * mcanis_ea1) / Ms_AFM.j;
		cuBReal d2Bnew = (M2[spin_index] * mcanis_ea2) / Ms_AFM.j;
		cuBReal d3Bnew = (M2[spin_index] * mcanis_ea3) / Ms_AFM.j;
		cuBReal d123Bnew = d1Bnew * d2Bnew*d3Bnew;

		cuBReal energyAnew = K1_AFM.i * (d1new*d1new*d2new*d2new + d1new * d1new*d3new*d3new + d2new * d2new*d3new*d3new) + K2_AFM.i * d123new*d123new;
		cuBReal energyBnew = K1_AFM.j * (d1Bnew*d1Bnew*d2Bnew*d2Bnew + d1Bnew * d1Bnew*d3Bnew*d3Bnew + d2Bnew * d2Bnew*d3Bnew*d3Bnew) + K2_AFM.j * d123Bnew*d123Bnew;

		return M.h.dim() * cuReal2(energyAnew - energyA, energyBnew - energyB);
	}
	else return M.h.dim() * cuReal2(energyA, energyB);
}

//AnisotropyBiaxialCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_AnisotropyBiaxialCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 K1_AFM = *pK1_AFM;
	cuReal2 K2_AFM = *pK2_AFM;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pMs_AFM, Ms_AFM, *pK1_AFM, K1_AFM, *pK2_AFM, K2_AFM, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1_A = (M[spin_index] * mcanis_ea1) / Ms_AFM.i;
	cuBReal u1_B = (M2[spin_index] * mcanis_ea1) / Ms_AFM.j;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1_A = (M[spin_index] * mcanis_ea2) / Ms_AFM.i;
	cuBReal b1_B = (M2[spin_index] * mcanis_ea2) / Ms_AFM.j;

	cuBReal b2_A = (M[spin_index] * mcanis_ea3) / Ms_AFM.i;
	cuBReal b2_B = (M2[spin_index] * mcanis_ea3) / Ms_AFM.j;

	cuBReal energyA = K1_AFM.i * (1 - u1_A * u1_A) + K2_AFM.i * b1_A*b1_A*b2_A*b2_A;
	cuBReal energyB = K1_AFM.j * (1 - u1_B * u1_B) + K2_AFM.j * b1_B*b1_B*b2_B*b2_B;

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		//calculate m.ea1 dot product (uniaxial contribution)
		cuBReal u1_Anew = (Mnew_A * mcanis_ea1) / Ms_AFM.i;
		cuBReal u1_Bnew = (Mnew_B * mcanis_ea1) / Ms_AFM.j;

		//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
		cuBReal b1_Anew = (Mnew_A * mcanis_ea2) / Ms_AFM.i;
		cuBReal b1_Bnew = (Mnew_B * mcanis_ea2) / Ms_AFM.j;

		cuBReal b2_Anew = (Mnew_A * mcanis_ea3) / Ms_AFM.i;
		cuBReal b2_Bnew = (Mnew_B * mcanis_ea3) / Ms_AFM.j;

		cuBReal energyAnew = K1_AFM.i * (1 - u1_Anew * u1_Anew) + K2_AFM.i * b1_Anew*b1_Anew*b2_Anew*b2_Anew;
		cuBReal energyBnew = K1_AFM.j * (1 - u1_Bnew * u1_Bnew) + K2_AFM.j * b1_Bnew*b1_Bnew*b2_Bnew*b2_Bnew;

		return M.h.dim() * cuReal2(energyAnew - energyA, energyBnew - energyB);
	}
	else return M.h.dim() * cuReal2(energyA, energyB);
}

//AnisotropyTensorialCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_AnisotropyTensorialCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	cuVEC_VC<cuReal3>& M = *pM;
	cuVEC_VC<cuReal3>& M2 = *pM2;

	cuReal2 Ms_AFM = *pMs_AFM;
	cuReal2 K1_AFM = *pK1_AFM;
	cuReal2 K2_AFM = *pK2_AFM;
	cuReal2 K3_AFM = *pK3_AFM;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pMs_AFM, Ms_AFM, *pK1_AFM, K1_AFM, *pK2_AFM, K2_AFM, *pK3_AFM, K3_AFM, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	auto Get_Energy = [&](cuBReal a, cuBReal b, cuBReal c, cuBReal a2, cuBReal b2, cuBReal c2) -> cuReal2 {

		cuBReal energyA = 0.0, energyB = 0.0;

		for (int tidx = 0; tidx < pKt->linear_size(); tidx++) {

			cuBReal coeff;
			int order = (*pKt)[tidx].j + (*pKt)[tidx].k + (*pKt)[tidx].l;
			if (order == 2) coeff = K1_AFM.i * (*pKt)[tidx].i;
			else if (order == 4) coeff = K2_AFM.i * (*pKt)[tidx].i;
			else if (order == 6) coeff = K3_AFM.i * (*pKt)[tidx].i;
			else coeff = (*pKt)[tidx].i;

			energyA += coeff * pow(a, (*pKt)[tidx].j)*pow(b, (*pKt)[tidx].k)*pow(c, (*pKt)[tidx].l);
		}

		for (int tidx = 0; tidx < pKt2->linear_size(); tidx++) {

			cuBReal coeff;
			int order = (*pKt2)[tidx].j + (*pKt2)[tidx].k + (*pKt2)[tidx].l;
			if (order == 2) coeff = K1_AFM.j * (*pKt2)[tidx].i;
			else if (order == 4) coeff = K2_AFM.j * (*pKt2)[tidx].i;
			else if (order == 6) coeff = K3_AFM.j * (*pKt2)[tidx].i;
			else coeff = (*pKt2)[tidx].i;

			energyB += coeff * pow(a2, (*pKt2)[tidx].j)*pow(b2, (*pKt2)[tidx].k)*pow(c2, (*pKt2)[tidx].l);
		}

		return cuReal2(energyA, energyB);
	};

	//calculate dot products
	cuBReal a = (M[spin_index] * mcanis_ea1) / Ms_AFM.i;
	cuBReal b = (M[spin_index] * mcanis_ea2) / Ms_AFM.i;
	cuBReal c = (M[spin_index] * mcanis_ea3) / Ms_AFM.i;

	cuBReal a2 = (M2[spin_index] * mcanis_ea1) / Ms_AFM.j;
	cuBReal b2 = (M2[spin_index] * mcanis_ea2) / Ms_AFM.j;
	cuBReal c2 = (M2[spin_index] * mcanis_ea3) / Ms_AFM.j;

	cuReal2 energy_ = Get_Energy(a, b, c, a2, b2, c2);

	if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

		cuBReal anew = Mnew_A * mcanis_ea1 / Ms_AFM.j;
		cuBReal bnew = Mnew_A * mcanis_ea2 / Ms_AFM.j;
		cuBReal cnew = Mnew_A * mcanis_ea3 / Ms_AFM.j;

		cuBReal a2new = Mnew_B * mcanis_ea1 / Ms_AFM.j;
		cuBReal b2new = Mnew_B * mcanis_ea2 / Ms_AFM.j;
		cuBReal c2new = Mnew_B * mcanis_ea3 / Ms_AFM.j;

		cuReal2 energynew_ = Get_Energy(anew, bnew, cnew, a2new, b2new, c2new);

		return M.h.dim() * (energynew_ - energy_);
	}
	else return M.h.dim() * energy_;
}

//RoughnessCUDA
__device__ cuReal2 ManagedMeshCUDA::Get_EnergyChange_AFM_RoughnessCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B)
{
	if (pFmul_rough && pFmul_rough->linear_size()) {

		cuVEC_VC<cuReal3>& M = *pM;
		cuVEC_VC<cuReal3>& M2 = *pM2;

		cuReal33 Fmat = cuReal33(
			cuReal3((*pFmul_rough)[spin_index].x, (*pFomul_rough)[spin_index].x, (*pFomul_rough)[spin_index].y),
			cuReal3((*pFomul_rough)[spin_index].x, (*pFmul_rough)[spin_index].y, (*pFomul_rough)[spin_index].z),
			cuReal3((*pFomul_rough)[spin_index].y, (*pFomul_rough)[spin_index].z, (*pFmul_rough)[spin_index].z));

		cuBReal energy_ = 0.0;

		cuReal3 Mval = (M[spin_index] + M2[spin_index]) / 2;
		cuReal3 Hrough = Fmat * Mval;

		if (Mnew_A != cuReal3() && Mnew_B != cuReal3()) {

			cuReal3 Mvalnew = (Mnew_A + Mnew_B) / 2;
			cuReal3 Hrough_new = Fmat * Mvalnew;

			energy_ = -M.h.dim() * MU0 * (Hrough_new * Mvalnew - Hrough * Mval);
		}
		else energy_ = -M.h.dim() * MU0 * Hrough * Mval;

		return cuReal2(energy_, energy_);
	}
	else return cuReal2();
}

#endif