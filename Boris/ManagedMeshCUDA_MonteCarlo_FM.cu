#include "ManagedMeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"
#include "MeshParamsControlCUDA.h"

#include "ModulesDefs.h"

//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

//Energy Deltas

//Ferromagnetic

//switch function which adds all assigned energy contributions in this mesh to calculate energy change from current spin to Mnew spin : return energy change as new - old
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

		case MOD_EXCHANGE:
			energy += Get_EnergyChange_FM_ExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_DMEXCHANGE:
			energy += Get_EnergyChange_FM_DMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_EnergyChange_FM_iDMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_SURFEXCHANGE:
			energy += Get_EnergyChange_FM_SurfExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_ZEEMAN:
			energy += Get_EnergyChange_FM_ZeemanCUDA(spin_index, Mnew, Ha);
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

	return -((cuBReal)MU0 / 2) * M.h.dim() * (
		(Mnew * cuReal3(-Nxy.x * Mnew.x, -Nxy.y * Mnew.y, -(1 - Nxy.x - Nxy.y) * Mnew.z)) -
		(M[spin_index] * cuReal3(-Nxy.x * M[spin_index].x, -Nxy.y * M[spin_index].y, -(1 - Nxy.x - Nxy.y) * M[spin_index].z)));
}

//Demag
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_DemagCUDA(int spin_index, cuReal3 Mnew)
{
	if (pDemag_Heff && pDemag_Heff->linear_size()) {

		cuVEC_VC<cuReal3>& M = *pM;

		return -(cuBReal)MU0 * M.h.dim() * (*pDemag_Heff)[M.cellidx_to_position(spin_index)] * (Mnew - M[spin_index]);
	}
	else return 0.0;
}

//Exch_6ngbr_Neu
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_ExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
	Mnew.renormalize(M[spin_index].norm());

	cuBReal Ms = *pMs;
	cuBReal A = *pA;
	update_parameters_mcoarse(spin_index, *pA, A, *pMs, Ms);

	cuReal3 Hexch = (2 * A / ((cuBReal)MU0*Ms*Ms)) * M.delsq_neu(spin_index);
	cuBReal energy_ = -(cuBReal)MU0 * M[spin_index] * Hexch / 2;

	cuReal3 Mold = M[spin_index];
	M[spin_index] = Mnew;
	Hexch = (2 * A / ((cuBReal)MU0*Ms*Ms)) * M.delsq_neu(spin_index);
	cuBReal energynew_ = -(cuBReal)MU0 * M[spin_index] * Hexch / 2;
	M[spin_index] = Mold;

	//multiply by 2 as we are not double-counting here (same as for Demag)
	return M.h.dim() * (energynew_ - energy_) * 2;
}

//DMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_DMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
	Mnew.renormalize(M[spin_index].norm());

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

	//new spin energy
	cuReal3 Mold = M[spin_index];
	M[spin_index] = Mnew;
	cuBReal energy_delta = Get_Energy();

	//minus old spin energy
	M[spin_index] = Mold;
	energy_delta -= Get_Energy();

	//do not divide by 2 as we are not double-counting here
	return -MU0 * M.h.dim() * energy_delta;
}

//iDMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_iDMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M = *pM;

	//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
	Mnew.renormalize(M[spin_index].norm());

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

	//new spin energy
	cuReal3 Mold = M[spin_index];
	M[spin_index] = Mnew;
	cuBReal energy_delta = Get_Energy();

	//minus old spin energy
	M[spin_index] = Mold;
	energy_delta -= Get_Energy();

	//do not divide by 2 as we are not double-counting here
	return -MU0 * M.h.dim() * energy_delta;
}

//SurfExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_SurfExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	//TO DO

	return 0.0;
}

//ZeemanCUDA
__device__ cuBReal ManagedMeshCUDA::Get_EnergyChange_FM_ZeemanCUDA(int spin_index, cuReal3 Mnew, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal cHA = *pcHA;
	update_parameters_mcoarse(spin_index, *pcHA, cHA);

	return -M.h.dim() * (Mnew - M[spin_index]) * (cuBReal)MU0 * (cHA * Ha);
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

	//calculate m.ea dot product
	cuBReal dotprod = (M[spin_index] * mcanis_ea1) / Ms;
	cuBReal dotprod_new = Mnew * mcanis_ea1 / Ms;
	cuBReal dpsq = dotprod * dotprod;
	cuBReal dpsq_new = dotprod_new * dotprod_new;

	return M.h.dim() * ((K1 + K2 * (1 - dpsq_new)) * (1 - dpsq_new) - (K1 + K2 * (1 - dpsq)) * (1 - dpsq));
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

	cuBReal d1_new = S_new * mcanis_ea1;
	cuBReal d2_new = S_new * mcanis_ea2;
	cuBReal d3_new = S_new * mcanis_ea3;
	cuBReal d123_new = d1_new * d2_new * d3_new;

	return M.h.dim() * (
		(K1 * (d1_new*d1_new*d2_new*d2_new + d1_new * d1_new*d3_new*d3_new + d2_new * d2_new*d3_new*d3_new) + K2 * d123_new*d123_new) -
		(K1 * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3) + K2 * d123*d123));
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

	//OLD
	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1 = M[spin_index] * mcanis_ea1 / Ms;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1 = M[spin_index] * mcanis_ea2 / Ms;
	cuBReal b2 = M[spin_index] * mcanis_ea3 / Ms;

	//NEW
	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1new = Mnew * mcanis_ea1 / Ms;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1new = Mnew * mcanis_ea2 / Ms;
	cuBReal b2new = Mnew * mcanis_ea3 / Ms;

	return M.h.dim() * ((K1 * (1 - u1new * u1new) + K2 * b1new*b1new*b2new*b2new) - (K1 * (1 - u1 * u1) + K2 * b1*b1*b2*b2));
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

	//calculate dot products
	cuBReal a = M[spin_index].normalized() * mcanis_ea1 / Ms;
	cuBReal b = M[spin_index].normalized() * mcanis_ea2 / Ms;
	cuBReal c = M[spin_index].normalized() * mcanis_ea3 / Ms;

	cuBReal anew = Mnew * mcanis_ea1 / Ms;
	cuBReal bnew = Mnew * mcanis_ea2 / Ms;
	cuBReal cnew = Mnew * mcanis_ea3 / Ms;

	cuBReal energy_ = 0.0, energynew_ = 0.0;

	for (int tidx = 0; tidx < pKt->linear_size(); tidx++) {

		cuBReal coeff;
		int order = (*pKt)[tidx].j + (*pKt)[tidx].k + (*pKt)[tidx].l;
		if (order == 2) coeff = K1 * (*pKt)[tidx].i;
		else if (order == 4) coeff = K2 * (*pKt)[tidx].i;
		else if (order == 6) coeff = K3 * (*pKt)[tidx].i;
		else coeff = (*pKt)[tidx].i;

		energy_ += coeff * pow(a, (*pKt)[tidx].j)*pow(b, (*pKt)[tidx].k)*pow(c, (*pKt)[tidx].l);
		energynew_ += coeff * pow(anew, (*pKt)[tidx].j)*pow(bnew, (*pKt)[tidx].k)*pow(cnew, (*pKt)[tidx].l);
	}

	return M.h.dim() * (energynew_ - energy_);
}

//Spin Energy

//Ferromagnetic

//switch function which adds all assigned energy contributions in this mesh to calculate energy change from current spin to Mnew spin : return energy change as new - old
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM(int spin_index, int*& cuModules, int numModules, cuReal3& Ha)
{
	cuBReal energy = 0.0;

	for (int idx = 0; idx < numModules; idx++) {

		switch (cuModules[idx]) {

		case MOD_DEMAG_N:
			energy += Get_Energy_FM_DemagNCUDA(spin_index);
			break;

		case MOD_DEMAG:
			energy += Get_Energy_FM_DemagCUDA(spin_index);
			break;

		case MOD_EXCHANGE:
			energy += Get_Energy_FM_ExchangeCUDA(spin_index);
			break;

		case MOD_DMEXCHANGE:
			energy += Get_Energy_FM_DMExchangeCUDA(spin_index);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_Energy_FM_iDMExchangeCUDA(spin_index);
			break;

		case MOD_SURFEXCHANGE:
			energy += Get_Energy_FM_SurfExchangeCUDA(spin_index);
			break;

		case MOD_ZEEMAN:
			energy += Get_Energy_FM_ZeemanCUDA(spin_index, Ha);
			break;

		case MOD_ANIUNI:
			energy += Get_Energy_FM_AnisotropyCUDA(spin_index);
			break;

		case MOD_ANICUBI:
			energy += Get_Energy_FM_AnisotropyCubiCUDA(spin_index);
			break;

		case MOD_ANIBI:
			energy += Get_Energy_FM_AnisotropyBiaxialCUDA(spin_index);
			break;

		case MOD_ANITENS:
			energy += Get_Energy_FM_AnisotropyTensorialCUDA(spin_index);
			break;

		default:
			break;
		}
	}

	return energy;
}

//Demag_N
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_DemagNCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M = *pM;

	//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
	cuReal2 Nxy = *pNxy;

	cuBReal Nz = (1 - Nxy.x - Nxy.y);

	return -((cuBReal)MU0 / 2) * M.h.dim() * (M[spin_index] * cuReal3(-Nxy.x * M[spin_index].x, -Nxy.y * M[spin_index].y, -(1 - Nxy.x - Nxy.y) * M[spin_index].z));
}

//Demag
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_DemagCUDA(int spin_index)
{
	if (pDemag_Heff && pDemag_Heff->linear_size()) {

		cuVEC_VC<cuReal3>& M = *pM;

		return -(cuBReal)MU0 * M.h.dim() * (*pDemag_Heff)[M.cellidx_to_position(spin_index)] * M[spin_index];
	}
	else return 0.0;
}

//Exch_6ngbr_Neu
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_ExchangeCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal Ms = *pMs;
	cuBReal A = *pA;
	update_parameters_mcoarse(spin_index, *pA, A, *pMs, Ms);
	
	cuReal3 Hexch = (2 * A / ((cuBReal)MU0*Ms*Ms)) * M.delsq_neu(spin_index);
	cuBReal energy_ = -(cuBReal)MU0 * M[spin_index] * Hexch / 2;

	//multiply by 2 as we are not double-counting here (same as for Demag)
	return M.h.dim() * energy_ * 2;
}

//DMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_DMExchangeCUDA(int spin_index)
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

	//do not divide by 2 as we are not double-counting here
	return -MU0 * M.h.dim() * Get_Energy();
}

//iDMExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_iDMExchangeCUDA(int spin_index)
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

	//do not divide by 2 as we are not double-counting here
	return -MU0 * M.h.dim() * Get_Energy();
}

//SurfExchangeCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_SurfExchangeCUDA(int spin_index)
{
	return 0.0;
}

//ZeemanCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_ZeemanCUDA(int spin_index, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal cHA = *pcHA;
	update_parameters_mcoarse(spin_index, *pcHA, cHA);

	return -M.h.dim() * M[spin_index] * (cuBReal)MU0 * (cHA * Ha);
}

//AnisotropyCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_AnisotropyCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M = *pM;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal Ms = *pMs;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	update_parameters_mcoarse(spin_index, *pMs, Ms, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1);

	//calculate m.ea dot product
	cuBReal dotprod = (M[spin_index] * mcanis_ea1) / Ms;
	cuBReal dpsq = dotprod * dotprod;

	return M.h.dim() * (K1 + K2 * (1 - dpsq)) * (1 - dpsq);
}

//AnisotropyCubiCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_AnisotropyCubiCUDA(int spin_index)
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

	//calculate m.ea1, m.ea2 and m.ea3 dot products
	cuBReal d1 = S * mcanis_ea1;
	cuBReal d2 = S * mcanis_ea2;
	cuBReal d3 = S * mcanis_ea3;
	cuBReal d123 = d1 * d2 * d3;

	return M.h.dim() * (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123);
}

//AnisotropyBiaxialCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_AnisotropyBiaxialCUDA(int spin_index)
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

	return M.h.dim() * (K1 * (1 - u1 * u1) + K2 * b1*b1*b2*b2);
}

//AnisotropyTensorialCUDA
__device__ cuBReal ManagedMeshCUDA::Get_Energy_FM_AnisotropyTensorialCUDA(int spin_index)
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

	//calculate dot products
	cuBReal a = M[spin_index].normalized() * mcanis_ea1 / Ms;
	cuBReal b = M[spin_index].normalized() * mcanis_ea2 / Ms;
	cuBReal c = M[spin_index].normalized() * mcanis_ea3 / Ms;

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

	return M.h.dim() * energy_;
}

#endif