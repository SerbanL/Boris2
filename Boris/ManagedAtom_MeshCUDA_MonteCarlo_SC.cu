#include "stdafx.h"
#include "ManagedAtom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"
#include "Atom_MeshParamsControlCUDA.h"

#include "ModulesDefs.h"

//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

//SIMPLE CUBIC

//switch function which adds all assigned energy contributions in this mesh
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC(int spin_index, cuReal3 Mnew, int*& cuaModules, int numModules, cuReal3& Ha)
{
	cuBReal energy = 0.0;

	for (int idx = 0; idx < numModules; idx++) {

		switch (cuaModules[idx]) {

		case MOD_DEMAG_N:
			energy += Get_Atomistic_EnergyChange_SC_DemagNCUDA(spin_index, Mnew);
			break;

		case MOD_EXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_ExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_DMEXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_DMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_iDMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_ZEEMAN:
			energy += Get_Atomistic_EnergyChange_SC_ZeemanCUDA(spin_index, Mnew, Ha);
			break;

		case MOD_ANIUNI:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyCUDA(spin_index, Mnew);
			break;

		case MOD_ANICUBI:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyCubiCUDA(spin_index, Mnew);
			break;

		case MOD_ANIBI:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyBiaxialCUDA(spin_index, Mnew);
			break;

		case MOD_ANITENS:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyTensorialCUDA(spin_index, Mnew);
			break;

		default:
			break;
		}
	}

	return energy;
}

//Atom_Demag_N
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_DemagNCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
	cuReal2 Nxy = *pNxy;

	cuBReal Nz = (1 - Nxy.x - Nxy.y);

	cuBReal r = (cuBReal)MUB / M1.h.dim();
	cuReal3 S = M1[spin_index];

	return ((cuBReal)MUB_MU0 / 2) * r * (Nxy.x * (Mnew.x*Mnew.x - S.x*S.x) + Nxy.y * (Mnew.y*Mnew.y - S.y*S.y) + Nz * (Mnew.z*Mnew.z - S.z*S.z));
}

//Atom_ExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_ExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	update_parameters_mcoarse(spin_index, *pJ, J);

	return -J * ((Mnew.normalized() - M1[spin_index].normalized()) * M1.ngbr_dirsum(spin_index));
}

//Atom_DMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_DMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D);

	//local spin energy given:
	//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
	//2) DM exchange: D * Sum_over_neighbors_j(rij . (Si x Sj))

	//Now anisotropic_ngbr_dirsum returns rij x Sj, and Si . (rij x Sj) = -Si. (Sj x rij) = -rij . (Si x Sj)
	//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.anisotropic_ngbr_dirsum(spin_index));

	return (Mnew.normalized() - M1[spin_index].normalized()) * (-J * M1.ngbr_dirsum(spin_index) - D * M1.anisotropic_ngbr_dirsum(spin_index));
}

//Atom_iDMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_iDMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D);

	//local spin energy given:
	//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
	//2) iDM exchange: D * Sum_over_neighbors_j((rij x z) . (Si x Sj))

	//Now zanisotropic_ngbr_dirsum returns (rij x z) x Sj, and Si . ((rij x z) x Sj) = -Si. (Sj x (rij x z)) = -(rij x z) . (Si x Sj)
	//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.zanisotropic_ngbr_dirsum(spin_index));

	return (Mnew.normalized() - M1[spin_index].normalized()) * (-J * M1.ngbr_dirsum(spin_index) - D * M1.zanisotropic_ngbr_dirsum(spin_index));
}

//Atom_ZeemanCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_ZeemanCUDA(int spin_index, cuReal3 Mnew, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal cHA = *pcHA;
	update_parameters_mcoarse(spin_index, *pcHA, cHA);

	return -(cuBReal)MUB * (Mnew - M1[spin_index]) * (cuBReal)MU0 * (cHA * Ha);
}

//Atom_AnisotropyCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	update_parameters_mcoarse(spin_index, *pK2, K2, *pmcanis_ea1, mcanis_ea1);

	//calculate m.ea dot product
	cuBReal dotprod = M1[spin_index].normalized() * mcanis_ea1;
	cuBReal dotprod_new = Mnew.normalized() * mcanis_ea1;
	cuBReal dpsq = dotprod * dotprod;
	cuBReal dpsq_new = dotprod_new * dotprod_new;

	//Hamiltonian contribution as -Ku * (S * ea)^2, where S is the local spin direction
	return -K1 * (dotprod_new * dotprod_new - dotprod * dotprod) - K2 * (dpsq_new * dpsq_new - dpsq * dpsq);
}

//Atom_AnisotropyCubiCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyCubiCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	cuReal3 S = M1[spin_index].normalized();
	cuReal3 Snew = Mnew.normalized();

	//calculate m.ea1, m.ea2 and m.ea3 dot products
	cuBReal d1 = S * mcanis_ea1;
	cuBReal d2 = S * mcanis_ea2;
	cuBReal d3 = S * mcanis_ea3;

	cuBReal d1_new = Snew * mcanis_ea1;
	cuBReal d2_new = Snew * mcanis_ea2;
	cuBReal d3_new = Snew * mcanis_ea3;

	//Hamiltonian contribution as K * (Sx^2*Sy^2 + Sx^2*Sz^2 + Sy^2*Sz^2), where S is the local spin direction (for easy axes coinciding with the xyz system)
	//This is equivalent to the form -K/2 * (Sx^4 + Sy^4 + Sz^4) - energy zero point differs but that's immaterial.
	//Also note the correct signs here for given easy axes (need to be careful, some publications have this wrong).
	return K1 * (d1_new*d1_new*d2_new*d2_new + d1_new*d1_new*d3_new*d3_new + d2_new*d2_new*d3_new*d3_new - d1*d1*d2*d2 - d1*d1*d3*d3 - d2*d2*d3*d3)
		+ K2 * (d1_new*d2_new*d3_new*d1_new*d2_new*d3_new - d1*d2*d3*d1*d2*d3);
}

//Atom_AnisotropyBiaxialCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyBiaxialCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	//OLD
	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1 = M1[spin_index].normalized() * mcanis_ea1;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1 = M1[spin_index].normalized() * mcanis_ea2;
	cuBReal b2 = M1[spin_index].normalized() * mcanis_ea3;

	//NEW
	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1new = Mnew.normalized() * mcanis_ea1;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1new = Mnew.normalized() * mcanis_ea2;
	cuBReal b2new = Mnew.normalized() * mcanis_ea3;

	return (K1 * (1 - u1new*u1new) + K2 * b1new*b1new*b2new*b2new) - (K1 * (1 - u1*u1) + K2 * b1*b1*b2*b2);
}

//Atom_AnisotropyTensorialCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyTensorialCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal K3 = *pK3;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pK1, K1, *pK2, K2, *pK3, K3, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	//calculate dot products
	cuBReal a = M1[spin_index].normalized() * mcanis_ea1;
	cuBReal b = M1[spin_index].normalized() * mcanis_ea2;
	cuBReal c = M1[spin_index].normalized() * mcanis_ea3;

	cuBReal anew = Mnew.normalized() * mcanis_ea1;
	cuBReal bnew = Mnew.normalized() * mcanis_ea2;
	cuBReal cnew = Mnew.normalized() * mcanis_ea3;

	cuBReal energy_ = 0.0, energynew_ = 0.0;

	for (int tidx = 0; tidx < pKt->linear_size(); tidx++) {

		//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
		//(-d / mu0 mu_s) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

		cuBReal coeff;
		int order = (*pKt)[tidx].j + (*pKt)[tidx].k + (*pKt)[tidx].l;
		if (order == 2) coeff = K1 * (*pKt)[tidx].i;
		else if (order == 4) coeff = K2 * (*pKt)[tidx].i;
		else if (order == 6) coeff = K3 * (*pKt)[tidx].i;
		else coeff = (*pKt)[tidx].i;

		energy_ += coeff * pow(a, (*pKt)[tidx].j)*pow(b, (*pKt)[tidx].k)*pow(c, (*pKt)[tidx].l);
		energynew_ += coeff * pow(anew, (*pKt)[tidx].j)*pow(bnew, (*pKt)[tidx].k)*pow(cnew, (*pKt)[tidx].l);
	}

	return energynew_ - energy_;
}

#endif