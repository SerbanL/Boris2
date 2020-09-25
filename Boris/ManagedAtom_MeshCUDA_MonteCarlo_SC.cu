#include "stdafx.h"
#include "ManagedAtom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"
#include "Atom_MeshParamsControlCUDA.h"

#include "ModulesDefs.h"

//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

//SIMPLE CUBIC

//switch function which adds all assigned energy contributions in this mesh
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC(int spin_index, int*& cuaModules, int numModules, cuReal3& Ha)
{
	cuBReal energy = 0.0;

	for (int idx = 0; idx < numModules; idx++) {

		switch (cuaModules[idx]) {

		case MOD_EXCHANGE:
			energy += Get_Atomistic_Energy_SC_ExchangeCUDA(spin_index);
			break;
		
		case MOD_DMEXCHANGE:
			energy += Get_Atomistic_Energy_SC_DMExchangeCUDA(spin_index);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_Atomistic_Energy_SC_iDMExchangeCUDA(spin_index);
			break;

		case MOD_ZEEMAN:
			energy += Get_Atomistic_Energy_SC_ZeemanCUDA(spin_index, Ha);
			break;

		case MOD_ANIUNI:
			energy += Get_Atomistic_Energy_SC_AnisotropyCUDA(spin_index);
			break;
			
		case MOD_ANICUBI:
			energy += Get_Atomistic_Energy_SC_AnisotropyCubiCUDA(spin_index);
			break;

		default:
			break;
		}
	}
	
	return energy;
}

//Atom_ExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC_ExchangeCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;
	
	cuBReal J = *pJ;
	update_parameters_mcoarse(spin_index, *pJ, J);

	return -J * (M1[spin_index].normalized() * M1.ngbr_dirsum(spin_index));
}

//Atom_DMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC_DMExchangeCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D);

	return -J * (M1[spin_index].normalized() * M1.ngbr_dirsum(spin_index));

	//local spin energy given:
	//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
	//2) DM exchange: D * Sum_over_neighbors_j(rij . (Si x Sj))

	//Now anisotropic_ngbr_dirsum returns rij x Sj, and Si . (rij x Sj) = -Si. (Sj x rij) = -rij . (Si x Sj)
	//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.anisotropic_ngbr_dirsum(spin_index));

	return M1[spin_index].normalized() * (-J * M1.ngbr_dirsum(spin_index) - D * M1.anisotropic_ngbr_dirsum(spin_index));
}

//Atom_iDMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC_iDMExchangeCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D);

	return -J * (M1[spin_index].normalized() * M1.ngbr_dirsum(spin_index));

	//local spin energy given:
		//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
		//2) iDM exchange: D * Sum_over_neighbors_j((rij x z) . (Si x Sj))

		//Now zanisotropic_ngbr_dirsum returns (rij x z) x Sj, and Si . ((rij x z) x Sj) = -Si. (Sj x (rij x z)) = -(rij x z) . (Si x Sj)
		//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.zanisotropic_ngbr_dirsum(spin_index));

	return M1[spin_index].normalized() * (-J * M1.ngbr_dirsum(spin_index) - D * M1.zanisotropic_ngbr_dirsum(spin_index));
}

//Atom_ZeemanCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC_ZeemanCUDA(int spin_index, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal cHA = *pcHA;
	update_parameters_mcoarse(spin_index, *pcHA, cHA);

	return -(cuBReal)MUB * M1[spin_index] * (cuBReal)MU0 * (cHA * Ha);
}

//Atom_AnisotropyCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC_AnisotropyCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K = *pK;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	update_parameters_mcoarse(spin_index, *pK, K, *pmcanis_ea1, mcanis_ea1);

	//calculate m.ea dot product
	cuBReal dotprod = M1[spin_index].normalized() * mcanis_ea1;

	//Hamiltonian contribution as -Ku * (S * ea)^2, where S is the local spin direction
	return -K * dotprod * dotprod;
}

//Atom_AnisotropyCubiCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_Energy_SC_AnisotropyCubiCUDA(int spin_index)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K = *pK;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	update_parameters_mcoarse(spin_index, *pK, K, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2);

	//vector product of ea1 and ea2 : the third orthogonal axis
	cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

	cuReal3 S = M1[spin_index].normalized();

	//calculate m.ea1, m.ea2 and m.ea3 dot products
	cuBReal d1 = S * mcanis_ea1;
	cuBReal d2 = S * mcanis_ea2;
	cuBReal d3 = S * mcanis_ea3;

	//Hamiltonian contribution as K * (Sx^2*Sy^2 + Sx^2*Sz^2 + Sy^2*Sz^2), where S is the local spin direction (for easy axes coinciding with the xyz system)
	//This is equivalent to the form -K/2 * (Sx^4 + Sy^4 + Sz^4) - energy zero point differs but that's immaterial.
	//Also note the correct signs here for given easy axes (need to be careful, some publications have this wrong).
	return K * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3);
}

#endif