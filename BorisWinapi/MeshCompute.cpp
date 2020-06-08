#include "stdafx.h"
#include "Mesh.h"

#include "Exchange.h"
#include "DMExchange.h"
#include "iDMExchange.h"

//compute exchange energy density spatial variation and have it available to display in Cust_S
void Mesh::Compute_Exchange(void)
{
	if (IsModuleSet(MOD_EXCHANGE6NGBR)) reinterpret_cast<Exch_6ngbr_Neu*>(pMod(MOD_EXCHANGE6NGBR))->Compute_Exchange(displayVEC_SCA);
	else if (IsModuleSet(MOD_DMEXCHANGE)) reinterpret_cast<DMExchange*>(pMod(MOD_DMEXCHANGE))->Compute_Exchange(displayVEC_SCA);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) reinterpret_cast<iDMExchange*>(pMod(MOD_IDMEXCHANGE))->Compute_Exchange(displayVEC_SCA);
}

//get exchange energy density over entire mesh
double Mesh::Get_Exchange_EnergyDensity(void)
{
	if (IsModuleSet(MOD_EXCHANGE6NGBR)) return pMod(MOD_EXCHANGE6NGBR)->GetEnergyDensity();
	else if (IsModuleSet(MOD_DMEXCHANGE)) return pMod(MOD_DMEXCHANGE)->GetEnergyDensity();
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return pMod(MOD_IDMEXCHANGE)->GetEnergyDensity();
	else return 0.0;
}

//get exchange energy density over specified rectangle only
double Mesh::Get_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_EXCHANGE6NGBR)) return reinterpret_cast<Exch_6ngbr_Neu*>(pMod(MOD_EXCHANGE6NGBR))->GetEnergyDensity(rectangle);
	else if (IsModuleSet(MOD_DMEXCHANGE)) return reinterpret_cast<DMExchange*>(pMod(MOD_DMEXCHANGE))->GetEnergyDensity(rectangle);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return reinterpret_cast<iDMExchange*>(pMod(MOD_IDMEXCHANGE))->GetEnergyDensity(rectangle);
	else return 0.0;
}

//get maximum exchange energy density modulus over specified rectangle
double Mesh::Get_Max_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_EXCHANGE6NGBR)) return reinterpret_cast<Exch_6ngbr_Neu*>(pMod(MOD_EXCHANGE6NGBR))->GetEnergy_Max(rectangle);
	else if (IsModuleSet(MOD_DMEXCHANGE)) return reinterpret_cast<DMExchange*>(pMod(MOD_DMEXCHANGE))->GetEnergy_Max(rectangle);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return reinterpret_cast<iDMExchange*>(pMod(MOD_IDMEXCHANGE))->GetEnergy_Max(rectangle);
	else return 0.0;
}

//compute topological charge density spatial dependence and have it available to display in Cust_S
//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
void Mesh::Compute_TopoChargeDensity(void)
{
	if (M.linear_size()) {

		displayVEC_SCA.resize(h, meshRect);

#if COMPILECUDA == 1
		if (pMeshCUDA) {

			//compute topological charge density spatial variation in MeshCUDA::aux_vec_sca
			pMeshCUDA->Compute_TopoChargeDensity();

			//now transfer it to displayVEC_SCA ready to display
			pMeshCUDA->copy_aux_vec_sca(displayVEC_SCA);
			return;
		}
#endif

#pragma omp parallel for
		for (int idx = 0; idx < M.linear_size(); idx++) {

			if (M.is_not_empty(idx)) {

				double M_mag = M[idx].norm();

				DBL33 M_grad = M.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / M_mag;
				DBL3 dm_dy = M_grad.y / M_mag;

				displayVEC_SCA[idx] = (M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y / (4 * PI);
			}
			else displayVEC_SCA[idx] = 0.0;
		}
	}
}