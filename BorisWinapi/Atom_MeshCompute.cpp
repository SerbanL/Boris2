#include "stdafx.h"
#include "Atom_Mesh.h"

#include "Atom_Exchange.h"
#include "Atom_DMExchange.h"
#include "Atom_iDMExchange.h"

//compute exchange energy density spatial variation and have it available to display in Cust_S
void Atom_Mesh::Compute_Exchange(void)
{
	if (IsModuleSet(MOD_EXCHANGE)) dynamic_cast<Atom_Exchange*>(pMod(MOD_EXCHANGE))->Compute_Exchange(displayVEC_SCA);
	else if (IsModuleSet(MOD_DMEXCHANGE)) dynamic_cast<Atom_DMExchange*>(pMod(MOD_DMEXCHANGE))->Compute_Exchange(displayVEC_SCA);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) dynamic_cast<Atom_iDMExchange*>(pMod(MOD_IDMEXCHANGE))->Compute_Exchange(displayVEC_SCA);
}

//get exchange energy density over entire mesh
double Atom_Mesh::Get_Exchange_EnergyDensity(void)
{
	if (IsModuleSet(MOD_EXCHANGE)) return pMod(MOD_EXCHANGE)->GetEnergyDensity();
	else if (IsModuleSet(MOD_DMEXCHANGE)) return pMod(MOD_DMEXCHANGE)->GetEnergyDensity();
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return pMod(MOD_IDMEXCHANGE)->GetEnergyDensity();
	else return 0.0;
}

//get exchange energy density over specified rectangle only
double Atom_Mesh::Get_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_EXCHANGE)) return dynamic_cast<Atom_Exchange*>(pMod(MOD_EXCHANGE))->GetEnergyDensity(rectangle);
	else if (IsModuleSet(MOD_DMEXCHANGE)) return dynamic_cast<Atom_DMExchange*>(pMod(MOD_DMEXCHANGE))->GetEnergyDensity(rectangle);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return dynamic_cast<Atom_iDMExchange*>(pMod(MOD_IDMEXCHANGE))->GetEnergyDensity(rectangle);
	else return 0.0;
}

//get maximum exchange energy density modulus over specified rectangle
double Atom_Mesh::Get_Max_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_EXCHANGE)) return dynamic_cast<Atom_Exchange*>(pMod(MOD_EXCHANGE))->GetEnergy_Max(rectangle);
	else if (IsModuleSet(MOD_DMEXCHANGE)) return dynamic_cast<Atom_DMExchange*>(pMod(MOD_DMEXCHANGE))->GetEnergy_Max(rectangle);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return dynamic_cast<Atom_iDMExchange*>(pMod(MOD_IDMEXCHANGE))->GetEnergy_Max(rectangle);
	else return 0.0;
}