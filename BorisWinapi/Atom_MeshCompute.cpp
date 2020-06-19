#include "stdafx.h"
#include "Atom_Mesh.h"

#include "Atom_Exchange.h"
#include "Atom_DMExchange.h"
#include "Atom_iDMExchange.h"

//TO DO: review use of these 2 methods, as they could be combined with the micromagnetic ones : same module ids used, so could just define Compute_Exchange and GetEnergy_Max in ExchangeBase, then call using base pointer if module active.

//compute exchange energy density spatial variation and have it available to display in Cust_S
void Atom_Mesh::Compute_Exchange(void)
{
	if (IsModuleSet(MOD_EXCHANGE)) dynamic_cast<Atom_Exchange*>(pMod(MOD_EXCHANGE))->Compute_Exchange(displayVEC_SCA);
	else if (IsModuleSet(MOD_DMEXCHANGE)) dynamic_cast<Atom_DMExchange*>(pMod(MOD_DMEXCHANGE))->Compute_Exchange(displayVEC_SCA);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) dynamic_cast<Atom_iDMExchange*>(pMod(MOD_IDMEXCHANGE))->Compute_Exchange(displayVEC_SCA);
}

//get maximum exchange energy density modulus over specified rectangle
double Atom_Mesh::Get_Max_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_EXCHANGE)) return dynamic_cast<Atom_Exchange*>(pMod(MOD_EXCHANGE))->GetEnergy_Max(rectangle);
	else if (IsModuleSet(MOD_DMEXCHANGE)) return dynamic_cast<Atom_DMExchange*>(pMod(MOD_DMEXCHANGE))->GetEnergy_Max(rectangle);
	else if (IsModuleSet(MOD_IDMEXCHANGE)) return dynamic_cast<Atom_iDMExchange*>(pMod(MOD_IDMEXCHANGE))->GetEnergy_Max(rectangle);
	else return 0.0;
}