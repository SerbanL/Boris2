#include "stdafx.h"
#include "Atom_Mesh.h"

#include "Atom_Exchange.h"

//compute exchange energy density spatial variation and have it available to display in Cust_S
void Atom_Mesh::Compute_Exchange(void)
{
	if (IsModuleSet(MOD_ATOM_EXCHANGE)) reinterpret_cast<Atom_Exchange*>(pMod(MOD_ATOM_EXCHANGE))->Compute_Exchange(displayVEC_SCA);
}

//get exchange energy density over entire mesh
double Atom_Mesh::Get_Exchange_EnergyDensity(void)
{
	if (IsModuleSet(MOD_ATOM_EXCHANGE)) return pMod(MOD_ATOM_EXCHANGE)->GetEnergyDensity();
	else return 0.0;
}

//get exchange energy density over specified rectangle only
double Atom_Mesh::Get_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_ATOM_EXCHANGE)) return reinterpret_cast<Atom_Exchange*>(pMod(MOD_ATOM_EXCHANGE))->GetEnergyDensity(rectangle);
	else return 0.0;
}

//get maximum exchange energy density modulus over specified rectangle
double Atom_Mesh::Get_Max_Exchange_EnergyDensity(Rect& rectangle)
{
	if (IsModuleSet(MOD_ATOM_EXCHANGE)) return reinterpret_cast<Atom_Exchange*>(pMod(MOD_ATOM_EXCHANGE))->GetEnergy_Max(rectangle);
	else return 0.0;
}