#include "stdafx.h"
#include "TMR.h"

#ifdef MODULE_COMPILATION_TMR

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TMRCUDA.h"
#endif

//-------------------Display Calculation Methods

//return x, y, or z component of spin current (component = 0, 1, or 2). Not used for TMR.
//VERIFIED - CORRECT
VEC<DBL3>& TMR::GetSpinCurrent(int component)
{
	return displayVEC;
}

//return spin torque computed from spin accumulation. Not used for TMR module.
//VERIFIED - CORRECT
VEC<DBL3>& TMR::GetSpinTorque(void)
{
	return displayVEC;
}

//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC. Not used for TMR module.
//VERIFIED - CORRECT
void TMR::CalculateDisplaySAInterfaceTorque(TransportBase* ptrans_sec, CMBNDInfo& contact)
{
}

#endif