#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

//----------------------------------- QUANTITY GETTERS

//returns M on the cpu, thus transfers M from gpu to cpu before returning if cuda enabled
VEC_VC<DBL3>& Mesh::Get_M(void)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
	}
#endif

	return M;
}

VEC_VC<DBL3>& Mesh::Get_M2(void)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

	return M2;
}

//returns charge current on the cpu, assuming transport module is enabled
VEC_VC<DBL3>& Mesh::Get_Jc(void)
{
	return dynamic_cast<Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent();
}

//returns bulk self-consistent spin torque on the cpu, assuming transport module is enabled and spin solver is enabled
VEC<DBL3>& Mesh::Get_SpinTorque(void)
{
	return dynamic_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque();
}

//returns interfacial self-consistent spin torque on the cpu, assuming transport module is enabled and spin solver is enabled
VEC<DBL3>& Mesh::Get_InterfacialSpinTorque(void)
{
	return dynamic_cast<STransport*>(pSMesh->GetSuperMeshModule(MODS_STRANSPORT))->GetInterfacialSpinTorque(dynamic_cast<Transport*>(pMod(MOD_TRANSPORT)));
}