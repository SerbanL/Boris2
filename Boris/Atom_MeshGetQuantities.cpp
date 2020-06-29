#include "stdafx.h"
#include "Atom_Mesh.h"

#include "SuperMesh.h"

//----------------------------------- QUANTITY GETTERS

//returns M on the cpu, thus transfers M from gpu to cpu before returning if cuda enabled
VEC_VC<DBL3>& Atom_Mesh::Get_M1(void)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		paMeshCUDA->M1()->copy_to_cpuvec(M1);
	}
#endif

	return M1;
}