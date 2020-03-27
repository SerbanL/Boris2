#include "stdafx.h"
#include "Mesh_DiamagneticCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"

DiaMeshCUDA::DiaMeshCUDA(DiaMesh* pMesh) :
	MeshCUDA(pMesh)
{
	pDiaMesh = pMesh;
}

DiaMeshCUDA::~DiaMeshCUDA()
{

}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
BError DiaMeshCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(FMeshCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		//resize arrays held in Mesh - if changing mesh to new size then use resize : this maps values from old mesh to new mesh size (keeping magnitudes). Otherwise assign magnetization along x in whole mesh.
		if (M()->size_cpu().dim()) {

			if (!M()->resize((cuReal3)h, (cuRect)meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else if (!M()->assign((cuReal3)h, (cuRect)meshRect, cuReal3(0, 0, 0))) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		if (!Heff()->assign((cuReal3)h, (cuRect)meshRect, cuReal3())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	return error;
}

#endif
#endif