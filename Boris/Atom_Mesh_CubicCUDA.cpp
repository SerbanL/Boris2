#include "stdafx.h"
#include "Atom_Mesh_CubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_Mesh_Cubic.h"

Atom_Mesh_CubicCUDA::Atom_Mesh_CubicCUDA(Atom_Mesh_Cubic* paMesh) :
	Atom_MeshCUDA(paMesh)
{
	paMeshCubic = paMesh;
}

Atom_Mesh_CubicCUDA::~Atom_Mesh_CubicCUDA()
{

}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError Atom_Mesh_CubicCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Mesh_CubicCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {

		//resize arrays held in Mesh - if changing mesh to new size then use resize : this maps values from old mesh to new mesh size (keeping magnitudes). Otherwise assign magnetization along x in whole mesh.
		if (M1()->size_cpu().dim()) {

			if (!M1()->resize((cuReal3)h, (cuRect)meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		}
		else if (!M1()->assign((cuReal3)h, (cuRect)meshRect, cuReal3(-(mu_s()->get_current_cpu()), 0, 0))) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		if (!Heff1()->assign((cuReal3)h, (cuRect)meshRect, cuReal3())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		//setup prng for MC methods
		if (prng()->initialize(GetSystemTickCount(), n.dim() / 128) != cudaSuccess) error(BERROR_OUTOFGPUMEMORY_NCRIT);

		//Setup cuaModules and cuaNumModules so they cann be used at runtime by Monte-Carlo methods
		if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {
			
			std::vector<int> modules_ids(paMeshCubic->pMod.size());

			for (int idx = 0; idx < modules_ids.size(); idx++) {

				modules_ids[idx] = paMeshCubic->pMod.get_ID_from_index(idx);
			}

			cuaModules.resize(modules_ids.size());
			cuaModules.copy_from_cpuvector(modules_ids);
			cuaNumModules.from_cpu((int)modules_ids.size());
		}
	}

	return error;
}

cuBReal Atom_Mesh_CubicCUDA::CheckMoveMesh(bool antisymmetric, double threshold)
{
	//move mesh algorithm applied to systems containing domain walls in order to simulate domain wall movement.
	//In this case the moments at the ends of the mesh must be fixed (ODE solver not applied to the ends).
	//Find the average moment projected along the absolute direction of the moment at the ends.
	//For a domain wall at the centre this value should be close to zero. If it exceeds a threshold then move mesh by a unit cellsize along +x or -x.

	if (!paMeshCubic->move_mesh_trigger) return 0;

	//the ends should not be completely empty, and must have a constant magnetisation direction
	int cells_fixed = (int)ceil_epsilon((double)n.x * MOVEMESH_ENDRATIO);

	cuReal3 M_end = M1()->average(n.dim(), cuBox(0, 0, 0, cells_fixed, n.y, n.z));
	cuReal3 M_av_left = M1()->average(n.dim(), cuBox(cells_fixed, 0, 0, n.x / 2, n.y, n.z));
	cuReal3 M_av_right = M1()->average(n.dim(), cuBox(n.x / 2, 0, 0, n.x - cells_fixed, n.y, n.z));

	cuBReal direction = (2 * cuBReal(antisymmetric) - 1);

	cuReal3 M_av = M_av_right + direction * M_av_left;

	if (GetMagnitude(M_end) * GetMagnitude(M_av)) {

		double Mratio = (M_end * M_av) / (GetMagnitude(M_end) * GetMagnitude(M_av));

		if (Mratio > threshold) return -h.x * direction;
		if (Mratio < -threshold) return h.x * direction;
	}

	return 0.0;
}

//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

//get exchange_couple_to_meshes status flag from the cpu version
bool Atom_Mesh_CubicCUDA::GetMeshExchangeCoupling(void)
{
	return paMeshCubic->GetMeshExchangeCoupling();
}

#endif
#endif