#include "stdafx.h"
#include "SuperMesh.h"

SuperMesh::SuperMesh(void) :
	ProgramStateNames(this,
		{
			VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(n_fm), VINFO(h_fm), VINFO(sMeshRect_fm),
			VINFO(n_e), VINFO(h_e), VINFO(sMeshRect_e),
			VINFO(odeSolver),
			VINFO(pMesh),
			VINFO(pSMod),
			VINFO(activeMeshName), VINFO(superMeshHandle),
			VINFO(scale_rects), VINFO(coupled_dipoles)
		}, 
		{
			//Mesh implementations
			IINFO(FMesh), IINFO(DipoleMesh), IINFO(MetalMesh), IINFO(InsulatorMesh), IINFO(AFMesh),
			//Super-mesh modules implementations (for pSMod)
			IINFO(SDemag), IINFO(StrayField), IINFO(STransport), IINFO(Oersted), IINFO(SHeat)
		})
{
	//default ODE settings
	SetODE(ODE_LLG, EVAL_RKF);

	//default starting state
	AddMesh(activeMeshName, MESH_FERROMAGNETIC, Rect(DBL3(80e-9, 80e-9, 10e-9)));
}

SuperMesh::~SuperMesh() 
{
	//clean-up meshes
	clear_vector(pMesh);
	clear_vector(pSMod);
}

BError SuperMesh::Error_On_Create(void)
{ 
	//go through meshes (and their modules)
	for (int idx = 0; idx < (int)pMesh.size(); idx++) {

		if(!error_on_create) error_on_create = pMesh[idx]->Error_On_Create();
	}

	//go through super-mesh modules
	for (int idx = 0; idx < (int)pSMod.size(); idx++) {

		if (!error_on_create) error_on_create = pSMod[idx]->Error_On_Create();
	}

	BError return_error = error_on_create;

	//clear error_on_create after reading it - the caller will restore program state
	error_on_create.clear();

	return return_error;
}

//----------------------------------- ODE SOLVER CONTROL

BError SuperMesh::SetODE(ODE_ setOde, EVAL_ evalMethod)
{ 
	BError error(__FUNCTION__); 

	if (setOde <= ODE_ERROR || evalMethod <= EVAL_ERROR) return error(BERROR_INCORRECTNAME);

	error = odeSolver.SetODE(setOde, evalMethod);

	error = UpdateConfiguration(UPDATECONFIG_ODE_SOLVER);

	return error; 
}

void SuperMesh::SetMoveMeshTrigger(bool status, string meshName)
{
	//if meshName not contained then we want meshId to be -1 : this means we'll be setting trigger on first ferromagnetic mesh
	int meshId = -1;

	if (contains(meshName)) meshId = pMesh[meshName]->get_id();

	//if calling with meshId = -1 then first ferromagnetic mesh will be used
	odeSolver.SetMoveMeshTrigger(status, meshId);

	UpdateConfiguration(UPDATECONFIG_ODE_MOVEMESH);
}

//--------------------------------------------------------- VALUE SETTERS

BError SuperMesh::SetFMSMeshCellsize(DBL3 h_fm_)
{ 
	BError error(__FUNCTION__);

	h_fm = h_fm_;
	
	error = UpdateConfiguration(UPDATECONFIG_SMESH_CELLSIZE);

	return error;
}

BError SuperMesh::SetESMeshCellsize(DBL3 h_e_)
{
	BError error(__FUNCTION__);
	
	h_e = h_e_;

	error = UpdateConfiguration(UPDATECONFIG_SMESH_CELLSIZE);

	return error;
}

//--------------------------------------------------------- GET PROPERTIES / VALUES

//get total volume energy density
double SuperMesh::GetTotalEnergy(void)
{
#if COMPILECUDA == 1

	//If CUDA enabled then obtain total energy density from all modules -> this involves potentially many individual cuBReal transfers from gpu to cpu so not ideal
	//This is fine since we don't need to get the total energy often (but will need optimizing if we need to use the total energy density every iteration in some solvers in the future).
	//If we've just switched CUDA on but didn't initialize yet we won't have the correct energy_density_weights (or even memory size!), so skip this for now.
	if (pSMeshCUDA && energy_density_weights.size() == pMesh.size()) {

		total_energy_density = 0.0;

		//collect energy density values from modules in individual meshes
		for (int idx = 0; idx < (int)pMesh.size(); idx++) {

			total_energy_density += (pMesh[idx]->GetEnergy(MOD_ALL) * energy_density_weights[idx]);
		}

		//collect energy density values from supermesh modules
		for (int idx = 0; idx < (int)pSMod.size(); idx++) {

			total_energy_density += pSMod[idx]->GetEnergy();
		}
	}
#endif

	return total_energy_density;
}