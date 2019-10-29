#include "stdafx.h"
#include "Zeeman.h"

#ifdef MODULE_ZEEMAN

#include "Mesh_Ferromagnetic.h"
#include "Mesh_AntiFerromagnetic.h"

#if COMPILECUDA == 1
#include "ZeemanCUDA.h"
#endif

Zeeman::Zeeman(Mesh *pMesh_) : 
	Modules(),
	ProgramStateNames(this, { VINFO(Ha) }, {})
{
	pMesh = pMesh_;

	Ha = DBL3(0,0,0);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Zeeman::~Zeeman() 
{
}

BError Zeeman::Initialize(void)
{
	BError error(CLASS_STR(Zeeman));

	initialized = true;

	return error;
}

BError Zeeman::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Zeeman));

	Uninitialize();

	Initialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Zeeman::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Zeeman));

#if COMPILECUDA == 1

		if (pMesh->pMeshCUDA) {

			//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
			//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
			pModuleCUDA = new ZeemanCUDA(pMesh->pMeshCUDA, this);
			error = pModuleCUDA->Error_On_Create();
		}

#endif

	return error;
}

double Zeeman::UpdateField(void) 
{
	if (IsZ(Ha.norm())) {

		this->energy = 0;
		return 0.0;
	}

	double energy = 0;
	
	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double cHA = pMesh->cHA;
			pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

			pMesh->Heff[idx] += (cHA * Ha);

			energy += pMesh->M[idx] * (cHA * Ha);
		}
	}
	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double cHA = pMesh->cHA;
			pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

			pMesh->Heff[idx] += (cHA * Ha);
			pMesh->Heff2[idx] += (cHA * Ha);

			energy += (pMesh->M[idx] + pMesh->M2[idx]) * (cHA * Ha) / 2;
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / pMesh->M.get_nonempty_cells();
	else energy = 0;

	this->energy = energy;

	return energy;
}

//----------------------------------------------- Others

void Zeeman::SetField(DBL3 Hxyz)
{
	Ha = Hxyz;

	//if atomic_moment is not zero then changing the applied field also changes the temperature dependence of me (normalised equilibrium magnetisation). This in turn affects a number of material parameters.
	//Note, if only using LLG then set atomic_moment= 0 as calling SetCurieTemperature every time the applied field changes can be costly (if the field changes very often)

	double atomic_moment_ub = pMesh->GetAtomicMoment();

	if (IsNZ(atomic_moment_ub) && IsNZ(pMesh->GetCurieTemperature())) {

		//calling SetCurieTemperature forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten
		pMesh->SetCurieTemperature(pMesh->GetCurieTemperature());
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) reinterpret_cast<ZeemanCUDA*>(pModuleCUDA)->SetField(Ha);
#endif
}

DBL3 Zeeman::GetField(void)
{ 
#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<ZeemanCUDA*>(pModuleCUDA)->GetField();
#endif

	return Ha;
}

#endif