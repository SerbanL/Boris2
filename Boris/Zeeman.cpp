#include "stdafx.h"
#include "Zeeman.h"

#ifdef MODULE_COMPILATION_ZEEMAN

#include "Mesh.h"
#include "MeshParamsControl.h"

#include "SuperMesh.h"

#if COMPILECUDA == 1
#include "ZeemanCUDA.h"
#endif

Zeeman::Zeeman(Mesh *pMesh_) : 
	Modules(),
	ZeemanBase(),
	ProgramStateNames(this, { VINFO(Ha), VINFO(H_equation) }, {})
{
	pMesh = pMesh_;

	pSMesh = pMesh->pSMesh;

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

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();

		//update mesh dimensions in equation constants
		if (H_equation.is_set()) {

			DBL3 meshDim = pMesh->GetMeshDimensions();

			H_equation.set_constant("Lx", meshDim.x, false);
			H_equation.set_constant("Ly", meshDim.y, false);
			H_equation.set_constant("Lz", meshDim.z, false);
			H_equation.remake_equation();
		}

		Initialize();
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void Zeeman::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		UpdateTEquationUserConstants(false);
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (H_equation.is_set()) H_equation.clear();
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
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
	/////////////////////////////////////////
	// Fixed set field
	/////////////////////////////////////////

	double energy = 0;

	if (!H_equation.is_set()) {

		if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
			for (int idx = 0; idx < pMesh->n.dim(); idx++) {

				double cHA = pMesh->cHA;
				pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

				pMesh->Heff[idx] = (cHA * Ha);
				pMesh->Heff2[idx] = (cHA * Ha);

				energy += (pMesh->M[idx] + pMesh->M2[idx]) * (cHA * Ha) / 2;
			}
		}

		else {

#pragma omp parallel for reduction(+:energy)
			for (int idx = 0; idx < pMesh->n.dim(); idx++) {

				double cHA = pMesh->cHA;
				pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

				pMesh->Heff[idx] = (cHA * Ha);

				energy += pMesh->M[idx] * (cHA * Ha);
			}
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		double time = pSMesh->GetStageTime();

		if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
			for (int j = 0; j < pMesh->n.y; j++) {
				for (int k = 0; k < pMesh->n.z; k++) {
					for (int i = 0; i < pMesh->n.x; i++) {

						int idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

						//on top of spatial dependence specified through an equation, also allow spatial dependence through the cHA parameter
						double cHA = pMesh->cHA;
						pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

						DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h;
						DBL3 H = H_equation.evaluate_vector(relpos.x, relpos.y, relpos.z, time);

						pMesh->Heff[idx] = (cHA * H);
						pMesh->Heff2[idx] = (cHA * H);

						energy += (pMesh->M[idx] + pMesh->M2[idx]) * (cHA * H) / 2;
					}
				}
			}
		}

		else {

#pragma omp parallel for reduction(+:energy)
			for (int j = 0; j < pMesh->n.y; j++) {
				for (int k = 0; k < pMesh->n.z; k++) {
					for (int i = 0; i < pMesh->n.x; i++) {

						int idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

						//on top of spatial dependence specified through an equation, also allow spatial dependence through the cHA parameter
						double cHA = pMesh->cHA;
						pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

						DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h;
						DBL3 H = H_equation.evaluate_vector(relpos.x, relpos.y, relpos.z, time);

						pMesh->Heff[idx] = (cHA * H);

						energy += pMesh->M[idx] * (cHA * H);
					}
				}
			}
		}
	}

	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / pMesh->M.get_nonempty_cells();
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//-------------------Energy density methods

double Zeeman::GetEnergyDensity(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->GetEnergyDensity(avRect);
#endif

	/////////////////////////////////////////
	// Fixed set field
	/////////////////////////////////////////

	double energy = 0;

	int num_points = 0;

	if (!H_equation.is_set()) {

		if (IsZ(Ha.norm())) {

			this->energy = 0;
			return 0.0;
		}

		if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy, num_points)
			for (int idx = 0; idx < pMesh->n.dim(); idx++) {

				//only average over values in given rectangle
				if (!avRect.contains(pMesh->M.cellidx_to_position(idx))) continue;

				double cHA = pMesh->cHA;
				pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

				energy += (pMesh->M[idx] + pMesh->M2[idx]) * (cHA * Ha) / 2;
				num_points++;
			}
		}

		else {

#pragma omp parallel for reduction(+:energy, num_points)
			for (int idx = 0; idx < pMesh->n.dim(); idx++) {

				//only average over values in given rectangle
				if (!avRect.contains(pMesh->M.cellidx_to_position(idx))) continue;

				double cHA = pMesh->cHA;
				pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

				energy += pMesh->M[idx] * (cHA * Ha);
				num_points++;
			}
		}
	}

	/////////////////////////////////////////
	// Field set from user equation
	/////////////////////////////////////////

	else {

		double time = pSMesh->GetStageTime();

		if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy, num_points)
			for (int j = 0; j < pMesh->n.y; j++) {
				for (int k = 0; k < pMesh->n.z; k++) {
					for (int i = 0; i < pMesh->n.x; i++) {

						int idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

						//on top of spatial dependence specified through an equation, also allow spatial dependence through the cHA parameter
						double cHA = pMesh->cHA;
						pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

						DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h;
						DBL3 H = H_equation.evaluate_vector(relpos.x, relpos.y, relpos.z, time);

						energy += (pMesh->M[idx] + pMesh->M2[idx]) * (cHA * H) / 2;
						num_points++;
					}
				}
			}
		}

		else {

#pragma omp parallel for reduction(+:energy, num_points)
			for (int j = 0; j < pMesh->n.y; j++) {
				for (int k = 0; k < pMesh->n.z; k++) {
					for (int i = 0; i < pMesh->n.x; i++) {

						int idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

						//on top of spatial dependence specified through an equation, also allow spatial dependence through the cHA parameter
						double cHA = pMesh->cHA;
						pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

						DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h;
						DBL3 H = H_equation.evaluate_vector(relpos.x, relpos.y, relpos.z, time);

						energy += pMesh->M[idx] * (cHA * H);
						num_points++;
					}
				}
			}
		}
	}

	if (num_points) energy *= -MU0 / num_points;
	else energy = 0;

	return energy;
}

//----------------------------------------------- Others

void Zeeman::SetField(DBL3 Hxyz)
{
	//fixed field is being set - remove any equation settings
	if (H_equation.is_set()) H_equation.clear();

	Ha = Hxyz;

	//if atomic_moment is not zero then changing the applied field also changes the temperature dependence of me (normalised equilibrium magnetization). This in turn affects a number of material parameters.
	//Note, if only using LLG then set atomic_moment= 0 as calling SetCurieTemperature every time the applied field changes can be costly (if the field changes very often)

	bool atomic_moment_zero = false;
	if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) atomic_moment_zero = IsZ(pMesh->GetAtomicMoment_AFM().i + pMesh->GetAtomicMoment_AFM().j);
	else if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) atomic_moment_zero = IsZ(pMesh->GetAtomicMoment());

	if (!atomic_moment_zero && IsNZ(pMesh->GetCurieTemperature())) {

		//calling SetCurieTemperature forces recalculation of affected material parameters temperature dependence - any custom dependence set will be overwritten
		pMesh->SetCurieTemperature(pMesh->GetCurieTemperature(), false);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->SetField(Ha);
#endif
}

DBL3 Zeeman::GetField(void) 
{ 
	if (H_equation.is_set()) {

		DBL3 meshDim = pMesh->GetMeshDimensions();

		return H_equation.evaluate_vector(meshDim.x / 2, meshDim.y / 2, meshDim.z / 2, pSMesh->GetStageTime());
	}
	else return Ha; 
}

BError Zeeman::SetFieldEquation(string equation_string, int step)
{
	BError error(CLASS_STR(Zeeman));

	DBL3 meshDim = pMesh->GetMeshDimensions();

	//set equation if not already set, or this is the first step in a stage
	if (!H_equation.is_set() || step == 0) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants(false);

		if (!H_equation.make_from_string(equation_string, { {"Lx", meshDim.x}, {"Ly", meshDim.y}, {"Lz", meshDim.z}, {"Tb", pMesh->base_temperature}, {"Ss", step} })) return error(BERROR_INCORRECTSTRING);
	}
	else {

		//equation set and not the first step : adjust Ss constant
		H_equation.set_constant("Ss", step);
		UpdateTEquationUserConstants(false);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) error = dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->SetFieldEquation(H_equation.get_vector_fspec());
#endif

	return error;
}

//Update TEquation object with user constants values
void Zeeman::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMesh->userConstants.size()) {

		vector<pair<string, double>> constants(pMesh->userConstants.size());
		for (int idx = 0; idx < pMesh->userConstants.size(); idx++) {

			constants[idx] = { pMesh->userConstants.get_key_from_index(idx), pMesh->userConstants[idx] };
		}

		H_equation.set_constants(constants);

		//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
		if (pModuleCUDA && makeCuda) dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->SetFieldEquation(H_equation.get_vector_fspec());
#endif
	}
}

//if base temperature changes we need to adjust Tb in H_equation if it's used.
void Zeeman::SetBaseTemperature(double Temperature)
{
	if (H_equation.is_set()) {

		H_equation.set_constant("Tb", Temperature);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->SetFieldEquation(H_equation.get_vector_fspec());
#endif
}

#endif
