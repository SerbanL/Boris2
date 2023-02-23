#include "stdafx.h"
#include "Zeeman.h"
#include "OVF2_Handlers.h"

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
	ProgramStateNames(this, { VINFO(Ha), VINFO(H_equation), VINFO(Havec) }, {})
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

//setup globalField transfer
BError Zeeman::InitializeGlobalField(void)
{
	BError error(__FUNCTION__);

	if (pSMesh->GetGlobalField().linear_size()) {

		//globalField here has to have same mesh rectangle and discretization as Heff.
		//we could initialize mesh transfer directly in Heff, but better not as it could be used for something else in the future
		if (!globalField.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_NCRIT);

		std::vector< VEC<DBL3>* > pVal_from;
		std::vector< VEC<DBL3>* > pVal_to;
		pVal_from.push_back(&(pSMesh->GetGlobalField()));

		if (!globalField.Initialize_MeshTransfer(pVal_from, pVal_to, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);
		globalField.transfer_in();
	}
	else globalField.clear();

	return error;
}

BError Zeeman::Initialize(void)
{
	BError error(CLASS_STR(Zeeman));

	//If using Havec make sure size and resolution matches M
	if (Havec.linear_size() && (Havec.size() != pMesh->M.size())) {
		if (!Havec.resize(pMesh->h, pMesh->meshRect)) {

			Havec.clear();
			error(BERROR_OUTOFMEMORY_NCRIT);
			initialized = false;
		}
	}

	//if using global field, then initialize mesh transfer if needed
	error = InitializeGlobalField();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_ZEEMAN || pMesh->IsOutputDataSet_withRect(DATA_E_ZEE),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_ZEEMAN || pMesh->IsOutputDataSet_withRect(DATA_E_ZEE),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError Zeeman::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Zeeman));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_SMESH_GLOBALFIELD)) {

		Uninitialize();

		//update mesh dimensions in equation constants
		if (H_equation.is_set()) {

			DBL3 meshDim = pMesh->GetMeshDimensions();

			H_equation.set_constant("Lx", meshDim.x, false);
			H_equation.set_constant("Ly", meshDim.y, false);
			H_equation.set_constant("Lz", meshDim.z, false);
			H_equation.remake_equation();
		}

		//if global field not set, then also clear it here
		if (!pSMesh->GetGlobalField().linear_size()) globalField.clear();
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
	double energy = 0;

	if (!H_equation.is_set()) {

		if (Havec.linear_size()) {

			/////////////////////////////////////////
			// Field VEC set
			/////////////////////////////////////////

			if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
				for (int idx = 0; idx < pMesh->n.dim(); idx++) {

					DBL3 Hext = Havec[idx];
					if (globalField.linear_size()) Hext += globalField[idx];

					pMesh->Heff[idx] = Hext;
					pMesh->Heff2[idx] = Hext;

					energy += (pMesh->M[idx] + pMesh->M2[idx]) * Hext / 2;

					if (Module_Heff.linear_size()) Module_Heff[idx] = Hext;
					if (Module_Heff2.linear_size()) Module_Heff2[idx] = Hext;
					if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * Hext;
					if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * pMesh->M2[idx] * Hext;
				}
			}

			else {

#pragma omp parallel for reduction(+:energy)
				for (int idx = 0; idx < pMesh->n.dim(); idx++) {

					DBL3 Hext = Havec[idx];
					if (globalField.linear_size()) Hext += globalField[idx];

					pMesh->Heff[idx] = Hext;

					energy += pMesh->M[idx] * Hext;

					if (Module_Heff.linear_size()) Module_Heff[idx] = Hext;
					if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * Hext;
				}
			}
		}
		else {

			/////////////////////////////////////////
			// Fixed set field
			/////////////////////////////////////////

			if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
				for (int idx = 0; idx < pMesh->n.dim(); idx++) {

					double cHA = pMesh->cHA;
					pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

					DBL3 Hext = cHA * Ha;
					if (globalField.linear_size()) Hext += globalField[idx];

					pMesh->Heff[idx] = Hext;
					pMesh->Heff2[idx] = Hext;

					energy += (pMesh->M[idx] + pMesh->M2[idx]) * Hext / 2;

					if (Module_Heff.linear_size()) Module_Heff[idx] = Hext;
					if (Module_Heff2.linear_size()) Module_Heff2[idx] = Hext;
					if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * Hext;
					if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * pMesh->M2[idx] * Hext;
				}
			}

			else {

#pragma omp parallel for reduction(+:energy)
				for (int idx = 0; idx < pMesh->n.dim(); idx++) {

					double cHA = pMesh->cHA;
					pMesh->update_parameters_mcoarse(idx, pMesh->cHA, cHA);

					DBL3 Hext = cHA * Ha;
					if (globalField.linear_size()) Hext += globalField[idx];

					pMesh->Heff[idx] = Hext;

					energy += pMesh->M[idx] * Hext;

					if (Module_Heff.linear_size()) Module_Heff[idx] = Hext;
					if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * Hext;
				}
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

						DBL3 Hext = cHA * H;
						if (globalField.linear_size()) Hext += globalField[idx];

						pMesh->Heff[idx] = Hext;
						pMesh->Heff2[idx] = Hext;

						energy += (pMesh->M[idx] + pMesh->M2[idx]) * Hext / 2;

						if (Module_Heff.linear_size()) Module_Heff[idx] = Hext;
						if (Module_Heff2.linear_size()) Module_Heff2[idx] = Hext;
						if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * Hext;
						if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * pMesh->M2[idx] * Hext;
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

						DBL3 Hext = cHA * H;
						if (globalField.linear_size()) Hext += globalField[idx];

						pMesh->Heff[idx] = Hext;

						energy += pMesh->M[idx] * Hext;

						if (Module_Heff.linear_size()) Module_Heff[idx] = Hext;
						if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * pMesh->M[idx] * Hext;
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

//-------------------Energy methods

//FM mesh
double Zeeman::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		if (!H_equation.is_set()) {

			if (Havec.linear_size()) {

				/////////////////////////////////////////
				// Field VEC set
				/////////////////////////////////////////

				DBL3 Hext = Havec[spin_index];
				if (globalField.linear_size()) Hext += globalField[spin_index];

				if (Mnew != DBL3()) return -pMesh->h.dim() * (Mnew - pMesh->M[spin_index]) * MU0 * Hext;
				else return -pMesh->h.dim() * pMesh->M[spin_index] * MU0 * Hext;
			}
			else {

				/////////////////////////////////////////
				// Fixed set field
				/////////////////////////////////////////

				if (IsZ(Ha.norm())) return 0.0;

				double cHA = pMesh->cHA;
				pMesh->update_parameters_mcoarse(spin_index, pMesh->cHA, cHA);

				DBL3 Hext = cHA * Ha;
				if (globalField.linear_size()) Hext += globalField[spin_index];

				if (Mnew != DBL3()) return -pMesh->h.dim() * (Mnew - pMesh->M[spin_index]) * MU0 * Hext;
				else return -pMesh->h.dim() * pMesh->M[spin_index] * MU0 * Hext;
			}
		}

		/////////////////////////////////////////
		// Field set from user equation
		/////////////////////////////////////////

		else {

			//on top of spatial dependence specified through an equation, also allow spatial dependence through the cHA parameter
			double cHA = pMesh->cHA;
			pMesh->update_parameters_mcoarse(spin_index, pMesh->cHA, cHA);

			DBL3 relpos = pMesh->M.cellidx_to_position(spin_index);
			DBL3 H = H_equation.evaluate_vector(relpos.x, relpos.y, relpos.z, pSMesh->GetStageTime());

			DBL3 Hext = cHA * H;
			if (globalField.linear_size()) Hext += globalField[spin_index];

			if (Mnew != DBL3()) return -pMesh->h.dim() * (Mnew - pMesh->M[spin_index]) * MU0 * Hext;
			else return -pMesh->h.dim() * pMesh->M[spin_index] * MU0 * Hext;
		}
	}
	else return 0.0;
}

//AFM mesh
DBL2 Zeeman::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		if (!H_equation.is_set()) {

			if (Havec.linear_size()) {

				/////////////////////////////////////////
				// Field VEC set
				/////////////////////////////////////////

				DBL3 Hext = Havec[spin_index];
				if (globalField.linear_size()) Hext += globalField[spin_index];

				if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

					return -pMesh->h.dim() * DBL2((Mnew_A - pMesh->M[spin_index]) * MU0 * Hext, (Mnew_B - pMesh->M2[spin_index]) * MU0 * Hext);
				}
				else return -pMesh->h.dim() * DBL2(pMesh->M[spin_index] * MU0 * Hext, pMesh->M2[spin_index] * MU0 * Hext);
			}
			else {

				/////////////////////////////////////////
				// Fixed set field
				/////////////////////////////////////////

				if (IsZ(Ha.norm())) return 0.0;

				double cHA = pMesh->cHA;
				pMesh->update_parameters_mcoarse(spin_index, pMesh->cHA, cHA);

				DBL3 Hext = cHA * Ha;
				if (globalField.linear_size()) Hext += globalField[spin_index];

				if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

					return -pMesh->h.dim() * DBL2((Mnew_A - pMesh->M[spin_index]) * MU0 * Hext, (Mnew_B - pMesh->M2[spin_index]) * MU0 * Hext);
				}
				else return -pMesh->h.dim() * DBL2(pMesh->M[spin_index] * MU0 * Hext, pMesh->M2[spin_index] * MU0 * Hext);
			}
		}

		/////////////////////////////////////////
		// Field set from user equation
		/////////////////////////////////////////

		else {

			//on top of spatial dependence specified through an equation, also allow spatial dependence through the cHA parameter
			double cHA = pMesh->cHA;
			pMesh->update_parameters_mcoarse(spin_index, pMesh->cHA, cHA);

			DBL3 relpos = pMesh->M.cellidx_to_position(spin_index);
			DBL3 H = H_equation.evaluate_vector(relpos.x, relpos.y, relpos.z, pSMesh->GetStageTime());

			DBL3 Hext = cHA * H;
			if (globalField.linear_size()) Hext += globalField[spin_index];

			if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

				return -pMesh->h.dim() * DBL2((Mnew_A - pMesh->M[spin_index]) * MU0 * Hext, (Mnew_B - pMesh->M2[spin_index]) * MU0 * Hext);
			}
			else return -pMesh->h.dim() * DBL2(pMesh->M[spin_index] * MU0 * Hext, pMesh->M2[spin_index] * MU0 * Hext);
		}
	}
	else return DBL2();
}

//----------------------------------------------- Others

void Zeeman::SetField(DBL3 Hxyz)
{
	//fixed field is being set - remove any equation settings
	if (H_equation.is_set()) H_equation.clear();
	//also release any memory in field VEC
	if (Havec.linear_size()) Havec.clear();

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

BError Zeeman::SetFieldEquation(std::string equation_string, int step)
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

	//also release any memory in field VEC
	if (Havec.linear_size()) Havec.clear();

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) error = dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->SetFieldEquation(H_equation.get_vector_fspec());
#endif

	return error;
}

BError Zeeman::SetFieldVEC_FromOVF2(std::string fileName)
{
	BError error(CLASS_STR(Zeeman));

	//Load data from file
	OVF2 ovf2;
	error = ovf2.Read_OVF2_VEC(fileName, Havec);
	if (error) {

		Havec.clear();
		return error;
	}

	//Make sure size and resolution matches M
	if (Havec.size() != pMesh->M.size()) {
		if (!Havec.resize(pMesh->h, pMesh->meshRect)) {

			Havec.clear();
			return error(BERROR_OUTOFMEMORY_NCRIT);
		}
	}

	//all good, clear any equation settings
	if (H_equation.is_set()) H_equation.clear();

	//if displaying module effective field also need to update these
	if (Module_Heff.linear_size()) Module_Heff.copy_values(Havec);
	if (Module_Heff2.linear_size()) Module_Heff2.copy_values(Havec);

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) error = dynamic_cast<ZeemanCUDA*>(pModuleCUDA)->SetFieldVEC(Havec);
#endif

	return error;
}

//Update TEquation object with user constants values
void Zeeman::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMesh->userConstants.size()) {

		std::vector<std::pair<std::string, double>> constants(pMesh->userConstants.size());
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

//-------------------Torque methods

DBL3 Zeeman::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif
