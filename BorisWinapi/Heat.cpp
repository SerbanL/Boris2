#include "stdafx.h"
#include "Heat.h"

#ifdef MODULE_HEAT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

Heat::Heat(Mesh *pMesh_) :
	Modules(),
	Q_equation({ "x", "y", "z", "t" }),
	ProgramStateNames(this, 
		{ 
			VINFO(T_ambient), VINFO(alpha_boundary), 
			VINFO(insulate_px), VINFO(insulate_nx), VINFO(insulate_py), VINFO(insulate_ny), VINFO(insulate_pz), VINFO(insulate_nz),
			VINFO(Q_equation)
		}, {})
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

Heat::~Heat()
{
	//free memory used for heat solver
	pMesh->Temp.clear();
}

//-------------------Abstract base class method implementations

BError Heat::Initialize(void)
{
	BError error(CLASS_STR(Heat));

	initialized = true;

	SetRobinBoundaryConditions();

	return error;
}

BError Heat::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Heat));

	Uninitialize();

	bool success = true;

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {

		//update mesh dimensions in equation constants
		if (Q_equation.is_set()) {

			DBL3 meshDim = pMesh->GetMeshDimensions();

			Q_equation.set_constant("Lx", meshDim.x, false);
			Q_equation.set_constant("Ly", meshDim.y, false);
			Q_equation.set_constant("Lz", meshDim.z, false);
			Q_equation.remake_equation();
		}

		//make sure the cellsize divides the mesh rectangle
		pMesh->n_t = round(pMesh->meshRect / pMesh->h_t);
		if (pMesh->n_t.x < 2) pMesh->n_t.x = 2;
		if (pMesh->n_t.y < 2) pMesh->n_t.y = 2;
		if (pMesh->n_t.z < 2) pMesh->n_t.z = 2;
		pMesh->h_t = pMesh->meshRect / pMesh->n_t;

		//make sure correct memory is assigned for heat solver quantities

		if (pMesh->M.linear_size()) {

			//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
			if (pMesh->Temp.linear_size()) {

				success = pMesh->Temp.resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMesh->Temp.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature, pMesh->M);
			}
		}
		else if (pMesh->elC.linear_size()) {

			//in a normal metal the electrical conductivity sets the shape

			//before doing this must make sure elC was itself set in the Transport module (it could be this module is being updated before the Transport module)
			//Transport module is guaranteed to be set otherwise elC would have zero size - it does mean Transport has UpdateConfiguration called twice but it doesn't matter.
			if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

				error = (*pMesh)(MOD_TRANSPORT)->UpdateConfiguration(cfgMessage);
				if (error) return error;
			}

			if (success) {

				if (pMesh->Temp.linear_size()) {

					success = pMesh->Temp.resize(pMesh->h_t, pMesh->meshRect, pMesh->elC);
				}
				else {

					success = pMesh->Temp.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature, pMesh->elC);
				}
			}
		}
		else {

			//in an insulator (or conductor with Transport module not set) the shape is set directly in Temp
			if (pMesh->Temp.linear_size()) {

				success = pMesh->Temp.resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMesh->Temp.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature);
			}
		}

		//allocate memory for the heatEq_RHS auxiliary vector
		if (success) {

			success = malloc_vector(heatEq_RHS, pMesh->n_t.dim(), 0.0);
		}
	}

	if (!success) return error(BERROR_OUTOFMEMORY_CRIT);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void Heat::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		UpdateTEquationUserConstants(false);
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (Q_equation.is_set()) Q_equation.clear();
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError Heat::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Heat));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new HeatCUDA(pMesh, pSMesh, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Heat::UpdateField(void)
{
	//no contribution to total energy density
	return 0.0;
}

//-------------------Calculation Methods

void Heat::IterateHeatEquation(double dT)
{
	//FTCS:

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

			if (!pMesh->Temp.is_not_empty(idx) || !pMesh->Temp.is_not_cmbnd(idx)) continue;

			double density = pMesh->density;
			double shc = pMesh->shc;
			double thermCond = pMesh->thermCond;
			pMesh->update_parameters_tcoarse(idx, pMesh->density, density, pMesh->shc, shc, pMesh->thermCond, thermCond);

			double cro = density * shc;
			double K = thermCond;

			//heat equation with Robin boundaries (based on Newton's law of cooling)
			heatEq_RHS[idx] = pMesh->Temp.delsq_robin(idx, K) * K / cro;

			//add Joule heating if set
			if (pMesh->E.linear_size()) {

				DBL3 position = pMesh->Temp.cellidx_to_position(idx);

				double elC_value = pMesh->elC.weighted_average(position, pMesh->Temp.h);
				DBL3 E_value = pMesh->E.weighted_average(position, pMesh->Temp.h);

				//add Joule heating source term
				heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro;
			}

			//add heat source contribution if set
			if (IsNZ(pMesh->Q.get0())) {

				double Q = pMesh->Q;
				pMesh->update_parameters_tcoarse(idx, pMesh->Q, Q);

				heatEq_RHS[idx] += Q / cro;
			}
		}
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		double time = pSMesh->GetStageTime();

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		for (int j = 0; j < pMesh->n_t.y; j++) {
			for (int k = 0; k < pMesh->n_t.z; k++) {
				for (int i = 0; i < pMesh->n_t.x; i++) {

					int idx = i + j * pMesh->n_t.x + k * pMesh->n_t.x*pMesh->n_t.y;

					if (!pMesh->Temp.is_not_empty(idx) || !pMesh->Temp.is_not_cmbnd(idx)) continue;

					double density = pMesh->density;
					double shc = pMesh->shc;
					double thermCond = pMesh->thermCond;
					pMesh->update_parameters_tcoarse(idx, pMesh->density, density, pMesh->shc, shc, pMesh->thermCond, thermCond);

					double cro = density * shc;
					double K = thermCond;

					//heat equation with Robin boundaries (based on Newton's law of cooling)
					heatEq_RHS[idx] = pMesh->Temp.delsq_robin(idx, K) * K / cro;

					//add Joule heating if set
					if (pMesh->E.linear_size()) {

						DBL3 position = pMesh->Temp.cellidx_to_position(idx);

						double elC_value = pMesh->elC.weighted_average(position, pMesh->Temp.h);
						DBL3 E_value = pMesh->E.weighted_average(position, pMesh->Temp.h);

						//add Joule heating source term
						heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro;
					}

					//add heat source contribution
					DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h_t;
					double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

					heatEq_RHS[idx] += Q / cro;
				}
			}
		}
	}

	//2. Now use forward time to advance by dT:
	#pragma omp parallel for
	for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

		pMesh->Temp[idx] += dT * heatEq_RHS[idx];
	}

	//Can 1 and 2 be combined into 1 step so heatEq_RHS is not needed? Seems it shouldn't work correctly - TO DO : investigate.
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of temperature and heat flux

//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
double Heat::afunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

double Heat::Heat::afunc_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

double Heat::bfunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	double thermCond = pMesh->thermCond;
	pMesh->update_parameters_atposition(relpos_m1, pMesh->thermCond, thermCond);

	return -1.0 * thermCond;
}

double Heat::bfunc_pri(int cell1_idx, int cell2_idx) const
{
	double thermCond = pMesh->thermCond;
	pMesh->update_parameters_tcoarse(cell1_idx, pMesh->thermCond, thermCond);

	return -1.0 * thermCond;
}

//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC
double Heat::diff2_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	double thermCond = pMesh->thermCond;

	if (pMesh->E.linear_size() || IsNZ(pMesh->Q.get0())) {

		pMesh->update_parameters_atposition(relpos_m1, pMesh->thermCond, thermCond);
	}
	else return 0.0;

	double value = 0.0;

	//Joule heating
	if (pMesh->E.linear_size()) {

		double elCval = pMesh->elC.weighted_average(relpos_m1, stencil);
		DBL3 Eval = pMesh->E.weighted_average(relpos_m1, stencil);
		value = -(elCval * Eval * Eval) / thermCond;
	}
	
	//heat source contribution if set
	if (IsNZ(pMesh->Q.get0())) {

		double Q = pMesh->Q;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->Q, Q);

		value -= Q / thermCond;
	}
	
	return value;
}

double Heat::diff2_pri(int cell1_idx, DBL3 shift) const
{
	double thermCond = pMesh->thermCond;

	if (pMesh->E.linear_size() || IsNZ(pMesh->Q.get0())) {

		pMesh->update_parameters_tcoarse(cell1_idx, pMesh->thermCond, thermCond);
	}
	else return 0.0;

	double value = 0.0;

	//Joule heating
	if (pMesh->E.linear_size()) {

		int idx1_E = pMesh->E.position_to_cellidx(pMesh->Temp.cellidx_to_position(cell1_idx));
		value = -(pMesh->elC[idx1_E] * pMesh->E[idx1_E] * pMesh->E[idx1_E]) / thermCond;
	}

	//heat source contribution if set
	if (IsNZ(pMesh->Q.get0())) {

		double Q = pMesh->Q;
		pMesh->update_parameters_tcoarse(cell1_idx, pMesh->Q, Q);

		value -= Q / thermCond;
	}

	return value;
}

//-------------------Setters

void Heat::SetAmbientTemperature(double T_ambient_)
{
	if (IsZoP(T_ambient))
		T_ambient = T_ambient_;

	SetRobinBoundaryConditions();
}

void Heat::SetAlphaBoundary(double alpha_boundary_)
{
	if (IsZoP(alpha_boundary_))
		alpha_boundary = alpha_boundary_;

	SetRobinBoundaryConditions();
}

//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
void Heat::SetBaseTemperature(double Temperature)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!pMesh->cT.is_sdep()) {

			reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetBaseTemperature(Temperature);
		}
		else {

			reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetBaseTemperature_Nonuniform(Temperature);
		}

		return;
	}
#endif

	if (!pMesh->cT.is_sdep()) {

		//uniform temperature
		pMesh->Temp.setnonempty(Temperature);
	}
	else {

		//non-uniform temperature setting
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

			if (pMesh->Temp.is_not_empty(idx)) {

				double cT = pMesh->cT;
				pMesh->update_parameters_tcoarse(idx, pMesh->cT, cT);

				pMesh->Temp[idx] = cT * Temperature;
			}
		}
	}
}

//Set Q_equation text equation object
BError Heat::SetQEquation(string equation_string, int step)
{
	BError error(CLASS_STR(Heat));

	DBL3 meshDim = pMesh->GetMeshDimensions();

	//set equation if not already set, or this is the first step in a stage
	if (!Q_equation.is_set() || step == 0) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants(false);

		if (!Q_equation.make_from_string(equation_string, { {"Lx", meshDim.x}, {"Ly", meshDim.y}, {"Lz", meshDim.z}, {"Ss", 0} })) return error(BERROR_INCORRECTSTRING);
	}
	else {

		//equation set and not the first step : adjust Ss constant
		Q_equation.set_constant("Ss", step);
		UpdateTEquationUserConstants(false);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) error = reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetQEquation(Q_equation.get_scalar_fspec());
#endif

	return error;
}

//Update TEquation object with user constants values
void Heat::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMesh->userConstants.size()) {

		vector<pair<string, double>> constants(pMesh->userConstants.size());
		for (int idx = 0; idx < pMesh->userConstants.size(); idx++) {

			constants[idx] = { pMesh->userConstants.get_key_from_index(idx), pMesh->userConstants[idx] };
		}

		Q_equation.set_constants(constants);

		//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
		if (pModuleCUDA && makeCuda) reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetQEquation(Q_equation.get_scalar_fspec());
#endif
	}
}

//set insulating mesh sides flags. str can be "x", "-x", "y", "-y", "z", "-z"
void Heat::SetInsulatingSides(string literal, bool status)
{
	if (literal == "x") insulate_px = status;
	else if (literal == "-x") insulate_nx = status;
	else if (literal == "y") insulate_py = status;
	else if (literal == "-y") insulate_ny = status;
	else if (literal == "z") insulate_pz = status;
	else if (literal == "-z") insulate_nz = status;

	SetRobinBoundaryConditions();
}

bool Heat::GetInsulatingSide(string literal)
{
	if (literal == "x") return insulate_px;
	else if (literal == "-x") return insulate_nx;
	else if (literal == "y") return insulate_py;
	else if (literal == "-y") return insulate_ny;
	else if (literal == "z") return insulate_pz;
	else if (literal == "-z") return insulate_nz;

	return false;
}

//-------------------Others

//called by MoveMesh method in this mesh - move relevant transport quantities
void Heat::MoveMesh_Heat(double x_shift)
{
	double mesh_end_size = pMesh->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(pMesh->meshRect.s + DBL3(mesh_end_size, 0, 0), pMesh->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->Temp()->shift_x(pMesh->Temp.linear_size(), x_shift, shift_rect);

		return;
	}
#endif

	//1. electrical conductivity
	pMesh->Temp.shift_x(x_shift, shift_rect);
}

void Heat::SetRobinBoundaryConditions(void)
{
	//set Robin boundary conditions
	DBL2 robin_values_px = DBL2(alpha_boundary, T_ambient) * (1 - insulate_px);
	DBL2 robin_values_nx = DBL2(alpha_boundary, T_ambient) * (1 - insulate_nx);
	DBL2 robin_values_py = DBL2(alpha_boundary, T_ambient) * (1 - insulate_py);
	DBL2 robin_values_ny = DBL2(alpha_boundary, T_ambient) * (1 - insulate_ny);
	DBL2 robin_values_pz = DBL2(alpha_boundary, T_ambient) * (1 - insulate_pz);
	DBL2 robin_values_nz = DBL2(alpha_boundary, T_ambient) * (1 - insulate_nz);

	pMesh->Temp.set_robin_conditions(DBL2(alpha_boundary, T_ambient), robin_values_px, robin_values_nx, robin_values_py, robin_values_ny, robin_values_pz, robin_values_nz);

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->Temp()->set_robin_conditions(cuReal2(alpha_boundary, T_ambient), robin_values_px, robin_values_nx, robin_values_py, robin_values_ny, robin_values_pz, robin_values_nz);
	}
#endif
}

#endif