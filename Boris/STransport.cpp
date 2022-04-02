#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "SuperMesh.h"

#include "Transport.h"
#include "Atom_Transport.h"
#include "TMR.h"

#include "STransport_Graph.h"

STransport::STransport(SuperMesh *pSMesh_) :
	Modules(),
	V_equation({ "t" }),
	I_equation({ "t" }),
	ProgramStateNames(this, { VINFO(electrode_rects), VINFO(electrode_potentials), 
							  VINFO(ground_electrode_index), VINFO(potential), VINFO(current), VINFO(net_current), VINFO(resistance), VINFO(constant_current_source), 
							  VINFO(errorMaxLaplace), VINFO(maxLaplaceIterations), VINFO(s_errorMax), VINFO(s_maxIterations), VINFO(SOR_damping),
							  VINFO(V_equation), VINFO(I_equation) }, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

//-------------------Abstract base class method implementations

BError STransport::Initialize(void)
{
	BError error(CLASS_STR(STransport));

	if (!initialized) {

		if (!pSMesh->disabled_transport_solver) {

			////////////////////////////////////////////////////////////////////////////
			//Calculate V, E and elC before starting
			////////////////////////////////////////////////////////////////////////////

			//initialize V with a linear slope between ground and another electrode (in most problems there are only 2 electrodes setup) - do this for all transport meshes
			initialize_potential_values();

			//set electric field and  electrical conductivity in individual transport modules (in this order!)
			for (int idx = 0; idx < (int)pTransport.size(); idx++) {

				pTransport[idx]->CalculateElectricField();
				pTransport[idx]->CalculateElectricalConductivity(true);
			}

			//solve only for charge current (V and Jc with continuous boundaries)
			if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
			//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
			else solve_spin_transport_sor();

			recalculate_transport = true;
			transport_recalculated = true;
		}

		initialized = true;
	}

	return error;
}

BError STransport::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(STransport));
	
	Uninitialize();
	
	if (ucfg::check_cfgflags(cfgMessage,
		UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, 
		UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED,
		UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, 
		UPDATECONFIG_TRANSPORT_ELECTRODE, UPDATECONFIG_TRANSPORT,
		UPDATECONFIG_ODE_SOLVER)) {
		
		////////////////////////////////////////////////////////////////////////////
		//check meshes to set transport boundary flags (NF_CMBND flags for V)
		////////////////////////////////////////////////////////////////////////////

		//clear everything then rebuild
		pTransport.clear();
		CMBNDcontacts.clear();
		pV.clear();
		pS.clear();
		
		//now build pTransport (and pV)
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->IsModuleSet(MOD_TRANSPORT) || (*pSMesh)[idx]->IsModuleSet(MOD_TMR)) {

				if ((*pSMesh)[idx]->IsModuleSet(MOD_TRANSPORT)) {

					Modules* pModule = (*pSMesh)[idx]->GetModule(MOD_TRANSPORT);
					if (dynamic_cast<Transport*>(pModule)) pTransport.push_back(dynamic_cast<Transport*>(pModule));
#if ATOMISTIC == 1
					else if (dynamic_cast<Atom_Transport*>(pModule)) pTransport.push_back(dynamic_cast<Atom_Transport*>(pModule));
#endif
				}
				else if ((*pSMesh)[idx]->IsModuleSet(MOD_TMR)) {

					Modules* pModule = (*pSMesh)[idx]->GetModule(MOD_TMR);
					pTransport.push_back(dynamic_cast<TMR*>(pModule));
				}

				pV.push_back(&(*pSMesh)[idx]->V);
				pS.push_back(&(*pSMesh)[idx]->S);
			}
		}

		//set fixed potential cells and cmbnd flags (also building contacts)
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			//1. Dirichlet conditions

			//make sure all fixed potential cells are marked
			pTransport[idx]->ClearFixedPotentialCells();

			for (int el_idx = 0; el_idx < electrode_rects.size(); el_idx++) {

				if (!pTransport[idx]->SetFixedPotentialCells(electrode_rects[el_idx], electrode_potentials[el_idx])) return error(BERROR_OUTOFMEMORY_NCRIT);
			}

			//2. CMBND conditions

			//build CMBND contacts and set flags for V
			CMBNDcontacts.push_back(pV[idx]->set_cmbnd_flags(idx, pV));

			//set flags for S also (same mesh dimensions as V so CMBNDcontacts are the same)
			if (pSMesh->SolveSpinCurrent()) pS[idx]->set_cmbnd_flags(idx, pS);
		}
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void STransport::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		UpdateTEquationUserConstants();
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (V_equation.is_set()) V_equation.clear();
		if (I_equation.is_set()) I_equation.clear();
	}
}

BError STransport::MakeCUDAModule(void)
{
	BError error(CLASS_STR(STransport));

#if COMPILECUDA == 1

	pModuleCUDA = new STransportCUDA(pSMesh, this);
	pSTransportCUDA = dynamic_cast<STransportCUDA*>(pModuleCUDA);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

//-------------------

double STransport::UpdateField(void)
{
	if (pSMesh->disabled_transport_solver) return 0.0;

	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be iterated only at the end of a step or stage
	//however, we still want to compute self-consistent spin torques if SolveSpinCurrent()

	//only need to update this after an entire magnetization equation time step is solved (but always update spin accumulation field if spin current solver enabled)
	if (pSMesh->CurrentTimeStepSolved() && !pSMesh->static_transport_solver) {

		//use V or I equation to set electrode potentials? time dependence only
		if (V_equation.is_set() || I_equation.is_set()) {

			if (V_equation.is_set()) SetPotential(V_equation.evaluate(pSMesh->GetStageTime()), false);
			else SetCurrent(I_equation.evaluate(pSMesh->GetStageTime()), false);

			recalculate_transport = true;
		}

		transport_recalculated = recalculate_transport;

		if (recalculate_transport) {

			recalculate_transport = false;

			//solve only for charge current (V and Jc with continuous boundaries)
			if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
			//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
			else solve_spin_transport_sor();

			//if constant current source is set then need to update potential to keep a constant current
			if (constant_current_source) {

				GetCurrent();

				//the electrode voltage values will have changed so should iterate to convergence threshold again
				//even though we've adjusted potential values these won't be quite correct
				//moreover the charge current density hasn't been recalculated
				//reiterating the transport solver to convergence will fix all this
				//Note : the electrode current will change again slightly so really you should be iterating to some electrode current convergence threshold
				//In normal running mode this won't be an issue as this is done every iteration; in static transport solver mode this could be a problem so must be tested

				double iters_to_conv_previous = iters_to_conv;

				//solve only for charge current (V and Jc with continuous boundaries)
				if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
				//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
				else solve_spin_transport_sor();

				//in constant current mode we spend more iterations so the user should be aware of this
				iters_to_conv += iters_to_conv_previous;
			}
		}
		else iters_to_conv = 0;
	}
	
	if (pSMesh->SolveSpinCurrent()) {
		
		//Calculate the spin accumulation field so a torque is generated when used in the LLG (or LLB) equation
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			pTransport[idx]->CalculateSAField();
		}

		//Calculate effective field from interface spin accumulation torque (in magnetic meshes for NF interfaces with G interface conductance set)
		CalculateSAInterfaceField();
	}
	
	//no contribution to total energy density
	return 0.0;
}

//-------------------

//set fixed SOR damping values (for V and S solvers)
void STransport::SetSORDamping(DBL2 _SOR_damping)
{
	SOR_damping = _SOR_damping;

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pSTransportCUDA->SOR_damping_V.from_cpu(SOR_damping.i);
		pSTransportCUDA->SOR_damping_S.from_cpu(SOR_damping.j);
	}
#endif
}

//-------------------

DBL2 STransport::GetCurrent(void)
{ 
	double save_current = current;

	current = 0.0;
	net_current = 0.0;

	//get total current flowing into the ground electrode (if set) as well as net current and total number of electrodes
	for (int el_idx = 0; el_idx < electrode_rects.size(); el_idx++) {
		
		double electrode_current = 0.0;

		for (int idx = 0; idx < (int)pTransport.size(); idx++) {
			
			electrode_current += pTransport[idx]->CalculateElectrodeCurrent(electrode_rects[el_idx]);
		}

		net_current += electrode_current;
		if (ground_electrode_index == el_idx) current = electrode_current;
	}

	//adjust the ground current for net current error (net current should really be zero)
	if (electrode_rects.size()) current = current - net_current / electrode_rects.size();
	
	if (constant_current_source) {

		//adjust potential in order to keep the current constant
		if (IsNZ(current)) {

			adjust_potential(fabs(potential / current) * save_current);
		}

		current = save_current;
	}
	
	return DBL2(current, net_current);
}

double STransport::GetResistance(void)
{ 
	//calculate new current value for the set potential. If a constant current source is used, the potential is adjusted instead in order to keep the set current value.
	GetCurrent();

	//calculate resistance value
	if (IsNZ(current)) resistance = fabs(potential / current);
	else resistance = 0;
	
	return resistance; 
}

void STransport::SetPotential(double potential_, bool clear_equation)
{
	if (clear_equation) V_equation.clear();

	constant_current_source = false;

	adjust_potential(potential_);
}

void STransport::SetCurrent(double current_, bool clear_equation)
{
	if (clear_equation) I_equation.clear();

	//before setting a constant current, make sure the currently set potential is not zero otherwise we won't be able to set a constant current source
	if (IsZ(potential)) SetPotential(1.0);

	constant_current_source = true;

	current = current_;

	//get_current will now also adjust the potential depending on the sample resistance
	GetCurrent();
}

//set text equation from std::string
BError STransport::SetPotentialEquation(std::string equation_string, int step)
{
	BError error(CLASS_STR(STransport));

	//set equation if not already set, or this is the first step in a stage
	if (!V_equation.is_set() || step == 0) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants();

		if (!V_equation.make_from_string(equation_string, { {"Ss", 0} })) return error(BERROR_INCORRECTSTRING);
	}
	else {

		//equation set and not the first step : adjust Ss constant
		V_equation.set_constant("Ss", step);
		UpdateTEquationUserConstants();
	}

	//don't use both V and I equations at the same time
	I_equation.clear();

	return error;
}

//set text equation from std::string
BError STransport::SetCurrentEquation(std::string equation_string, int step)
{
	BError error(CLASS_STR(STransport));

	//set equation if not already set, or this is the first step in a stage
	if (!I_equation.is_set() || step == 0) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants();

		if (!I_equation.make_from_string(equation_string, { {"Ss", 0} })) return error(BERROR_INCORRECTSTRING);
	}
	else {

		//equation set and not the first step : adjust Ss constant
		I_equation.set_constant("Ss", step);
		UpdateTEquationUserConstants();
	}

	//don't use both V and I equations at the same time
	V_equation.clear();

	return error;
}

//Update TEquation object with user constants values
void STransport::UpdateTEquationUserConstants(void)
{
	if (pSMesh->userConstants.size()) {

		std::vector<std::pair<std::string, double>> constants(pSMesh->userConstants.size());
		for (int idx = 0; idx < pSMesh->userConstants.size(); idx++) {

			constants[idx] = { pSMesh->userConstants.get_key_from_index(idx), pSMesh->userConstants[idx] };
		}

		V_equation.set_constants(constants);
		I_equation.set_constants(constants);
	}
}

void STransport::adjust_potential(double potential_)
{
	double electrode_voltage = 0, ground_voltage = 0;

	for (int el_idx = 0; el_idx < (int)electrode_rects.size(); el_idx++) {

		if (el_idx == ground_electrode_index) {

			if (SetElectrodePotential(el_idx , -potential_ / 2)) {

				//found ground electrode and could set its value
				ground_voltage = -potential_ / 2;
			}
		}
		else if (SetElectrodePotential(el_idx , +potential_ / 2)) {

			//found at least another electrode and could set its value
			electrode_voltage = +potential_ / 2;
		}
	}
	
	//this is now the old potential value
	potential_ = potential;

	//the new potential value
	potential = electrode_voltage - ground_voltage;
	
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (IsNZ(potential_)) {

			pSTransportCUDA->scale_potential_values(potential / potential_);
		}
		else pSTransportCUDA->initialize_potential_values();

		return;
	}
#endif
	
	//if previous potential value was not zero then scale all V values - saves on transport solver computations
	if (IsNZ(potential_)) {

		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			pV[idx]->scale_values(potential / potential_);
		}
	}
	else initialize_potential_values();
}

//set potential values using a slope between the potential values of ground and another electrode (if set)
void STransport::initialize_potential_values(void)
{
	//Note, it's possible V already has values, e.g. we've just loaded a simulation file with V saved.
	//We don't want to re-initialize the V values as this will force the transport solver to iterate many times to get back the correct V values - which we already have!
	//Then, only apply the default V initialization if the voltage values are zero - if the average V is exactly zero (averaged over all meshes) then it's highly probable V is zero everywhere.
	//It could be that V has a perfectly anti-symmetrical set of values, in which case the average will also be zero. But in this case there's also no point to re-initialize the values.
	double V_average = 0;

	for (int idx = 0; idx < pV.size(); idx++) {

		V_average += pV[idx]->average_nonempty();
	}
	
	if (IsZ(V_average)) set_linear_potential_drops();
}

//auxiliary used by initialize_potential_values for setting potential drops in all required meshes
void STransport::set_linear_potential_drops(void)
{
	//need a ground electrode and at least 2 electrodes overall
	if (ground_electrode_index < 0 || electrode_rects.size() < 2) return;

	//mRects : available mesh rectangles; eRects : available electrode rectangles, other than ground
	std::vector<Rect> mRects(pTransport.size()), eRects;
	//ground electrode rectangle
	Rect gndRect = electrode_rects[ground_electrode_index];

	//ground and elecrode potentials
	std::vector<double> eV;
	double gndV = electrode_potentials[ground_electrode_index];

	//mesh conductivities
	std::vector<double> mCond(pTransport.size());

	//collect all relevant mesh rectangles, and make sure each is initialized and get average conductivities
	for (int idx = 0; idx < pTransport.size(); idx++) {

		//make sure modules are initialized before obtaining conductivities (as in some, e.g. TMR this needs to be computed)
		if (pTransport[idx]->pMeshBase->GetMeshType() == MESH_INSULATOR) {

#if COMPILECUDA == 1
			if (!dynamic_cast<TMR*>(pTransport[idx])->IsInitialized()) dynamic_cast<TMR*>(pTransport[idx])->InitializeCUDA();
#else
			if (!dynamic_cast<TMR*>(pTransport[idx])->IsInitialized()) dynamic_cast<TMR*>(pTransport[idx])->Initialize();
#endif
		}
		else {

#if COMPILECUDA == 1
			if (!dynamic_cast<Transport*>(pTransport[idx])->IsInitialized()) dynamic_cast<Transport*>(pTransport[idx])->InitializeCUDA();
#else
			if (!dynamic_cast<Transport*>(pTransport[idx])->IsInitialized()) dynamic_cast<Transport*>(pTransport[idx])->Initialize();
#endif
		}

		mRects[idx] = pTransport[idx]->pMeshBase->meshRect;
		mCond[idx] = pTransport[idx]->pMeshBase->GetAverageElectricalConductivity();
	}

	//collect electrode rectangles and potentials
	for (int idx = 0; idx < electrode_rects.size(); idx++) {

		if (idx != ground_electrode_index) {

			eRects.push_back(electrode_rects[idx]);
			eV.push_back(electrode_potentials[idx]);
		}
	}

	//make graph and calculate potentials for each node in each path
	Graph graph(mRects, eRects, gndRect);
	graph.calculate_potentials(mCond, eV, gndV);

	//traverse graph and set potentials
	for (int pidx = 0; pidx < graph.num_paths(); pidx++) {
		for (int nidx = 0; nidx < graph.path_size(pidx); nidx++) {

			int idx = graph.get_node_mesh_idx(pidx, nidx);
			graph.set_potential_drop<TransportBase>(pidx, nidx, &TransportBase::Set_Linear_PotentialDrop, *pTransport[idx]);
		}
	}
}

void STransport::SetConvergenceError(double errorMaxLaplace_, int maxLaplaceIterations_)
{
	errorMaxLaplace = errorMaxLaplace_;
	if (maxLaplaceIterations_ > 0) maxLaplaceIterations = maxLaplaceIterations_;
	else maxLaplaceIterations = 1000;

	recalculate_transport = true;
}

void STransport::SetSConvergenceError(double s_errorMax_, int s_maxIterations_)
{
	s_errorMax = s_errorMax_;
	if (s_maxIterations_ > 0) s_maxIterations = s_maxIterations_;
	else s_maxIterations = 1000;

	recalculate_transport = true;
}

//-------------------Electrodes methods

void STransport::AddElectrode(double electrode_potential, Rect electrode_rect)
{
	electrode_rects.push_back(electrode_rect);
	electrode_potentials.push_back(electrode_potential);

	UpdateConfiguration(UPDATECONFIG_TRANSPORT_ELECTRODE);
}

bool STransport::DelElectrode(int index)
{
	if (GoodIdx(electrode_rects.size(), index)) {

		electrode_rects.erase(index);
		electrode_potentials.erase(electrode_potentials.begin() + index);

		//have we deleted the ground electrode?
		if (ground_electrode_index == index) ground_electrode_index = -1;

		UpdateConfiguration(UPDATECONFIG_TRANSPORT_ELECTRODE);

		return true;
	}
	else return false;
}

void STransport::ClearElectrodes(void)
{
	electrode_rects.clear();
	electrode_potentials.clear();

	UpdateConfiguration(UPDATECONFIG_TRANSPORT_ELECTRODE);
}

void STransport::DesignateGroundElectrode(int index)
{
	if (GoodIdx(electrode_rects.size(), index)) {

		ground_electrode_index = index;
	}
}

bool STransport::SetElectrodeRect(int index, Rect electrode_rect)
{
	if (GoodIdx(electrode_rects.size(), index)) {

		electrode_rects[index] = electrode_rect;

		UpdateConfiguration(UPDATECONFIG_TRANSPORT_ELECTRODE);

		return true;
	}
	else return false;
}

bool STransport::SetElectrodePotential(int index, double electrode_potential)
{
	//NOTE : all potential value changes (including adjustments due to constant current source) eventually call this method. Set flag to recalculate transport.

	if (GoodIdx(electrode_rects.size(), index)) {

		electrode_potentials[index] = electrode_potential;

		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			pTransport[idx]->SetFixedPotentialCells(electrode_rects[index], electrode_potentials[index]);
		}

		recalculate_transport = true;

		return true;
	}
	else return false;
}

std::pair<Rect, double> STransport::GetElectrodeInfo(int index)
{
	std::pair<Rect, double> elInfo;

	if (GoodIdx(electrode_rects.size(), index)) {

		elInfo.first = electrode_rects[index];
		elInfo.second = electrode_potentials[index];
	}

	return elInfo;
}

int STransport::GetElectrodeid(int index)
{
	if (GoodIdx(electrode_rects.size(), index)) {

		return electrode_rects.get_id_from_index(index).minor;
	}

	return -1;
}

#endif
