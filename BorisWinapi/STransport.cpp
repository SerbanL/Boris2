#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_TRANSPORT

#include "SuperMesh.h"

STransport::STransport(SuperMesh *pSMesh_) :
	Modules(),
	ProgramStateNames(this, { VINFO(electrode_rects), VINFO(electrode_potentials), 
							  VINFO(ground_electrode_index), VINFO(potential), VINFO(current), VINFO(net_current), VINFO(resistance), VINFO(constant_current_source), 
							  VINFO(errorMaxLaplace), VINFO(maxLaplaceIterations), VINFO(s_errorMax), VINFO(s_maxIterations), VINFO(fixed_SOR_damping), VINFO(SOR_damping) }, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration();

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

		//Calculate V and Jc before starting

		//initialize V with a linear slope between ground and another electrode (in most problems there are only 2 electrodes setup) - do this for all transport meshes
		initialize_potential_values();

		//solve only for charge current (V and Jc with continuous boundaries)
		if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
		//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
		else solve_spin_transport_sor();

		recalculate_transport = true;
		transport_recalculated = true;

		initialized = true;
	}

	return error;
}

BError STransport::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(STransport));

	Uninitialize();
	
	//check meshes to set transport boundary flags (NF_CMBND flags for V)

	//clear everything then rebuild
	pTransport.clear();
	CMBNDcontacts.clear();
	pV.clear();
	pS.clear();

	//now build pTransport (and pV)
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->IsModuleSet(MOD_TRANSPORT)) {

			pTransport.push_back(dynamic_cast<Transport*>((*pSMesh)[idx]->GetModule(MOD_TRANSPORT)));
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

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
	}
#endif

	return error;
}

BError STransport::MakeCUDAModule(void)
{
	BError error(CLASS_STR(STransport));

#if COMPILECUDA == 1

	pModuleCUDA = new STransportCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

//-------------------

double STransport::UpdateField(void)
{
	//only need to update this after an entire magnetisation equation time step is solved (but always update spin accumulation field if spin current solver enabled)
	if (pSMesh->CurrentTimeStepSolved()) {

		transport_recalculated = recalculate_transport;

		if (recalculate_transport) {

			recalculate_transport = false;

			//solve only for charge current (V and Jc with continuous boundaries)
			if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
			//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
			else solve_spin_transport_sor();

			//if constant current source is set then need to update potential to keep a constant current
			if (constant_current_source) GetCurrent();
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

		reinterpret_cast<STransportCUDA*>(pModuleCUDA)->SOR_damping_V.from_cpu(SOR_damping.i);
		reinterpret_cast<STransportCUDA*>(pModuleCUDA)->SOR_damping_S.from_cpu(SOR_damping.j);
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

void STransport::SetPotential(double potential_)
{
	constant_current_source = false;

	adjust_potential(potential_);
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

		if (IsNZ(potential_) && initialized) {

			reinterpret_cast<STransportCUDA*>(pModuleCUDA)->scale_potential_values(potential / potential_);

			return;
		}
		else reinterpret_cast<STransportCUDA*>(pModuleCUDA)->initialize_potential_values();
	}
#endif

	//if previous potential value was not zero then scale all V values - saves on transport solver computations
	if (IsNZ(potential_) && initialized) {

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
	
	if (IsZ(V_average)) {

		//initialize V with a linear slope between ground and another electrode (in most problems there are only 2 electrodes setup) - do this for all transport meshes
		if (ground_electrode_index >= 0 && electrode_rects.size() >= 2) {

			DBL3 ground_electrode_center = electrode_rects[ground_electrode_index].get_c();
			double ground_potential = electrode_potentials[ground_electrode_index];

			//pick another electrode that is not the ground electrode
			int electrode_idx = (ground_electrode_index < electrode_rects.size() - 1 ? electrode_rects.size() - 1 : electrode_rects.size() - 2);

			//not get its center and potential
			DBL3 electrode_center = electrode_rects[electrode_idx].get_c();
			double electrode_potential = electrode_potentials[electrode_idx];

			for (int idx = 0; idx < pTransport.size(); idx++) {

				pTransport[idx]->pMesh->V.set_linear(ground_electrode_center, ground_potential, electrode_center, electrode_potential);
			}
		}
	}
}

void STransport::SetCurrent(double current_)
{
	constant_current_source = true;

	current = current_;
	
	//get_current will now also adjust the potential depending on the sample resistance
	GetCurrent();
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

	UpdateConfiguration();
}

bool STransport::DelElectrode(int index)
{
	if (GoodIdx(electrode_rects.size(), index)) {

		electrode_rects.erase(index);
		electrode_potentials.erase(electrode_potentials.begin() + index);

		//have we deleted the ground electrode?
		if (ground_electrode_index == index) ground_electrode_index = -1;

		UpdateConfiguration();

		return true;
	}
	else return false;
}

void STransport::ClearElectrodes(void)
{
	electrode_rects.clear();
	electrode_potentials.clear();

	UpdateConfiguration();
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

		UpdateConfiguration();

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

pair<Rect, double> STransport::GetElectrodeInfo(int index)
{
	pair<Rect, double> elInfo;

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