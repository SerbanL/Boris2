#include "stdafx.h"
#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "STransport.h"
#include "SuperMesh.h"

STransportCUDA::STransportCUDA(SuperMesh* pSMesh_, STransport* pSTrans_) :
	ModulesCUDA()
{
	pSMesh = pSMesh_;
	pSTrans = pSTrans_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

STransportCUDA::~STransportCUDA()
{

}

//-------------------Abstract base class method implementations

BError STransportCUDA::Initialize(void)
{
	BError error(CLASS_STR(STransportCUDA));

	//no energy density contribution here
	ZeroEnergy();

	if (!initialized) {

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

		pSTrans->recalculate_transport = true;
		pSTrans->transport_recalculated = true;

		initialized = true;
	}

	return error;
}

BError STransportCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(STransportCUDA));

	Uninitialize();
	
	if (ucfg::check_cfgflags(cfgMessage,
		UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE,
		UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED,
		UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED,
		UPDATECONFIG_TRANSPORT_ELECTRODE,
		UPDATECONFIG_ODE_SOLVER)) {
		
		////////////////////////////////////////////////////////////////////////////
		//check meshes to set transport boundary flags (NF_CMBND flags for V)
		////////////////////////////////////////////////////////////////////////////
		
		//clear everything then rebuild
		pTransport.clear();
		CMBNDcontactsCUDA.clear();
		CMBNDcontacts.clear();
		pV.clear();
		pS.clear();
		
		//now build pTransport (and pV)
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->IsModuleSet(MOD_TRANSPORT) && !(*pSMesh)[idx]->is_atomistic()) {

				pTransport.push_back(dynamic_cast<TransportCUDA*>((*pSMesh)[idx]->GetCUDAModule(MOD_TRANSPORT)));
				pV.push_back(&dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->V);
				pS.push_back(&dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA->S);
			}
		}
		
		//set fixed potential cells and cmbnd flags
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			//it's easier to just copy the flags entirely from the cpu versions.
			//Notes :
			//1. By design the cpu versions are required to keep size and flags up to date (but not mesh values)
			//2. pTransport in STransport has exactly the same size and order
			//3. STransport UpdateConfiguration was called just before, which called this CUDA version at the end.

			if (!(*pV[idx])()->copyflags_from_cpuvec(*pSTrans->pV[idx])) error(BERROR_GPUERROR_CRIT);

			if (pSMesh->SolveSpinCurrent()) {

				if (!(*pS[idx])()->copyflags_from_cpuvec(*pSTrans->pS[idx])) error(BERROR_GPUERROR_CRIT);
			}
		}

		for (int idx = 0; idx < pSTrans->CMBNDcontacts.size(); idx++) {

			vector<cu_obj<CMBNDInfoCUDA>> mesh_contacts;
			vector<CMBNDInfoCUDA> mesh_contacts_cpu;

			for (int idx_contact = 0; idx_contact < pSTrans->CMBNDcontacts[idx].size(); idx_contact++) {

				cu_obj<CMBNDInfoCUDA> contact;

				contact()->copy_from_CMBNDInfo<CMBNDInfo>(pSTrans->CMBNDcontacts[idx][idx_contact]);

				mesh_contacts.push_back(contact);

				mesh_contacts_cpu.push_back(pSTrans->CMBNDcontacts[idx][idx_contact]);
			}

			CMBNDcontactsCUDA.push_back(mesh_contacts);
			CMBNDcontacts.push_back(mesh_contacts_cpu);
		}

		//copy fixed SOR damping from STransport
		SOR_damping_V.from_cpu(pSTrans->SOR_damping.i);
		SOR_damping_S.from_cpu(pSTrans->SOR_damping.j);
	}

	return error;
}

//scale all potential values in all V cuVECs by given scaling value
void STransportCUDA::scale_potential_values(cuBReal scaling)
{
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		(*pV[idx])()->scale_values(scaling);
	}
}

//set potential values using a slope between the potential values of ground and another electrode (if set)
void STransportCUDA::initialize_potential_values(void)
{
	//Note, it's possible V already has values, e.g. we've just loaded a simulation file with V saved.
	//We don't want to re-initialize the V values as this will force the transport solver to iterate many times to get back the correct V values - which we already have!
	//Then, only apply the default V initialization if the voltage values are zero - if the average V is exactly zero (averaged over all meshes) then it's highly probable V is zero everywhere.
	//It could be that V has a perfectly anti-symmetrical set of values, in which case the average will also be zero. But in this case there's also no point to re-initialize the values.
	double V_average = 0;

	for (int idx = 0; idx < pV.size(); idx++) {

		V_average += (*pV[idx])()->average_nonempty((*pV[idx])()->size_cpu().dim());
	}

	if (IsZ(V_average)) {

		if (pSTrans->ground_electrode_index >= 0 && pSTrans->electrode_rects.size() >= 2) {

			DBL3 ground_electrode_center = pSTrans->electrode_rects[pSTrans->ground_electrode_index].get_c();
			double ground_potential = pSTrans->electrode_potentials[pSTrans->ground_electrode_index];

			//pick another electrode that is not the ground electrode
			int electrode_idx = (pSTrans->ground_electrode_index < pSTrans->electrode_rects.size() - 1 ? pSTrans->electrode_rects.size() - 1 : pSTrans->electrode_rects.size() - 2);

			//not get its center and potential
			DBL3 electrode_center = pSTrans->electrode_rects[electrode_idx].get_c();
			double electrode_potential = pSTrans->electrode_potentials[electrode_idx];

			for (int idx = 0; idx < pTransport.size(); idx++) {

				pTransport[idx]->pMeshCUDA->V()->set_linear(ground_electrode_center, ground_potential, electrode_center, electrode_potential);
			}
		}
	}
}

void STransportCUDA::UpdateField(void)
{
	//skip any transport solver computations if static_transport_solver is enabled : transport solver will be interated only at the end of a step or stage
	if (pSMesh->static_transport_solver) return;

	//only need to update this after an entire magnetization equation time step is solved (but always update spin accumulation field if spin current solver enabled)
	if (pSMesh->CurrentTimeStepSolved()) {

		//use V or I equation to set electrode potentials? time dependence only
		if (pSTrans->V_equation.is_set() || pSTrans->I_equation.is_set()) {

			if (pSTrans->V_equation.is_set()) pSTrans->SetPotential(pSTrans->V_equation.evaluate(pSMesh->GetStageTime()), false);
			else pSTrans->SetCurrent(pSTrans->I_equation.evaluate(pSMesh->GetStageTime()), false);

			pSTrans->recalculate_transport = true;
		}

		pSTrans->transport_recalculated = pSTrans->recalculate_transport;

		if (pSTrans->recalculate_transport) {

			pSTrans->recalculate_transport = false;

			//solve only for charge current (V and Jc with continuous boundaries)
			if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
			//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
			else solve_spin_transport_sor();

			//if constant current source is set then need to update potential to keep a constant current
			if (pSTrans->constant_current_source) {

				pSTrans->GetCurrent();

				//the electrode voltage values will have changed so should iterate to convergence threshold again
				//even though we've adjusted potential values these won't be quite correct
				//moreover the charge current density hasn't been recalculated
				//reiterating the transport solver to convergence will fix all this
				//Note : the electrode current will change again slightly so really you should be iterating to some electrode current convergence threshold
				//In normal running mode this won't be an issue as this is done every iteration; in static transport solver mode this could be a problem so must be tested

				double iters_to_conv_previous = pSTrans->iters_to_conv;

				//solve only for charge current (V and Jc with continuous boundaries)
				if (!pSMesh->SolveSpinCurrent()) solve_charge_transport_sor();
				//solve both spin and charge currents (V, Jc, S with appropriate boundaries : continuous, except between N and F layers where interface conductivities are specified)
				else solve_spin_transport_sor();

				//in constant current mode we spend more iterations so the user should be aware of this
				pSTrans->iters_to_conv += iters_to_conv_previous;
			}
		}
		else pSTrans->iters_to_conv = 0;
	}

	if (pSMesh->SolveSpinCurrent()) {

		//Calculate the spin accumulation field so a torque is generated when used in the LLG (or LLB) equation
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			pTransport[idx]->CalculateSAField();
		}

		//Calculate effective field from interface spin accumulation torque (in magnetic meshes for NF interfaces with G interface conductance set)
		CalculateSAInterfaceField();
	}
}

//-------------------Setters

#endif

#endif

