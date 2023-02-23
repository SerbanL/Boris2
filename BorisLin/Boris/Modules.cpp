#include "stdafx.h"

#include "Modules.h"

//-------------------------- CUDA Switch

//switch CUDA state on/off
BError Modules::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(Modules));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		if (!pModuleCUDA) {

			error = MakeCUDAModule();
		}
	}
	else {

		//cuda switched off so delete cuda module object
		if (pModuleCUDA) delete pModuleCUDA;
		pModuleCUDA = nullptr;
	}

#endif

	return error;
}

//-------------------------- Effective field and energy VECs

//Make sure memory is allocated correctly for display data if used, else free memory
BError Modules::Update_Module_Display_VECs(DBL3 h, Rect meshRect, bool Module_Heff_used, bool Module_Energy_used, bool twosublattice)
{
	BError error(CLASS_STR(Modules));
	
	//1. Heff - sub-lattice A

	if (Module_Heff.linear_size()) {

		if (Module_Heff_used && !Module_Heff.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		else if(!Module_Heff_used) Module_Heff.clear();
	}
	else if (Module_Heff_used) {

		if (!Module_Heff.assign(h, meshRect, DBL3())) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//2. Heff - sub-lattice B

	if (Module_Heff2.linear_size()) {

		if (twosublattice && Module_Heff_used && !Module_Heff2.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		else Module_Heff2.clear();
	}
	else if (twosublattice && Module_Heff_used) {

		if (!Module_Heff2.assign(h, meshRect, DBL3())) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//3. Energy Density - sub-lattice A

	if (Module_energy.linear_size()) {

		if (Module_Energy_used && !Module_energy.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		else if (!Module_Energy_used) Module_energy.clear();
	}
	else if (Module_Energy_used) {

		if (!Module_energy.assign(h, meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//4. Energy Density - sub-lattice B

	if (Module_energy2.linear_size()) {

		if (twosublattice && Module_Energy_used && !Module_energy2.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		else Module_energy2.clear();
	}
	else if (twosublattice && Module_Energy_used) {

		if (!Module_energy2.assign(h, meshRect, 0.0)) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	return error;
}

//zero the VECs if not empty
void Modules::ZeroModuleVECs(void)
{
	if (Module_Heff.linear_size()) Module_Heff.set();
	if (Module_Heff2.linear_size()) Module_Heff2.set();

	if (Module_energy.linear_size()) Module_energy.set();
	if (Module_energy2.linear_size()) Module_energy2.set();
}

//return cross product of M with Module_Heff, averaged in given rect (relative)
DBL3 Modules::CalculateTorque(VEC_VC<DBL3>& M, Rect& avRect)
{
	if (Module_Heff.size() != M.size()) return DBL3();

	Box box = M.box_from_rect_max(avRect + M.rect.s);

	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * M.n.x + k * M.n.x*M.n.y;
				if (M.is_not_empty(idx)) reduction.reduce_average(M[idx] ^ Module_Heff[idx]);
			}
		}
	}

	return reduction.average();
}

//-------------------------- Energy density calculation

//Get energy density averaged over the entire mesh during the UpdateField call
double Modules::GetEnergyDensity(void)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetEnergyDensity();
#endif

	return energy;
}

//Calculate the energy density in the given rect only
double Modules::GetEnergyDensity(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetEnergyDensity(avRect);
#endif

	double energy1 = 0.0, energy2 = 0.0;

	if (Module_energy.linear_size()) {

		energy1 = Module_energy.average_nonempty_omp(avRect);
	}

	if (Module_energy2.linear_size()) {

		energy2 = Module_energy2.average_nonempty_omp(avRect);
		return (energy1 + energy2) / 2;
	}
	else return energy1;
}

