#include "ModulesCUDA.h"

#if COMPILECUDA == 1

__global__ void ZeroEnergy_kernel(cuBReal& energy, size_t& points_count)
{
	if (threadIdx.x == 0) energy = 0.0;
	if (threadIdx.x == 1) points_count = 0;
}

void ModulesCUDA::ZeroEnergy(void)
{
	ZeroEnergy_kernel <<< 1, CUDATHREADS >>> (energy, points_count);
}

__global__ void ZeroModuleVECs_kernel(size_t size, cuVEC<cuReal3>& Module_Heff, cuVEC<cuReal3>& Module_Heff2, cuVEC<cuBReal>& Module_energy, cuVEC<cuBReal>& Module_energy2)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		if (Module_Heff.linear_size()) Module_Heff[idx] = cuReal3();
		if (Module_Heff2.linear_size()) Module_Heff2[idx] = cuReal3();
		if (Module_energy.linear_size()) Module_energy[idx] = 0.0;
		if (Module_energy2.linear_size()) Module_energy2[idx] = 0.0;
	}
}

void ModulesCUDA::ZeroModuleVECs(void)
{
	size_t size = (Module_Heff_size > Module_Heff2_size ? Module_Heff_size : Module_Heff2_size);
	size = (size > Module_energy_size ? size : Module_energy_size);
	size = (size > Module_energy2_size ? size : Module_energy2_size);

	ZeroModuleVECs_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, Module_Heff, Module_Heff2, Module_energy, Module_energy2);
}

//-------------------------- Effective field and energy VECs

//Make sure memory is allocated correctly for display data if used, else free memory
BError ModulesCUDA::Update_Module_Display_VECs(cuReal3 h, cuRect meshRect, bool Module_Heff_used, bool Module_Energy_used, bool twosublattice)
{
	BError error(CLASS_STR(ModulesCUDA));

	//1. Heff - sub-lattice A

	if (Module_Heff()->size_cpu().dim()) {

		if (Module_Heff_used && !Module_Heff()->resize(h, meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (!Module_Heff_used) Module_Heff()->clear();	
	}
	else if (Module_Heff_used) {

		if (!Module_Heff()->assign(h, meshRect, cuReal3())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	//2. Heff - sub-lattice B

	if (Module_Heff2()->size_cpu().dim()) {

		if (twosublattice && Module_Heff_used && !Module_Heff2()->resize(h, meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else Module_Heff2()->clear();
	}
	else if (twosublattice && Module_Heff_used) {

		if (!Module_Heff2()->assign(h, meshRect, cuReal3())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	//3. Energy Density - sub-lattice A

	if (Module_energy()->size_cpu().dim()) {

		if (Module_Energy_used && !Module_energy()->resize(h, meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (!Module_Energy_used) Module_energy()->clear();
	}
	else if (Module_Energy_used) {

		if (!Module_energy()->assign(h, meshRect, 0.0)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	//4. Energy Density - sub-lattice B

	if (Module_energy2()->size_cpu().dim()) {

		if (twosublattice && Module_Energy_used && !Module_energy2()->resize(h, meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else Module_energy2()->clear();
	}
	else if (twosublattice && Module_Energy_used) {

		if (!Module_energy2()->assign(h, meshRect, 0.0)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	}

	Module_Heff_size = Module_Heff()->linear_size_cpu();
	Module_Heff2_size = Module_Heff2()->linear_size_cpu();
	Module_energy_size = Module_energy()->linear_size_cpu();
	Module_energy2_size = Module_energy2()->linear_size_cpu();

	return error;
}

#endif