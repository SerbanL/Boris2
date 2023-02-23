#include "stdafx.h"

#include "ModulesCUDA.h"
#if COMPILECUDA == 1

//Calculate the energy density in the given rect only
cuBReal ModulesCUDA::GetEnergyDensity(cuRect avRect)
{
	cuBReal energy1 = 0.0, energy2 = 0.0;

	if (Module_energy_size) {

		energy1 = Module_energy()->average_nonempty(Module_energy_size, avRect);
	}

	if (Module_energy2_size) {

		energy2 = Module_energy2()->average_nonempty(Module_energy2_size, avRect);
		return (energy1 + energy2) / 2;
	}
	else return energy1;
}

#endif