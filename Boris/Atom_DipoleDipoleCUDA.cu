#include "Atom_DipoleDipoleCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

__global__ void Energy_to_EnergyDensity_Kernel(cuBReal& energy, cuVEC<cuReal3>& V)
{
	if (threadIdx.x == 0) energy *= (cuBReal)MUB / V.h.dim();
}

//convert value in energy to energy density by dividing by cellsize volume of V
void Atom_DipoleDipoleCUDA::Energy_to_EnergyDensity(cu_obj<cuVEC<cuReal3>>& V)
{
	Energy_to_EnergyDensity_Kernel <<< 1, CUDATHREADS >>> (energy, V);
}


#endif

#endif