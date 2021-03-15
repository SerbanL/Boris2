#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

class ModulesCUDA;

class ManagedModulesCUDA {

public:

	//energy value for this effective field term
	cuBReal* penergy;

	//auxiliary for obtaining average energy in a custom rectangle : count non-zero points for average.
	size_t* ppoints_count;

	//effective field in this module : if sized then module should be updating this if appropriate, else skip saving effective field
	//2nd VEC used for 2-sublattice modules
	cuVEC<cuReal3>* pModule_Heff;
	cuVEC<cuReal3>* pModule_Heff2;

	//energy (density) spatial variation : if sized then module should be updating this is appropriate, else skip saving energy density
	//2nd VEC used for 2-sublattice modules
	cuVEC<cuBReal>* pModule_energy;
	cuVEC<cuBReal>* pModule_energy2;

private:

public:

	void construct_cu_obj(void) {}

	void destruct_cu_obj(void) {}

	BError set_pointers(ModulesCUDA* pModulesCUDA);
};

#endif
