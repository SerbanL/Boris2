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