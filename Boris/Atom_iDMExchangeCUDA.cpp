#include "stdafx.h"
#include "Atom_iDMExchangeCUDA.h"

#if defined(MODULE_COMPILATION_IDMEXCHANGE) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisLib.h"

Atom_iDMExchangeCUDA::Atom_iDMExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_iDMExchangeCUDA::~Atom_iDMExchangeCUDA()
{}

BError Atom_iDMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_iDMExchangeCUDA));

	initialized = true;

	return error;
}

BError Atom_iDMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DMExchangeCUDA));

	Uninitialize();

	return error;
}

void Atom_iDMExchangeCUDA::Compute_Exchange(VEC<double>& displayVEC)
{
	Compute_ExchangeCUDA();

	exchange_displayVEC()->copy_to_cpuvec(displayVEC);
}

#endif

#endif

