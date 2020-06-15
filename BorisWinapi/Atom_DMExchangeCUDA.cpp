#include "stdafx.h"
#include "Atom_DMExchangeCUDA.h"

#if defined(MODULE_DMEXCHANGE) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisLib.h"

Atom_DMExchangeCUDA::Atom_DMExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_DMExchangeCUDA::~Atom_DMExchangeCUDA()
{}

BError Atom_DMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_DMExchangeCUDA));

	initialized = true;

	return error;
}

BError Atom_DMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DMExchangeCUDA));

	Uninitialize();

	return error;
}

void Atom_DMExchangeCUDA::Compute_Exchange(VEC<double>& displayVEC)
{
	Compute_ExchangeCUDA();

	exchange_displayVEC()->copy_to_cpuvec(displayVEC);
}

#endif

#endif

