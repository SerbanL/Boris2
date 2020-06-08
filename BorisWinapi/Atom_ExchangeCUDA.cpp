#include "stdafx.h"
#include "Atom_ExchangeCUDA.h"

#ifdef MODULE_ATOM_EXCHANGE

#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisLib.h"

Atom_ExchangeCUDA::Atom_ExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_ExchangeCUDA::~Atom_ExchangeCUDA()
{}

BError Atom_ExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_ExchangeCUDA));

	initialized = true;

	return error;
}

BError Atom_ExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_ExchangeCUDA));

	Uninitialize();

	return error;
}

void Atom_ExchangeCUDA::Compute_Exchange(VEC<double>& displayVEC)
{
	Compute_ExchangeCUDA();

	exchange_displayVEC()->copy_to_cpuvec(displayVEC);
}

#endif

#endif

