#include "stdafx.h"
#include "MeshBase.h"

//Set Temp from given data VEC -> stretch data to mesh dimensions if needed.
void MeshBase::SetTempFromData(VEC<double>& data, const Rect& dstRect)
{
#if COMPILECUDA == 1
	//refresh Temp from gpu memory
	if (pMeshBaseCUDA) pMeshBaseCUDA->Temp()->copy_to_cpuvec(Temp);
#endif

	if (!Temp.linear_size()) return;

	//make sure dimensions match
	if (!data.resize(Temp.n)) return;

	//copy values in data, as well as shape
	Temp.copy_values(data, dstRect);

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshBaseCUDA) pMeshBaseCUDA->Temp()->copy_from_cpuvec(Temp);
#endif
}

//Set E from current density data, by dividing by electrical conductivity (elC). Set V and S to zero. This is meant for computations with fixed Jc and transport solver iteration disabled.
void MeshBase::SetEFromJcData(VEC<DBL3>& data)
{
#if COMPILECUDA == 1
	//refresh Temp from gpu memory
	if (pMeshBaseCUDA) {

		pMeshBaseCUDA->E()->copy_to_cpuvec(E);
		pMeshBaseCUDA->elC()->copy_to_cpuvec(elC);
	}
#endif

	if (!E.linear_size()) return;

	//make sure dimensions match
	if (!data.resize(E.n)) return;

	//copy values in data, as well as shape
	E.copy_values(data);

	//now divide by elC so the correct electric field results for the required Jc (Jc = sigma * E)
#pragma omp parallel for
	for (int idx = 0; idx < E.linear_size(); idx++) {

		if (elC.is_not_empty(idx)) {

			E[idx] /= elC[idx];
		}
		else {

			E[idx] = 0.0;
		}

		V[idx] = 0.0;
		if (S.linear_size()) S[idx] = 0.0;
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshBaseCUDA) {

		pMeshBaseCUDA->E()->copy_from_cpuvec(E);
		pMeshBaseCUDA->V()->copy_from_cpuvec(V);
		if (S.linear_size()) pMeshBaseCUDA->S()->copy_from_cpuvec(S);
	}
#endif
}

//set electric field VEC from a constant Jc value
void MeshBase::SetEFromJcValue(DBL3 Jcvalue)
{
#if COMPILECUDA == 1
	//refresh Temp from gpu memory
	if (pMeshBaseCUDA) {

		pMeshBaseCUDA->E()->copy_to_cpuvec(E);
		pMeshBaseCUDA->elC()->copy_to_cpuvec(elC);
	}
#endif

	if (!E.linear_size()) return;

	//now divide Jc by elC so the correct electric field results
#pragma omp parallel for
	for (int idx = 0; idx < E.linear_size(); idx++) {

		if (elC.is_not_empty(idx)) {

			E[idx] = Jcvalue / elC[idx];
		}
		else {

			E[idx] = 0.0;
		}

		V[idx] = 0.0;
		if (S.linear_size()) S[idx] = 0.0;
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshBaseCUDA) {

		pMeshBaseCUDA->E()->copy_from_cpuvec(E);
		pMeshBaseCUDA->V()->copy_from_cpuvec(V);
		if (S.linear_size()) pMeshBaseCUDA->S()->copy_from_cpuvec(S);
	}
#endif
}