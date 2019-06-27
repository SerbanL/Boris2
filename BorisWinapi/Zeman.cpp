#include "stdafx.h"
#include "Zeeman.h"

Zeeman::Zeeman(Mesh *pMesh) {
	
	Ha = DBL3(0,0,0);

	this->pMesh = pMesh;

	MeshSizeChange();
}

Zeeman::~Zeeman() {

}

bool Zeeman::Initialize(void) {

	initialized = true;

	return initialized;
}

void Zeeman::MeshSizeChange(void) {

	n = pMesh->n;

	Initialize();
}

void Zeeman::UpdateField(void) {

	double energy = 0;

	#pragma omp parallel for reduction(+:energy)
	for(int idx = 0; idx < n.dim(); idx++) {

		pMesh->Heff[idx] += Ha;
		
		energy += pMesh->M[idx]*Ha;
	}
	
	energy *= -MU0 / pMesh->magpoints;
	this->energy = energy;
}