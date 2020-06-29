#include "stdafx.h"
#include "Mesh_Ferromagnetic.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

//----------------------------------- ODE CONTROL METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_..._ODEControl.cpp

//get rate of change of magnetization (overloaded by Ferromagnetic meshes)
DBL3 FMesh::dMdt(int idx)
{
	return meshODE.dMdt(idx);
}

#endif