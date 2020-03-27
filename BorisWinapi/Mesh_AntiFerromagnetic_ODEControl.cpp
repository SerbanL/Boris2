#include "stdafx.h"
#include "Mesh_AntiFerromagnetic.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

//----------------------------------- ODE CONTROL METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_..._ODEControl.cpp

//get rate of change of magnetization (overloaded by Ferromagnetic meshes)
DBL3 AFMesh::dMdt(int idx)
{
	return meshODE.dMdt(idx);
}

#endif