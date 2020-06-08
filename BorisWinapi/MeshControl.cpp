#include "stdafx.h"
#include "Mesh.h"

#include "Transport.h"
#include "Heat.h"

//----------------------------------- MESH QUANTITIES CONTROL

void Mesh::MoveMesh(double x_shift)
{
	double mesh_end_size = meshRect.size().x * MOVEMESH_ENDRATIO;

	//the rectangle in which to shift mesh values
	Rect shift_rect = Rect(meshRect.s + DBL3(mesh_end_size, 0, 0), meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		//moving mesh for gpu memory quantities
		
		//1. shift M
		if (M.linear_size()) pMeshCUDA->M()->shift_x(M.linear_size(), x_shift, shift_rect);
		if (M2.linear_size()) pMeshCUDA->M2()->shift_x(M2.linear_size(), x_shift, shift_rect);

		//2. shift elC
		if (elC.linear_size()) CallModuleMethod(&Transport::MoveMesh_Transport, x_shift);

		//3. shift Temp
		if (Temp.linear_size()) CallModuleMethod(&Heat::MoveMesh_Heat, x_shift);

		return;
	}
#endif

	//moving mesh for cpu memory quantities

	//1. shift M
	if (M.linear_size()) M.shift_x(x_shift, shift_rect);
	if (M2.linear_size()) M2.shift_x(x_shift, shift_rect);

	//2. shift elC
	if (elC.linear_size()) CallModuleMethod(&Transport::MoveMesh_Transport, x_shift);

	//3. shift Temp
	if (Temp.linear_size()) CallModuleMethod(&Heat::MoveMesh_Heat, x_shift);
}

//set PBC for required VECs : should only be called from a demag module
BError Mesh::Set_Magnetic_PBC(INT3 pbc_images)
{
	BError error(__FUNCTION__);

	if (M.linear_size()) M.set_pbc(pbc_images.x, pbc_images.y, pbc_images.z);
	if (M2.linear_size()) M2.set_pbc(pbc_images.x, pbc_images.y, pbc_images.z);

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M.linear_size() && !pMeshCUDA->M()->copyflags_from_cpuvec(M)) return error(BERROR_GPUERROR_CRIT);
		if (M2.linear_size() && !pMeshCUDA->M2()->copyflags_from_cpuvec(M2)) return error(BERROR_GPUERROR_CRIT);
	}
#endif

	return error;
}