#include "stdafx.h"
#include "Atom_Mesh.h"

//----------------------------------- MESH QUANTITIES CONTROL

void Atom_Mesh::MoveMesh(double x_shift)
{
	double mesh_end_size = meshRect.size().x * MOVEMESH_ENDRATIO;

	//the rectangle in which to shift mesh values
	Rect shift_rect = Rect(meshRect.s + DBL3(mesh_end_size, 0, 0), meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (paMeshCUDA) {

		//moving mesh for gpu memory quantities
		
		//1. shift M1
		if (M1.linear_size()) paMeshCUDA->M1()->shift_x(M1.linear_size(), x_shift, shift_rect);

		//TO DO
		//2. shift elC
		//if (elC.linear_size()) CallModuleMethod(&Transport::MoveMesh_Transport, x_shift);

		//3. shift Temp
		//if (Temp.linear_size()) CallModuleMethod(&Heat::MoveMesh_Heat, x_shift);

		return;
	}
#endif

	//moving mesh for cpu memory quantities

	//1. shift M1
	if (M1.linear_size()) M1.shift_x(x_shift, shift_rect);

	//TO DO
	//2. shift elC
	//if (elC.linear_size()) CallModuleMethod(&Transport::MoveMesh_Transport, x_shift);

	//3. shift Temp
	//if (Temp.linear_size()) CallModuleMethod(&Heat::MoveMesh_Heat, x_shift);
}

//set PBC for required VECs : should only be called from a demag module
BError Atom_Mesh::Set_Magnetic_PBC(INT3 pbc_images)
{
	BError error(__FUNCTION__);

	if (M1.linear_size()) M1.set_pbc(pbc_images.x, pbc_images.y, pbc_images.z);

#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (M1.linear_size() && !paMeshCUDA->M1()->copyflags_from_cpuvec(M1)) return error(BERROR_GPUERROR_CRIT);
	}
#endif

	return error;
}