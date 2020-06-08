#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

//----------------------------------- PARAMETERS TEMPERATURE

//set mesh base temperature
void Mesh::SetBaseTemperature(double Temperature, bool clear_equation)
{
	if (clear_equation) T_equation.clear();

	//new base temperature and adjust parameter output values for new base temperature
	base_temperature = Temperature;

#if COMPILECUDA == 1
	if (pMeshCUDA) pMeshCUDA->base_temperature.from_cpu(base_temperature);
#endif

	//update any text equations used in mesh parameters (base temperature used as a constant - Tb)
	//need to update equations before parameters : if a parameter has a temperature dependence set using an equation by Temp.size() is zero, then updating parameters sets values based on the mesh base temperature
	update_all_meshparam_equations();

	//update parameter current values : if they have a temperature dependence set the base temperature will change their values
	update_parameters();

	//1a. reset Temp VEC to base temperature
	CallModuleMethod(&Heat::SetBaseTemperature, Temperature);

	//1b. Zeeman module might set field using a custom user equation where the base temperature is a parameter
	CallModuleMethod(&Zeeman::SetBaseTemperature, Temperature);

	//NOTE : do not call UpdateConfiguration here - this is to allow Temperature sequences in simulation stages
	//Instead deal with any adjustments required on a module by module basis (e.g. StrayField)

	//2. If this Mesh is the base for a DipoleMesh then need to force recalculateStrayField to be set to true, in case StrayField super-mesh module is set : dipole magnetisation value could have changed now.
	if (dynamic_cast<DipoleMesh*>(this)) dynamic_cast<DipoleMesh*>(this)->Reset_Mdipole();

	//3. electrical conductivity might also need updating so force it here - if Transport module not set then nothing happens (note elC will have zero size in this case)
	CallModuleMethod(&Transport::CalculateElectricalConductivity, true);
}

//----------------------------------- OTHERS

//copy all parameters from another Mesh
BError Mesh::copy_mesh_parameters(MeshBase& copy_this)
{
	BError error(__FUNCTION__);

	if (meshType != copy_this.GetMeshType()) return error(BERROR_INCORRECTVALUE);
	else copy_parameters(*dynamic_cast<Mesh*>(&copy_this));

	//make sure to update the spatial variation as any copied spatial variation may not have the correct cellsize now
	update_meshparam_var();

	return error;
}