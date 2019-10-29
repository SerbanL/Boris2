#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

int Mesh::meshIdCounter = 0;

Mesh::Mesh(MESH_ meshType, SuperMesh *pSMesh_) :
	MeshParams(params_for_meshtype(meshType)),
	pSMesh(pSMesh_)
{
	this->meshType = meshType;

	OmpThreads = omp_get_num_procs();
}

Mesh::~Mesh() 
{
	//delete all allocated Modules
	clear_vector(pMod);

#if COMPILECUDA == 1
	//free cuda memory by deleting allocated pMeshCUDA
	if (pMeshCUDA) {
		
		//mark implementation of Mesh as destroyed so the CUDA mesh version doesn't attempt to use its data in destructor
		pMeshCUDA->Holder_Mesh_Destroyed();

		delete pMeshCUDA;
		pMeshCUDA = nullptr;
	}
#endif
}

BError Mesh::Error_On_Create(void)
{ 
	for (int idx = 0; idx < pMod.size(); idx++) {

		if(!error_on_create) error_on_create = pMod[idx]->Error_On_Create();
	}

	return error_on_create; 
}

//----------------------------------- VALUE GETTERS

//get energy value for given module or one of its exclusive modules (if none active return 0); call it with MOD_ALL to return total energy density in this mesh.
double Mesh::GetEnergy(MOD_ moduleType)
{
	if (moduleType == MOD_ALL) {

		double energy = 0.0;

		//get total energy density in currently set modules
		for (int idx = 0; idx < pMod.size(); idx++) {

			energy += pMod[idx]->GetEnergy();
		}

		return energy;
	}
	else {

		for (int idx = 0; idx < (int)exclusiveModules[moduleType].size(); idx++) {

			MOD_ module = exclusiveModules[moduleType][idx];
			if (IsModuleSet(module)) { return pMod(module)->GetEnergy(); }
		}

		return 0.0;
	}
}

//get average magnetisation in given rectangle (entire mesh if none specified)
DBL3 Mesh::GetAverageMagnetisation(Rect rectangle) 
{ 
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetAverageMagnetisation(rectangle);
#endif

	if (M.linear_size()) return M.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

//get average magnetisation in given rectangle (entire mesh if none specified); sub-lattice B
DBL3 Mesh::GetAverageMagnetisation2(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetAverageMagnetisation2(rectangle);
#endif

	if (M2.linear_size()) return M2.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

DBL3 Mesh::GetAverageChargeCurrentDensity(Rect rectangle)
{ 
	return CallModuleMethod(&Transport::GetAverageChargeCurrent, rectangle);
}

DBL3 Mesh::GetAverageSpinCurrentX(Rect rectangle)
{
	return CallModuleMethod(&Transport::GetAverageSpinCurrent, 0, rectangle);
}

DBL3 Mesh::GetAverageSpinCurrentY(Rect rectangle)
{
	return CallModuleMethod(&Transport::GetAverageSpinCurrent, 1, rectangle);
}

DBL3 Mesh::GetAverageSpinCurrentZ(Rect rectangle)
{
	return CallModuleMethod(&Transport::GetAverageSpinCurrent, 2, rectangle);
}

double Mesh::GetAverageElectricalPotential(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetAverageElectricalPotential(rectangle);
#endif

	if (V.linear_size()) return V.average_nonempty_omp(rectangle);
	else return 0.0;
}

DBL3 Mesh::GetAverageSpinAccumulation(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetAverageSpinAccumulation(rectangle);
#endif

	if (S.linear_size()) return S.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

double Mesh::GetAverageElectricalConductivity(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetAverageElectricalConductivity(rectangle);
#endif

	if (elC.linear_size()) return elC.average_nonempty_omp(rectangle);
	else return 0.0;
}

double Mesh::GetAverageTemperature(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetAverageTemperature(rectangle);
#endif

	if (Temp.linear_size()) return Temp.average_nonempty_omp(rectangle); 
	else return base_temperature; 
}

//----------------------------------- QUANTITY GETTERS

//returns M on the cpu, thus transfers M from gpu to cpu before returning if cuda enabled
VEC_VC<DBL3>& Mesh::Get_M(void)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
	}
#endif

	return M;
}

//returns charge current on the cpu, assuming transport module is enabled
VEC_VC<DBL3>& Mesh::Get_Jc(void)
{
	return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent();
}

//returns bulk self-consistent spin torque on the cpu, assuming transport module is enabled and spin solver is enabled
VEC<DBL3>& Mesh::Get_SpinTorque(void)
{
	return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque();
}

//returns interfacial self-consistent spin torque on the cpu, assuming transport module is enabled and spin solver is enabled
VEC<DBL3>& Mesh::Get_InterfacialSpinTorque(void)
{
	return reinterpret_cast<STransport*>(pSMesh->GetSuperMeshModule(MODS_STRANSPORT))->GetInterfacialSpinTorque(reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT)));
}

//----------------------------------- MESH INFO AND SIZE GET/SET METHODS

BError Mesh::SetMeshRect(Rect meshRect_)
{
	BError error(__FUNCTION__);

	//cannot set a plane, must be a proper 3D rect
	if (meshRect_.IsPlane() || meshRect_.e <= meshRect_.s) return error(BERROR_INCORRECTVALUE);

	meshRect = meshRect_;

	auto adjust_default_cellsize = [](Rect meshRect, DBL3& cellsize) {

		if (!cellsize.dim()) cellsize = DBL3(DEFAULTCELLSIZE);
		INT3 cells = round(meshRect / cellsize);

		if (cells.x > MAXSTARTINGCELLS_X) cells.x = MAXSTARTINGCELLS_X;
		if (cells.y > MAXSTARTINGCELLS_Y) cells.y = MAXSTARTINGCELLS_Y;
		if (cells.z > MAXSTARTINGCELLS_Z) cells.z = MAXSTARTINGCELLS_Z;

		//adjusted cellsize which will result in an integer number of cells with upper limits set
		cellsize = meshRect / cells;
	};

	//adjust cellsizes from their current values so they result in an integer number of cells with upper limits on number of cells in each dimension. The cellsizes can be manually adjusted later if needed.
	adjust_default_cellsize(meshRect, h);
	adjust_default_cellsize(meshRect, h_e);
	adjust_default_cellsize(meshRect, h_t);

	error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHCHANGE);

	return error;
}

//ferromagnetic properties
BError Mesh::SetMeshCellsize(DBL3 h_)
{ 
	BError error(__FUNCTION__);

	h = h_;
	error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHCHANGE);

	return error;
}

//electrical conduction properties
BError Mesh::SetMeshECellsize(DBL3 h_e_)
{ 
	BError error(__FUNCTION__);

	h_e = h_e_;
	error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHCHANGE);

	return error;
}

//thermal conduction properties
BError Mesh::SetMeshTCellsize(DBL3 h_t_)
{ 
	BError error(__FUNCTION__);

	h_t = h_t_;
	error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHCHANGE);

	return error;
}

//----------------------------------- OTHER MESH SHAPE CONTROL

BError Mesh::copy_mesh_data(Mesh& copy_this)
{
	BError error(__FUNCTION__);

	if (meshType != copy_this.GetMeshType()) return error(BERROR_INCORRECTVALUE);

#if COMPILECUDA == 1
	//if CUDA on then first transfer data to cpu as there may be a mismatch
	if (pMeshCUDA) {

		error = pMeshCUDA->copy_shapes_to_cpu();
	}
#endif

	//1a. shape magnetization
	if (M.linear_size() && copy_this.M.linear_size()) {

		//if Roughness module is enabled then apply shape via the Roughness module instead
		if (IsModuleSet(MOD_ROUGHNESS) && copy_this.IsModuleSet(MOD_ROUGHNESS)) {

			error = reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->copy_roughness(reinterpret_cast<Roughness*>(copy_this.pMod(MOD_ROUGHNESS)));

			if (error) return error;
		}
		else {

			M.copy_values(copy_this.M);

			//1b. shape magnetization in AF meshes
			if (M2.linear_size() && copy_this.M2.linear_size()) {

				M2.copy_values(copy_this.M2);
			}
		}
	}

	//2. shape electrical conductivity
	if (elC.linear_size() && copy_this.elC.linear_size()) elC.copy_values(copy_this.elC);

	//3. shape temperature
	if (Temp.linear_size() && copy_this.Temp.linear_size()) Temp.copy_values(copy_this.Temp);

#if COMPILECUDA == 1
	//if CUDA on then load back to gpu
	if (pMeshCUDA) {

		error = pMeshCUDA->copy_shapes_from_cpu();
	}
#endif

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}

	return error;
}

//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
BError Mesh::applymask(double zDepth_m, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, double zDepth_m, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader) -> BError {

		BError error;

		INT3 cells = VEC_quantity.n;

		if(!VEC_quantity.apply_bitmap_mask(bitmap_loader(fileName, INT2(cells.x, cells.y)), zDepth_m))
			return error(BERROR_COULDNOTLOADFILE);
		
		return error;
	};

	error = change_mesh_shape(run_this, zDepth_m, fileName, bitmap_loader);

	return error;
}

//set cells to empty in given box (delete by setting entries to zero)
BError Mesh::delrect(Rect rectangle)
{ 
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, Rect& rectangle) -> BError {

		VEC_quantity.delrect(rectangle);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	return error;
}

//set cells to non-empty in given box
BError Mesh::setrect(Rect rectangle)
{ 
	BError error(__FUNCTION__);

	auto run_this = [](auto& VEC_quantity, auto default_value, Rect& rectangle) -> BError {

		VEC_quantity.setrect(rectangle, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	return error;
}

//roughen mesh sides perpendicular to a named axis (axis = "x", "y", "z") to given depth (same units as h) with prng instantiated with given seed.
BError Mesh::RoughenMeshSides(string axis, double depth, unsigned seed)
{
	BError error(__FUNCTION__);

	if ((axis != "x" && axis != "y" && axis != "z") || depth <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, string axis, double depth, unsigned seed) -> BError {

		bool success = true;

		string side_m = "-z", side_p = "z";

		if (axis == "x") {

			side_m = "-x";
			side_p = "x";
		}

		if (axis == "y") {

			side_m = "-y";
			side_p = "y";
		}	

		success &= VEC_VC_quantity.generate_roughside(side_m, depth, seed);
		success &= VEC_VC_quantity.generate_roughside(side_p, depth, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, axis, depth, seed);

	return error;
}

//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
//Rough both top and bottom if sides is empty, else it should be either -z or z.
BError Mesh::RoughenMeshSurfaces_Jagged(double depth, double spacing, unsigned seed, string sides)
{
	BError error(__FUNCTION__);

	if (depth <= 0 || spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double depth, double spacing, unsigned seed, string sides) -> BError {

		bool success = true;

		success &= VEC_VC_quantity.generate_jagged_surfaces(depth, spacing, seed, sides);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, depth, spacing, seed, sides);

	return error;
}

//Generate Voronoi 2D grains in xy plane (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
BError Mesh::GenerateGrains2D(double spacing, unsigned seed)
{
	BError error(__FUNCTION__);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double spacing, unsigned seed) -> BError {

		bool success = true;

		success = VEC_VC_quantity.generate_Voronoi2D(spacing, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, spacing, seed);

	return error;
}

//Generate Voronoi 3D grains (boundaries between Voronoi cells set to empty) at given average spacing with prng instantiated with given seed.
BError Mesh::GenerateGrains3D(double spacing, unsigned seed)
{
	BError error(__FUNCTION__);

	if (spacing <= 0 || seed < 1) return error(BERROR_INCORRECTVALUE);

	auto run_this = [](auto& VEC_VC_quantity, auto default_value, double spacing, unsigned seed) -> BError {

		bool success = true;

		success = VEC_VC_quantity.generate_Voronoi3D(spacing, seed);

		if (success) return BError();
		else {

			BError error;
			return error(BERROR_INCORRECTVALUE);
		}
	};

	error = change_mesh_shape(run_this, spacing, seed);

	return error;
}