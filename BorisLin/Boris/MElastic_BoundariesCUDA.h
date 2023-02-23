#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.h"

class MElastic_BoundaryCUDA {

public:

	//surf_rect not required here, only the derived cells_box and surf_orientation

	//cellsize of mesh for which this boundary applies
	cuReal3 h_m;

	//box of mesh cells corresponding to the surface rect
	//start and end coordinates are cell index coordinates in the mesh to which this boundary applies
	cuBox cells_box;

	//what is the surface rectangle orientation (surface normal)?
	//1 : +x, -1 : -x
	//2 : +y, -2 : -y
	//3 : +z, -3 : -z
	int surf_orientation;

	//equation for external surface force (pressure - N/m^2) applicable on surface_rect
	//depends on x, y, t, where the x, y coordinates are relative to surf_rect. t is the simulation time.
	//when evaluating the equation with x y coordinates, the origin of surf_rect is the bottom-left corner when looking at surf_rect from outside the mesh (i.e. against surface normal direction)
	//if nullptr then equation not assigned
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal> *pFext_equation_x;
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal> *pFext_equation_y;
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal> *pFext_equation_z;

	//external surface force (pressure - N/m^2), used if Fext_equation not assigned
	cuReal3 Fext_constant;

	//use Fext_constant instead of Fext_equation
	bool use_constant_force;

public:

	//---------------------------------

	__host__ void construct_cu_obj(void) 
	{
		nullgpuptr(pFext_equation_x);
		nullgpuptr(pFext_equation_y);
		nullgpuptr(pFext_equation_z);

		set_gpu_value(Fext_constant, cuReal3());
		set_gpu_value(use_constant_force, true);
	}

	__host__ void destruct_cu_obj(void) {}

	//copy constructor needed so we can store cu_obj<MElastic_BoundaryCUDA> in a std::vector
	//it doesn't need to copy anything though, as empty cu_obj<MElastic_BoundaryCUDA> will be created when vector is resized
	__host__ void construct_cu_obj(const MElastic_BoundaryCUDA& copy_this)
	{
	}

	////////////////////////////////////// SETUP

	__host__ void setup_surface(cuBox cells_box_, cuReal3 h_m_, int surf_orientation_)
	{
		set_gpu_value(cells_box, cells_box_);
		set_gpu_value(h_m, h_m_);
		set_gpu_value(surf_orientation, surf_orientation_);
	}

	////////////////////////////////////// EQUATION STIMULUS

	//from Fext_equation make the CUDA version
	__host__ bool make_cuda_equation(TEquationCUDA<cuBReal, cuBReal, cuBReal>& Fext_equationCUDA, std::vector<std::vector< std::vector<EqComp::FSPEC> >> fspec)
	{
		bool success = Fext_equationCUDA.make_vector(fspec);

		if (success) {

			set_gpu_value(pFext_equation_x, Fext_equationCUDA.get_pcu_obj_x()->get_managed_object());
			set_gpu_value(pFext_equation_y, Fext_equationCUDA.get_pcu_obj_y()->get_managed_object());
			set_gpu_value(pFext_equation_z, Fext_equationCUDA.get_pcu_obj_z()->get_managed_object());

			set_gpu_value(use_constant_force, false);
		}
		else {

			nullgpuptr(pFext_equation_x);
			nullgpuptr(pFext_equation_y);
			nullgpuptr(pFext_equation_z);
		}

		return success;
	}

	////////////////////////////////////// FIXED STIMULUS

	__host__ void setup_fixed_stimulus(cuReal3 Fext_constant_)
	{
		set_gpu_value(Fext_constant, Fext_constant_);
		set_gpu_value(use_constant_force, true);

		nullgpuptr(pFext_equation_x);
		nullgpuptr(pFext_equation_y);
		nullgpuptr(pFext_equation_z);
	}

	////////////////////////////////////// SEARCH (runtime)

	//does cells_box contain cell with given i, j, k indexes? (indexed on the same mesh resolution, i.e. at the h_m, n_m mesh resolution)
	//return 0 if it doesn't, else return the corresponding surface orientiation
	__device__ int contains(const cuINT3& ijk)
	{
		if (cells_box.Contains(ijk)) return surf_orientation;
		else return 0;
	}

	////////////////////////////////////// GET (runtime)

	//get external force (Fx, Fy, Fz) at vertex with given index
	//this is used for diagonal components only when solving for stress
	//(sxx = +/-Fx at x normal face, syy = +/-Fy at y normal face, szz = +/-Fz at z normal face)
	//external force along surface normal returned (specified by user in xyz system)
	__device__ cuBReal get_ext_force_vertices(const cuINT3& ijk_vertex, const cuBReal& time)
	{
		if (use_constant_force) {

			//Fixed external surface force
			switch (surf_orientation) {

			case -1:
				return -Fext_constant.x;
				break;

			case 1:
				return Fext_constant.x;
				break;

			case -2:
				return -Fext_constant.y;
				break;

			case 2:
				return Fext_constant.y;
				break;

			case -3:
				return -Fext_constant.z;
				break;

			case 3:
				return Fext_constant.z;
				break;

			default:
				return 0.0;
				break;
			}
		}
		else {

			//External surface force from text equation

			//x, y vertex position relative to surf_rect
			cuBReal x, y;

			switch (surf_orientation) {

				//-x (so surface is in yz plane). origin at +y, -z
			case -1:
				x = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return -pFext_equation_x->evaluate(x, y, time);
				break;

				//+x (so surface is in yz plane). origin at -y, -z
			case +1:
				x = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return pFext_equation_x->evaluate(x, y, time);
				break;

				//-y (so surface is in xz plane). origin at -x, -z
			case -2:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return -pFext_equation_y->evaluate(x, y, time);
				break;

				//+y (so surface is in xz plane). origin at +x, -z
			case +2:
				x = (cells_box.e.i - ijk_vertex.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return pFext_equation_y->evaluate(x, y, time);
				break;

				//-z (so surface is in xy plane). origin at -x, +y
			case -3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				return -pFext_equation_z->evaluate(x, y, time);
				break;

				//+z (so surface is in xy plane). origin at -x, -y
			case +3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				return pFext_equation_z->evaluate(x, y, time);
				break;

			default:
				return 0.0;
				break;
			}
		}
	}

	//get external force (Fx, Fy, Fz) at 2 edge centers, where edges meet at given vertex index
	//this is used for off-diagonal components of stress only when solving for velocity
	//x normal face:
	//sxy = +/-Fy, sxz = +/-Fz, so return (0.0, Fy, Fz)
	//y normal face
	//sxy = +/-Fx, syz = +/-Fz, so return (Fx, 0.0, Fz)
	//z normal face
	//sxz = +/-Fx, syz = +/-Fy, so return (Fx, Fy, 0.0)
	//In all cases the forces Fx, Fy, Fz are evaluated at centre of corresponding +x, +y, +z axis edges originating from given ijk_vertex
	__device__ cuReal3 get_ext_force_edges(const cuINT3& ijk_vertex, const cuBReal& time)
	{
		if (use_constant_force) {

			//Fixed external surface force
			switch (surf_orientation) {

			case -1:
				return cuReal3(0.0, -Fext_constant.y, -Fext_constant.z);
				break;

			case 1:
				return cuReal3(0.0, Fext_constant.y, Fext_constant.z);
				break;

			case -2:
				return cuReal3(-Fext_constant.x, 0.0, -Fext_constant.z);
				break;

			case 2:
				return cuReal3(Fext_constant.x, 0.0, Fext_constant.z);
				break;

			case -3:
				return cuReal3(-Fext_constant.x, -Fext_constant.y, 0.0);
				break;

			case 3:
				return cuReal3(Fext_constant.x, Fext_constant.y, 0.0);
				break;

			default:
				return cuReal3();
				break;
			}
		}
		else {

			//External surface force from text equation

			//x, y vertex position relative to surf_rect
			cuBReal x, y;

			switch (surf_orientation) {

				//-x (so surface is in yz plane). origin at +y, -z
			case -1:
				x = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return cuReal3(0.0, -pFext_equation_y->evaluate(x + h_m.y / 2, y, time), -pFext_equation_z->evaluate(x, y + h_m.z / 2, time));
				break;

				//+x (so surface is in yz plane). origin at -y, -z
			case +1:
				x = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return cuReal3(0.0, pFext_equation_y->evaluate(x + h_m.y / 2, y, time), pFext_equation_z->evaluate(x, y + h_m.z / 2, time));
				break;

				//-y (so surface is in xz plane). origin at -x, -z
			case -2:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return cuReal3(-pFext_equation_x->evaluate(x + h_m.x / 2, y, time), 0.0, -pFext_equation_z->evaluate(x, y + h_m.z / 2, time));
				break;

				//+y (so surface is in xz plane). origin at +x, -z
			case +2:
				x = (cells_box.e.i - ijk_vertex.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return cuReal3(pFext_equation_x->evaluate(x + h_m.x / 2, y, time), 0.0, pFext_equation_z->evaluate(x, y + h_m.z / 2, time));
				break;

				//-z (so surface is in xy plane). origin at -x, +y
			case -3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				return cuReal3(-pFext_equation_x->evaluate(x + h_m.x / 2, y, time), -pFext_equation_y->evaluate(x, y + h_m.y / 2, time), 0.0);
				break;

				//+z (so surface is in xy plane). origin at -x, -y
			case +3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				return cuReal3(pFext_equation_x->evaluate(x + h_m.x / 2, y, time), pFext_equation_y->evaluate(x, y + h_m.y / 2, time), 0.0);
				break;

			default:
				return cuReal3();
				break;
			}
		}
	}
};

#endif

#endif
