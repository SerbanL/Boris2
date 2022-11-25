#pragma once

#include "Boris_Enums_Defs.h"
#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisLib.h"

class MElastic_Boundary
{

private:

	//surface rectangle for which this boundary applies (absolute coordinates)
	Rect surf_rect;

	//cellsize of mesh for which this boundary applies
	DBL3 h_m;

	//box of mesh cells corresponding to the surface rect
	//start and end coordinates are cell index coordinates in the mesh to which this boundary applies
	Box cells_box;

	//what is the surface rectangle orientation (surface normal)?
	//1 : +x, -1 : -x
	//2 : +y, -2 : -y
	//3 : +z, -3 : -z
	int surf_orientation = 0;

	//equation for external surface force (pressure - N/m^2) applicable on surface_rect
	//depends on x, y, t, where the x, y coordinates are relative to surf_rect. t is the simulation time.
	//when evaluating the equation with x y coordinates, the origin of surf_rect is the bottom-left corner when looking at surf_rect from outside the mesh (i.e. against surface normal direction)
	//if nullptr then equation not assigned
	TEquation<double, double, double> Fext_equation;

	//external surface force (pressure - N/m^2), used if Fext_equation not assigned
	DBL3 Fext_constant;

	//use Fext_constant instead of Fext_equation
	bool use_constant_force = true;

public:

	//////////////////////////////////////

	MElastic_Boundary(void) :
		Fext_equation({"x", "y", "t"})
	{}
	
	~MElastic_Boundary()
	{
	}

	////////////////////////////////////// SETUP

	void setup_surface(const VEC_VC<DBL3>& u_disp, Rect surf_rect_)
	{
		surf_rect = surf_rect_; 
		h_m = u_disp.h;

		Rect intersection = u_disp.rect.get_intersection(surf_rect);

		//y-z plane
		if (IsZ(intersection.s.x - intersection.e.x)) {

			//on lower x side
			if (IsZ(u_disp.rect.s.x - intersection.s.x)) {

				intersection.e.x += h_m.x;
				surf_orientation = -1;
			}
			//on upper x side
			else if (IsZ(u_disp.rect.e.x - intersection.s.x)) {

				intersection.s.x -= h_m.x;
				surf_orientation = +1;
			}
		}
		//x-z plane
		else if (IsZ(intersection.s.y - intersection.e.y)) {

			//on lower y side
			if (IsZ(u_disp.rect.s.y - intersection.s.y)) {

				intersection.e.y += h_m.y;
				surf_orientation = -2;
			}
			//on upper y side
			else if (IsZ(u_disp.rect.e.y - intersection.s.y)) {

				intersection.s.y -= h_m.y;
				surf_orientation = +2;
			}
		}
		//x-y plane
		else if (IsZ(intersection.s.z - intersection.e.z)) {

			//on lower z side
			if (IsZ(u_disp.rect.s.z - intersection.s.z)) {

				intersection.e.z += h_m.z;
				surf_orientation = -3;
			}
			//on upper z side
			else if (IsZ(u_disp.rect.e.z - intersection.s.z)) {

				intersection.s.z -= h_m.z;
				surf_orientation = +3;
			}
		}

		cells_box = u_disp.box_from_rect_max(intersection);
	}

	////////////////////////////////////// EQUATION STIMULUS

	//setup stimulus using a user supplied text equation
	bool setup_equation_stimulus(std::string eq_text, std::vector<std::pair<std::string, double>> constants)
	{
		bool success = Fext_equation.make_from_string(eq_text, constants);
		if (success) use_constant_force = false;

		return success;
	}

	void set_constants(std::vector<std::pair<std::string, double>> constants)
	{
		Fext_equation.set_constants(constants);
	}

	////////////////////////////////////// FIXED STIMULUS

	//setup a fixed stimulus
	void setup_fixed_stimulus(DBL3 Fext) 
	{
		Fext_constant = Fext;
		use_constant_force = true;
		
		Fext_equation.clear();
	}

	////////////////////////////////////// SEARCH (runtime)

	//does cells_box contain cell with given i, j, k indexes? (indexed on the same mesh resolution, i.e. at the h_m, n_m mesh resolution)
	//return 0 if it doesn't, else return the corresponding surface orientiation
	int contains(const INT3& ijk_cell_centre) 
	{ 
		if (cells_box.Contains(ijk_cell_centre)) return surf_orientation;
		else return 0;
	}

	////////////////////////////////////// GET (runtime)

	//get external force (Fx, Fy, Fz) at vertex with given index
	//this is used for diagonal components only when solving for stress
	//(sxx = +/-Fx at x normal face, syy = +/-Fy at y normal face, szz = +/-Fz at z normal face)
	//external force along surface normal returned (specified by user in xyz system)
	double get_ext_force_vertices(const INT3& ijk_vertex, const double& time)
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
			double x, y;

			switch (surf_orientation) {

				//-x (so surface is in yz plane). origin at +y, -z
			case -1:
				x = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return -Fext_equation.evaluate_vector_x(x, y, time);
				break;

				//+x (so surface is in yz plane). origin at -y, -z
			case +1:
				x = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return Fext_equation.evaluate_vector_x(x, y, time);
				break;

				//-y (so surface is in xz plane). origin at -x, -z
			case -2:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return -Fext_equation.evaluate_vector_y(x, y, time);
				break;

				//+y (so surface is in xz plane). origin at +x, -z
			case +2:
				x = (cells_box.e.i - ijk_vertex.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return Fext_equation.evaluate_vector_y(x, y, time);
				break;

				//-z (so surface is in xy plane). origin at -x, +y
			case -3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				return -Fext_equation.evaluate_vector_z(x, y, time);
				break;

				//+z (so surface is in xy plane). origin at -x, -y
			case +3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				return Fext_equation.evaluate_vector_z(x, y, time);
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
	DBL3 get_ext_force_edges(const INT3& ijk_vertex, const double& time)
	{
		if (use_constant_force) {

			//Fixed external surface force
			switch (surf_orientation) {

			case -1:
				return DBL3(0.0, -Fext_constant.y, -Fext_constant.z);
				break;

			case 1:
				return DBL3(0.0, Fext_constant.y, Fext_constant.z);
				break;

			case -2:
				return DBL3(-Fext_constant.x, 0.0, -Fext_constant.z);
				break;

			case 2:
				return DBL3(Fext_constant.x, 0.0, Fext_constant.z);
				break;

			case -3:
				return DBL3(-Fext_constant.x, -Fext_constant.y, 0.0);
				break;

			case 3:
				return DBL3(Fext_constant.x, Fext_constant.y, 0.0);
				break;

			default:
				return DBL3();
				break;
			}
		}
		else {

			//External surface force from text equation

			//x, y vertex position relative to surf_rect
			double x, y;

			switch (surf_orientation) {

				//-x (so surface is in yz plane). origin at +y, -z
			case -1:
				x = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return DBL3(0.0, -Fext_equation.evaluate_vector_y(x + h_m.y / 2, y, time), -Fext_equation.evaluate_vector_z(x, y + h_m.z / 2, time));
				break;

				//+x (so surface is in yz plane). origin at -y, -z
			case +1:
				x = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return DBL3(0.0, Fext_equation.evaluate_vector_y(x + h_m.y / 2, y, time), Fext_equation.evaluate_vector_z(x, y + h_m.z / 2, time));
				break;

				//-y (so surface is in xz plane). origin at -x, -z
			case -2:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return DBL3(-Fext_equation.evaluate_vector_x(x + h_m.x / 2, y, time), 0.0, -Fext_equation.evaluate_vector_z(x, y + h_m.z / 2, time));
				break;

				//+y (so surface is in xz plane). origin at +x, -z
			case +2:
				x = (cells_box.e.i - ijk_vertex.i) * h_m.x;
				y = (ijk_vertex.k - cells_box.s.k) * h_m.z;
				return DBL3(Fext_equation.evaluate_vector_x(x + h_m.x / 2, y, time), 0.0, Fext_equation.evaluate_vector_z(x, y + h_m.z / 2, time));
				break;

				//-z (so surface is in xy plane). origin at -x, +y
			case -3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (cells_box.e.j - ijk_vertex.j) * h_m.y;
				return DBL3(-Fext_equation.evaluate_vector_x(x + h_m.x / 2, y, time), -Fext_equation.evaluate_vector_y(x, y + h_m.y / 2, time), 0.0);
				break;

				//+z (so surface is in xy plane). origin at -x, -y
			case +3:
				x = (ijk_vertex.i - cells_box.s.i) * h_m.x;
				y = (ijk_vertex.j - cells_box.s.j) * h_m.y;
				return DBL3(Fext_equation.evaluate_vector_x(x + h_m.x / 2, y, time), Fext_equation.evaluate_vector_y(x, y + h_m.y / 2, time), 0.0);
				break;

			default:
				return DBL3();
				break;
			}
		}
	}

	////////////////////////////////////// INFO

	const Rect& get_surface(void) const { return surf_rect; }
	
	const Box& get_box(void) const { return cells_box; }
	const DBL3& get_cellsize(void) const { return h_m; }
	const int& get_orientation(void) const { return surf_orientation; }
	
	bool is_constant_force(void) const { return use_constant_force; }
	bool is_equation_set(void) const { return Fext_equation.is_set(); }
	DBL3 get_constant_force(void) const { return Fext_constant; }
	std::vector<std::vector< std::vector<EqComp::FSPEC> >> get_equation_fspec(void) { return Fext_equation.get_vector_fspec(); }
};

#endif