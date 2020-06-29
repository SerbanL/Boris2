#pragma once

#include "VEC.h"
#include "BLib_prng.h"

//--------------------------------------------VEC GENERATORS

//generate custom values from grayscale bitmap : black = 0, white = 1. Apply scaling and offset also.
//bitmap size must match n.x * n.y obtained from new_h and new_rect
template <>
inline bool VEC<double>::generate_custom_2D(SZ3 new_n, Rect new_rect, double offset, double scale, const std::vector<unsigned char>& bitmap)
{
	DBL3 new_h = new_rect / new_n;

	if (!resize(new_h, new_rect)) return false;

	//bitmap must have the right size (i.e. have n.x * n.y pixels, remembering each pixel has 4 bytes as B-G-R-A)
	if (bitmap.size() != n.x*n.y * 4) return false;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			//In the file the 0, 0 point is in the left-top corner, but in the mesh the 0, 0 point is in the left-bottom corner
			int indexBitMap = i + (n.y - 1 - j)*n.x;
			int indexMesh = i + j * n.x;

			double B = (double)bitmap[indexBitMap * 4] / 255;
			//double G = (double)bitmap[indexBitMap * 4 + 1] / 255;
			//double R = (double)bitmap[indexBitMap * 4 + 2] / 255;

			quantity[indexMesh] = offset + scale * B;
		}
	}

	return true;
}

template <>
inline void VEC<double>::set_linear(DBL3 position1, double value1, DBL3 position2, double value2)
{
	if (position1 == position2) return;

	DBL3 del_p = position2 - position1;
	double del_p_sq = del_p * del_p;

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		//the absolute position of this cell center for which we are setting the value
		DBL3 P = cellidx_to_position(idx) + rect.s;

		//p is the point on the line containing position1 and position2, such that the segment (P - p) is orthogonal to it
		DBL3 p = (del_p * (P - position1)) * del_p / del_p_sq + position1;

		//use linear interpolation along the line containing position1 and position2 to set value for this cell
		quantity[idx] = (value1 * (position2 - p).norm() + value2 * (p - position1).norm()) / sqrt(del_p_sq);
	}
}

//linear : use interpolation to set values in this VEC based on projected distance between position1 and position2 and given fixed end values.
template <>
inline bool VEC<double>::generate_linear(DBL3 new_h, Rect new_rect, DBL3 position1, double value1, DBL3 position2, double value2)
{
	if (!resize(new_h, new_rect)) return false;

	set_linear(position1, value1, position2, value2);

	return true;
}

//random: set VEC dimensions and generate random values in given range (prng instantiated with given seed)
template <typename VType>
bool VEC<VType>::generate_random(DBL3 new_h, Rect new_rect, DBL2 range, unsigned seed)
{
	if (!resize(new_h, new_rect)) return false;

	BorisRand prng(seed);

	//Don't use parallel loop here, otherwise the random values will be different on different computers with differing number of cores
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] = VType(prng.rand() * (range.j - range.i) + range.i);
	}

	return true;
}

//jagged: set VEC dimensions (force 2D in xy plane) and generate random values in given range (prng instantiated with given seed) at a given spacing. 
//In between these random values use bi-linear interpolation. The random values are spaced in the xy plane at equal distances along x or y using the spacing value (same units as the VEC rect)
template <>
inline bool VEC<double>::generate_jagged(DBL3 new_h, Rect new_rect, DBL2 range, double spacing, unsigned seed)
{
	//force 2D in xy plane
	new_h.k = new_rect.size().k;

	if (!resize(new_h, new_rect) || !spacing) return false;

	BorisRand prng(seed);

	int cells_spacing_x = round(spacing / h.x);
	int cells_spacing_y = round(spacing / h.y);

	//first generate random peak values at given square spacing
	//Don't use parallel loop here, otherwise the random values will be different on different computers with differing number of cores
	for (int j = 0; j < n.y; j += cells_spacing_y) {

		for (int i = 0; i < n.x; i += cells_spacing_x) {

			int idx = i + j * n.x;

			quantity[idx] = prng.rand() * (range.j - range.i) + range.i;

			//last point on this column must also have a value generated (only generate it once, e.g. when j == 0 easiest)
			if (j == 0) quantity[i + (n.y - 1) * n.x] = prng.rand() * (range.j - range.i) + range.i;
		}

		//last point on this row must also have a value generated
		quantity[n.x - 1 + j*n.x] = prng.rand() * (range.j - range.i) + range.i;
	}

	//upper-right point
	quantity[n.x - 1 + (n.y - 1) * n.x] = prng.rand() * (range.j - range.i) + range.i;

	//now fill in the rest using bi-linear interpolation
#pragma omp parallel for
	for (int j = 0; j < n.y; j += cells_spacing_y) {

		//the row indexes and positions where we've generated random points
		int j1 = j;
		int j2 = (j + cells_spacing_y < n.y ? j + cells_spacing_y : n.y - 1);

		double y1 = j1 * h.y;
		double y2 = j2 * h.y;

		for (int i = 0; i < n.x; i += cells_spacing_x) {

			//the column indexes and positions where we've generated random points
			int i1 = i;
			int i2 = (i + cells_spacing_x < n.x ? i + cells_spacing_x : n.x - 1);

			double x1 = i1 * h.x;
			double x2 = i2 * h.x;

			//now loop over all points in between the 4 outlined above
			for (int jj = j1; jj <= j2; jj++) {
				
				//position
				double y = jj * h.y;

				for (int ii = i1; ii <= i2; ii++) {

					if ((ii == i1 && jj == j1) || (ii == i1 && jj == j2) || (ii == i2 && jj == j1) || (ii == i2 && jj == j2)) continue;

					//position
					double x = ii * h.x;

					//generate using bilinear interpolation
					quantity[ii + jj * n.x] = interpolate_bilinear(x1, y1, x2, y2, quantity[i1 + j1 * n.x], quantity[i1 + j2 * n.x], quantity[i2 + j1 * n.x], quantity[i2 + j2 * n.x], DBL2(x, y));
				}
			}
		}
	}

	return true;
}

//defects: generate circular defects with a tanh radial profile with values in the given range, diameter range and average spacing (prng instantiated with given seed). The defect positioning is random. 
template <>
inline bool VEC<double>::generate_defects(DBL3 new_h, Rect new_rect, DBL2 range, double base_value, DBL2 diameter_range, double spacing, unsigned seed)
{
	//force 2D in xy plane
	new_h.k = new_rect.size().k;

	if (!resize(new_h, new_rect) || !spacing || diameter_range.i <= 0 || diameter_range.j <= 0) return false;

	BorisRand prng(seed);

	//set base value everywhere before generating defects
	set(base_value);

	//generate a circular defect at position (relative) with tanh radial profile (ranges from outside value of 1 to given center value) for given radius
	auto set_defect = [&](double center_value, double radius, DBL2 center_position) -> void {

		//relative rectangle in which to generate defect
		Rect defect_rect = Rect(DBL3(center_position.x - radius, center_position.y - radius, 0), DBL3(center_position.x + radius, center_position.y + radius, rect.e.z-rect.s.z));

		//get start and end cell indexes
		INT3 start = round(defect_rect.s / h);
		INT3 end = round(defect_rect.e / h);

		//go over all cells in the defect box which are contained in the VEC
		for (int i = (start.i >= 0 ? start.i : 0); i < end.i && i < n.x; i++) {
			for (int j = (start.j >= 0 ? start.j : 0); j < end.j && j < n.y; j++) {

				int idx = i + j * n.x;

				//distance from origin
				DBL2 position = DBL2(((double)i + 0.5) * h.x, ((double)j + 0.5) * h.y);
				double distance = (center_position - position).norm();

				double TANH_STRETCH = 3;

				double current_value = quantity[idx];
				double new_value = tanh(((distance / radius) - 0.5) * TANH_STRETCH) * (1 - center_value) / 2 + (center_value + 1) / 2;

				//blend in new value since defects may overlap
				quantity[idx] = current_value * (distance / radius) + new_value * (1 - (distance / radius));
			}
		}
	};

	//spacing determines defect density
	int num_defects = round((rect.e.x - rect.s.x) * (rect.e.y - rect.s.y) / (spacing * spacing));

	for (int idx_defect = 0; idx_defect < num_defects; idx_defect++) {

		double center_value = prng.rand() * (range.j - range.i) + range.i;
		double radius = (prng.rand() * (diameter_range.j - diameter_range.i) + diameter_range.i) / 2;
		
		DBL2 center_position = DBL2(prng.rand() * (rect.e.x - rect.s.x), prng.rand() * (rect.e.y - rect.s.y));

		set_defect(center_value, radius, center_position);
	}

	return true;
}

//faults: set VEC dimensions (force 2D in xy plane) and generate line faults in the given range length, orientation length (degrees azimuthal) and average spacing (prng instantiated with given seed).
template <>
inline bool VEC<double>::generate_faults(DBL3 new_h, Rect new_rect, DBL2 range, double base_value, DBL2 length_range, DBL2 orientation_range, double spacing, unsigned seed)
{
	//force 2D in xy plane
	new_h.k = new_rect.size().k;

	if (!resize(new_h, new_rect) || !spacing || length_range <= 0) return false;

	BorisRand prng(seed);

	//set base value everywhere before generating defects
	set(base_value);

	auto set_fault = [&](DBL3 s, DBL3 e, double value) -> void {

		double w = maximum(h.x, h.y, (e - s).norm()*0.05);

		DBL2 s_rhs = DBL2(s.x, s.y) + DBL2(w, 0);
		DBL2 s_lhs = DBL2(s.x, s.y) + DBL2(0, w);
		DBL2 e_rhs = DBL2(e.x, e.y) - DBL2(0, w);
		DBL2 e_lhs = DBL2(e.x, e.y) - DBL2(w, 0);

		//get start and end cell indexes
		INT3 start = round(s / h);
		INT3 end = round(e / h);

		//go over all cells in the defect box which are contained in the VEC
		for (int i = (start.i >= 0 ? start.i : 0); i < end.i && i < n.x; i++) {
			for (int j = (start.j >= 0 ? start.j : 0); j < end.j && j < n.y; j++) {

				int idx = i + j * n.x;

				DBL2 position = DBL2(((double)i + 0.5) * h.x, ((double)j + 0.5) * h.y);

				if (point_on_rhs_of_line(s_lhs, e_lhs, position) && !point_on_rhs_of_line(s_rhs, e_rhs, position)) {

					quantity[idx] = value;
				}
			}
		}
	};

	//spacing determines defect density
	int num_faults = round((rect.e.x - rect.s.x) * (rect.e.y - rect.s.y) / (spacing * spacing));

	for (int idx_fault = 0; idx_fault < num_faults; idx_fault++) {

		double value = prng.rand() * (range.j - range.i) + range.i;
		double length = prng.rand() * (length_range.j - length_range.i) + length_range.i;
		double orientation_rads = (prng.rand() * (orientation_range.j - orientation_range.i) + orientation_range.i) * PI / 180.0;
		DBL3 center_position = DBL3(prng.rand() * (rect.e.x - rect.s.x), prng.rand() * (rect.e.y - rect.s.y), 0);

		//build unit vector in direction of orientation
		DBL3 n = DBL3(cos(orientation_rads), sin(orientation_rads), 0);
		
		//now generate end points for the fault line
		DBL3 point_plus = center_position + n * length/2;
		DBL3 point_minus = center_position - n * length/2;

		DBL3 start, end;
		
		if (point_plus >= point_minus) {

			start = point_minus;
			end = point_plus;
		}
		else {

			start = point_plus;
			end = point_minus;
		}

		set_fault(start, end, value);
	}

	return true;
}