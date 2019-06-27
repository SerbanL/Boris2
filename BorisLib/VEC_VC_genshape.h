#pragma once

#include "VEC_VC.h"
#include "prng.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE GENERATORS : VEC_VC_genshape.h, VEC_VC_Voronoi.h

//roughen a mesh side (side = "-x", "x", "-y", "y", "-z", "z") to given depth (same units as h) with prng instantiated with given seed
template <typename VType>
bool VEC_VC<VType>::generate_roughside(std::string side, double depth, unsigned seed)
{
	BorisRand prng(seed);

	Rect side_rect_abs;

	if (side == "x") side_rect_abs = rect.get_face_px(h.x);
	else if (side == "-x") side_rect_abs = rect.get_face_mx(h.x);
	else if (side == "y") side_rect_abs = rect.get_face_py(h.y);
	else if (side == "-y") side_rect_abs = rect.get_face_my(h.y);
	else if (side == "z") side_rect_abs = rect.get_face_pz(h.z);
	else if (side == "-z") side_rect_abs = rect.get_face_mz(h.z);

	Box cells_box = box_from_rect_min(side_rect_abs); 

	for (int i = cells_box.s.x; i < cells_box.e.x; i++) {
		for (int j = cells_box.s.y; j < cells_box.e.y; j++) {
			for (int k = cells_box.s.z; k < cells_box.e.z; k++) {

				if (side == "x") {

					int cells_cut = round(prng.rand() * (depth / h.x));

					for (int idx = 0; idx < (cells_cut < n.x ? cells_cut : n.x); idx++) {

						mark_empty((i - idx) + j * n.x + k * n.x*n.y);
					}
				}

				else if (side == "-x") {

					int cells_cut = round(prng.rand() * (depth / h.x));

					for (int idx = 0; idx < (cells_cut < n.x ? cells_cut : n.x); idx++) {

						mark_empty((i + idx) + j * n.x + k * n.x*n.y);
					}
				}

				else if (side == "y") {

					int cells_cut = round(prng.rand() * (depth / h.y));

					for (int idx = 0; idx < (cells_cut < n.y ? cells_cut : n.y); idx++) {

						mark_empty(i + (j - idx) * n.x + k * n.x*n.y);
					}
				}

				else if (side == "-y") {

					int cells_cut = round(prng.rand() * (depth / h.y));

					for (int idx = 0; idx < (cells_cut < n.y ? cells_cut : n.y); idx++) {

						mark_empty(i + (j + idx) * n.x + k * n.x*n.y);
					}
				}

				else if (side == "z") {

					int cells_cut = round(prng.rand() * (depth / h.z));

					for (int idx = 0; idx < (cells_cut < n.z ? cells_cut : n.z); idx++) {

						mark_empty(i + j * n.x + (k - idx) * n.x*n.y);
					}
				}

				else if (side == "-z") {

					int cells_cut = round(prng.rand() * (depth / h.z));

					for (int idx = 0; idx < (cells_cut < n.z ? cells_cut : n.z); idx++) {

						mark_empty(i + j * n.x + (k + idx) * n.x*n.y);
					}
				}
			}
		}
	}

	set_ngbrFlags();

	return true;
}

//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
//Rough both top and bottom if sides is empty, else it should be either -z or z.
template <typename VType>
bool VEC_VC<VType>::generate_jagged_surfaces(double depth, double spacing, unsigned seed, std::string sides)
{
	//this VEC may already have a shape and we want to roughen the shape as it is
	//for this routine to work the jagged surface must be generated on a full rectangle shape so:
	//1. save current shape
	//2. reset shape to full rectangle
	//3. generate jagged surface
	//4. mask with original shape
	//5. recalculate flags

	//cut top if sides is "z" or empty
	bool cut_top = (!sides.length() || sides == "z");

	//cut bottom if sides is "-z" or empty
	bool cut_bot = (!sides.length() || sides == "-z");

	//1.
	vector<int> save_ngbrFlags = ngbrFlags;

	//2.
	ngbrFlags.assign(n.dim(), NF_NOTEMPTY);

	//3.
	BorisRand prng(seed);

	//cut given number of cells from starting index, either from -z face (cells_cut > 0) or from +z face (cells_cut < 0)
	auto cut_cells = [&](int cells_cut, INT3 starting_idx) -> void {

		for (int cut_idx = 0; cut_idx < (mod(cells_cut) < n.z ? mod(cells_cut) : n.z); cut_idx++) {

			INT3 cut_int3_idx = starting_idx + get_sign(cells_cut) * INT3(0, 0, cut_idx);

			mark_empty(cut_int3_idx.x + cut_int3_idx.y * n.x + cut_int3_idx.z * n.x*n.y);
		}
	};

	//cout number of empty cells along the z aixs from starting index, either from -z face (direction > 0) or from +z face (direction < 0)
	auto count_cells_cut = [&](INT3 starting_idx, int direction) -> int {

		int cells_cut = 0;

		for (int count_idx = 0; count_idx < n.z; count_idx++) {

			if (is_empty(starting_idx + get_sign(direction) * INT3(0, 0, count_idx))) cells_cut++;
			else return cells_cut;
		}

		return cells_cut;
	};

	int cells_spacing_x = round(spacing / h.x);
	int cells_spacing_y = round(spacing / h.y);

	//first generate random cut values at given square spacing
#pragma omp parallel for
	for (int j = 0; j < n.y; j += cells_spacing_y) {

		for (int i = 0; i < n.x; i += cells_spacing_x) {
			
			if (cut_top) cut_cells(-round(prng.rand() * (depth / h.z)), INT3(i, j, n.z - 1));
			if (cut_bot) cut_cells(+round(prng.rand() * (depth / h.z)), INT3(i, j, 0));

			//last point on this column must also have a value generated (only generate it once, e.g. when j == 0 easiest)
			if (j == 0) {

				if (cut_top) cut_cells(-round(prng.rand() * (depth / h.z)), INT3(i, (n.y - 1), n.z - 1));
				if (cut_bot) cut_cells(+round(prng.rand() * (depth / h.z)), INT3(i, (n.y - 1), 0));
			}
		}

		//last point on this row must also have a value generated
		if (cut_top) cut_cells(-round(prng.rand() * (depth / h.z)), INT3(n.x - 1, j, n.z - 1));
		if (cut_bot) cut_cells(+round(prng.rand() * (depth / h.z)), INT3(n.x - 1, j, 0));
	}

	//upper-right point
	if (cut_top) cut_cells(-round(prng.rand() * (depth / h.z)), INT3(n.x - 1, (n.y - 1), n.z - 1));
	if (cut_bot) cut_cells(+round(prng.rand() * (depth / h.z)), INT3(n.x - 1, (n.y - 1), 0));

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
					int cells_cut11_top = count_cells_cut(INT3(i1, j1, n.z - 1), -1);
					int cells_cut12_top = count_cells_cut(INT3(i1, j2, n.z - 1), -1);
					int cells_cut21_top = count_cells_cut(INT3(i2, j1, n.z - 1), -1);
					int cells_cut22_top = count_cells_cut(INT3(i2, j2, n.z - 1), -1);
					int cells_cut11_bot = count_cells_cut(INT3(i1, j1, 0), +1);
					int cells_cut12_bot = count_cells_cut(INT3(i1, j2, 0), +1);
					int cells_cut21_bot = count_cells_cut(INT3(i2, j1, 0), +1);
					int cells_cut22_bot = count_cells_cut(INT3(i2, j2, 0), +1);

					int cells_cut = interpolate_bilinear(x1, y1, x2, y2, cells_cut11_top, cells_cut12_top, cells_cut21_top, cells_cut22_top, DBL2(x, y));

					if (cut_top) cut_cells(-cells_cut, INT3(ii, jj, n.z - 1));

					cells_cut = interpolate_bilinear(x1, y1, x2, y2, cells_cut11_bot, cells_cut12_bot, cells_cut21_bot, cells_cut22_bot, DBL2(x, y));

					if (cut_bot) cut_cells(+cells_cut, INT3(ii, jj, 0));
				}
			}
		}
	}

	//4.
#pragma omp parallel for
	for (int idx = 0; idx < save_ngbrFlags.size(); idx++) {

		if (!(save_ngbrFlags[idx] & NF_NOTEMPTY)) mark_empty(idx);
	}

	//5.
	set_ngbrFlags();

	return true;
}