#pragma once

#include "VEC_VC.h"

//CMBNDInfo describes a composite media boundary contact between 2 meshes of same type, used to calculate values at CMBND cells using boundary conditions
struct CMBNDInfo {

	//index of contacting meshes as : (secondary mesh, primary mesh)
	INT2 mesh_idx;

	//primary cellsize shift perpendicular to contact (from first primary cell-center position add half this DBL3 to reach the contact interface)
	DBL3 hshift_primary;

	//secondary cellsize shift perpendicular to contact (from contact interface position add half this DBL3 to reach the first position to read value from on the contacting mesh, then add the full DBL3 to reach the second position)
	DBL3 hshift_secondary;

	//weights depending on primary and secondary cellsizes perpendicular to the interface; in particular : i: max{h_primary, h_secondary} / h_secondary and j: max{h_primary, h_secondary} / h_primary
	DBL2 weights;

	//integer cell shift on primary side (add full INT3 to reach second cell on primary side to read value from) -> if positive then primary side is on +ve side, else on -ve side.
	INT3 cell_shift;

	//Box containing all cells on primary side of contact (though not all may be CMBND cells so check flag before setting boundary conditions)
	Box cells_box;

	//check if primary is on +ve or -ve side of contact (e.g. top or bottom)
	bool IsPrimaryTop(void) { return (cell_shift >= INT3(0)); }
};

//-------------------------------- CALCULATE CONTACTS and SET FLAGS

//set cmbnd flags by identifying contacts with other vecs
template <typename VType>
std::vector<CMBNDInfo> VEC_VC<VType>::set_cmbnd_flags(int primary_mesh_idx, std::vector<VEC_VC<VType>*> &pVECs, bool check_neighbors)
{
	std::vector<CMBNDInfo> contacts;

	CMBNDInfo contact;

	//first clear all NF_CMBND flags before recounting them
	clear_cmbnd_flags();

	//check each mesh against this
	for (int check_mesh = 0; check_mesh < (int)pVECs.size(); check_mesh++) {

		Rect check_rect = pVECs[check_mesh]->rect;

		Rect mesh_intersection = VEC<VType>::rect.get_intersection(check_rect);

		//intersection must be exactly a plane
		if (mesh_intersection.IsPlane()) {

			int flag_value;

			//store contacting meshes indexes in pVECs
			contact.mesh_idx = INT2(check_mesh, primary_mesh_idx);

			//y-z plane : x is the perpendicular direction
			if (IsZ(mesh_intersection.s.x - mesh_intersection.e.x)) {

				//is the selected mesh on the positive or negative side of the boundary? Set flag to use. Also adjust mesh_intersection so it contains all the cells to set flags for.
				if (IsZ(VEC<VType>::rect.s.x - mesh_intersection.s.x)) {

					flag_value = NF_CMBNDPX;
					mesh_intersection.e.x += VEC<VType>::h.x;
					contact.hshift_primary = DBL3(-VEC<VType>::h.x, 0, 0);
					contact.hshift_secondary = DBL3(-pVECs[check_mesh]->h.x, 0, 0);
					contact.cell_shift = INT3(1, 0, 0);
				}
				else {

					flag_value = NF_CMBNDNX;
					mesh_intersection.s.x -= VEC<VType>::h.x;
					contact.hshift_primary = DBL3(+VEC<VType>::h.x, 0, 0);
					contact.hshift_secondary = DBL3(+pVECs[check_mesh]->h.x, 0, 0);
					contact.cell_shift = INT3(-1, 0, 0);
				}

				contact.weights.i = (VEC<VType>::h.x > pVECs[check_mesh]->h.x ? VEC<VType>::h.x : pVECs[check_mesh]->h.x) / pVECs[check_mesh]->h.x;		//L or secondary mesh weight
				contact.weights.j = (VEC<VType>::h.x > pVECs[check_mesh]->h.x ? VEC<VType>::h.x : pVECs[check_mesh]->h.x) / VEC<VType>::h.x;			//R or primary mesh weight
			}
			//x-z plane : y is the perpendicular direction
			else if (IsZ(mesh_intersection.s.y - mesh_intersection.e.y)) {

				//is the selected mesh on the positive or negative side of the boundary? Set flag to use. Also adjust mesh_intersection so it contains all the cells to set flags for.
				if (IsZ(VEC<VType>::rect.s.y - mesh_intersection.s.y)) {

					flag_value = NF_CMBNDPY;
					mesh_intersection.e.y += VEC<VType>::h.y;
					contact.hshift_primary = DBL3(0, -VEC<VType>::h.y, 0);
					contact.hshift_secondary = DBL3(0, -pVECs[check_mesh]->h.y, 0);
					contact.cell_shift = INT3(0, 1, 0);
				}
				else {

					flag_value = NF_CMBNDNY;
					mesh_intersection.s.y -= VEC<VType>::h.y;
					contact.hshift_primary = DBL3(0, +VEC<VType>::h.y, 0);
					contact.hshift_secondary = DBL3(0, +pVECs[check_mesh]->h.y, 0);
					contact.cell_shift = INT3(0, -1, 0);
				}

				contact.weights.i = (VEC<VType>::h.y > pVECs[check_mesh]->h.y ? VEC<VType>::h.y : pVECs[check_mesh]->h.y) / pVECs[check_mesh]->h.y;
				contact.weights.j = (VEC<VType>::h.y > pVECs[check_mesh]->h.y ? VEC<VType>::h.y : pVECs[check_mesh]->h.y) / VEC<VType>::h.y;
			}
			//x-y plane : z is the perpendicular direction
			else if (IsZ(mesh_intersection.s.z - mesh_intersection.e.z)) {

				//is the selected mesh on the positive or negative side of the boundary? Set flag to use. Also adjust mesh_intersection so it contains all the cells to set flags for.
				if (IsZ(VEC<VType>::rect.s.z - mesh_intersection.s.z)) {

					flag_value = NF_CMBNDPZ;
					mesh_intersection.e.z += VEC<VType>::h.z;
					contact.hshift_primary = DBL3(0, 0, -VEC<VType>::h.z);
					contact.hshift_secondary = DBL3(0, 0, -pVECs[check_mesh]->h.z);
					contact.cell_shift = INT3(0, 0, 1);
				}
				else {

					flag_value = NF_CMBNDNZ;
					mesh_intersection.s.z -= VEC<VType>::h.z;
					contact.hshift_primary = DBL3(0, 0, +VEC<VType>::h.z);
					contact.hshift_secondary = DBL3(0, 0, +pVECs[check_mesh]->h.z);
					contact.cell_shift = INT3(0, 0, -1);
				}

				contact.weights.i = (VEC<VType>::h.z > pVECs[check_mesh]->h.z ? VEC<VType>::h.z : pVECs[check_mesh]->h.z) / pVECs[check_mesh]->h.z;
				contact.weights.j = (VEC<VType>::h.z > pVECs[check_mesh]->h.z ? VEC<VType>::h.z : pVECs[check_mesh]->h.z) / VEC<VType>::h.z;
			}

			//box of cells in this primary mesh completely included in the intersection - not all will be marked as CMBND cells, only marked if:
			// 1. not empty
			// 2. cell next to it in primary mesh along boundary normal is also not empty (unless check_neighbors = false)
			// 3. space for cell on secondary mesh is completely contained in non-empty cells; 
			// "space for cell" is the rectangle with same footprint on the boundary as h, but thickness (i.e. along boundary normal) given by cell thickness from secondary mesh
			// 4. same as in 3. for space for cell one further cell thickness along the boundary normal (unless check_neighbors = false)

			contact.cells_box = VEC<VType>::box_from_rect_min(mesh_intersection);

			//new contact done
			contacts.push_back(contact);

			//now set flags
			INT3 ijk_start = contact.cells_box.s;
			INT3 ijk_end = contact.cells_box.e;

			for (int i = ijk_start.i; i < ijk_end.i; i++) {
				for (int j = ijk_start.j; j < ijk_end.j; j++) {
					for (int k = ijk_start.k; k < ijk_end.k; k++) {

						int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

						//check cell 1 on primary
						if (!(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

						//check cell 2 on primary (just next to cell 1 along interface perpendicular)
						if (check_neighbors) {

							int idx2 = idx + contact.cell_shift.x + contact.cell_shift.y*VEC<VType>::n.x + contact.cell_shift.z*VEC<VType>::n.x*VEC<VType>::n.y;
							if (idx2 >= VEC<VType>::n.dim() || !(ngbrFlags[idx2] & NF_NOTEMPTY)) continue;
						}

						//check cell 1 space on secondary
						DBL3 abspos = VEC<VType>::cellidx_to_position(idx) + VEC<VType>::rect.s + (contact.hshift_primary + contact.hshift_secondary) / 2;
						//imageStencil has the footprint of h but thickness set by the secondary mesh
						DBL3 imageStencil = VEC<VType>::h - mod(contact.hshift_primary) + mod(contact.hshift_secondary);
						//abscellRect is the space for cell 1 on secondary mesh (first cell space from the boundary)
						Rect abscellRect = Rect(abspos - imageStencil / 2, abspos + imageStencil / 2);
						if (!pVECs[check_mesh]->rect.contains(abscellRect) || !(pVECs[check_mesh]->is_not_empty(abscellRect))) continue;

						if (check_neighbors) {

							//check cell 2 space on secondary
							abscellRect = abscellRect + contact.hshift_secondary;
							if (!pVECs[check_mesh]->rect.contains(abscellRect) || !(pVECs[check_mesh]->is_not_empty(abscellRect))) continue;
						}

						ngbrFlags[idx] |= flag_value;
					}
				}
			}
		}
	}

	return contacts;
}

//-------------------------------- CONTINUOUS FLUX and POTENTIAL

template <typename VType>
template <typename Owner>
void VEC_VC<VType>::set_cmbnd_continuous(
	VEC_VC<VType> &V_sec, CMBNDInfo& contact,
	std::function<VType(const Owner&, DBL3, DBL3, DBL3)> a_func_sec, std::function<VType(const Owner&, int, int, DBL3)> a_func_pri,
	std::function<double(const Owner&, DBL3, DBL3, DBL3)> b_func_sec, std::function<double(const Owner&, int, int)> b_func_pri,
	std::function<VType(const Owner&, DBL3, DBL3, DBL3)> diff2_sec, std::function<VType(const Owner&, int, DBL3)> diff2_pri,
	Owner& instance_sec, Owner& instance_pri)
{
	INT3 box_sizes = contact.cells_box.size();

	//cellsizes perpendicular to interface
	double hL = contact.hshift_secondary.norm();
	double hR = contact.hshift_primary.norm();
	double hmax = (hL > hR ? hL : hR);

	//primary cells in this contact
#pragma omp parallel for
	for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		int cell1_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

		if (is_empty(cell1_idx) || is_not_cmbnd(cell1_idx)) continue;

		//calculate second primary cell index
		int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * VEC<VType>::n.x + (k + contact.cell_shift.k) * VEC<VType>::n.x*VEC<VType>::n.y;

		//cell values either side of the boundary: V_m2 V_m1 | V_1 V_2; positions as : -2 -1 | 1 2
		//V_m2 and V_2 are known. We need to set values V_m1 and V_1. Here we only set V_1. For this primary, secondary mesh pair there will be another one with order reversed, and there our V_m1 value will be set - so don't worry about it here!
		//NOTE : meshes at composite media boundaries must always be at least 2 non-empty cells thick in directions perpendicular to interface !!! -> this is actually checked when the list of contacts is made

		//relative position of cell -1 in secondary mesh
		DBL3 relpos_m1 = VEC<VType>::rect.s - V_sec.rect.s + ((DBL3(i, j, k) + DBL3(0.5)) & VEC<VType>::h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

		//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
		DBL3 stencil = VEC<VType>::h - mod(contact.hshift_primary) + mod(contact.hshift_secondary);

		//potential values at cells -2 and 2
		VType V_2 = VEC<VType>::quantity[cell2_idx];
		VType V_m2 = V_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

		//obtain a and b values used to define the flux as f(V) = a + b V', both on primary and secondary
		//a values
		VType a_val_sec = a_func_sec(instance_sec, relpos_m1, contact.hshift_secondary, stencil);
		VType a_val_pri = a_func_pri(instance_pri, cell1_idx, cell2_idx, contact.hshift_secondary);

		//b values adjusted with weights
		double b_val_sec = b_func_sec(instance_sec, relpos_m1, contact.hshift_secondary, stencil) * contact.weights.i;
		double b_val_pri = b_func_pri(instance_pri, cell1_idx, cell2_idx) * contact.weights.j;

		//V'' values at cell positions -1 and 1
		VType Vdiff2_sec = diff2_sec(instance_sec, relpos_m1, stencil, contact.hshift_secondary);
		VType Vdiff2_pri = diff2_pri(instance_pri, cell1_idx, contact.hshift_secondary);

		//Formula for V1
		VEC<VType>::quantity[cell1_idx] = (V_m2 * 2 * b_val_sec / 3 + V_2 * (b_val_pri + b_val_sec / 3)
			- Vdiff2_sec * b_val_sec * hL * hL - Vdiff2_pri * b_val_pri * hR * hR
			+ (a_val_pri - a_val_sec) * hmax) / (b_val_sec + b_val_pri);
	}
}

//-------------------------------- CONTINUOUS FLUX ONLY

//calculate cmbnd values based on continuity of flux only. The potential is allowed to drop across the interface as:
//f_sec(V) = f_pri(V) = A + B * delV, where f_sec and f_pri are the fluxes on the secondary and primary sides of the interface, and delV = V_pri - V_sec, the drop in potential across the interface.
//Thus in addition to the functions in set_cmbnd_continuous we need two extra sets of functions A_func_sec, A_func_pri returning VType and B_func_sec, B_func_pri returning double, s.t. A = (A_func_pri + A_func_sec), B = (B_func_pri + B_func_sec)
template <typename VType>
template <typename Owner, typename SOwner>
void VEC_VC<VType>::set_cmbnd_continuousflux(VEC_VC<VType> &V_sec, CMBNDInfo& contact,
	std::function<VType(const Owner&, DBL3, DBL3, DBL3)> a_func_sec, std::function<VType(const Owner&, int, int, DBL3)> a_func_pri,
	std::function<double(const Owner&, DBL3, DBL3, DBL3)> b_func_sec, std::function<double(const Owner&, int, int)> b_func_pri,
	std::function<VType(const Owner&, DBL3, DBL3, DBL3)> diff2_sec, std::function<VType(const Owner&, int, DBL3)> diff2_pri,
	std::function<VType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> A_func, std::function<double(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> B_func,
	Owner& instance_sec, Owner& instance_pri, SOwner& instance_s)
{
	INT3 box_sizes = contact.cells_box.size();

	//cellsizes perpendicular to interface
	double hL = contact.hshift_secondary.norm();
	double hR = contact.hshift_primary.norm();
	double hmax = (hL > hR ? hL : hR);

	//primary cells in this contact
#pragma omp parallel for
	for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		int cell1_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

		if (is_empty(cell1_idx) || is_not_cmbnd(cell1_idx)) continue;

		//calculate second primary cell index
		int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * VEC<VType>::n.x + (k + contact.cell_shift.k) * VEC<VType>::n.x*VEC<VType>::n.y;

		//cell values either side of the boundary: V_m2 V_m1 | V_1 V_2; positions as : -2 -1 | 1 2
		//V_m2 and V_2 are known. We need to set values V_m1 and V_1. Here we only set V_1. For this primary, secondary mesh pair there will be another one with order reversed, and there our V_m1 value will be set - so don't worry about it here!
		//NOTE : meshes at composite media boundaries must always be at least 2 non-empty cells thick in directions perpendicular to interface !!! -> this is actually checked when the list of contacts is made

		//relative position of cell -1 in secondary mesh
		DBL3 relpos_m1 = VEC<VType>::rect.s - V_sec.rect.s + ((DBL3(i, j, k) + DBL3(0.5)) & VEC<VType>::h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

		//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
		DBL3 stencil = VEC<VType>::h - mod(contact.hshift_primary) + mod(contact.hshift_secondary);

		//potential values at cells -2 and 2
		VType V_2 = VEC<VType>::quantity[cell2_idx];
		VType V_m2 = V_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

		//obtain a and b values used to define the flux as f(V) = a + b V', both on primary and secondary
		//a values
		VType a_val_sec = a_func_sec(instance_sec, relpos_m1, contact.hshift_secondary, stencil);
		VType a_val_pri = a_func_pri(instance_pri, cell1_idx, cell2_idx, contact.hshift_secondary);

		//b values adjusted with weights
		double b_val_sec = b_func_sec(instance_sec, relpos_m1, contact.hshift_secondary, stencil) * contact.weights.i;
		double b_val_pri = b_func_pri(instance_pri, cell1_idx, cell2_idx) * contact.weights.j;

		//V'' values at cell positions -1 and 1
		VType Vdiff2_sec = diff2_sec(instance_sec, relpos_m1, stencil, contact.hshift_secondary);
		VType Vdiff2_pri = diff2_pri(instance_pri, cell1_idx, contact.hshift_secondary);

		//A value
		VType A_val = A_func(instance_s, cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, instance_sec, instance_pri);

		//B value, then use it to find G value
		double B_val = B_func(instance_s, cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, instance_sec, instance_pri);
		
		double G_val = 0.0;
		if (B_val) G_val = 2 * b_val_pri * b_val_sec / (3 * B_val * hmax);

		//Formula for V1 - note, this reduces to the continuous case for G_val = 0 (or B_val tends to infinity) and A_val = VType(0)
		VEC<VType>::quantity[cell1_idx] = (V_m2 * 2 * b_val_sec / 3 + V_2 * (b_val_pri + b_val_sec / 3 + G_val)
			- Vdiff2_sec * b_val_sec * hL * hL - Vdiff2_pri * (b_val_pri + G_val) * hR * hR
			+ (a_val_pri * (1 + G_val / b_val_pri) - a_val_sec) * hmax
			- (G_val / b_val_pri) * A_val * hmax) / (b_val_sec + b_val_pri + G_val);
	}
}

//-------------------------------- DISCONTINUOUS FLUX and POTENTIAL

//most general case of composite media boundary conditions
//calculate cmbnd values based on set boundary flux values; both the flux and potential is allowed to be discontinuous across the interface.
//Fluxes at the interface are specified as: f_sec(V) = A_sec + B_sec * delVs, f_pri(V) = A_pri + B_pri * delVs, with directions from secondary to primary
//B functions may return a double, or a DBL33 (3x3 matrix) if VType is DBL3 (Cartesian vector).
//delVs = c_pri * V_pri - c_sec * V_sec, c are double values specified by given functions
template <typename VType>
template <typename Owner, typename SOwner, typename BType,
	std::enable_if_t<(std::is_same<VType, double>::value && std::is_same<BType, double>::value) ||
	(std::is_same<VType, DBL3>::value && (std::is_same<BType, double>::value || std::is_same<BType, DBL33>::value))>*>
void VEC_VC<VType>::set_cmbnd_discontinuous(
	VEC_VC<VType> &V_sec, CMBNDInfo& contact,
	std::function<VType(const Owner&, DBL3, DBL3, DBL3)> a_func_sec, std::function<VType(const Owner&, int, int, DBL3)> a_func_pri,
	std::function<double(const Owner&, DBL3, DBL3, DBL3)> b_func_sec, std::function<double(const Owner&, int, int)> b_func_pri,
	std::function<VType(const Owner&, DBL3, DBL3, DBL3)> diff2_sec, std::function<VType(const Owner&, int, DBL3)> diff2_pri,
	std::function<VType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> A_func_sec, std::function<BType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> B_func_sec,
	std::function<VType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> A_func_pri, std::function<BType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> B_func_pri,
	std::function<double(const Owner&, DBL3, DBL3)> c_func_sec, std::function<double(const Owner&, int)> c_func_pri,
	Owner& instance_sec, Owner& instance_pri, SOwner& instance_s)
{
	INT3 box_sizes = contact.cells_box.size();

	//cellsizes perpendicular to interface
	double hL = contact.hshift_secondary.norm();
	double hR = contact.hshift_primary.norm();
	double hmax = (hL > hR ? hL : hR);

	//primary cells in this contact
#pragma omp parallel for
	for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		int cell1_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

		if (is_empty(cell1_idx) || is_not_cmbnd(cell1_idx)) continue;

		//calculate second primary cell index
		int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * VEC<VType>::n.x + (k + contact.cell_shift.k) * VEC<VType>::n.x*VEC<VType>::n.y;

		//cell values either side of the boundary: V_m2 V_m1 | V_1 V_2; positions as : -2 -1 | 1 2
		//V_m2 and V_2 are known. We need to set values V_m1 and V_1. Here we only set V_1. For this primary, secondary mesh pair there will be another one with order reversed, and there our V_m1 value will be set - so don't worry about it here!
		//NOTE : meshes at composite media boundaries must always be at least 2 non-empty cells thick in directions perpendicular to interface !!! -> this is actually checked when the list of contacts is made

		//relative position of cell -1 in secondary mesh
		DBL3 relpos_m1 = VEC<VType>::rect.s - V_sec.rect.s + ((DBL3(i, j, k) + DBL3(0.5)) & VEC<VType>::h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

		//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
		DBL3 stencil = VEC<VType>::h - mod(contact.hshift_primary) + mod(contact.hshift_secondary);

		//potential values at cells -2 and 2
		VType V_2 = VEC<VType>::quantity[cell2_idx];
		VType V_m2 = V_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

		//obtain a and b values used to define the flux as f(V) = a + b V', both on primary and secondary
		//a values
		VType a_val_sec = a_func_sec(instance_sec, relpos_m1, contact.hshift_secondary, stencil);
		VType a_val_pri = a_func_pri(instance_pri, cell1_idx, cell2_idx, contact.hshift_secondary);

		//b values adjusted with weights
		double b_val_sec = b_func_sec(instance_sec, relpos_m1, contact.hshift_secondary, stencil) * contact.weights.i;
		double b_val_pri = b_func_pri(instance_pri, cell1_idx, cell2_idx) * contact.weights.j;

		//V'' values at cell positions -1 and 1
		VType Vdiff2_sec = diff2_sec(instance_sec, relpos_m1, stencil, contact.hshift_secondary);
		VType Vdiff2_pri = diff2_pri(instance_pri, cell1_idx, contact.hshift_secondary);

		//A values
		VType A_val_sec = A_func_sec(instance_s, cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, instance_sec, instance_pri);
		VType A_val_pri = A_func_pri(instance_s, cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, instance_sec, instance_pri);

		//B values
		BType B_val_sec = B_func_sec(instance_s, cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, instance_sec, instance_pri);
		BType B_val_pri = B_func_pri(instance_s, cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, instance_sec, instance_pri);

		//c values
		double c_m1 = c_func_sec(instance_sec, relpos_m1, stencil);
		double c_m2 = c_func_sec(instance_sec, relpos_m1 + contact.hshift_secondary, stencil);
		double c_1 = c_func_pri(instance_pri, cell1_idx);
		double c_2 = c_func_pri(instance_pri, cell2_idx);
		
		//Form N and N3 values:
		BType N = hmax * c_m2 * B_val_sec + 2.0 * b_val_sec * ident<BType>();
		BType N3 = 3.0 * hmax * c_m1 * B_val_sec + 2.0 * b_val_sec * ident<BType>();

		//Form G value:
		BType G = (B_val_pri * inverse<BType>(N3)) * hmax;

		//Form inverse D value:
		BType D_inv = inverse<BType>(9.0 * hmax * c_m1 * c_1 * G * B_val_sec - 3.0 * hmax * c_1 * B_val_pri - 2.0 * b_val_pri * ident<BType>());

		//Formula for V1
		VEC<VType>::quantity[cell1_idx] = D_inv *
			((hmax * c_m2 * B_val_pri - 3.0 * c_m1 * G * N) * V_m2
			+ (3.0 * hmax * c_m1 * c_2 * G * B_val_sec - hmax * c_2 * B_val_pri - 2.0 * b_val_pri * ident<BType>()) * V_2
			+ 6.0 * c_m1 * b_val_sec * hL * hL * G * Vdiff2_sec + 2.0 * b_val_pri * hR * hR * Vdiff2_pri
			+ (2.0 * (A_val_pri - a_val_pri) - 6.0 * c_m1 * G * (A_val_sec - a_val_sec)) * hmax);
	}
}