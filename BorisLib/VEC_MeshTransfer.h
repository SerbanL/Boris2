#pragma once

#include "VEC.h"

//mesh transfer type to use for transfer in direction : 
//for fully covered super-mesh cells always use weighted average.
//for partially covered super-mesh cells you can use:
//MESHTRANSFERTYPE_WEIGHTED : partially covered cells are weighted by the covered volume ratio (applicable to transfer in only)
//MESHTRANSFERTYPE_CLIPPED : ignore partially covered super-mesh cells.
//MESHTRANSFERTYPE_ENLARGED : as with MESHTRANSFERTYPE_WEIGHTED except the value is not multiplied by the covered volume ratio, so effectively the mesh appears as enlarged.
//for the transfer out direction always use weighted average
enum MESHTRANSFERTYPE_ { MESHTRANSFERTYPE_WEIGHTED = 0, MESHTRANSFERTYPE_CLIPPED, MESHTRANSFERTYPE_ENLARGED };

//list of input mesh cells with pre-calculated weights
//these cells will contribute to some super-mesh cell, and are all the contributions to that cell
struct InMeshCellsWeights {

private:

	//INT2 for (mesh index, cell index); double for the weight
	std::vector< std::pair<INT2, double> > cells;

public:

	InMeshCellsWeights(void) {}

	//copy constructor
	InMeshCellsWeights(const InMeshCellsWeights& copy_weights) { *this = copy_weights; }

	//assignment operator
	InMeshCellsWeights& operator=(const InMeshCellsWeights& copy_weights)
	{
		cells = copy_weights.cells;
		return *this;
	}

	std::pair<INT2, double>& operator[](int idx) { return cells[idx]; }

	void clear(void) { cells.clear(); }

	size_t size(void) const { return cells.size(); }

	//new contribution entry
	void push_back(const std::pair<INT2, double> &newCell) { cells.push_back(newCell); }

	//adjust all stored weights by multiplying with the given value
	//this is used when you first obtain all the required contributions, then at the end apply a common factor multiplication (e.g. multiply by 1 over the total reciprocal distance)
	void multiply_weights(double value)
	{
		for (int idx = 0; idx < (int)cells.size(); idx++)
			cells[idx].second *= value;
	}

	//multiply weights with a constant on the right and return new weights as a new InMeshCellsWeights object
	InMeshCellsWeights operator*(double multiplier) const 
	{ 
		InMeshCellsWeights new_weights = *this; 
		new_weights.multiply_weights(multiplier);

		return new_weights;
	}

	//use accumulation operator to add new contributions to already calculated contributions.
	//New entries are added, or if an entry for (mesh index, cell index) already exists then add to the total weight for that cell.
	void operator+=(const InMeshCellsWeights& new_contributions)
	{
		//most of the time there are none already existing contributions
		if (!cells.size()) {

			cells = new_contributions.cells;
			return;
		}

		//check if this object already has the given INT2 entry (mesh index, cell index) -> if so return index in cells, else -1
		//not the most efficient way to add new entries but the number of contributions are small so it's fine.
		auto has_entry = [&](INT2 entry) -> int {

			for (int idx = 0; idx < cells.size(); idx++) {

				if (cells[idx].first == entry) return idx;
			}

			return -1;
		};

		for (int contrib_idx = 0; contrib_idx < new_contributions.size(); contrib_idx++) {
			
			int existing_idx = has_entry(new_contributions.cells[contrib_idx].first);

			if (existing_idx >= 0) {

				//existing entry : add to weight
				cells[existing_idx].second += new_contributions.cells[contrib_idx].second;
			}
			else {

				//new entry
				cells.push_back(new_contributions.cells[contrib_idx]);
			}
		}
	}
};

//list of super-mesh cells with pre-calculated weights
//these super-mesh cells will contribute to some out mesh cell, and are all the contributions to that cell
struct SuperMeshCellsWeights {

private:

	std::vector< std::pair<int, double> > cells;

public:

	std::pair<int, double>& operator[](int idx) { return cells[idx]; }

	void clear(void) { cells.clear(); }

	void push_back(const std::pair<int, double> &newCell) { cells.push_back(newCell); }

	void multiply_weights(double value)
	{
		for (int idx = 0; idx < (int)cells.size(); idx++)
			cells[idx].second *= value;
	}

	//find weight for given super-mesh cell index idx in cells
	double find_weight(int idx)
	{
		for (int cidx = 0; cidx < cells.size(); cidx++) {

			if (cells[cidx].first == idx) return cells[cidx].second;
		}

		//not found : zero weight
		return 0.0;
	}

	size_t size(void) { return cells.size(); }
};

//Transfer class held by a VEC
template <typename VType>
class Transfer {

private:

	//This mesh, or supermesh : means the VEC holding this Transfer object
	//External mesh : means the VEC we'll transfer values from into this mesh using transfer_info 

	//the VEC holding this Transfer object
	VEC<VType>* pVEC;

	//mesh_in and mesh_out VECs have exactly the same rectangle and cellsize for each index, but may differ in value stored (e.g. magnetisation and effective field) - they could also be exactly the same VEC
	std::vector<VEC<VType>*> mesh_in, mesh_out;

	//input mesh list of contributing cells and weights - transfer_in_info has size pVEC->linear_size()
	//for each super-mesh cell, we have a list of contributing cells from the in meshes together with pre-calculated weights - InMeshCellsWeights
	std::vector< InMeshCellsWeights > transfer_in_info;

	//transfer_out_info has size mesh_out.size() : number of output meshes; for each transfer_out_info entry, we have a vector the size of that out mesh.
	//For each out mesh cell we have a list of contributing super-mesh cells with pre-calculated weights - SuperMeshCellsWeights
	std::vector< std::vector< SuperMeshCellsWeights > > transfer_out_info;

	//total number of transfers from input meshes (i.e. in the flattened transfer info)
	size_t transfer_in_info_size = 0;

	//total number of transfers to output meshes (i.e. in the flattened transfer info)
	size_t transfer_out_info_size = 0;

private:
	
	//----------------------------------- INITIALIZATION HELPERS

	//build InMeshCellsWeights, used to transfer values from mesh_in cells to a super-mesh cell with rectangle rect_c. This calculates the contribution from cells in mesh_in[mesh_idx]
	//Return total reciprocal distance
	double build_meshcells_weights(InMeshCellsWeights &cellsWeights, Rect rect_c, int mesh_idx);

	double build_supermeshcells_weights(SuperMeshCellsWeights &cellsWeights, Rect rect_mc);

	//MESHTRANSFERTYPE_WEIGHTED
	bool initialize_transfer_in_weighted(const std::vector< VEC<VType>* >& mesh_in_);

	//MESHTRANSFERTYPE_CLIPPED
	bool initialize_transfer_in_clipped(const std::vector< VEC<VType>* >& mesh_in_);

	//MESHTRANSFERTYPE_ENLARGED
	bool initialize_transfer_in_enlarged(const std::vector< VEC<VType>* >& mesh_in_);

	bool initialize_transfer_out(const std::vector< VEC<VType>* >& mesh_out_);

public:

	//----------------------------------- CONSTRUCTOR

	//initialize with owner VEC
	Transfer(VEC<VType>* pVEC_) :
		pVEC(pVEC_)
	{}

	//----------------------------------- RUN-TIME TRANSFER METHODS

	//transfer values from the external meshes (mesh_in) into supermesh
	void transfer_from_external_meshes(bool clear = true);

	//transfer values to the external meshes (mesh_out) from the supermesh
	void transfer_to_external_meshes(bool clear = true);

	//----------------------------------- CONFIGURATION

	//clear all transfer information stored - called every time the VEC is resized
	void clear(void);

	//----------------------------------- INITIALIZE TRANSFER

	bool initialize_transfer(const std::vector< VEC<VType>* >& mesh_in_, const std::vector< VEC<VType>* >& mesh_out_, int correction_type);
	
	//----------------------------------- INFO

	size_t size_transfer_in(void) { return transfer_in_info_size; }
	size_t size_transfer_out(void) { return transfer_out_info_size; }

	//----------------------------------- FLATTENED TRANSFER INFO

	//from transfer_in_info and transfer_out_info build flatted transfer_info and pass it on (note vector perfect forwarding makes this ok - build the vector inside this function and return it, the caller can then use it)
	std::vector<std::pair<INT3, double>> get_flattened_transfer_in_info(void);
	std::vector<std::pair<INT3, double>> get_flattened_transfer_out_info(void);

};

//----------------------------------- INITIALIZATION HELPERS

template <typename VType>
double Transfer<VType>::build_meshcells_weights(InMeshCellsWeights &cellsWeights, Rect rect_c, int mesh_idx)
{
	//use mesh_in : this cannot be empty
	VEC<VType>& mesh = *mesh_in[mesh_idx];

	//mesh cellsize
	DBL3 h = mesh.h;
	//mesh size
	INT3 n = mesh.n;

	//center of rect_c relative to mesh
	DBL3 coord = (rect_c.e + rect_c.s) / 2 - mesh.rect.s;

	//positions of lower-left and uper-right corners of cell relative to containing mesh
	DBL3 pos_ll = rect_c.s - mesh.rect.s;
	DBL3 pos_ur = rect_c.e - mesh.rect.s;

	//indexes of cells containing the lower-left and uper-right coordinates
	INT3 idx_ll = floor(pos_ll / h);
	INT3 idx_ur = ceil(pos_ur / h);

	//maximum distance
	double d_max = GetMagnitude(rect_c.size() / 2 + h / 2);
	//total reciprocal distance
	double d_recip_total = 0;

	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//the mesh cell index
				int cell_idx = ii + jj * n.x + kk * n.x*n.y;

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				double d_recip = d_max - get_distance(coord, DBL3((ii + 0.5)*h.x, (jj + 0.5)*h.y, (kk + 0.5)*h.z));
				d_recip_total += d_recip;

				//store indexes and reciprocal distances. Will need to divided later by total reciprocal distance, thus accumulate it and return (d_recip_total).
				cellsWeights.push_back(std::pair<INT2, double>(INT2(mesh_idx, cell_idx), d_recip));
			}
		}
	}

	return d_recip_total;
}

template <typename VType>
double Transfer<VType>::build_supermeshcells_weights(SuperMeshCellsWeights &cellsWeights, Rect rect_mc)
{
	//supermesh cellsize
	DBL3 h = pVEC->h;
	//supermesh size
	INT3 n = pVEC->n;

	//center of rect_mc relative to supermesh
	DBL3 coord = (rect_mc.e + rect_mc.s) / 2 - pVEC->rect.s;

	//positions of lower-left and uper-right corners of cell relative to supermesh
	DBL3 pos_ll = rect_mc.s - pVEC->rect.s;
	DBL3 pos_ur = rect_mc.e - pVEC->rect.s;

	//indexes of cells containing the lower-left and uper-right coordinates
	INT3 idx_ll = floor(pos_ll / h);
	INT3 idx_ur = ceil(pos_ur / h);

	//maximum distance
	double d_max = GetMagnitude(rect_mc.size() / 2 + h / 2);
	//total reciprocal distance
	double d_recip_total = 0;

	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				double d_recip = d_max - get_distance(coord, DBL3((ii + 0.5)*h.x, (jj + 0.5)*h.y, (kk + 0.5)*h.z));
				d_recip_total += d_recip;

				int cell_idx = ii + jj * n.x + kk * n.x*n.y;
				
				//store indexes and reciprocal distances. Will need to divided later by total reciprocal distance, thus accumulate it and return (d_recip_total).
				cellsWeights.push_back(std::pair<int, double>(cell_idx, d_recip));
			}
		}
	}

	return d_recip_total;
}

//----------------------------------- RUN-TIME TRANSFER METHODS

//transfer values from the external meshes (mesh_in) into supermesh
template <typename VType>
void Transfer<VType>::transfer_from_external_meshes(bool clear)
{
	//go through all super-mesh cells
#pragma omp parallel for
	for (int idx = 0; idx < pVEC->linear_size(); idx++) {

		if (transfer_in_info[idx].size()) {

			//first contribution to cell idx : set or add depending on clear flag
			INT2 full_index = transfer_in_info[idx][0].first;
			double weight = transfer_in_info[idx][0].second;

			//obtain weighted value from external mesh
			VType weighted_value = (*mesh_in[full_index.i])[full_index.j] * weight;

			//stored contribution in supermesh
			if (clear) (*pVEC)[idx] = weighted_value;
			else (*pVEC)[idx] += weighted_value;

			//go through all other contributions to cell idx
			for (int cidx = 1; cidx < transfer_in_info[idx].size(); cidx++) {

				INT2 full_index = transfer_in_info[idx][cidx].first;
				double weight = transfer_in_info[idx][cidx].second;

				//obtain weighted value from external mesh
				VType weighted_value = (*mesh_in[full_index.i])[full_index.j] * weight;

				//stored contribution in supermesh
				(*pVEC)[idx] += weighted_value;
			}
		}
		else if(clear) (*pVEC)[idx] = VType();
	}
}

//transfer values to the external meshes (mesh_out) from the supermesh
template <typename VType>
void Transfer<VType>::transfer_to_external_meshes(bool clear)
{	
	//go through all out meshes
	for (int meshIdx = 0; meshIdx < mesh_out.size(); meshIdx++) {

		//for each out mesh go through all its cells to build mesh_entry
#pragma omp parallel for
		for (int idx = 0; idx < mesh_out[meshIdx]->linear_size(); idx++) {

			if (transfer_out_info[meshIdx][idx].size()) {

				//first contribution to cell idx : set or add depending on clear flag
				int index = transfer_out_info[meshIdx][idx][0].first;
				double weight = transfer_out_info[meshIdx][idx][0].second;

				if (clear) (*mesh_out[meshIdx])[idx] = (*pVEC)[index] * weight;
				else (*mesh_out[meshIdx])[idx] += (*pVEC)[index] * weight;

				//go through all other contributions to cell idx
				for (int cidx = 1; cidx < transfer_out_info[meshIdx][idx].size(); cidx++) {

					int index = transfer_out_info[meshIdx][idx][cidx].first;
					double weight = transfer_out_info[meshIdx][idx][cidx].second;

					(*mesh_out[meshIdx])[idx] += (*pVEC)[index] * weight;
				}
			}
			else if (clear) (*mesh_out[meshIdx])[idx] = VType();
		}
	}
}

//----------------------------------- CONFIGURATION

//clear all transfer information stored - called every time the VEC is resized
template <typename VType>
void Transfer<VType>::clear(void)
{
	mesh_in.clear();
	mesh_out.clear();

	transfer_in_info.clear();
	transfer_in_info.shrink_to_fit();
	transfer_out_info.clear();
	transfer_out_info.shrink_to_fit();

	transfer_in_info_size = 0;
	transfer_out_info_size = 0;
}

//----------------------------------- INITIALIZE TRANSFER

template <typename VType>
bool Transfer<VType>::initialize_transfer(const std::vector< VEC<VType>* >& mesh_in_, const std::vector< VEC<VType>* >& mesh_out_, int correction_type)
{
	//-------------------------------------------------------------- Build transfer_in_info

	switch (correction_type) {

	case MESHTRANSFERTYPE_WEIGHTED:
		if (!initialize_transfer_in_weighted(mesh_in_)) return false;
		break;

	case MESHTRANSFERTYPE_CLIPPED:
		if (!initialize_transfer_in_clipped(mesh_in_)) return false;
		break;

	case MESHTRANSFERTYPE_ENLARGED:
		if (!initialize_transfer_in_enlarged(mesh_in_)) return false;
		break;
	};

	//-------------------------------------------------------------- Build transfer_out_info

	if (!initialize_transfer_out(mesh_out_)) return false;

	return true;
}

//MESHTRANSFERTYPE_WEIGHTED
template <typename VType>
bool Transfer<VType>::initialize_transfer_in_weighted(const std::vector< VEC<VType>* >& mesh_in_)
{
	mesh_in = mesh_in_;

	//-------------------------------------------------------------- Build transfer_in_info

	//set size for transfer_in_info : number of cells in the super-mesh
	if (!malloc_vector(transfer_in_info, pVEC->linear_size())) return false;

	transfer_in_info_size = 0;

	//go through all super-mesh cells
	for (int idx = 0; idx < pVEC->linear_size(); idx++) {

		//super-mesh cell rectangle (absolute)
		Rect rect_c = pVEC->get_cellrect(idx);

		//total reciprocal distance
		double d_recip_total = 0;
		double covered_volume_ratio = 0.0;

		//list of contributions generated from considering intersection of this super-mesh cell with in meshes;
		InMeshCellsWeights mesh_cellsWeights;

		//check all input meshes
		for (int mesh_idx = 0; mesh_idx < mesh_in.size(); mesh_idx++) {

			//mesh rectangle
			Rect rect_mesh = mesh_in[mesh_idx]->rect;

			Rect rect_i = rect_mesh.get_intersection(rect_c);

			//does this mesh intersect with the supermesh cell, such that the intersection has a non-zero volume?
			//Don't check the volume directly!! Cells can be on the nm (or even smaller) scale, which means the volume will be tiny, making floating point comparisons tricky.
			//Best to check if it intersects, and the intersection is neither a plane nor a point.
			if (rect_mesh.intersects(rect_c) && !rect_i.IsPlane() && !rect_i.IsPoint()) {

				covered_volume_ratio += (rect_i.size().x / rect_c.size().x) * (rect_i.size().y / rect_c.size().y) * (rect_i.size().z / rect_c.size().z);

				//yes it does. Build intersecting cells weights, keeping track of total reciprocal distance.
				d_recip_total += build_meshcells_weights(mesh_cellsWeights, rect_c, mesh_idx);
			}
		}

		//cellsWeights now contains info about all mesh cells intersecting with rect_c. Divide by total reciprocal distance to obtain proper weights.
		if (d_recip_total > 0) mesh_cellsWeights.multiply_weights(covered_volume_ratio / d_recip_total);

		//store calculated contributions for this super-mesh cell.
		transfer_in_info[idx] = mesh_cellsWeights;

		//calculate the total number of contributing cell transfers : sum of all in meshes transfers to each super-mesh cell. This is the flattened total number of transfers.
		transfer_in_info_size += transfer_in_info[idx].size();
	}

	return true;
}

//MESHTRANSFERTYPE_CLIPPED
template <typename VType>
bool Transfer<VType>::initialize_transfer_in_clipped(const std::vector< VEC<VType>* >& mesh_in_)
{
	mesh_in = mesh_in_;

	//-------------------------------------------------------------- Build transfer_in_info

	//set size for transfer_in_info : number of cells in the super-mesh
	if (!malloc_vector(transfer_in_info, pVEC->linear_size())) return false;

	transfer_in_info_size = 0;

	//go through all super-mesh cells
	for (int idx = 0; idx < pVEC->linear_size(); idx++) {

		//super-mesh cell rectangle (absolute)
		Rect rect_c = pVEC->get_cellrect(idx);

		//total reciprocal distance
		double d_recip_total = 0;

		//list of contributions generated from considering intersection of this super-mesh cell with in meshes;
		InMeshCellsWeights mesh_cellsWeights;

		//check all input meshes
		for (int mesh_idx = 0; mesh_idx < mesh_in.size(); mesh_idx++) {

			//mesh rectangle
			Rect rect_mesh = mesh_in[mesh_idx]->rect;

			//get intersection rectangle
			Rect rect_i = rect_mesh.get_intersection(rect_c);

			//use only fully covered super-mesh cells
			if (rect_i == rect_c) {

				//yes it does. Build intersecting cells weights, keeping track of total reciprocal distance.
				d_recip_total = build_meshcells_weights(mesh_cellsWeights, rect_c, mesh_idx);

				//break as no other meshes can intersect with rect_c
				break;
			}
		}

		//cellsWeights now contains info about all mesh cells intersecting with rect_c. Divide by total reciprocal distance to obtain proper weights.
		if (d_recip_total > 0) mesh_cellsWeights.multiply_weights(1.0 / d_recip_total);

		//store calculated contributions for this super-mesh cell.
		transfer_in_info[idx] = mesh_cellsWeights;

		//calculate the total number of contributing cell transfers : sum of all in meshes transfers to each super-mesh cell. This is the flattened total number of transfers.
		transfer_in_info_size += transfer_in_info[idx].size();
	}

	return true;
}

//MESHTRANSFERTYPE_ENLARGED
template <typename VType>
bool Transfer<VType>::initialize_transfer_in_enlarged(const std::vector< VEC<VType>* >& mesh_in_)
{
	mesh_in = mesh_in_;

	//-------------------------------------------------------------- Build transfer_in_info

	//set size for transfer_in_info : number of cells in the super-mesh
	if (!malloc_vector(transfer_in_info, pVEC->linear_size())) return false;

	transfer_in_info_size = 0;

	//go through all super-mesh cells
	for (int idx = 0; idx < pVEC->linear_size(); idx++) {

		//super-mesh cell rectangle (absolute)
		Rect rect_c = pVEC->get_cellrect(idx);

		//total reciprocal distance
		double d_recip_total = 0;

		//list of contributions generated from considering intersection of this super-mesh cell with in meshes;
		InMeshCellsWeights mesh_cellsWeights;

		//check all input meshes
		for (int mesh_idx = 0; mesh_idx < mesh_in.size(); mesh_idx++) {

			//mesh rectangle
			Rect rect_mesh = mesh_in[mesh_idx]->rect;

			//does this mesh intersect with the supermesh cell, such that the intersection has a non-zero volume?
			//Don't check the volume directly!! Cells can be on the nm (or even smaller) scale, which means the volume will be tiny, making floating point comparisons tricky.
			//Best to check if it intersects, and the intersection is neither a plane nor a point.
			if (rect_mesh.intersects(rect_c) && !rect_mesh.get_intersection(rect_c).IsPlane() && !rect_mesh.get_intersection(rect_c).IsPoint()) {

				//yes it does. Build intersecting cells weights, keeping track of total reciprocal distance.
				d_recip_total += build_meshcells_weights(mesh_cellsWeights, rect_c, mesh_idx);
			}
		}

		//cellsWeights now contains info about all mesh cells intersecting with rect_c. Divide by total reciprocal distance to obtain proper weights.
		if (d_recip_total > 0) mesh_cellsWeights.multiply_weights(1.0 / d_recip_total);

		//store calculated contributions for this super-mesh cell.
		transfer_in_info[idx] = mesh_cellsWeights;

		//calculate the total number of contributing cell transfers : sum of all in meshes transfers to each super-mesh cell. This is the flattened total number of transfers.
		transfer_in_info_size += transfer_in_info[idx].size();
	}

	return true;
}

template <typename VType>
bool Transfer<VType>::initialize_transfer_out(const std::vector< VEC<VType>* >& mesh_out_)
{
	mesh_out = mesh_out_;

	//-------------------------------------------------------------- Build transfer_out_info

	transfer_out_info.clear();
	transfer_out_info.shrink_to_fit();

	transfer_out_info_size = 0;

	//go through all out meshes
	for (int meshIdx = 0; meshIdx < mesh_out.size(); meshIdx++) {

		//build the entry for out mesh with meshIdx
		std::vector< SuperMeshCellsWeights > mesh_entry;

		//for each out mesh go through all its cells to build mesh_entry
		for (int idx = 0; idx < mesh_out[meshIdx]->linear_size(); idx++) {

			//mesh cell rectangle (absolute)
			Rect rect_mc = mesh_out[meshIdx]->get_cellrect(idx);

			//list of all supermesh cells intersecting with this mesh cell
			SuperMeshCellsWeights supermesh_cellsWeights;

			//total reciprocal distance
			double d_recip_total = build_supermeshcells_weights(supermesh_cellsWeights, rect_mc);

			if (d_recip_total > 0) supermesh_cellsWeights.multiply_weights(1.0 / d_recip_total);

			//store in mesh_entry : new entry for cell idx
			mesh_entry.push_back(supermesh_cellsWeights);

			//for this mesh and mesh cell supermesh_cellsWeights.size() has the number of contributing super-mesh cells transfers : add them to flattened total number of transfers
			transfer_out_info_size += supermesh_cellsWeights.size();
		}

		try {

			transfer_out_info.push_back(mesh_entry);
		}
		catch (...) {

			return false;
		}
	}

	return true;
}

//----------------------------------- FLATTENED TRANSFER INFO

//this is used to pass transfer information to a cuVEC for copying to gpu memory : for gpu computations we use "flattened" transfers so it can be parallelized better
//return type: vector of transfers, where INT3 contains : i - input mesh index, j - input mesh cell index, k - super-mesh cell index. the double entry is the weight for the value contribution
template <typename VType>
std::vector<std::pair<INT3, double>> Transfer<VType>::get_flattened_transfer_in_info(void)
{
	std::vector<std::pair<INT3, double>> flattened_transfer_info;

	if (!malloc_vector(flattened_transfer_info, transfer_in_info_size)) return flattened_transfer_info;

	//go through all super-mesh cells (from transfer_in_info)
	int store_index = 0;

	for (int smcIdx = 0; smcIdx < transfer_in_info.size(); smcIdx++) {

		//go through all contributions to this cell
		for (int cidx = 0; cidx < transfer_in_info[smcIdx].size(); cidx++) {

			//in mesh and contributing cell index
			INT2 full_index = transfer_in_info[smcIdx][cidx].first;

			//weight for in transfer for this super-mesh cell and in mesh cell
			double in_weight = transfer_in_info[smcIdx][cidx].second;

			//store flattened info
			flattened_transfer_info[store_index++] = std::pair<INT3, double>(INT3(full_index.i, full_index.j, smcIdx), in_weight);
		}
	}

	return flattened_transfer_info;
}

//this is used to pass transfer information to a cuVEC for copying to gpu memory : for gpu computations we use "flattened" transfers so it can be parallelized better
//return type: vector of transfers, where INT3 contains : i - output mesh index, j - output mesh cell index, k - super-mesh cell index. the double entry is the weight for the value contribution
template <typename VType>
std::vector<std::pair<INT3, double>> Transfer<VType>::get_flattened_transfer_out_info(void)
{
	std::vector<std::pair<INT3, double>> flattened_transfer_info;

	if (!malloc_vector(flattened_transfer_info, transfer_out_info_size)) return flattened_transfer_info;

	//go through all output meshes cells (from transfer_out_info)
	int store_index = 0;

	//parse output meshes
	for (int meshIdx = 0; meshIdx < transfer_out_info.size(); meshIdx++) {

		//parse all cells in each output mesh
		for (int cellIdx = 0; cellIdx < transfer_out_info[meshIdx].size(); cellIdx++) {

			//go through all super-mesh cells contributions to this mesh cell
			for (int idx = 0; idx < transfer_out_info[meshIdx][cellIdx].size(); idx++) {

				INT3 full_index = INT3(meshIdx, cellIdx, transfer_out_info[meshIdx][cellIdx][idx].first);
				double out_weight = transfer_out_info[meshIdx][cellIdx][idx].second;

				//store flattened info
				flattened_transfer_info[store_index++] = std::pair<INT3, double>(full_index, out_weight);
			}
		}
	}

	return flattened_transfer_info;
}