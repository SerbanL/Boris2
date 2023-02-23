#pragma once

#include "BorisLib.h"

#include "Boris_Enums_Defs.h"

#ifdef MODULE_COMPILATION_TRANSPORT

class Graph
{
	//////////////////////////////////////////////////////////

	struct Graph_Node
	{
		//mesh or electrode index : start from none set (-1), but one of them will be set to 0 or greater
		int mesh_idx = -1;
		int electrode_idx = -1;

		//settings for applying potential drops for this node
		Rect start_contact, end_contact;
		//potential drop across node
		double start_potential = 0.0, end_potential = 0.0;
		//node resistance along potential drop
		double Resistance = 0.0;

		//some nodes may be degenerate, i.e. used in multiple paths.
		//in this case count the degeneracy index (first), and total degeneracy
		DBL2 degeneracy;

		//------------------------------------

		Graph_Node(int mesh_idx_, int electrode_idx_ = -1) :
			mesh_idx(mesh_idx_), electrode_idx(electrode_idx_)
		{}

		//copy constructor
		Graph_Node(const Graph_Node &copyThis) { *this = copyThis; }

		//------------------------------------

		//comparison
		bool operator==(const Graph_Node &rhs) const { return mesh_idx == rhs.mesh_idx && electrode_idx == rhs.electrode_idx; }

		//assignment operator
		Graph_Node& operator=(const Graph_Node &rhs)
		{
			mesh_idx = rhs.mesh_idx;
			electrode_idx = rhs.electrode_idx;
			return *this;
		}

		//------------------------------------

		bool is_mesh(void) const { return mesh_idx != -1; }

		void set_contacts(Rect start_contact_, Rect end_contact_)
		{
			start_contact = start_contact_;
			end_contact = end_contact_;
		}
	};

	//////////////////////////////////////////////////////////

	struct Graph_Path
	{
		mutable std::vector<Graph_Node> path;

		//total path resistance (sum of node Resistance values)
		double Path_Resistance = 0.0;

		//------------------------------------

		//constructors
		Graph_Path(void) {}
		Graph_Path(const Graph_Node& node) { path = { node }; }

		//------------------------------------

		//indexing
		Graph_Node& operator[](int idx) const { return path[idx]; }

		Graph_Node front(void) { return path.front(); }
		Graph_Node back(void) { return path.back(); }

		//vector method
		size_t size(void) const { return path.size(); }
		void push_back(const Graph_Node& node) { path.push_back(node); };

		//iterators
		Graph_Node* begin(void) { return &path[0]; }
		Graph_Node* data(void) { return path.data(); }
		Graph_Node* end(void) { return &path[size()]; }

		//------------------------------------

		//check if path contains node
		bool path_contains(const Graph_Node& node) const { return vector_contains(path, node); }

		//check if path intersects with given node
		bool path_intersects(const Graph_Node& node, const std::vector<Rect>& mRects, int end_skip = 1) const
		{
			for (int idx = 0; idx < path.size() - end_skip; idx++) {

				if (mRects[path[idx].mesh_idx].get_intersection(mRects[node.mesh_idx]).IsPlane()) return true;
			}
			return false;
		}

		//------------------------------------

		//Comparison used for sorting a vector of paths
		bool operator<(const Graph_Path &rhs) const { return (!path.back().is_mesh() && rhs.path.back().is_mesh()) || size() < rhs.size(); }
	};

	//////////////////////////////////////////////////////////

	//the constructed graph with paths from ground to electrodes
	std::vector<Graph_Path> paths;

	//mesh rectangles
	std::vector<Rect> mRects;

	//electrode rectangles
	std::vector<Rect> eRects;
	//ground electrode rectangle
	Rect gndRect;

public:

	////////////////////////// CTOR

	//input : mesh rectangles, electrode rectangles (other than ground), ground rectangle
	//result : all graph paths from ground to any other electrode, stored in paths
	Graph(const std::vector<Rect>& mRects_, const std::vector<Rect>& eRects_, const Rect& gndRect_) :
		mRects(mRects_), eRects(eRects_), gndRect(gndRect_)
	{
		//build a tree of current paths starting at the ground electrode and ending at any other electrode
		//all meshes in a path (called nodes) must contact at least 1 other mesh, but loops are forbidden; each node in the path touches the previous and next node
		//(apart from ends : starting node touches ground, ending node touches another electrode, or none, in which case all nodes become dangling nodes)
		//start at ground electrode and identify all nodes - meshes in contact with ground electrode
		//take each node in turn and add it to the current path
		//from current node, identify the next set of nodes and keep going recursively until terminating at an electrode (other than ground)
		//if at any point a loop is detected then stop building current path
		//if at any point no more nodes are available, but path is not terminated yet, then trace back (but keep nodes as dangling nodes) until an alternative node is found to continue path
		//if path reaches a dead end (i.e. tracing back through dangling nodes all the way to ground electrode) then terminate (all nodes will be marked as dangling)

		//return all nodes which intersect the given rectangle (and do not coincide with it)
		auto find_nodes = [&](const Rect& rect) -> Graph_Path
		{
			Graph_Path nodes;
			for (int idx = 0; idx < mRects.size(); idx++) {

				if (mRects[idx] != rect && mRects[idx].get_intersection(rect).IsPlane()) nodes.push_back(Graph_Node(idx));
			}
			return nodes;
		};

		//check all electrodes for intersection with given rect : return index if found, else -1
		auto check_electrodes = [&](const Rect& rect) -> int
		{
			for (int idx = 0; idx < eRects.size(); idx++) {

				if (eRects[idx].get_intersection(rect).IsPlane()) return idx;
			}
			return -1;
		};

		//find all starting nodes
		Graph_Path nodes = find_nodes(gndRect);
		for (auto node : nodes) paths.push_back(Graph_Path(node));

		std::function<void(void)> build_paths = [&](void) -> void
		{
			bool new_nodes_found = false;
			size_t num_paths = paths.size();
			for (int pidx = 0; pidx < num_paths; pidx++) {

				//1. if last node touches an electrode, then terminate path immediately
				int eidx = check_electrodes(mRects[paths[pidx].back().mesh_idx]);
				if (eidx >= 0) {
					//found electrode, so terminate path
					paths[pidx].push_back(Graph_Node(-1, eidx));
					continue;
				}

				//2. find all nodes which intersect with last node in path.
				//Note : if no new nodes found, this is a dangling path, so just let it be - code below will not do anything. It will be eliminated when simplifying the graph if redundant.
				Graph_Path nodes = find_nodes(mRects[paths[pidx].back().mesh_idx]);

				//3. for each node found, then add it to current or new path if more than 1 new node found, if following conditions satisfied:
				//a) not already in current path
				//b) doesn't intersect with ground electrode
				//c) doesn't intersect with another node already in this path (other than the last one).

				Graph_Path path = paths[pidx];

				bool first_node_added = false;
				for (auto node : nodes) {

					//a) not already in current path
					if (path.path_contains(node)) continue;

					//b) doesn't intersect with ground electrode
					if (mRects[node.mesh_idx].get_intersection(gndRect).IsPlane()) continue;

					//c) doesn't intersect with another node already in this path (other than the last one)
					if (path.path_intersects(node, mRects)) continue;

					//passed all conditions for viable new nodes
					new_nodes_found = true;
					//on first node add to current path
					if (!first_node_added) { paths[pidx].push_back(node); first_node_added = true; }
					else {

						//if there are more than 1 nodes found, then must create a new path for each
						Graph_Path newpath = path;
						newpath.push_back(node);
						paths.push_back(newpath);
					}
				}
			}

			//keep going until no new nodes found
			if (new_nodes_found) build_paths();
		};

		build_paths();

		//now sort paths by number of nodes, such that paths terminating in an electrode always go before ones without, and starting from the end check if any path has all its nodes in other paths : delete all such redundant paths
		//thus from the end we first check any paths without electrode termination, as we want to get rid of them first, then get rid of longest redundant paths.
		std::sort(paths.begin(), paths.end());

		for (int idx = paths.size() - 1; idx > 0; idx--) {

			//check if this path is redundant (all its nodes contained in other paths)
			bool redundant = true;
			for (auto node : paths[idx]) {

				//check all its nodes to see if contained in other paths
				bool node_contained = false;
				for (int sidx = idx - 1; sidx >= 0; sidx--) {

					if (paths[sidx].path_contains(node)) {

						//yes, this one is. now check the remaining ones
						node_contained = true;
						break;
					}
				}
				if (node_contained) continue;
				//found a node not contained in other paths, so this is not a redundant path and can stop checking it
				redundant = false;
				break;
			}

			if (redundant) paths.pop_back();
		}

		//finally count node degeneracy
		for (int idx = 0; idx < mRects.size(); idx++) {

			//total degeneracy of this mesh
			int degeneracy = 0;
			//store path and node indexes for this mesh, if degenerate
			std::vector<std::pair<int, int>> degeneracy_info;

			for (int pidx = 0; pidx < num_paths(); pidx++) {
				for (int nidx = 0; nidx < path_size(pidx); nidx++) {
					
					if (paths[pidx][nidx].mesh_idx == idx) {

						degeneracy++;
						degeneracy_info.push_back(std::pair<int, int>(pidx, nidx));
						//found node in this path, so stop searching here : a path can only contain a node once
						break;
					}
				}
			}

			//set degeneracy if found
			if (degeneracy) {

				for (int idx = 0; idx < degeneracy_info.size(); idx++) {

					paths[degeneracy_info[idx].first][degeneracy_info[idx].second].degeneracy = DBL2(idx, degeneracy);
				}
			}
		}
	}

	////////////////////////// POTENTIALS IN GRAPH

	//compute potential drops in constructed graph paths
	//input : mesh average conductivities (mCond), electrode potentials (eV), ground potential (gndV)
	//result : all nodes in paths will had start and end coordinates, start and end potentials, ready to be set
	void calculate_potentials(const std::vector<double>& mCond, const std::vector<double>& eV, double gndV)
	{
		//1. find contact start and end coordinates for each node in each path
		for (int pidx = 0; pidx < num_paths(); pidx++) {
			for (int nidx = 0; nidx < path_size(pidx); nidx++) {

				Rect start_contact, end_contact;

				Graph_Node& node = paths[pidx][nidx];

				//cannot setup potential drops if not a mesh node
				if (!node.is_mesh()) continue;

				//find start contact rectangle
				if (nidx == 0) {

					//first node is in contact with ground electrode
					start_contact = mRects[node.mesh_idx].get_intersection(gndRect);
				}
				else {
					//all other nodes are in contact with the previous node, which must be a mesh
					Graph_Node& pnode = paths[pidx][nidx - 1];
					start_contact = mRects[node.mesh_idx].get_intersection(mRects[pnode.mesh_idx]);
				}

				//find end contact rectangle
				if (nidx == path_size(pidx) - 1) {

					//if last node is a mesh, then this is a dangling path - we will be setting all node potentials to ground potential
					end_contact = start_contact;
				}
				else {
					//all other nodes are in contact with the next node, noting the last node could be an electrode
					Graph_Node& nnode = paths[pidx][nidx + 1];

					if (nnode.is_mesh()) {

						end_contact = mRects[node.mesh_idx].get_intersection(mRects[nnode.mesh_idx]);
					}
					else {

						end_contact = mRects[node.mesh_idx].get_intersection(eRects[nnode.electrode_idx]);
					}
				}

				//now set contact start and end points (contact centres)
				node.set_contacts(start_contact, end_contact);
			}
		}

		//2. find resistance for each node in each path
		for (int pidx = 0; pidx < num_paths(); pidx++) {

			//don't bother setting resistances for a dangling path (identified by termination not being an electrode)
			if (paths[pidx].back().is_mesh()) continue;

			//total path resistance
			double Path_Resistance = 0.0;

			for (int nidx = 0; nidx < path_size(pidx); nidx++) {

				Graph_Node& node = paths[pidx][nidx];

				//cannot setup resistance if not a mesh node
				if (!node.is_mesh()) continue;

				//which of x (0), y (1), z (2) component do we want to orient potential drop along? 
				//Choose component perpendicular to largest area contact
				DBL3 normal = (
					node.start_contact.max_area() > node.end_contact.max_area() ?
					node.start_contact.get_normal() : node.end_contact.get_normal());
				//component selection x (0), y (1), z (2)
				int component = (normal.x == 1.0 ? 0 : (normal.y == 1.0 ? 1 : 2));

				double length = 0.0, area = 0.0;
				Rect meshRect = mRects[node.mesh_idx];
				double elecCond = mCond[node.mesh_idx];

				switch (component) {

					//x direction
				case 0:
					length = meshRect.length();
					area = meshRect.width() * meshRect.height();
					break;

					//y direction
				case 1:
					length = meshRect.width();
					area = meshRect.length() * meshRect.height();
					break;

					//z direction
				case 2:
					length = meshRect.height();
					area = meshRect.length() * meshRect.width();
					break;
				}

				node.Resistance = length / (elecCond * area);
				Path_Resistance += node.Resistance;
			}

			paths[pidx].Path_Resistance = Path_Resistance;
		}

		//3. set potential start and end values for each node in each path
		for (int pidx = 0; pidx < num_paths(); pidx++) {

			//for a dangling path (identified by termination not being an electrode) set potentials to ground potential
			if (paths[pidx].back().is_mesh() || IsZ(paths[pidx].Path_Resistance)) {

				for (int nidx = 0; nidx < path_size(pidx); nidx++) {

					paths[pidx][nidx].start_potential = gndV;
					paths[pidx][nidx].end_potential = gndV;
				}
			}
			else {

				//final electrode potential for this path - potential drops from gndV to electrodeV for this path
				double electrodeV = eV[paths[pidx].back().electrode_idx];

				//potential value at start of node
				double V = gndV;

				for (int nidx = 0; nidx < path_size(pidx); nidx++) {

					Graph_Node& node = paths[pidx][nidx];

					//cannot setup potentials if not a mesh node
					if (!node.is_mesh()) continue;

					//potential drop across this node
					double Vdrop = (node.Resistance / paths[pidx].Path_Resistance) * (gndV - electrodeV);

					//start and end potential values for this node
					node.start_potential = V;
					node.end_potential = V - Vdrop;

					//next
					V -= Vdrop;
				}
			}
		}
	}

	////////////////////////// GRAPH PARSING

	//typical structure : 
	//1) loop over paths
	//2) inner loop over nodes in path
	//3) get mesh index for each node
	//4) call set_potential_drop by passing in the required method for the given mesh index for that node

	size_t num_paths(void) const { return paths.size(); }
	size_t path_size(int path_index) const { return paths[path_index].size(); }
	int get_node_mesh_idx(int path_index, int node_index) const { return paths[path_index][node_index].mesh_idx; }

	//set configured node potential drop using method for that node passed in here
	template <typename Owner>
	void set_potential_drop(
		int path_index, int node_index,
		std::function<void(Owner&, Rect, double, Rect, double, DBL2)> Set_Linear_PotentialDrop, Owner& Transport)
	{
		Graph_Node& node = paths[path_index][node_index];
		if (node.is_mesh())
			Transport.Set_Linear_PotentialDrop(node.start_contact, node.start_potential, node.end_contact, node.end_potential, node.degeneracy);
	}
};

#endif