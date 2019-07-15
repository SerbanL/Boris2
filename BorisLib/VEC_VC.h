#pragma once

#include "VEC.h"
#include "ProgramState.h"

////////////////////////////////////////////////////////////////////////////////////////////////// VEC_VC<VType>
//
// extends VEC with vector calculus operations

/////////////////////////////////////////////////////////////////////

template <typename VType> class CGSolve;

struct CMBNDInfo;

template <typename VType>
class VEC_VC : 
	public VEC<VType>,
	public ProgramState<VEC_VC<VType>, 
	std::tuple<SZ3, DBL3, Rect, std::vector<VType>, std::vector<int>, int,
	std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>,
	DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL3, double, int, int>,
	std::tuple<>>
{

	friend CGSolve<VType>;

//the following are used as masks for ngbrFlags. 32 bits in total (4 bytes for an int)

//magnetic neighbor existence masks (+x, -x, +y, -y, +z, -z). Bits 0, 1, 2, 3, 4, 5
#define NF_NPX	1
#define NF_NNX	2
#define NF_NPY	4
#define NF_NNY	8
#define NF_NPZ	16
#define NF_NNZ	32

//Existance of at least one neighbor along a given axis : use as masks (test using &).
#define NF_NGBRX	3	//test bits 0 and 1
#define NF_NGBRY	12	//test bits 2 and 3
#define NF_NGBRZ	48  //test bits 4 and 5

//existence of both neighbors along axes x, y, z. Bits 6, 7, 8
#define NF_BOTHX	64
#define NF_BOTHY	128
#define NF_BOTHZ	256

//mask to check for cell with zero value set : bit 9
#define NF_NOTEMPTY	512

//this is not necessarily an empty cell, but mark them to be skipped during computations for some algorithms (e.g. moving mesh algorithm where the ends of the magnetic mesh must not be updated by the ODE solver) : bit 10
#define NF_SKIPCELL	1024

//these are used in conjunction with dirichlet vectors to indicate dirichlet boundary conditions should be used
//cell on +x side of boundary : bit 11
#define NF_DIRICHLETPX	2048
//cell on -x side of boundary : bit 12
#define NF_DIRICHLETNX	4096
//cell on +y side of boundary : bit 13
#define NF_DIRICHLETPY	8192
//cell on -y side of boundary : bit 14
#define NF_DIRICHLETNY	16384
//cell on +z side of boundary : bit 15
#define NF_DIRICHLETPZ	32768
//cell on -z side of boundary : bit 16
#define NF_DIRICHLETNZ	65536

//masks for x, y, z directions for Dirichlet cells
#define NF_DIRICHLETX (NF_DIRICHLETPX + NF_DIRICHLETNX)
#define NF_DIRICHLETY (NF_DIRICHLETPY + NF_DIRICHLETNY)
#define NF_DIRICHLETZ (NF_DIRICHLETPZ + NF_DIRICHLETNZ)
#define NF_DIRICHLET (NF_DIRICHLETX + NF_DIRICHLETY + NF_DIRICHLETZ)

//composite media boundary cells (used to flag cells where boundary conditions must be applied). These flags are not set using set_ngbrFlags, but must be externally set
//cell on positive x side of boundary : bit 17 -> use dirichlet_nx values
#define NF_CMBNDPX	131072
//cell on negative x side of boundary : bit 18 -> use dirichlet_px values, etc.
#define NF_CMBNDNX	262144
//cell on positive y side of boundary : bit 19
#define NF_CMBNDPY	524288
//cell on negative y side of boundary : bit 20
#define NF_CMBNDNY	1048576
//cell on positive z side of boundary : bit 21
#define NF_CMBNDPZ	2097152
//cell on negative z side of boundary : bit 22
#define NF_CMBNDNZ	4194304

//mask for all cmbnd flags
#define NF_CMBND	(NF_CMBNDPX + NF_CMBNDNX + NF_CMBNDPY + NF_CMBNDNY + NF_CMBNDPZ + NF_CMBNDNZ)
//masks for cmbnd flags along the x, y, or z axes
#define NF_CMBNDX	(NF_CMBNDPX + NF_CMBNDNX)
#define NF_CMBNDY	(NF_CMBNDPY + NF_CMBNDNY)
#define NF_CMBNDZ	(NF_CMBNDPZ + NF_CMBNDNZ)

//Robin boundary conditions flags
//cell on positive x side of boundary : bit 23 -> use robin_nx values
#define NF_ROBINPX	8388608
//cell on negative x side of boundary : bit 24 -> use robin_px values, etc.
#define NF_ROBINNX	16777216
//cell on positive y side of boundary : bit 25
#define NF_ROBINPY	33554432
//cell on negative y side of boundary : bit 26
#define NF_ROBINNY	67108864
//cell on positive z side of boundary : bit 27
#define NF_ROBINPZ	134217728
//cell on negative z side of boundary : bit 28
#define NF_ROBINNZ	268435456
//flag Robin boundary with a void cell (use robin_v values) : bit 29
#define NF_ROBINV	536870912

//mask for all Robin flags
#define NF_ROBIN	(NF_ROBINPX + NF_ROBINNX + NF_ROBINPY + NF_ROBINNY + NF_ROBINPZ + NF_ROBINNZ + NF_ROBINV)
//masks for Robin flags along the x, y, or z axes
#define NF_ROBINX	(NF_ROBINPX + NF_ROBINNX)
#define NF_ROBINY	(NF_ROBINPY + NF_ROBINNY)
#define NF_ROBINZ	(NF_ROBINPZ + NF_ROBINNZ)

//periodic boundary condition along x. Set at x sides only if there is a neighbor present at the other side - this is how we know which side to use, +x or -x : bit 30
#define NF_PBCX	1073741824

//periodic boundary condition along y. Set at y sides only if there is a neighbor present at the other side - this is how we know which side to use, +y or -y : bit 31 (last bit)
#define NF_PBCY	2147483648

//mask for pbc flags, either x or y
#define NF_PBC (NF_PBCX + NF_PBCY)

private:
	
	//used for reductions by magnitude : i.e. for a value of VType reduce (e.g. maximum) for GetMagnitude(value)
	OmpReduction<decltype(GetMagnitude(std::declval<VType>()))> magnitude_reduction, magnitude_reduction2;

	//mark cells with various flags to indicate properties of neighboring cells
	std::vector<int> ngbrFlags;

	int nonempty_cells = 0;

	//store dirichlet boundary conditions at mesh sides - only allocate memory as required
	//these vectors are of sizes equal to 1 cell deep at each respective side. dirichlet_nx are the dirichlet values at the -x side of the mesh, etc.
	std::vector<VType> dirichlet_nx, dirichlet_px, dirichlet_ny, dirichlet_py, dirichlet_nz, dirichlet_pz;

	//Robin boundary conditions values : diff_norm(u) = alpha * (u - h), where diff_norm means differential along surface normal (e.g. positive sign at +x boundary, negative sign at -x boundary).
	//alpha is a positive constant : robins_nx.i. Note, if this is zero then homogeneous Neumann boundary condition results.
	//h is a value (e.g. ambient temperature for heat equation): robins_nx.j
	//nx, px, etc... for mesh boundaries - use these when flagged with NF_ROBINNX etc. and not flagged with NF_ROBINV
	DBL2 robin_px, robin_nx;
	DBL2 robin_py, robin_ny;
	DBL2 robin_pz, robin_nz;
	//robin_v applies for boundary conditions at void cells - more precisely for cells next to a void cell. Use this when flagged with NF_ROBINNX etc. and also flagged with NF_ROBINV
	DBL2 robin_v;

	//when used with moving mesh algorithms calls to shift... functions may be used. If the shift requested is smaller than the cellsize then we cannot perform the shift. 
	//Add it to shift_debt and on next shift call we might be able to shift the mesh values.
	DBL3 shift_debt;

	//adaptive SOR algorithm damping value
	double aSOR_damping = 1.0;
	//save last aSOR error
	double aSOR_lasterror = 0.0;
	//save last aSOR ln(error) gradient
	double aSOR_lastgrad = 0.0;

	//Periodic boundary conditions for evaluating differential operators. If these are set then neighbor flags are calculated accordingly, and applied when evaluating operators.
	//Only implemented x and/or y pbc, not along z.
	int pbc_x = 0;
	int pbc_y = 0;

private:

	//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS : VEC_VC_flags.h

	//set size of ngbrFlags to new_n also mapping shape from current size to new size (if current zero size set solid shape). Memory must be reserved in ngbrFlags to guarantee success. Also n should still have the old value : call this before changing it.
	void resize_ngbrFlags(SZ3 new_n);

	//initialization method for neighbor flags : set flags at size n, counting neighbors etc.
	//Set empty cell values using information in linked_vec (keep same shape) - this must have same rectangle
	template <typename LVType>
	void set_ngbrFlags(VEC_VC<LVType> &linked_vec);

	//initialization method for neighbor flags : set flags at size n, counting neighbors etc. Use current shape in ngbrFlags
	void set_ngbrFlags(void);

	//from NF_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
	VType get_dirichlet_value(int dirichlet_flag, int cell_idx) const;

	//set robin flags from robin values and shape. Doesn't affect any other flags. Call from set_ngbrFlags after counting neighbors, and after setting robin values
	void set_robin_flags(void);

	//set pbc flags depending on set conditions and currently calculated flags - ngbrFlags must already be calculated before using this
	void set_pbc_flags(void);

	//mark cell as not empty / empty : internal use only; routines that use these must finish with recalculating ngbrflags as neighbours will have changed
	void mark_not_empty(int index) { ngbrFlags[index] |= NF_NOTEMPTY; }
	void mark_empty(int index) { ngbrFlags[index] &= ~NF_NOTEMPTY; quantity[index] = VType(); }

	//--------------------------------------------aSOR helper method : VEC_VC_solve.h

	//adjust aSOR damping based on current gradient of ln(error), if error > err_limit.
	void adjust_aSOR_damping(double grad_lnerror, double error, double err_limit);

public:

	//--------------------------------------------CONSTRUCTORS : VEC_VC_mng.h

	VEC_VC(void);

	VEC_VC(const SZ3& n_);

	VEC_VC(const DBL3& h_, const Rect& rect_);

	VEC_VC(const DBL3& h_, const Rect& rect_, VType value);

	~VEC_VC() {}

	//implement ProgramState method
	void RepairObjectState() 
	{ 
		//any mesh transfer info will have to be remade
		transfer.clear(); 
	}

	//--------------------------------------------SPECIAL DATA ACCESS (typically used for copy to/from cuVECs)

	std::vector<int>& ngbrFlags_ref(void) { return ngbrFlags; }

	std::vector<VType>& dirichlet_px_ref(void) { return dirichlet_px; }
	std::vector<VType>& dirichlet_nx_ref(void) { return dirichlet_nx; }
	std::vector<VType>& dirichlet_py_ref(void) { return dirichlet_py; }
	std::vector<VType>& dirichlet_ny_ref(void) { return dirichlet_ny; }
	std::vector<VType>& dirichlet_pz_ref(void) { return dirichlet_pz; }
	std::vector<VType>& dirichlet_nz_ref(void) { return dirichlet_nz; }

	int& nonempty_cells_ref(void) { return nonempty_cells; }

	DBL2& robin_px_ref(void) { return robin_px; }
	DBL2& robin_nx_ref(void) { return robin_nx; }
	DBL2& robin_py_ref(void) { return robin_py; }
	DBL2& robin_ny_ref(void) { return robin_ny; }
	DBL2& robin_pz_ref(void) { return robin_pz; }
	DBL2& robin_nz_ref(void) { return robin_nz; }
	DBL2& robin_v_ref(void) { return robin_v; }

	DBL3& shift_debt_ref(void) { return shift_debt; }

	double& aSOR_damping_ref(void) { return aSOR_damping; }

	int& pbc_x_ref(void) { return pbc_x; }
	int& pbc_y_ref(void) { return pbc_y; }

	//--------------------------------------------SIZING : VEC_VC_mng.h

	//sizing methods return true or false (failed to resize) - if failed then no changes made

	//resize and set shape using linked vec
	template <typename LVType>
	bool resize(const SZ3& new_n, VEC_VC<LVType> &linked_vec);
	//resize but keep shape
	bool resize(const SZ3& new_n);

	//resize and set shape using linked vec
	template <typename LVType>
	bool resize(const DBL3& new_h, const Rect& new_rect, VEC_VC<LVType> &linked_vec);
	//resize but keep shape
	bool resize(const DBL3& new_h, const Rect& new_rect);

	//set value and shape from linked vec
	template <typename LVType>
	bool assign(const SZ3& new_n, VType value, VEC_VC<LVType> &linked_vec);
	//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
	bool assign(const SZ3& new_n, VType value);

	//set value and shape from linked vec
	template <typename LVType>
	bool assign(const DBL3& new_h, const Rect& new_rect, VType value, VEC_VC<LVType> &linked_vec);
	//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
	bool assign(const DBL3& new_h, const Rect& new_rect, VType value);

	void clear(void);

	void shrink_to_fit(void);

	//--------------------------------------------FLAG CHECKING : VEC_VC_flags.h

	int get_nonempty_cells(void) const { return nonempty_cells; }

	bool is_not_empty(int index) const { return (ngbrFlags[index] & NF_NOTEMPTY); }
	bool is_not_empty(const INT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_NOTEMPTY); }
	bool is_not_empty(const DBL3& rel_pos) const { return (ngbrFlags[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y] & NF_NOTEMPTY); }

	bool is_empty(int index) const { return !(ngbrFlags[index] & NF_NOTEMPTY); }
	bool is_empty(const INT3& ijk) const { return !(ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_NOTEMPTY); }
	bool is_empty(const DBL3& rel_pos) const { return !(ngbrFlags[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y] & NF_NOTEMPTY); }

	//check if all cells intersecting the rectangle (absolute coordinates) are empty
	bool is_empty(const Rect& rectangle) const;
	//check if all cells intersecting the rectangle (absolute coordinates) are not empty
	bool is_not_empty(const Rect& rectangle) const;

	bool is_not_cmbnd(int index) const { return !(ngbrFlags[index] & NF_CMBND); }
	bool is_not_cmbnd(const INT3& ijk) const { return !(ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_CMBND); }

	bool is_cmbnd(int index) const { return (ngbrFlags[index] & NF_CMBND); }
	bool is_cmbnd(const INT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_CMBND); }

	bool is_skipcell(int index) const { return (ngbrFlags[index] & NF_SKIPCELL); }

	//--------------------------------------------SET CELL FLAGS - EXTERNAL USE : VEC_VC_flags.h

	//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
	//return false on memory allocation failure only, otherwise return true even if surface_rect was not valid
	bool set_dirichlet_conditions(const Rect& surface_rect, VType value);

	//clear all dirichlet flags and vectors
	void clear_dirichlet_flags(void);

	//clear all pbc flags
	void clear_pbc_flags(void);

	//clear only pbc flags for x direction
	void clear_pbc_x(void);

	//clear only pbc flags for y direction
	void clear_pbc_y(void);

	//set pbc for both x and y
	void set_pbc(void);

	//set pbc for x direction only
	void set_pbc_x(void);

	//set pbc for y direction only
	void set_pbc_y(void);

	//clear all composite media boundary flags
	void clear_cmbnd_flags(void);

	//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations (if status true, else clear the skip cells flags in this rectangle)
	void set_skipcells(const Rect& rectangle, bool status = true);

	//clear all skip cell flags
	void clear_skipcells(void);

	void set_robin_conditions(DBL2 robin_v_, DBL2 robin_px_, DBL2 robin_nx_, DBL2 robin_py_, DBL2 robin_ny_, DBL2 robin_pz_, DBL2 robin_nz_);

	//clear all Robin boundary conditions and values
	void clear_robin_conditions(void);

	//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : VEC_VC_cmbnd.h

	//set cmbnd flags by identifying contacts with other vecs (listed in pVECs); this primary mesh index in that vector is given here as it needs to be stored in CMBNDInfo
	std::vector<CMBNDInfo> set_cmbnd_flags(int primary_mesh_idx, std::vector<VEC_VC<VType>*> &pVECs);

	//calculate and set values at composite media boundary cells in this mesh for a given contacting mesh (the "secondary" V_sec) and given contact description (previously calculated using set_cmbnd_flags)
	//The values are calculated based on the continuity of a potential and flux normal to the interface. The potential is this VEC_VC, call it V, and the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
	//The second order differential of V perpendicular to the interface, V'', is also used and specified using the methods diff2.
	//Boundary with labelled cells either side, on the left is the secondary mesh, on the right is the primary mesh (i.e. the mesh for which we are setting boundary values using this method now): -2 -1 | 1 2 
	//a_func_sec is for the secondary side and takes a position for cell -1, a position shift to add to position to reach cell -2 and finally a stencil to use when obtaining values at cells -1 and -2 (use weighted_average)
	//b_func_sec similar
	//diff2_sec only takes a position and stencil (since only need cell -1)
	//a_func_pri takes indexes for cells 1 and 2. It also takes a position shift vector perpendicular to the interface and pointing from primary to secondary.
	//b_func_pri takes indexes for cells 1 and 2
	//diff2_pri only takes index for cell 1
	//also need instances for the secondary and primary objects whose classes contain the above methods
	template <typename Owner>
	void set_cmbnd_continuous(VEC_VC<VType> &V_sec, CMBNDInfo& contact,
		std::function<VType(const Owner&, DBL3, DBL3, DBL3)> a_func_sec, std::function<VType(const Owner&, int, int, DBL3)> a_func_pri,
		std::function<double(const Owner&, DBL3, DBL3, DBL3)> b_func_sec, std::function<double(const Owner&, int, int)> b_func_pri,
		std::function<VType(const Owner&, DBL3, DBL3)> diff2_sec, std::function<VType(const Owner&, int)> diff2_pri,
		Owner& instance_sec, Owner& instance_pri);
	
	//calculate cmbnd values based on continuity of flux only. The potential is allowed to drop across the interface as:
	//f_sec(V) = f_pri(V) = A + B * delV, where f_sec and f_pri are the fluxes on the secondary and primary sides of the interface, and delV = V_pri - V_sec, the drop in potential across the interface.
	//Thus in addition to the functions in set_cmbnd_continuous we need two extra functions A, B.
	//Note, this reduces to the fully continuous case for B tending to infinity and A = VType(0)
	template <typename Owner, typename SOwner>
	void set_cmbnd_continuousflux(VEC_VC<VType> &V_sec, CMBNDInfo& contact,
		std::function<VType(const Owner&, DBL3, DBL3, DBL3)> a_func_sec, std::function<VType(const Owner&, int, int, DBL3)> a_func_pri,
		std::function<double(const Owner&, DBL3, DBL3, DBL3)> b_func_sec, std::function<double(const Owner&, int, int)> b_func_pri,
		std::function<VType(const Owner&, DBL3, DBL3)> diff2_sec, std::function<VType(const Owner&, int)> diff2_pri,
		std::function<VType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> A_func, std::function<double(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> B_func,
		Owner& instance_sec, Owner& instance_pri, SOwner& instance_s);

	//most general case of composite media boundary conditions
	//calculate cmbnd values based on set boundary flux values; both the flux and potential is allowed to be discontinuous across the interface.
	//Fluxes at the interface are specified as: f_sec(V) = A_sec + B_sec * delVs, f_pri(V) = A_pri + B_pri * delVs, with directions from secondary to primary
	//B functions may return a double, or a DBL33 (3x3 matrix) if VType is DBL3 (Cartesian vector).
	//delVs = c_pri * V_pri - c_sec * V_sec, c are double values specified by given functions
	template <typename Owner, typename SOwner, typename BType,
		std::enable_if_t<(std::is_same<VType, double>::value && std::is_same<BType, double>::value) ||
		(std::is_same<VType, DBL3>::value && (std::is_same<BType, double>::value || std::is_same<BType, DBL33>::value))>* = nullptr>
	void set_cmbnd_discontinuous(
		VEC_VC<VType> &V_sec, CMBNDInfo& contact,
		std::function<VType(const Owner&, DBL3, DBL3, DBL3)> a_func_sec, std::function<VType(const Owner&, int, int, DBL3)> a_func_pri,
		std::function<double(const Owner&, DBL3, DBL3, DBL3)> b_func_sec, std::function<double(const Owner&, int, int)> b_func_pri,
		std::function<VType(const Owner&, DBL3, DBL3)> diff2_sec, std::function<VType(const Owner&, int)> diff2_pri,
		std::function<VType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> A_func_sec, std::function<BType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> B_func_sec,
		std::function<VType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> A_func_pri, std::function<BType(const SOwner&, int, int, DBL3, DBL3, DBL3, Owner&, Owner&)> B_func_pri,
		std::function<double(const Owner&, DBL3, DBL3)> c_func_sec, std::function<double(const Owner&, int)> c_func_pri,
		Owner& instance_sec, Owner& instance_pri, SOwner& instance_s);


	//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS : VEC_VC_shape.h

	//set value in box (i.e. in cells entirely included in box) - all cells become non-empty cells irrespective of value set
	void setbox(const Box& box, VType value = VType());

	//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle - all cells become non-empty cells irrespective of value set
	void setrect(const Rect& rectangle, VType value = VType());

	//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells
	void delrect(const Rect& rectangle);

	//mask values in cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth depending on grayscale value (all if default 0 value).
	bool apply_bitmap_mask(std::vector<BYTE>& bitmap, double zDepth = 0.0);

	//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE GENERATORS : VEC_VC_genshape.h, VEC_VC_Voronoi.h

	//roughen a mesh side (side = "-x", "x", "-y", "y", "-z", "z") to given depth (same units as h) with prng instantiated with given seed
	bool generate_roughside(std::string side, double depth, unsigned seed);

	//Roughen mesh top and bottom surfaces using a jagged pattern to given depth and peak spacing (same units as h) with prng instantiated with given seed.
	//Rough both top and bottom if sides is empty, else it should be either -z or z.
	bool generate_jagged_surfaces(double depth, double spacing, unsigned seed, std::string sides);

	//Generate 2D Voronoi cells with boundaries between cells set to empty
	bool generate_Voronoi2D(double spacing, unsigned seed);

	//Generate 3D Voronoi cells with boundaries between cells set to empty
	bool generate_Voronoi3D(double spacing, unsigned seed);

	//--------------------------------------------MULTIPLE ENTRIES SETTERS - OTHERS : VEC_VC_oper.h

	//exactly the same as assign value - do not use assign as it is slow (sets flags)
	void setnonempty(VType value = VType());

	//set value in non-empty cells only in given rectangle (relative coordinates)
	void setrectnonempty(const Rect& rectangle, VType value = VType());

	//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
	template <typename PType = decltype(GetMagnitude(std::declval<VType>()))>
	void renormalize(PType new_norm);

	//copy values from copy_this but keep current dimensions - if necessary map values from copy_this to local dimensions; from flags only copy the shape but not the boundary condition values or anything else - these are reset
	void copy_values(const VEC_VC<VType>& copy_this);

	//scale all stored values by the given constant
	void scale_values(double constant);

	//shift all the values in this VEC by the given delta (units same as h). Shift values in given shift_rect (absolute coordinates). Not optimised, not parallel code.
	void shift(const DBL3& delta, const Rect& shift_rect);
	
	//shift along an axis - use this for moving mesh algorithms. Fast parallel code.
	void shift_x(double delta, const Rect& shift_rect);

	//shift along an axis - use this for moving mesh algorithms. Fast parallel code.
	void shift_y(double delta, const Rect& shift_rect);

	//shift all the values in this VEC by the given delta (units same as h). Shift values in given shift_rect (absolute coordinates).
	//Also keep magnitude in each cell (e.g. use for vectorial quantities, such as magnetization, to shift only the direction). Not optimised, not parallel code.
	void shift_keepmag(const DBL3& delta, const Rect& shift_rect);

	//--------------------------------------------OPERATIONS : VEC_VC_oper.h

	//overload VEC method : use NF_NOTEMPTY flags instead here
	VType average_nonempty(const Box& box) const;
	//average over non-empty cells over given rectangle (relative to this VEC's rect)
	VType average_nonempty(const Rect& rectangle = Rect()) const;

	//parallel processing versions - do not call from parallel code!!!
	VType average_nonempty_omp(const Box& box) const;
	VType average_nonempty_omp(const Rect& rectangle = Rect()) const;

	//--------------------------------------------OPERATORS and ALGORITHMS

	//----LAPLACE OPERATOR : VEC_VC_del.h

	//calculate Laplace operator at cell with given index. Use Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	VType delsq_neu(int idx) const;

	//calculate Laplace operator at cell with given index. Use non-homogeneous Neumann boundary conditions with the specified boundary differential.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Returns zero at composite media boundary cells.
	VType delsq_nneu(int idx, VAL3<VType>& bdiff) const;

	//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	VType delsq_diri(int idx) const;

	//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Returns zero at composite media boundary cells.
	VType delsq_diri_nneu(int idx, VAL3<VType>& bdiff) const;

	//calculate Laplace operator at cell with given index. Use Robin boundary conditions (defaulting to Neumann if not set).
	//Returns zero at composite media boundary cells.
	//The K constant is used in Robin boundary condition calculations, where -K*diff_norm(T) = alpha*(Tboundary - Tambient) is the flux normal to the boundary - K is the thermal conductivity in the heat equation
	VType delsq_robin(int idx, double K) const;

	//----GRADIENT OPERATOR : VEC_VC_grad.h

	//gradient operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	VAL3<VType> grad_neu(int idx) const;

	//gradient operator. Use non-homogeneous Neumann boundary conditions.
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	VAL3<VType> grad_nneu(int idx, VAL3<VType>& bdiff) const;

	//gradient operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	VAL3<VType> grad_diri(int idx) const;

	//gradient operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Can be used at composite media boundaries where sided differentials will be used instead.
	VAL3<VType> grad_diri_nneu(int idx, VAL3<VType>& bdiff) const;

	//gradient operator. Use sided differentials at boundaries (including at composite media boundaries)
	VAL3<VType> grad_sided(int idx) const;

	//----DIVERGENCE OPERATOR : VEC_VC_div.h

	//divergence operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	double div_neu(int idx) const;
	
	//divergence operator. Use non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	double div_nneu(int idx, VAL3<VType>& bdiff) const;

	//divergence operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	double div_diri(int idx) const;

	//divergence operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	double div_diri_nneu(int idx, VAL3<VType>& bdiff) const;

	//divergence operator. Use sided differentials (also at composite media boundaries)
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	double div_sided(int idx) const;
	
	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	VType diveps3_neu(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	VType diveps3_diri(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use sided differentials (also at composite media boundaries)
	VType diveps3_sided(int idx) const;

	//----CURL OPERATOR : VEC_VC_curl.h

	//curl operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	VType curl_neu(int idx) const;

	//curl operator. Use non-homogeneous Neumann boundary conditions.
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//can only be applied if VType is a VAL3
	VType curl_nneu(int idx, VAL3<VType>& bdiff) const;

	//curl operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	VType curl_diri(int idx) const;

	//curl operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	VType curl_diri_nneu(int idx, VAL3<VType>& bdiff) const;

	//curl operator. Use sided differentials at boundaries (including at composite media boundaries)
	//can only be applied if VType is a VAL3
	VType curl_sided(int idx) const;

	//----LAPLACE / POISSON EQUATION : VEC_VC_solve.h

	//Take one SOR iteration for Laplace equation on this VEC. Return error (maximum change in quantity from one iteration to the next)
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateLaplace_SOR(double relaxation_param = 1.9);

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = F : use Poisson_RHS
	//F must be a member const method of Owner taking an index value (the index ranges over this VEC) and returning a double value : F(index) evaluated at the index-th cell.
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	template <typename Owner>
	DBL2 IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, Owner& instance, double relaxation_param = 1.9);

	//This solves delsq V = F + M * V : For M use Tensor_RHS (For VType double M returns type double, For VType DBL3 M returns DBL33)
	//For Poisson equation we need a function to specify the RHS of the equation delsq V = F : use Poisson_RHS
	//F must be a member const method of Owner taking an index value (the index ranges over this VEC) and returning a double value : F(index) evaluated at the index-th cell.
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	template <typename Owner, typename MType>
	DBL2 IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<MType(const Owner&, int)> Tensor_RHS, Owner& instance, double relaxation_param = 1.9);

	//Poisson equation solved using adaptive SOR algorithm, using homogeneous Neumann boundary condition
	//start_iters must flag the start of a sequence of relaxation iterations (i.e. no equations parameters change during a relaxation sequence)
	//err_limit is the converegence error
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	template <typename Owner>
	DBL2 IteratePoisson_aSOR(std::function<VType(const Owner&, int)> Poisson_RHS, Owner& instance, bool start_iters, double err_limit);

	//Poisson equation solved using SOR, but using non-homogeneous Neumann boundary condition - this is evaluated using the bdiff call-back method.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	template <typename Owner>
	DBL2 IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<VAL3<VType>(const Owner&, int)> bdiff, Owner& instance, double relaxation_param = 1.9);

	//This solves delsq V = F + M * V : For M use Tensor_RHS (For VType double M returns type double, For VType DBL3 M returns DBL33)
	//Poisson equation solved using SOR, but using non-homogeneous Neumann boundary condition - this is evaluated using the bdiff call-back method.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	template <typename Owner, typename MType>
	DBL2 IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<MType(const Owner&, int)> Tensor_RHS, std::function<VAL3<VType>(const Owner&, int)> bdiff, Owner& instance, double relaxation_param = 1.9);

	//Poisson equation solved using adaptive SOR algorithm, using non-homogeneous Neumann boundary condition - this is evaluated using the bdiff call-back method.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	//start_iters must flag the start of a sequence of relaxation iterations (i.e. no equations parameters change during a relaxation sequence)
	//err_limit is the converegence error
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	template <typename Owner>
	DBL2 IteratePoisson_aSOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<VAL3<VType>(const Owner&, int)> bdiff, Owner& instance, bool start_iters, double err_limit);

	//return current aSOR damping value
	double aSOR_get_damping(void) { return aSOR_damping; }
};