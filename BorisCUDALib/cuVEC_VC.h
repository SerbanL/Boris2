#pragma once

#include <vector>

#include "cuVEC.h"

////////////////////////////////////////////////////////////////////////////////////////////////// VEC_VC<VType>
//
// extends VEC with vector calculus operations

/////////////////////////////////////////////////////////////////////

//CMBNDInfo describes a composite media boundary contact between 2 meshes of same type, used to calculate values at CMBND cells using boundary conditions
struct CMBNDInfoCUDA;

template <typename VType>
class cuVEC_VC :
	public cuVEC<VType>
{

	//the following are used as masks for ngbrFlags. 32 bits in total (4 bytes for an int)

	//magnetic neighbor existence masks (+x, -x, +y, -y, +z, -z). Bits 0, 1, 2, 3, 4, 5
#define NF_NPX	1
#define NF_NNX	2
#define NF_NPY	4
#define NF_NNY	8
#define NF_NPZ	16
#define NF_NNZ	32

//Existance of at least one neighbor along a given axis : use as masks (test using &).
#define NF_NGBRX	(NF_NPX + NF_NNX)	//test bits 0 and 1
#define NF_NGBRY	(NF_NPY + NF_NNY)	//test bits 2 and 3
#define NF_NGBRZ	(NF_NPZ + NF_NNZ)   //test bits 4 and 5

//existence of both neighbors along axes x, y, z
//the use this check mask with value and see if the result is the same as the mask, e.g. if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) { both neighbors are present }
#define NF_BOTHX	(NF_NPX + NF_NNX)
#define NF_BOTHY	(NF_NPY + NF_NNY)
#define NF_BOTHZ	(NF_NPZ + NF_NNZ)

//periodic boundary condition along x. Set at x sides only if there is a neighbor present at the other side - this is how we know which side to use, +x or -x : bit 6
#define NF_PBCX	64

//periodic boundary condition along y. Set at y sides only if there is a neighbor present at the other side - this is how we know which side to use, +y or -y : bit 7
#define NF_PBCY	128

//periodic boundary condition along z. Set at z sides only if there is a neighbor present at the other side - this is how we know which side to use, +z or -z : bit 8
#define NF_PBCZ	256

//mask for pbc flags, either x or y
#define NF_PBC (NF_PBCX + NF_PBCY + NF_PBCZ)

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

private:

	//mark cells with various flags to indicate properties of neighboring cells
	int* ngbrFlags;
	//allocated size for ngbrFlags (value stored on gpu)
	size_t ngbrFlags_size;

	//number of non-empty cells (i.e. cells  not marked with NF_NOTEMPTY)
	int nonempty_cells;

	//store dirichlet boundary conditions at mesh sides - only allocate memory as required
	//these vectors are of sizes equal to 1 cell deep at each respective side. dirichlet_nx are the dirichlet values at the -x side of the mesh, etc.
	VType* dirichlet_px;
	VType* dirichlet_nx;
	VType* dirichlet_py;
	VType* dirichlet_ny;
	VType* dirichlet_pz;
	VType* dirichlet_nz;

	//allocated sizes for dirichlet values (size values stored on gpu)
	size_t dirichlet_px_size;
	size_t dirichlet_nx_size;
	size_t dirichlet_py_size;
	size_t dirichlet_ny_size;
	size_t dirichlet_pz_size;
	size_t dirichlet_nz_size;

	//Robin boundary conditions values : diff_norm(u) = alpha * (u - h), where diff_norm means differential along surface normal (e.g. positive sign at +x boundary, negative sign at -x boundary).
	//alpha is a positive constant : robins_nx.i. Note, if this is zero then homogeneous Neumann boundary condition results.
	//h is a value (e.g. ambient temperature for heat equation): robins_nx.j
	//nx, px, etc... for mesh boundaries - use these when flagged with NF_ROBINNX etc. and not flagged with NF_ROBINV
	cuReal2 robin_px, robin_nx;
	cuReal2 robin_py, robin_ny;
	cuReal2 robin_pz, robin_nz;
	//robin_v applies for boundary conditions at void cells - more precisely for cells next to a void cell. Use this when flagged with NF_ROBINNX etc. and also flagged with NF_ROBINV
	cuReal2 robin_v;

	//when used with moving mesh algorithms calls to shift... functions may be used. If the shift requested is smaller than the cellsize then we cannot perform the shift. 
	//Add it to shift_debt and on next shift call we might be able to shift the mesh values.
	cuReal3 shift_debt;

public:

	//adaptive SOR algorithm damping value
	cuReal aSOR_damping;
	
private:

	//save last aSOR error
	cuReal aSOR_lasterror;
	//save last aSOR ln(error) gradient
	cuReal aSOR_lastgrad;

	//Periodic boundary conditions for evaluating differential operators. If these are set then neighbor flags are calculated accordingly, and applied when evaluating operators.
	//Only implemented x and/or y pbc, not along z.
	int pbc_x;
	int pbc_y;
	int pbc_z;

private:

	//--------------------------------------------MEMORY MANAGEMENT HELPER METHODS : in cuVEC_VC_mng.h

	//memory allocation for objects and initialize to default - only call at start of managed constructors
	__host__ void alloc_initialize_data(void);

	//set size of ngbrFlags, allocating memory
	__host__ cudaError_t set_ngbrFlags_size(size_t size);

	//get ngbrFlags_size value in cpu memory
	__host__ size_t get_ngbrFlags_size(void) const;

	//set robin values using a robin id
	__host__ void set_robin(cuReal2 robin_value, int robin_id);

	//set dirichlet array size, allocating memory as required
	__host__ cudaError_t set_dirichlet_size(size_t size, int dirichlet_id);

	//get dirichlet array size
	__host__ size_t get_dirichlet_size(int dirichlet_id) const;

	//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS : in cuVEC_VC_flags.cuh and cuVEC_VC_flags.h
	
	//set size of ngbrFlags to new_n also mapping shape from current size to new size (if current zero size set solid shape). Memory must be reserved in ngbrFlags to guarantee success. Also n should still have the old value : call this before changing it.
	__host__ cudaError_t resize_ngbrFlags(cuSZ3 new_n);

	//initialization method for neighbor flags : set flags at size n, counting neighbors etc.
	//Set empty cell values using information in linked_vec (keep same shape) - this must have same rectangle
	__host__ void set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);

	//initialization method for neighbor flags : set flags at size n, counting neighbors etc. Use current shape in ngbrFlags
	__host__ void set_ngbrFlags(void);
	
	//from NF_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
	__device__ VType get_dirichlet_value(int dirichlet_flag, int idx) const;
	__device__ VType get_dirichlet_value(int dirichlet_flag, const cuINT3& ijk) const;
	
	//set robin flags from robin values and shape. Doesn't affect any other flags. Call from set_ngbrFlags after counting neighbors, and after setting robin values
	__host__ void set_robin_flags(void);

	//set pbc flags depending on set conditions and currently calculated flags - ngbrFlags must already be calculated before using this
	__host__ void set_pbc_flags(void);

	//mark cell as not empty / empty : internal use only; routines that use these must finish with recalculating ngbrflags as neighbours will have changed
	__device__ void mark_not_empty(int index) { ngbrFlags[index] |= NF_NOTEMPTY; }
	__device__ void mark_empty(int index) { ngbrFlags[index] &= ~NF_NOTEMPTY; quantity[index] = VType(); }

public:

	//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only. Real constructors are never called since you should never make a real instance of a cuVEC_VC. : cuVEC_VC_mng.h

	//void constructor
	__host__ void construct_cu_obj(void);

	//construct to given number of cells : n_ is in cpu memory
	__host__ void construct_cu_obj(const cuSZ3& n_);

	//construct to given dimensions : h_ and rect_ are in cpu memory
	__host__ void construct_cu_obj(const cuReal3& h_, const cuRect& rect_);

	//construct to given dimensions and initialize to given value : h_ and rect_ are in cpu memory
	__host__ void construct_cu_obj(const cuReal3& h_, const cuRect& rect_, VType value);

	//copy constructor
	__host__ void construct_cu_obj(const cuVEC_VC& copyThis);

	//assignment operator
	__host__ void assign_cu_obj(const cuVEC_VC& copyThis);

	//destructor
	__host__ void destruct_cu_obj(void);

	//--------------------------------------------SPECIAL ACCESS

	__host__ int*& ngbrFlags_ref(void) { return ngbrFlags; }

	//--------------------------------------------COPY TO / FROM VEC_VC

	//copy everything from a VEC_VC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC_VC>
	__host__ bool set_from_cpuvec(cpuVEC_VC& vec_vc);

	//copy everything to a VEC_VC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC_VC>
	__host__ bool set_cpuvec(cpuVEC_VC& vec_vc);

	//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC : only quantity is copied.
	template <typename cpuVEC_VC>
	__host__ bool copy_from_cpuvec(cpuVEC_VC& vec_vc);

	//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC_VC : only quantity and ngbrFlags are copied.
	template <typename cpuVEC_VC>
	__host__ bool copy_to_cpuvec(cpuVEC_VC& vec_vc);

	//copy flags only from vec_vc, where sizes must match
	template <typename cpuVEC_VC>
	__host__ bool copyflags_from_cpuvec(cpuVEC_VC& vec_vc);

	//--------------------------------------------COPY TO ANOTHER cuVEC :  cuVEC_VC_mng.cuh

	//extract values from this and place them in cuvec : both must have same rectangle, but can differ in h - cuvec.h <= this->h needed (and hence n, where cuvec.n.dim() = size); e.g. this method allows extraction of a coarser cuvec.
	__host__ void extract_cuvec(size_t size, cuVEC_VC<VType>& cuvec);

	//--------------------------------------------SIZING : cuVEC_VC_shape.h and cuVEC_VC_shape.cuh

	//sizing methods return true or false (failed to resize) - if failed then no changes made

	//resize and set shape using linked vec
	template <typename LVType>
	__host__ bool resize(cuSZ3 new_n, cuVEC_VC<LVType>& linked_vec);
	//resize but keep shape
	__host__ bool resize(cuSZ3 new_n);
	
	//resize and set shape using linked vec
	template <typename LVType>
	__host__ bool resize(cuReal3 new_h, cuRect new_rect, cuVEC_VC<LVType>& linked_vec);
	//resize but keep shape
	__host__ bool resize(cuReal3 new_h, cuRect new_rect);

	//set value and shape from linked vec
	template <typename LVType>
	__host__ bool assign(cuSZ3 new_n, VType value, cuVEC_VC<LVType>& linked_vec);
	//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
	__host__ bool assign(cuSZ3 new_n, VType value);

	//set value and shape from linked vec
	template <typename LVType>
	__host__ bool assign(cuReal3 new_h, cuRect new_rect, VType value, cuVEC_VC<LVType>& linked_vec);
	//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
	__host__ bool assign(cuReal3 new_h, cuRect new_rect, VType value);

	void clear(void);

	//--------------------------------------------FLAG CHECKING : cuVEC_VC_flags.h

	__device__ int get_nonempty_cells(void) const { return nonempty_cells; }
	__host__ int get_nonempty_cells_cpu(void);
	
	__device__ bool is_not_empty(int index) const { return (ngbrFlags[index] & NF_NOTEMPTY); }
	__device__ bool is_not_empty(const cuINT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_NOTEMPTY); }
	__device__ bool is_not_empty(const cuReal3& rel_pos) const { return (ngbrFlags[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y] & NF_NOTEMPTY); }

	__device__ bool is_empty(int index) const { return !(ngbrFlags[index] & NF_NOTEMPTY); }
	__device__ bool is_empty(const cuINT3& ijk) const { return !(ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_NOTEMPTY); }
	__device__ bool is_empty(const cuReal3& rel_pos) const { return !(ngbrFlags[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y] & NF_NOTEMPTY); }

	__device__ bool is_empty(const cuRect& rectangle) const;

	__device__ bool is_not_cmbnd(int index) const { return !(ngbrFlags[index] & NF_CMBND); }
	__device__ bool is_not_cmbnd(const cuINT3& ijk) const { return !(ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_CMBND); }

	__device__ bool is_cmbnd(int index) const { return (ngbrFlags[index] & NF_CMBND); }
	__device__ bool is_cmbnd(const cuINT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_CMBND); }

	__device__ bool is_skipcell(int index) const { return (ngbrFlags[index] & NF_SKIPCELL); }

	//are all neighbors available? (for 2D don't check the z neighbors)
	__device__ bool is_interior(int index) const { return (((ngbrFlags[index] & NF_BOTHX) == NF_BOTHX) && ((ngbrFlags[index] & NF_BOTHY) == NF_BOTHY) && (n.z == 1 || ((ngbrFlags[index] & NF_BOTHZ) == NF_BOTHZ))); }
	//are all neighbors in the xy plane available?
	__device__ bool is_plane_interior(int index) const { return (((ngbrFlags[index] & NF_BOTHX) == NF_BOTHX) && ((ngbrFlags[index] & NF_BOTHY) == NF_BOTHY)); }

	//--------------------------------------------SET CELL FLAGS - EXTERNAL USE : cuVEC_VC_flags.h and cuVEC_VC_flags.cuh
	
	//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
	//return false on memory allocation failure only, otherwise return true even if surface_rect was not valid
	__host__ bool set_dirichlet_conditions(cuRect surface_rect, VType value);

	//clear all dirichlet flags and vectors
	__host__ void clear_dirichlet_flags(void);

	//set pbc conditions : setting any to false clears flags
	__host__ void set_pbc(bool pbc_x_, bool pbc_y_, bool pbc_z_);

	//clear all pbc flags : can also be achieved setting all flags to false in set_pbc but this one is more readable
	__host__ void clear_pbc(void);

	//clear all composite media boundary flags
	__host__ void clear_cmbnd_flags(void);

	//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations (if status true, else clear the skip cells flags in this rectangle)
	__host__ void set_skipcells(cuRect rectangle, bool status = true);

	//clear all skip cell flags
	__host__ void clear_skipcells(void);

	__host__ void set_robin_conditions(cuReal2 robin_v_, cuReal2 robin_px_, cuReal2 robin_nx_, cuReal2 robin_py_, cuReal2 robin_ny_, cuReal2 robin_pz_, cuReal2 robin_nz_);

	//clear all Robin boundary conditions and values
	__host__ void clear_robin_conditions(void);
	
	//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS : cuVEC_VC_shape.h and cuVEC_VC_shape.cuh

	//set value in box (i.e. in cells entirely included in box) - all cells become non-empty cells irrespective of value set
	__host__ void setbox(cuBox box, VType value = VType());

	//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle - all cells become non-empty cells irrespective of value set
	__host__ void setrect(cuRect rectangle, VType value = VType());

	//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells
	__host__ void delrect(cuRect rectangle);

	//mask values in cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
	__host__ bool apply_bitmap_mask(std::vector<BYTE>& bitmap, int zDepth = 0);
	
	//--------------------------------------------MULTIPLE ENTRIES SETTERS - OTHERS : cuVEC_oper.h and cuVEC_oper.cuh

	//exactly the same as assign value - do not use assign as it is slow (sets flags)
	__host__ void setnonempty(VType value = VType());

	//set value in non-empty cells only in given rectangle (relative coordinates)
	__host__ void setrectnonempty(const cuRect& rectangle, VType value = VType());

	//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ void renormalize(size_t arr_size, PType new_norm);

	//scale all stored values by the given constant
	__host__ void scale_values(size_t size, cuReal constant);
	__host__ void scale_values(cuReal constant) { scale_values(get_gpu_value(n).dim(), constant); }

	//shift all the values in this cuVEC by the given delta (units same as h). Shift values in given shift_rect (absolute coordinates).
	__host__ void shift_x(size_t size, cuReal delta, cuRect shift_rect);

	//--------------------------------------------OPERATIONS : cuVEC_VC_oper.h and cuVEC_VC_oper.cuh

	//overload VEC method : use NF_NOTEMPTY flags instead here
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	__host__ VType average_nonempty(size_t arr_size, cuBox box);
	//average over non-empty cells over given rectangle (relative to this VEC's rect)
	__host__ VType average_nonempty(size_t arr_size, cuRect rectangle = cuRect());
	
	//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : cuVEC_VC_cmbnd.h

	//calculate and set values at composite media boundary cells in this mesh for a given contacting mesh (the "secondary" V_sec) and given contact description (previously calculated using set_cmbnd_flags)
	//The values are calculated based on the continuity of a potential and flux normal to the interface. The potential is this VEC_VC, call it V, and the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
	//The second order differential of V perpendicular to the interface, V'', is also used and specified using the methods diff2 in Class_CMBND
	//Boundary with labelled cells either side, on the left is the secondary mesh, on the right is the primary mesh (i.e. the mesh for which we are setting boundary values using this method now): -2 -1 | 1 2 
	//Functions which must be defined in Class_CMBND:
	//a_func_sec is for the secondary side and takes a position for cell -1, a position shift to add to position to reach cell -2 and finally a stencil to use when obtaining values at cells -1 and -2 (use weighted_average)
	//b_func_sec similar
	//diff2_sec only takes a position and stencil (since only need cell -1)
	//a_func_pri takes indexes for cells 1 and 2. It also takes a position shift vector perpendicular to the interface and pointing from primary to secondary.
	//b_func_pri takes indexes for cells 1 and 2
	//diff2_pri only takes index for cell 1
	//also need instances for the secondary and primary objects whose classes contain the above methods
	template <typename Class_CMBND>
	__host__ void set_cmbnd_continuous(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, CMBNDInfoCUDA& contact);

	//calculate cmbnd values based on continuity of flux only. The potential is allowed to drop across the interface as:
	//f_sec(V) = f_pri(V) = A + B * delV, where f_sec and f_pri are the fluxes on the secondary and primary sides of the interface, and delV = V_pri - V_sec, the drop in potential across the interface.
	//Thus in addition to the functions in set_cmbnd_continuous we need two extra functions A_func, B_func, defined in Class_CMBND_S
	//Note, this reduces to the fully continuous case for B tending to infinity and A = VType(0)
	template <typename Class_CMBND, typename Class_CMBND_S>
	__host__ void set_cmbnd_continuousflux(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, Class_CMBND_S& cmbndFuncs_s, CMBNDInfoCUDA& contact);

	//most general case of composite media boundary conditions
	//calculate cmbnd values based on set boundary flux values; both the flux and potential is allowed to be discontinuous across the interface.
	//Fluxes at the interface are specified as: f_sec(V) = A_sec + B_sec * delVs, f_pri(V) = A_pri + B_pri * delVs, with directions from secondary to primary
	//B functions may return a double, or a DBL33 (3x3 matrix) if VType is DBL3 (Cartesian vector).
	//delVs = c_pri * V_pri - c_sec * V_sec, c are double values specified by given functions
	template <typename Class_CMBND, typename Class_CMBND_S>
	__host__ void set_cmbnd_discontinuous(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, Class_CMBND_S& cmbndFuncs_s, CMBNDInfoCUDA& contact);

	//--------------------------------------------OPERATORS and ALGORITHMS

	//----LAPLACE OPERATOR : cuVEC_VC_del.h

	//calculate Laplace operator at cell with given index. Use Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	__device__ VType delsq_neu(int idx) const;

	//calculate Laplace operator at cell with given index. Use non-homogeneous Neumann boundary conditions with the specified boundary differential.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Returns zero at composite media boundary cells.
	template <typename Class_BDiff>
	__device__ VType delsq_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType delsq_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	__device__ VType delsq_diri(int idx) const;

	//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Returns zero at composite media boundary cells.
	template <typename Class_BDiff>
	__device__ VType delsq_diri_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType delsq_diri_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//calculate Laplace operator at cell with given index. Use Robin boundary conditions (defaulting to Neumann if not set).
	//Returns zero at composite media boundary cells.
	//The K constant is used in Robin boundary condition calculations, where -K*diff_norm(T) = alpha*(Tboundary - Tambient) is the flux normal to the boundary - K is the thermal conductivity in the heat equation
	__device__ VType delsq_robin(int idx, cuReal K) const;

	//----GRADIENT OPERATOR : cuVEC_VC_grad.h

	//gradient operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ cuVAL3<VType> grad_neu(int idx) const;

	//gradient operator. Use non-homogeneous Neumann boundary conditions.
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	template <typename Class_BDiff>
	__device__ cuVAL3<VType> grad_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuVAL3<VType> grad_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//gradient operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ cuVAL3<VType> grad_diri(int idx) const;

	//gradient operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	template <typename Class_BDiff>
	__device__ cuVAL3<VType> grad_diri_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuVAL3<VType> grad_diri_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//gradient operator. Use sided differentials at boundaries (including at composite media boundaries)
	__device__ cuVAL3<VType> grad_sided(int idx) const;

	//----DIVERGENCE OPERATOR : cuVEC_VC_div.h

	//divergence operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	__device__ cuReal div_neu(int idx) const;

	//divergence operator. Use non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	template <typename Class_BDiff>
	__device__ cuReal div_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuReal div_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//divergence operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	__device__ cuReal div_diri(int idx) const;

	//divergence operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	template <typename Class_BDiff>
	__device__ cuReal div_diri_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuReal div_diri_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//divergence operator. Use sided differentials (also at composite media boundaries)
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	__device__ cuReal div_sided(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ VType diveps3_neu(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ VType diveps3_diri(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use sided differentials (also at composite media boundaries)
	__device__ VType diveps3_sided(int idx) const;

	//----CURL OPERATOR : cuVEC_VC_curl.h

	//curl operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	__device__ VType curl_neu(int idx) const;

	//curl operator. Use non-homogeneous Neumann boundary conditions.
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//can only be applied if VType is a VAL3
	template <typename Class_BDiff>
	__device__ VType curl_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType curl_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//curl operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	__device__ VType curl_diri(int idx) const;

	//curl operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	template <typename Class_BDiff>
	__device__ VType curl_diri_nneu(int idx, Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType curl_diri_nneu(int idx, cuVAL3<VType>& bdiff) const;

	//curl operator. Use sided differentials at boundaries (including at composite media boundaries)
	//can only be applied if VType is a VAL3
	__device__ VType curl_sided(int idx) const;

	//----LAPLACE / POISSON EQUATION : cuVEC_VC_solve.cuh
	
	//LAPLACE

	//Take one SOR iteration for Laplace equation on this VEC. Return error (maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a cu_obj managed cuReal.
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	__host__ void IterateLaplace_SOR(size_t arr_size, cuReal& damping, cuReal& max_error, cuReal& max_val);

	//red and black SOR passes launched from kernels : the flow is : 
	//IterateLaplace_SOR -> 1. global red pass kernel with *this as parameter -> call IterateLaplace_SOR_red, 
	//2. global black pass kernel with *this as parameter -> call IterateLaplace_SOR_black
	__device__ void IterateLaplace_SOR_red(cuReal damping);
	__device__ void IterateLaplace_SOR_black(cuReal damping, cuReal& max_error, cuReal& max_val);

	//POISSON with homogeneous Neumann boundaries

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = Poisson_RHS
	//Poisson_RHS must be a member const method of Class_Poisson_RHS taking an index value (the index ranges over this VEC) and returning a cuReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//obj of type Class_Poisson_RHS must be cu_obj managed so it is entirely in gpu memory.
	//Return error(maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a cu_obj managed cuReal.
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	template <typename Class_Poisson_RHS>
	__host__ void IteratePoisson_SOR(size_t arr_size, Class_Poisson_RHS& obj, cuReal& damping, cuReal& max_error, cuReal& max_val);

	//red and black SOR passes launched from kernels : the flow is : 
	//IteratePoisson_SOR -> 1. global red pass kernel with *this as parameter -> call IteratePoisson_SOR_red, 
	//2. global black pass kernel with *this as parameter -> call IteratePoisson_SOR_black
	template <typename Class_Poisson_RHS>
	__device__ void IteratePoisson_SOR_red(Class_Poisson_RHS& obj, cuReal damping);
	template <typename Class_Poisson_RHS>
	__device__ void IteratePoisson_SOR_black(Class_Poisson_RHS& obj, cuReal damping, cuReal& max_error, cuReal& max_val);

	//POISSON with homogeneous Neumann boundaries and adaptive SOR algorithm

	//Poisson equation solved using adaptive SOR algorithm, using homogeneous Neumann boundary condition
	//The algorithm also requires error limits : min error below which the solution is assumed to have converged, max error above which damping must be immediately set to smallest value to avoid solution divergence.
	//These errors must be specified un-normalized, noting the return error is also un-normalized. 
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	template <typename Class_Poisson_RHS>
	__host__ void IteratePoisson_aSOR(size_t arr_size, Class_Poisson_RHS& obj, bool start_iters, cuReal err_limit, cuReal& max_error, cuReal& max_val);

	//adjust aSOR parameters depending on passed error value (and also set the error value to the current aSOR error as required by the algorithm)
	__device__ void adjust_aSOR_damping(bool start_iters, cuReal err_limit, cuReal& error, cuReal& max_val);

	//POISSON with non-homogeneous Neumann boundaries

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = Poisson_RHS
	//Poisson_RHS must be a member const method of Class_Poisson_NNeu taking an index value (the index ranges over this VEC) and returning a cuReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//Class_Poisson_NNeu must also define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index) - this is the non-homogeneous Neumann boundary condition at that cell
	//obj of type Class_Poisson_NNeu must be cu_obj managed so it is entirely in gpu memory.
	//Return error(maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a cu_obj managed cuReal.
	//Dirichlet boundary conditions used, defaulting to non-homogeneous Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	template <typename Class_Poisson_NNeu>
	__host__ void IteratePoisson_NNeu_SOR(size_t arr_size, Class_Poisson_NNeu& obj, cuReal& damping, cuReal& max_error, cuReal& max_val);

	//red and black SOR passes launched from kernels : the flow is : 
	//IteratePoisson_NNeu_SOR -> 1. global red pass kernel with *this as parameter -> call IteratePoisson_NNeu_SOR_red, 
	//2. global black pass kernel with *this as parameter -> call IteratePoisson_NNeu_SOR_black
	template <typename Class_Poisson_NNeu>
	__device__ void IteratePoisson_NNeu_SOR_red(Class_Poisson_NNeu& obj, cuReal damping);
	template <typename Class_Poisson_NNeu>
	__device__ void IteratePoisson_NNeu_SOR_black(Class_Poisson_NNeu& obj, cuReal damping, cuReal& max_error, cuReal& max_val);

	//POISSON with non-homogeneous Neumann boundaries and adaptive SOR algorithm

	//Poisson equation solved using adaptive SOR algorithm, using homogeneous Neumann boundary condition
	//The algorithm also requires error limits : min error below which the solution is assumed to have converged, max error above which damping must be immediately set to smallest value to avoid solution divergence.
	//These errors must be specified un-normalized, noting the return error is also un-normalized. 
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	//Poisson_RHS must be a member const method of Class_Poisson_NNeu taking an index value (the index ranges over this VEC) and returning a cuReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//Class_Poisson_NNeu must also define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index) - this is the non-homogeneous Neumann boundary condition at that cell
	template <typename Class_Poisson_NNeu>
	__host__ void IteratePoisson_NNeu_aSOR(size_t arr_size, Class_Poisson_NNeu& obj, bool start_iters, cuReal err_limit, cuReal& max_error, cuReal& max_val);

	//return current aSOR damping value
	cuReal aSOR_get_damping_cpu(void) { return get_gpu_value(aSOR_damping); }
};
