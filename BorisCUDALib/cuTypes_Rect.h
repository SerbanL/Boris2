#pragma once

#include <cuda_runtime.h>

#include "cuFuncs_Aux.h"
#include "cuFuncs_Math.h"
#include "cuTypes_VAL.h"

#include "alloc_cpy.h"

#include "../BorisLib/Types_Rect.h"

////////////////////////////////////////////////////////////////////////////////////////////////// cuBox
//
// cuBox (or a Plane, or a Line, or a Point!, but probably only useful as an actual cuBox or Plane) - Discrete coordinates.

template <typename Type = void>
struct __cuBox {

	//----------------------------- DATA

	//start and end coordinates of discrete box (integer coordinates).
	cuINT3 s, e;

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(s, cuINT3());
		set_gpu_value(e, cuINT3());
	}

	__host__ void construct_cu_obj(const __cuBox& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const __cuBox& copyThis)
	{
		gpu_to_gpu(s, copyThis.s);
		gpu_to_gpu(e, copyThis.e);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ __cuBox(void) {}

	__host__ __device__ __cuBox(const cuINT3& e) { this->e = e; }
	__host__ __device__ __cuBox(int ex, int ey, int ez) { s = cuINT3(); e = cuINT3(ex, ey, ez); }

	__host__ __device__ __cuBox(const cuINT3& s, const cuINT3& e) { this->s = s; this->e = e; }
	__host__ __device__ __cuBox(int sx, int sy, int sz, int ex, int ey, int ez) { s = cuINT3(sx, sy, sz); e = cuINT3(ex, ey, ez); }

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__host__ __device__ __cuBox(const __cuBox &copyThis) { s = copyThis.s; e = copyThis.e; }

	//assignment operator
	__host__ __device__ __cuBox& operator=(const __cuBox &rhs) { s = rhs.s; e = rhs.e; return *this; }

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION

	__host__ operator Box() const
	{
		return Box((INT3)s, (INT3)e);
	}

	__host__ __cuBox& operator=(const Box& rhs)
	{
		s = (cuINT3)rhs.s; e = (cuINT3)rhs.e;
		return *this;
	}

	__host__ __cuBox(const Box& rhs)
	{
		s = (cuINT3)rhs.s; e = (cuINT3)rhs.e;
	}

	//----------------------------- COMPARISON OPERATORS

	//comparison operator
	__host__ __device__ bool operator==(const __cuBox &rhs) const { return (s == rhs.s && e == rhs.e); }
	__host__ __device__ bool operator!=(const __cuBox &rhs) const { return (s != rhs.s || e != rhs.e); }

	//comparisons with a cuINT3 : check start and end points against a point
	__host__ __device__ bool operator>=(const cuINT3 &rhs) const { return (s >= rhs); }
	__host__ __device__ bool operator>(const cuINT3 &rhs) const { return (s > rhs); }
	__host__ __device__ bool operator<=(const cuINT3 &rhs) const { return (e <= rhs); }
	__host__ __device__ bool operator<(const cuINT3 &rhs) const { return (e < rhs); }

	//comparisons with another cuBox : e.g. >= means start point must be greater or equal, and size must be greater or equal, etc..
	__host__ __device__ bool operator>=(const __cuBox &rhs) const { return (s >= rhs.s && size() >= rhs.size()); }
	__host__ __device__ bool operator>(const __cuBox &rhs) const { return (s > rhs.s && size() >= rhs.size()); }
	__host__ __device__ bool operator<=(const __cuBox &rhs) const { return (e <= rhs.e && size() <= rhs.size()); }
	__host__ __device__ bool operator<(const __cuBox &rhs) const { return (e < rhs.e && size() < rhs.size()); }

	//----------------------------- OTHERS

	//is this box contained in (or equal to) box rect ?
	__host__ __device__ bool IsContainedIn(const __cuBox& rect) const 
	{ 
		return (rect.s <= rect.e && s.x >= rect.s.x && s.y >= rect.s.y && s.z >= rect.s.z && e.x <= rect.e.x && e.y <= rect.e.y && e.z <= rect.e.z); 
	}

	__host__ __device__ bool IsPlane(void) const
	{
		return ((s.x == e.x && s.y != e.y && s.z != e.z) ||
			(s.y == e.y && s.x != e.x && s.z != e.z) ||
			(s.z == e.z && s.x != e.x && s.y != e.y));
	}

	__host__ __device__ bool IsLine(void) const { return ((cu_mod(e.x - s.x) + cu_mod(e.y - s.y) == 0) || (cu_mod(e.x - s.x) + cu_mod(e.z - s.z) == 0) || (cu_mod(e.z - s.z) + cu_mod(e.y - s.y) == 0)); }

	//is this a null box? (i.e. start and end coordinates are 0)
	__host__ __device__ bool IsNull(void) const { return (s == cuINT3(0) && e == cuINT3(0)); }

	__host__ __device__ bool Contains(int i, int j, int k) const { return Contains(cuINT3(i, j, k)); }

	//is the given coordinate inside this cuBox? Convention is : starting coordinates are contained, ending coordinates are outside.
	__host__ __device__ bool Contains(const cuINT3& coord) const
	{
		if (!IsPlane())
			//This is a box, so coord must be smaller than end point
			return (coord.x >= s.x && coord.y >= s.y && coord.z >= s.z && coord.x < e.x && coord.y < e.y && coord.z < e.z);
		else
			//Not a box, so we need to take into account dimensions of zero width
			return (coord.x >= s.x && coord.x <= e.x - (e.x != s.x) &&
				coord.y >= s.y && coord.y <= e.y - (e.y != s.y) &&
				coord.z >= s.z && coord.z <= e.z - (e.z != s.z));
	}

	//get box size
	__host__ __device__ cuINT3 size(void) const { return (e - s); }

	//resize box according to change in mesh dimensions (mesh which contains box) from n_old to n_new
	__host__ __device__ void resize(const cuINT3& n_old, const cuINT3& n_new)
	{
		bool x_flat, y_flat, z_flat;
		if (s.x != e.x) x_flat = false; else x_flat = true;
		if (s.y != e.y) y_flat = false; else y_flat = true;
		if (s.z != e.z) z_flat = false; else z_flat = true;

		s = (s & (cuReal3(n_new) / cuReal3(n_old)));
		e = (e & (cuReal3(n_new) / cuReal3(n_old)));

		if (s.x == e.x && !x_flat) e.x = s.x + 1;
		if (s.y == e.y && !y_flat) e.y = s.y + 1;
		if (s.z == e.z && !z_flat) e.z = s.z + 1;
	}

	//get intersection box
	__host__ __device__ __cuBox get_intersection(const __cuBox& box) const
	{
		//construct start point using largest coordinate values from the 2 rect start points
		cuINT3 start = cuINT3((box.s.x > s.x ? box.s.x : s.x), (box.s.y > s.y ? box.s.y : s.y), (box.s.z > s.z ? box.s.z : s.z));

		//construct end point using smallest coordinate values from the 2 rect end points
		cuINT3 end = cuINT3((box.e.x < e.x ? box.e.x : e.x), (box.e.y < e.y ? box.e.y : e.y), (box.e.z < e.z ? box.e.z : e.z));

		if (end >= start) return __cuBox(start, end);

		return __cuBox();
	}

	//get union box
	__host__ __device__ __cuBox get_union(const __cuBox& box) const
	{
		if (IsNull()) return box;
		if (box.IsNull()) return *this;

		//construct start point using smaller coordinate values from the 2 rect start points
		cuINT3 start = cuINT3((box.s.x < s.x ? box.s.x : s.x), (box.s.y < s.y ? box.s.y : s.y), (box.s.z < s.z ? box.s.z : s.z));

		//construct end point using largest coordinate values from the 2 rect end points
		cuINT3 end = cuINT3((box.e.x > e.x ? box.e.x : e.x), (box.e.y > e.y ? box.e.y : e.y), (box.e.z > e.z ? box.e.z : e.z));

		return __cuBox(start, end);
	}
};

typedef __cuBox<void> cuBox;

////////////////////////////////////////////////////////////////////////////////////////////////// cuRect
//
//similar to a cuBox but uses floating-point start and end points

template <typename Type = void>
struct __cuRect {

	//----------------------------- DATA

	cuReal3 s, e;

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(s, cuReal3());
		set_gpu_value(e, cuReal3());
	}

	__host__ void construct_cu_obj(const __cuRect& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const __cuRect& copyThis)
	{
		gpu_to_gpu(s, copyThis.s);
		gpu_to_gpu(e, copyThis.e);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ __cuRect(void) { s = cuReal3(); e = cuReal3(); }

	__host__ __device__ __cuRect(const cuReal3& e) { s = cuReal3(); this->e = e; }
	__host__ __device__ __cuRect(cuBReal ex, cuBReal ey, cuBReal ez) { s = cuReal3(); e = cuReal3(ex, ey, ez); }

	__host__ __device__ __cuRect(const cuReal3& s, const cuReal3& e) { this->s = s; this->e = e; }
	__host__ __device__ __cuRect(cuBReal sx, cuBReal sy, cuBReal sz, cuBReal ex, cuBReal ey, cuBReal ez) { s = cuReal3(sx, sy, sz); e = cuReal3(ex, ey, ez); }

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__host__ __device__ __cuRect(const __cuRect &copyThis) { s = copyThis.s; e = copyThis.e; }

	//assignment operator
	__host__ __device__ __cuRect& operator=(const __cuRect &rhs) { s = rhs.s; e = rhs.e; return *this; }

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION

	__host__ operator Rect() const
	{
		return Rect((DBL3)s, (DBL3)e);
	}

	__host__ __cuRect& operator=(const Rect& rhs)
	{
		s = (cuReal3)rhs.s; e = (cuReal3)rhs.e;
		return *this;
	}

	__host__ __cuRect(const Rect& rhs)
	{
		s = (cuReal3)rhs.s; e = (cuReal3)rhs.e;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//division by a cuVAL3
	template <typename VType> 
	__host__ __device__ cuReal3 operator/(const cuVAL3<VType> &rhs) const { return cuReal3((e.x - s.x) / rhs.x, (e.y - s.y) / rhs.y, (e.z - s.z) / rhs.z); }

	//products with a constant (must be fundamental type) on the LHS
	template <class MVType, std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr>
	__host__ __device__ friend __cuRect operator*(const MVType &mult, const __cuRect &rhs) { return __cuRect(rhs.s * mult, rhs.e * mult); }

	//products with a constant (must be fundamental type) on the RHS
	template <class MVType, std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr>
	__host__ __device__ __cuRect operator*(const MVType &mult) const { return __cuRect(s * mult, e * mult); }

	//add a rectangle : enlarges rectangle to contain both (but don't contain rects with all dimensions zero)
	__host__ __device__ void operator+=(const __cuRect &rhs)
	{
		if (rhs.s == rhs.e) return;

		if (s == e) {

			s = rhs.s;
			e = rhs.e;
			return;
		}

		if (s.x > rhs.s.x) s.x = rhs.s.x;
		if (s.y > rhs.s.y) s.y = rhs.s.y;
		if (s.z > rhs.s.z) s.z = rhs.s.z;

		if (e.x < rhs.e.x) e.x = rhs.e.x;
		if (e.y < rhs.e.y) e.y = rhs.e.y;
		if (e.z < rhs.e.z) e.z = rhs.e.z;
	}

	//shift the rectangle : add a cuReal3
	__host__ __device__ void operator+=(const cuReal3 &shift) { s += shift; e += shift; }

	//sum and difference with a cuReal3 : return shifted rectangle
	__host__ __device__ __cuRect operator+(const cuReal3 &rhs) const { return __cuRect(s + rhs, e + rhs); }
	__host__ __device__ __cuRect operator-(const cuReal3 &rhs) const { return __cuRect(s - rhs, e - rhs); }

	//sum and difference with a rectangle
	__host__ __device__ __cuRect operator+(const __cuRect& rhs) const { return __cuRect(s + rhs.s, e + rhs.e); }
	__host__ __device__ __cuRect operator-(const __cuRect& rhs) const { return __cuRect(s - rhs.s, e - rhs.e); }

	//----------------------------- COMPARISON OPERATORS

	//comparison operator
	__host__ __device__ bool operator==(const __cuRect &rhs) const { return (s == rhs.s && e == rhs.e); }
	__host__ __device__ bool operator!=(const __cuRect &rhs) const { return (s != rhs.s || e != rhs.e); }

	//comparisons with a cuReal3 : check start and end points against a point
	__host__ __device__ bool operator>=(const cuReal3 &rhs) const { return (s >= rhs); }
	__host__ __device__ bool operator>(const cuReal3 &rhs) const { return (s > rhs); }
	__host__ __device__ bool operator<=(const cuReal3 &rhs) const { return (e <= rhs); }
	__host__ __device__ bool operator<(const cuReal3 &rhs) const { return (e < rhs); }

	//comparisons with another cuRect : e.g. >= means start point must be greater or equal, and size must be greater or equal, etc..
	__host__ __device__ bool operator>=(const __cuRect &rhs) const { return (s >= rhs.s && size() >= rhs.size()); }
	__host__ __device__ bool operator>(const __cuRect &rhs) const { return (s > rhs.s && size() >= rhs.size()); }
	__host__ __device__ bool operator<=(const __cuRect &rhs) const { return (e <= rhs.e && size() <= rhs.size()); }
	__host__ __device__ bool operator<(const __cuRect &rhs) const { return (e < rhs.e && size() < rhs.size()); }

	//----------------------------- OTHERS

	//check if the rectangle contains the given coordinate
	__host__ __device__ bool contains(const cuReal3& position) const { return (position >= s && position <= e); }
	//check if the rectangle contains the given rectangle
	__host__ __device__ bool contains(const __cuRect& rectangle) const { return (rectangle.s >= s && rectangle.e <= e); }

	__host__ __device__ bool IsNull(void) const { return (s == cuReal3() && e == cuReal3()); }

	__host__ __device__ bool IsPlane(void) const
	{
		return ((cuIsZ(s.x - e.x) && cuIsNZ(s.y - e.y) && cuIsNZ(s.z - e.z)) ||
			(cuIsZ(s.y - e.y) && cuIsNZ(s.x - e.x) && cuIsNZ(s.z - e.z)) ||
			(cuIsZ(s.z - e.z) && cuIsNZ(s.x - e.x) && cuIsNZ(s.y - e.y)));
	}

	//check if the rectangle intersects with the given rectangle
	__host__ __device__ bool intersects(const __cuRect& rect) const { return (get_intersection(rect) != __cuRect()); }

	//get intersection rectangle
	__host__ __device__ __cuRect get_intersection(const __cuRect& rect) const
	{
		//construct start point using largest coordinate values from the 2 rect start points
		cuReal3 start = cuReal3((rect.s.x > s.x ? rect.s.x : s.x), (rect.s.y > s.y ? rect.s.y : s.y), (rect.s.z > s.z ? rect.s.z : s.z));

		//construct end point using smallest coordinate values from the 2 rect end points
		cuReal3 end = cuReal3((rect.e.x < e.x ? rect.e.x : e.x), (rect.e.y < e.y ? rect.e.y : e.y), (rect.e.z < e.z ? rect.e.z : e.z));

		if (end >= start) return __cuRect(start, end);

		return __cuRect();
	}

	//get union rectangle
	__host__ __device__ __cuRect get_union(const __cuRect& rect) const
	{
		if (IsNull()) return rect;
		if (rect.IsNull()) return *this;

		//construct start point using smaller coordinate values from the 2 rect start points
		cuReal3 start = cuReal3((rect.s.x < s.x ? rect.s.x : s.x), (rect.s.y < s.y ? rect.s.y : s.y), (rect.s.z < s.z ? rect.s.z : s.z));

		//construct end point using largest coordinate values from the 2 rect end points
		cuReal3 end = cuReal3((rect.e.x > e.x ? rect.e.x : e.x), (rect.e.y > e.y ? rect.e.y : e.y), (rect.e.z > e.z ? rect.e.z : e.z));

		return __cuRect(start, end);
	}

	//return the intersection volume
	__host__ __device__ cuBReal intersection_volume(const __cuRect &rect) const { return get_intersection(rect).volume(); }

	//get bottom-left quadrant as seen in xy-plane
	__host__ __device__ __cuRect get_quadrant_bl(void) const
	{
		return __cuRect(
			s,
			cuReal3((s.x + e.x) / 2, (s.y + e.y) / 2, e.z));
	}

	//get bottom-right quadrant as seen in xy-plane
	__host__ __device__ __cuRect get_quadrant_br(void) const
	{
		return __cuRect(
			cuReal3((s.x + e.x) / 2, s.y, s.z),
			cuReal3(e.x, (s.y + e.y) / 2, e.z));
	}

	//get top-left quadrant as seen in xy-plane
	__host__ __device__ __cuRect get_quadrant_tl(void) const
	{
		return __cuRect(
			cuReal3(s.x, (s.y + e.y) / 2, s.z),
			cuReal3((s.x + e.x) / 2, e.y, e.z));
	}

	//get top-right quadrant as seen in xy-plane
	__host__ __device__ __cuRect get_quadrant_tr(void) const
	{
		return __cuRect(
			cuReal3((s.x + e.x) / 2, (s.y + e.y) / 2, s.z),
			e);
	}

	//get rect face: -x face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__host__ __device__ __cuRect get_face_mx(cuBReal thickness = 0.0)
	{
		return __cuRect(s, cuReal3(s.x + thickness, e.y, e.z));
	}

	//get rect face: +x face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__host__ __device__ __cuRect get_face_px(cuBReal thickness = 0.0)
	{
		return __cuRect(cuReal3(e.x - thickness, s.y, s.z), e);
	}

	//get rect face: -y face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__host__ __device__ __cuRect get_face_my(cuBReal thickness = 0.0)
	{
		return __cuRect(s, cuReal3(e.x, s.y + thickness, e.z));
	}

	//get rect face: +y face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__host__ __device__ __cuRect get_face_py(cuBReal thickness = 0.0)
	{
		return __cuRect(cuReal3(s.x, e.y - thickness, s.z), e);
	}

	//get rect face: -z face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__host__ __device__ __cuRect get_face_mz(cuBReal thickness = 0.0)
	{
		return __cuRect(s, cuReal3(e.x, e.y, s.z + thickness));
	}

	//get rect face: +z face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__host__ __device__ __cuRect get_face_pz(cuBReal thickness = 0.0)
	{
		return __cuRect(cuReal3(s.x, s.y, e.z - thickness), e);
	}

	//are all the rectangle dimensions integer divisible by the given cuReal3 (direction by direction)?
	__host__ __device__ bool divisible(const cuReal3& rhs) const
	{
		return cuIsNZ(rhs.x) && cuIsNZ(rhs.y) && cuIsNZ(rhs.z) &&
			cuIsZ((e.x - s.x) / rhs.x - (cuBReal)round((e.x - s.x) / rhs.x)) &&
			cuIsZ((e.y - s.y) / rhs.y - (cuBReal)round((e.y - s.y) / rhs.y)) &&
			cuIsZ((e.z - s.z) / rhs.z - (cuBReal)round((e.z - s.z) / rhs.z));
	}

	//test intersection with a line from start to end and get intersection point closest to start
	__host__ __device__ bool intersection_test(const cuReal3& start, const cuReal3& end, cuReal3* pIntersection) const
	{
		//Note, only up to 3 faces need to be tested, hence the if checks below
		//also only one intersection is possible with the checked faces : the other intersection point (which is further from start) is on the "shadowed" faces and these are not checked.

		//test start x face
		if (start.x <= s.x) {

			(*pIntersection).x = s.x;
			if (cu_solve_line_equation_fixed_x(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test end x face
		if (start.x >= e.x) {

			(*pIntersection).x = e.x;
			if (cu_solve_line_equation_fixed_x(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test start y face
		if (start.y <= s.y) {

			(*pIntersection).y = s.y;
			if (cu_solve_line_equation_fixed_y(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test end y face
		if (start.y >= e.y) {

			(*pIntersection).y = e.y;
			if (cu_solve_line_equation_fixed_y(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test start z face
		if (start.z <= s.z) {

			(*pIntersection).z = s.z;
			if (cu_solve_line_equation_fixed_z(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test end z face
		if (start.z >= e.z) {

			(*pIntersection).z = e.z;
			if (cu_solve_line_equation_fixed_z(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		return false;
	}

	//get rectangle size
	__host__ __device__ cuReal3 size(void) const { return (e - s); }

	//these are used to get an lvalue for s and e
	__host__ __device__ cuReal3 get_s(void) const { return s; }
	__host__ __device__ cuReal3 get_e(void) const { return e; }

	//get the center point of the rect
	__host__ __device__ cuReal3 get_c(void) const { return (s + e) / 2; }

	__host__ __device__ cuBReal length(void) const { return (e.x - s.x); }
	__host__ __device__ cuBReal width(void) const { return (e.y - s.y); }
	__host__ __device__ cuBReal height(void) const { return (e.z - s.z); }
	__host__ __device__ cuBReal volume(void) const { return length()*width()*height(); }

	__host__ __device__ cuBReal maxDimension(void) const
	{
		cuBReal maxDim = length();

		maxDim = (width() > maxDim ? width() : maxDim);
		maxDim = (height() > maxDim ? height() : maxDim);

		return maxDim;
	}

	//resize rect according to change in rectangle sizes from old to new (i.e. the cuRect(s,e) is relative to the external rectangle which is changing from old to new)
	__host__ __device__ void resize(const cuReal3& size_old, const cuReal3& size_new)
	{
		cuReal3 scaling = size_new / size_old;
		s = s & scaling;
		e = e & scaling;
	}

	//snap coordinates to closest multiple of snap_unit
	__device__ void snap(cuBReal snap_unit)
	{
		s = cu_round(s / snap_unit) * snap_unit;
		e = cu_round(e / snap_unit) * snap_unit;
	}
};

typedef __cuRect<void> cuRect;
