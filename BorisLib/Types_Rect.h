#pragma once

#include "Types_VAL.h"
#include "Funcs_Math.h"

////////////////////////////////////////////////////////////////////////////////////////////////// Box
//
// Box (or a Plane, or a Line, or a Point!, but probably only useful as an actual Box or Plane) - Discrete coordinates.

template <typename Type = void>
struct __Box {

	//----------------------------- DATA

	//start and end coordinates of discrete box (integer coordinates).
	INT3 s, e;

	//----------------------------- VALUE CONSTRUCTORS

	__Box(void) {}
	
	__Box(const INT3& e) { this->e = e; }
	__Box(int ex, int ey, int ez) { s = INT3(); e = INT3(ex, ey, ez); }
	
	__Box(const INT3& s, const INT3& e) { this->s = s; this->e = e; }
	__Box(int sx, int sy, int sz, int ex, int ey, int ez) { s = INT3(sx, sy, sz); e = INT3(ex, ey, ez); }

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__Box(const __Box &copyThis) { s = copyThis.s; e = copyThis.e; }

	//assignment operator
	__Box& operator=(const __Box &rhs) { s = rhs.s; e = rhs.e; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const __Box &rhs) { os << ToString(rhs.s) << "; " << ToString(rhs.e); return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const __Box &rhs) { lhs << rhs.s << std::string("; ") << rhs.e; return lhs; }

	//allows conversions from std::string to Box
	friend __Box& operator>>(const std::stringstream &ss, __Box &rhs)
	{ 
		//normally the values in std::string representation are as: "sx, sy, sz; ex, ey, ez"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "sx sy sz ex ey ez"
		if(components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings as would be done in the available constructors if the values were passed directly
		switch (components.size()) {

		case 1:			//Box(INT3 e) ctor
			rhs.e = ToNum(trimspaces(components[0]), "");
			rhs.s = INT3();
			break;
		case 2:			//Box(INT3 s, INT3 e) ctor
			rhs.s = ToNum(trimspaces(components[0]), "");
			rhs.e = ToNum(trimspaces(components[1]), "");
			break;
		case 3:			//Box(int ex, int ey, int ez) ctor
			rhs.e.x = ToNum(trimspaces(components[0]), "");
			rhs.e.y = ToNum(trimspaces(components[1]), "");
			rhs.e.z = ToNum(trimspaces(components[2]), "");
			rhs.s = INT3();
			break;
		case 6:			//Box(int sx, int sy, int sz, int ex, int ey, int ez) ctor
			rhs.s.x = ToNum(trimspaces(components[0]), "");
			rhs.s.y = ToNum(trimspaces(components[1]), "");
			rhs.s.z = ToNum(trimspaces(components[2]), "");
			rhs.e.x = ToNum(trimspaces(components[3]), "");
			rhs.e.y = ToNum(trimspaces(components[4]), "");
			rhs.e.z = ToNum(trimspaces(components[5]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- COMPARISON OPERATORS

	//comparison operator
	bool operator==(const __Box &rhs) const { return (s == rhs.s && e == rhs.e); }
	bool operator!=(const __Box &rhs) const { return (s != rhs.s || e != rhs.e); }

	//comparisons with a INT3 : check start and end points against a point
	bool operator>=(const INT3 &rhs) const { return (s >= rhs); }
	bool operator>(const INT3 &rhs) const { return (s > rhs); }
	bool operator<=(const INT3 &rhs) const { return (e <= rhs); }
	bool operator<(const INT3 &rhs) const { return (e < rhs); }

	//comparisons with another Box : e.g. >= means start point must be greater or equal, and size must be greater or equal, etc..
	bool operator>=(const __Box &rhs) const { return (s >= rhs.s && size() >= rhs.size()); }
	bool operator>(const __Box &rhs) const { return (s > rhs.s && size() >= rhs.size()); }
	bool operator<=(const __Box &rhs) const { return (e <= rhs.e && size() <= rhs.size()); }
	bool operator<(const __Box &rhs) const { return (e < rhs.e && size() < rhs.size()); }

	//----------------------------- GEOMETRY-RELATED CHECKS

	//is this box contained in (or equal to) box rect ?
	bool IsContainedIn(const __Box& rect) const { return (rect.s <= rect.e && s.x >= rect.s.x && s.y >= rect.s.y && s.z >= rect.s.z && e.x <= rect.e.x && e.y <= rect.e.y && e.z <= rect.e.z); }

	//is the given coordinate inside this Box? Convention is : starting coordinates are contained, ending coordinates are outside.
	bool Contains(const INT3& coord) const
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

	bool Contains(int i, int j, int k) const { return Contains(INT3(i, j, k)); }

	//is it a plane (2D) box?
	bool IsPlane(void) const
	{
		return ((s.x == e.x && s.y != e.y && s.z != e.z) ||
				(s.y == e.y && s.x != e.x && s.z != e.z) ||
				(s.z == e.z && s.x != e.x && s.y != e.y));
	}

	bool IsLine(void) const { return ((mod(e.x-s.x) + mod(e.y-s.y) == 0) || (mod(e.x-s.x) + mod(e.z-s.z) == 0) || (mod(e.z-s.z) + mod(e.y-s.y) == 0)); }

	bool IsPoint(void) const { return (e == s); }

	//is this a null box? (i.e. start and end coordinates are 0)
	bool IsNull(void) const { return (s == INT3(0) && e == INT3(0)); }

	//----------------------------- GET PROPERTIES

	//get box size
	INT3 size(void) const { return (e - s); }

	//----------------------------- GET PROPERTIES INVOLVING ANOTHER BOX

	//get intersection box
	__Box get_intersection(const __Box& box) const
	{
		//construct start point using largest coordinate values from the 2 rect start points
		INT3 start = INT3((box.s.x > s.x ? box.s.x : s.x), (box.s.y > s.y ? box.s.y : s.y), (box.s.z > s.z ? box.s.z : s.z));

		//construct end point using smallest coordinate values from the 2 rect end points
		INT3 end = INT3((box.e.x < e.x ? box.e.x : e.x), (box.e.y < e.y ? box.e.y : e.y), (box.e.z < e.z ? box.e.z : e.z));

		if (end >= start) return __Box(start, end);

		return __Box();
	}

	//get union box
	__Box get_union(const __Box& box) const
	{
		if (IsNull()) return box;
		if (box.IsNull()) return *this;

		//construct start point using smaller coordinate values from the 2 rect start points
		INT3 start = INT3((box.s.x < s.x ? box.s.x : s.x), (box.s.y < s.y ? box.s.y : s.y), (box.s.z < s.z ? box.s.z : s.z));

		//construct end point using largest coordinate values from the 2 rect end points
		INT3 end = INT3((box.e.x > e.x ? box.e.x : e.x), (box.e.y > e.y ? box.e.y : e.y), (box.e.z > e.z ? box.e.z : e.z));

		return __Box(start, end);
	}

	//----------------------------- MODIFIERS

	//resize box according to change in mesh dimensions (mesh which contains box) from n_old to n_new
	void resize(const INT3& n_old, const INT3& n_new)
	{
		bool x_flat, y_flat, z_flat;
		if (s.x != e.x) x_flat = false; else x_flat = true;
		if (s.y != e.y) y_flat = false; else y_flat = true;
		if (s.z != e.z) z_flat = false; else z_flat = true;

		s = (s & (DBL3(n_new) / DBL3(n_old)));
		e = (e & (DBL3(n_new) / DBL3(n_old)));

		if (s.x == e.x && !x_flat) e.x = s.x + 1;
		if (s.y == e.y && !y_flat) e.y = s.y + 1;
		if (s.z == e.z && !z_flat) e.z = s.z + 1;
	}
};

typedef __Box<void> Box;

////////////////////////////////////////////////////////////////////////////////////////////////// Rect
//
//similar to a Box but uses floating-point start and end points

template <typename Type = void>
struct __Rect {

	//----------------------------- DATA

	DBL3 s, e;

	//----------------------------- VALUE CONSTRUCTORS

	__Rect(void) { s = DBL3(); e = DBL3(); }
	
	__Rect(const DBL3& e) { s = DBL3(); this->e = e; }
	__Rect(double ex, double ey, double ez) { s = DBL3(); e = DBL3(ex, ey, ez); }

	__Rect(const DBL3& s, const DBL3& e) { this->s = s; this-> e = e; }
	__Rect(double sx, double sy, double sz, double ex, double ey, double ez) { s = DBL3(sx, sy, sz); e = DBL3(ex, ey, ez); }
	
	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__Rect(const __Rect &copyThis) { s = copyThis.s; e = copyThis.e; }

	//assignment operator
	__Rect& operator=(const __Rect &rhs) { s = rhs.s; e = rhs.e; return *this; }

	//----------------------------- STREAM OPERATORS
	
	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const __Rect &rhs) { os << ToString(rhs.s) << "; " << ToString(rhs.e); return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const __Rect &rhs) { lhs << rhs.s << std::string("; ") << rhs.e; return lhs; }

	//allows conversions from std::string to Rect
	friend __Rect& operator>>(const std::stringstream &ss, __Rect &rhs)
	{ 
		//normally the values in std::string representation are as: "sx, sy, sz; ex, ey, ez"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "sx sy sz ex ey ez"
		if(components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings as would be done in the available constructors if the values were passed directly
		switch (components.size()) {

		case 1:			//Rect(INT3 e) ctor
			rhs.e = ToNum(trimspaces(components[0]), "");
			rhs.s = DBL3();
			break;
		case 2:			//Rect(INT3 s, INT3 e) ctor
			rhs.s = ToNum(trimspaces(components[0]), "");
			rhs.e = ToNum(trimspaces(components[1]), "");
			break;
		case 3:			//Rect(int ex, int ey, int ez) ctor
			rhs.e.x = ToNum(trimspaces(components[0]), "");
			rhs.e.y = ToNum(trimspaces(components[1]), "");
			rhs.e.z = ToNum(trimspaces(components[2]), "");
			rhs.s = DBL3();
			break;
		case 6:			//Rect(int sx, int sy, int sz, int ex, int ey, int ez) ctor
			rhs.s.x = ToNum(trimspaces(components[0]), "");
			rhs.s.y = ToNum(trimspaces(components[1]), "");
			rhs.s.z = ToNum(trimspaces(components[2]), "");
			rhs.e.x = ToNum(trimspaces(components[3]), "");
			rhs.e.y = ToNum(trimspaces(components[4]), "");
			rhs.e.z = ToNum(trimspaces(components[5]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//division by a VAL3
	template <typename VType> DBL3 operator/(const VAL3<VType> &rhs) const { return DBL3((e.x - s.x) / rhs.x, (e.y - s.y) / rhs.y, (e.z - s.z) / rhs.z); }

	//products with a constant (must be fundamental type) on the LHS
	template <class MVType, std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr>
	friend __Rect operator*(const MVType &mult, const __Rect &rhs) { return __Rect(rhs.s * mult, rhs.e * mult); }

	//products with a constant (must be fundamental type) on the RHS
	template <class MVType, std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr>
	__Rect operator*(const MVType &mult) const { return __Rect(s * mult, e * mult); }

	//add a rectangle : enlarges rectangle to contain both (but don't contain rects with all dimensions zero)
	void operator+=(const __Rect &rhs)
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

	//shift the rectangle : add a DBL3
	void operator+=(const DBL3 &shift) { s += shift; e += shift; }

	//sum and difference with a DBL3 : return shifted rectangle
	__Rect operator+(const DBL3 &rhs) const { return __Rect(s + rhs, e + rhs); }
	__Rect operator-(const DBL3 &rhs) const { return __Rect(s - rhs, e - rhs); }

	//sum and difference with a rectangle
	__Rect operator+(const __Rect& rhs) const { return __Rect(s + rhs.s, e + rhs.e); }
	__Rect operator-(const __Rect& rhs) const { return __Rect(s - rhs.s, e - rhs.e); }

	//----------------------------- COMPARISON OPERATORS
	
	//comparison operator
	bool operator==(const __Rect &rhs) const { return (s == rhs.s && e == rhs.e); }
	bool operator!=(const __Rect &rhs) const { return (s != rhs.s || e != rhs.e); }

	//comparisons with a DBL3 : check start and end points against a point
	bool operator>=(const DBL3 &rhs) const { return (s >= rhs); }
	bool operator>(const DBL3 &rhs) const { return (s > rhs); }
	bool operator<=(const DBL3 &rhs) const { return (e <= rhs); }
	bool operator<(const DBL3 &rhs) const { return (e < rhs); }

	//comparisons with another Rect : e.g. >= means start point must be greater or equal, and size must be greater or equal, etc..
	bool operator>=(const __Rect &rhs) const { return (s >= rhs.s && size() >= rhs.size()); }
	bool operator>(const __Rect &rhs) const { return (s > rhs.s && size() >= rhs.size()); }
	bool operator<=(const __Rect &rhs) const { return (e <= rhs.e && size() <= rhs.size()); }
	bool operator<(const __Rect &rhs) const { return (e < rhs.e && size() < rhs.size()); }

	//----------------------------- GEOMETRY-RELATED CHECKS

	//check if the rectangle contains the given coordinate
	bool contains(const DBL3& position) const { return (position >= s && position <= e); }
	
	//check if the rectangle contains the given rectangle
	bool contains(const __Rect& rectangle) const { return (rectangle.s >= s && rectangle.e <= e); }

	bool IsPlane(void) const
	{
		return ((IsZ(s.x - e.x) && IsNZ(s.y - e.y) && IsNZ(s.z - e.z)) ||
				(IsZ(s.y - e.y) && IsNZ(s.x - e.x) && IsNZ(s.z - e.z)) ||
				(IsZ(s.z - e.z) && IsNZ(s.x - e.x) && IsNZ(s.y - e.y)));
	}

	bool IsLine(void) const { return (IsZ(mod(e.x - s.x) + mod(e.y - s.y)) || IsZ(mod(e.x - s.x) + mod(e.z - s.z)) || IsZ(mod(e.z - s.z) + mod(e.y - s.y))); }

	bool IsPoint(void) const { return (s == e); }

	bool IsNull(void) const { return (s == DBL3() && e == DBL3()); }

	//check if the rectangle intersects with the given rectangle
	bool intersects(const __Rect& rect) const { return (get_intersection(rect) != __Rect()); }

	//are all the rectangle dimensions integer divisible by the given DBL3 (direction by direction)?
	bool divisible(const DBL3& rhs) const
	{
		return IsNZ(rhs.x) && IsNZ(rhs.y) && IsNZ(rhs.z) &&
			IsZ((e.x - s.x) / rhs.x - (double)round((e.x - s.x) / rhs.x)) &&
			IsZ((e.y - s.y) / rhs.y - (double)round((e.y - s.y) / rhs.y)) &&
			IsZ((e.z - s.z) / rhs.z - (double)round((e.z - s.z) / rhs.z));
	}

	//----------------------------- GET PROPERTIES

	//get rectangle size
	DBL3 size(void) const { return (e - s); }

	//these are used to get an lvalue for s and e
	DBL3 get_s(void) const { return s; }
	DBL3 get_e(void) const { return e; }

	//get the center point of the rect
	DBL3 get_c(void) const { return (s + e) / 2; }

	double length(void) const { return (e.x - s.x); }
	double width(void) const { return (e.y - s.y); }
	double height(void) const { return (e.z - s.z); }
	double volume(void) const { return length()*width()*height(); }

	//get maximum dimension (from length, width, height)
	double maxDimension(void) const
	{
		double maxDim = length();

		maxDim = (width() > maxDim ? width() : maxDim);
		maxDim = (height() > maxDim ? height() : maxDim);

		return maxDim;
	}

	//get bottom-left quadrant as seen in xy-plane
	__Rect get_quadrant_bl(void) const
	{
		return __Rect(
			s,
			DBL3((s.x + e.x) / 2, (s.y + e.y) / 2, e.z));
	}

	//get bottom-right quadrant as seen in xy-plane
	__Rect get_quadrant_br(void) const
	{
		return __Rect(
			DBL3((s.x + e.x) / 2, s.y, s.z),
			DBL3(e.x, (s.y + e.y) / 2, e.z));
	}

	//get top-left quadrant as seen in xy-plane
	__Rect get_quadrant_tl(void) const
	{
		return __Rect(
			DBL3(s.x, (s.y + e.y) / 2, s.z),
			DBL3((s.x + e.x) / 2, e.y, e.z));
	}

	//get top-right quadrant as seen in xy-plane
	__Rect get_quadrant_tr(void) const
	{
		return __Rect(
			DBL3((s.x + e.x) / 2, (s.y + e.y) / 2, s.z),
			e);
	}

	//get rect face: -x face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__Rect get_face_mx(double thickness = 0.0) const
	{
		return __Rect(s, DBL3(s.x + thickness, e.y, e.z));
	}

	//get rect face: +x face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__Rect get_face_px(double thickness = 0.0) const
	{
		return __Rect(DBL3(e.x - thickness, s.y, s.z), e);
	}

	//get rect face: -y face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__Rect get_face_my(double thickness = 0.0) const
	{
		return __Rect(s, DBL3(e.x, s.y + thickness, e.z));
	}

	//get rect face: +y face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__Rect get_face_py(double thickness = 0.0) const
	{
		return __Rect(DBL3(s.x, e.y - thickness, s.z), e);
	}

	//get rect face: -z face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__Rect get_face_mz(double thickness = 0.0) const
	{
		return __Rect(s, DBL3(e.x, e.y, s.z + thickness));
	}

	//get rect face: +z face; This is a plane rectangle, but it can be thickened inside the rect using the thickness value.
	__Rect get_face_pz(double thickness = 0.0) const
	{
		return __Rect(DBL3(s.x, s.y, e.z - thickness), e);
	}

	//return rectangle of given z layer, where layer_idx must range from 0 to (e.z - s.z) / thickness - 1; if outside of this range then return the entire rectangle.
	__Rect get_zlayer(int layer_idx, double thickness) const
	{
		if (layer_idx * thickness < 0 || (layer_idx + 1) * thickness > (e.z - s.z)) return *this;
		else return __Rect(DBL3(s.x, s.y, s.z + layer_idx * thickness), DBL3(e.x, e.y, s.z + (layer_idx + 1) * thickness));
	}

	//----------------------------- GET PROPERTIES INVOLVING ANOTHER RECT

	//get intersection rectangle
	__Rect get_intersection(const __Rect& rect) const
	{
		//construct start point using largest coordinate values from the 2 rect start points
		DBL3 start = DBL3((rect.s.x > s.x ? rect.s.x : s.x), (rect.s.y > s.y ? rect.s.y : s.y), (rect.s.z > s.z ? rect.s.z : s.z));

		//construct end point using smallest coordinate values from the 2 rect end points
		DBL3 end = DBL3((rect.e.x < e.x ? rect.e.x : e.x), (rect.e.y < e.y ? rect.e.y : e.y), (rect.e.z < e.z ? rect.e.z : e.z));

		if(end >= start) return __Rect(start, end);
		
		return __Rect();
	}

	//get union rectangle
	__Rect get_union(const __Rect& rect) const
	{
		if (IsNull()) return rect; 
		if (rect.IsNull()) return *this;

		//construct start point using smaller coordinate values from the 2 rect start points
		DBL3 start = DBL3((rect.s.x < s.x ? rect.s.x : s.x), (rect.s.y < s.y ? rect.s.y : s.y), (rect.s.z < s.z ? rect.s.z : s.z));

		//construct end point using largest coordinate values from the 2 rect end points
		DBL3 end = DBL3((rect.e.x > e.x ? rect.e.x : e.x), (rect.e.y > e.y ? rect.e.y : e.y), (rect.e.z > e.z ? rect.e.z : e.z));

		return __Rect(start, end);
	}

	//return the intersection volume
	double intersection_volume(const __Rect &rect) const { return get_intersection(rect).volume(); }

	//from this __Rect get a face to be fully contained in the given rect. If none found return empty __Rect.
	//If this __Rect is not fully contained in the given rect then only one face can be fully contained - this method is meant for these type of checks
	//if any found, indicate which face it is in the second pair argument as 0: -x, 1: +x, 2: -y, 3: +y, 4: -z, 5: +z
	std::pair<__Rect, int> get_contained_face(const __Rect &rect) const
	{
		//the intersection of the two rects
		__Rect rect_i = get_intersection(rect);

		if (get_face_mx() == rect_i.get_face_mx()) return std::pair<__Rect, int>(get_face_mx(), 0);
		if (get_face_px() == rect_i.get_face_px()) return std::pair<__Rect, int>(get_face_px(), 1);

		if (get_face_my() == rect_i.get_face_my()) return std::pair<__Rect, int>(get_face_my(), 2);
		if (get_face_py() == rect_i.get_face_py()) return std::pair<__Rect, int>(get_face_py(), 3);

		if (get_face_mz() == rect_i.get_face_mz()) return std::pair<__Rect, int>(get_face_mz(), 4);
		if (get_face_pz() == rect_i.get_face_pz()) return std::pair<__Rect, int>(get_face_pz(), 5);

		//none found
		return std::pair<__Rect, int>(__Rect(), -1);
	}

	//check if the given rect fully contains a face of this __Rect - use in conjuction with the above method if needed
	bool has_contained_face(const __Rect &rect) const { return (get_contained_face(rect).second >= 0); }

	//----------------------------- GET OTHER PROPERTIES

	//test intersection with a line from start to end and get intersection point closest to start
	bool intersection_test(const DBL3& start, const DBL3& end, DBL3* pIntersection) const
	{
		//Note, only up to 3 faces need to be tested, hence the if checks below
		//also only one intersection is possible with the checked faces : the other intersection point (which is further from start) is on the "shadowed" faces and these are not checked.
		
		//test start x face
		if (start.x <= s.x) {

			(*pIntersection).x = s.x;
			if (solve_line_equation_fixed_x(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test end x face
		if (start.x >= e.x) {

			(*pIntersection).x = e.x;
			if (solve_line_equation_fixed_x(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test start y face
		if (start.y <= s.y) {

			(*pIntersection).y = s.y;
			if (solve_line_equation_fixed_y(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test end y face
		if (start.y >= e.y) {

			(*pIntersection).y = e.y;
			if (solve_line_equation_fixed_y(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test start z face
		if (start.z <= s.z) {

			(*pIntersection).z = s.z;
			if (solve_line_equation_fixed_z(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		//test end z face
		if (start.z >= e.z) {

			(*pIntersection).z = e.z;
			if (solve_line_equation_fixed_z(start, end, pIntersection) && *pIntersection >= s && *pIntersection <= e) return true;
		}

		return false;
	}

	//----------------------------- MODIFIERS

	//resize rect according to change in rectangle sizes from old to new (i.e. the Rect(s,e) is relative to the external rectangle which is changing from old to new)
	void resize(const DBL3& size_old, const DBL3& size_new)
	{
		DBL3 scaling = size_new / size_old;
		s = s & scaling;
		e = e & scaling;
	}

	//snap coordinates to closest multiple of snap_unit
	void snap(double snap_unit)
	{
		s = round(s / snap_unit) * snap_unit;
		e = round(e / snap_unit) * snap_unit;
	}

	//given a containing rect which has a grid specified by the cellsize h, snap the coordinates of this __Rect to grid by increasing size
	void snap_to_grid_up(const __Rect& rect, const DBL3& h)
	{
		//rect must contain this __Rect
		if (!rect.contains(*this)) return;

		//number of cells in each direction of containing rect
		DBL3 ratio = rect.size() / h;
		INT3 cells = INT3(round(ratio.x), round(ratio.y), round(ratio.z));

		//move start point to grid - enlarge so use floor
		s = (floor((s - rect.s) / h) & h) + rect.s;

		//move end point to grid - enlarge so use ceil
		e = (ceil((e - rect.s) / h) & h) + rect.s;
	}

	//given a containing rect which has a grid specified by the cellsize h, snap the coordinates of this __Rect to grid by decreasing size
	void snap_to_grid_down(const __Rect& rect, const DBL3& h)
	{
		//rect must contain this __Rect
		if (!rect.contains(*this)) return;

		//number of cells in each direction of containing rect
		DBL3 ratio = rect.size() / h;
		INT3 cells = INT3(round(ratio.x), round(ratio.y), round(ratio.z));

		//move start point to grid - decrease so use ceil
		s = (ceil((s - rect.s) / h) & h) + rect.s;

		//move end point to grid - decrease so use floor
		e = (floor((e - rect.s) / h) & h) + rect.s;
	}

};

typedef __Rect<void> Rect;
