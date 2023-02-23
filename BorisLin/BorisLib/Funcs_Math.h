#pragma once

#include "Funcs_Math_base.h"
#include "Types_VAL.h"
#include "Funcs_Algorithms.h"
#include "Funcs_Aux_base.h"

///////////////////////////////////////////////////////////////////////////////

//get magnitude for fundamental types - need this to make VEC work for both fundamental and composite types
template <typename Type, std::enable_if_t<std::is_integral<Type>::value>* = nullptr>
Type GetMagnitude(const Type& V)
{
	return V * (2 * (V >= 0) - 1);
}

//get magnitude for fundamental types - need this to make VEC work for both fundamental and composite types
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type GetMagnitude(const Type& V)
{
	return fabs(V);
}

//get magnitude for fundamental types using 2 components
template <typename Type, std::enable_if_t<std::is_fundamental<Type>::value>* = nullptr>
Type GetMagnitude(const Type& Vx, const Type& Vy)
{
	return sqrt(Vx*Vx + Vy * Vy);
}

//get magnitude for VAL2 types (a VAL2 type can be identified by checking if an INT2 can be converted to it.
template <typename VAL2Type, std::enable_if_t<std::is_convertible<INT2, VAL2Type>::value>* = nullptr>
auto GetMagnitude(const VAL2Type& V) -> decltype(std::declval<VAL2Type>().y)
{
	return sqrt(V.x*V.x + V.y*V.y);
}

//get magnitude for fundamental types using 3 components
template <typename Type, std::enable_if_t<std::is_fundamental<Type>::value>* = nullptr>
Type GetMagnitude(const Type& Vx, const Type& Vy, const Type& Vz)
{ 
	return sqrt(Vx*Vx + Vy*Vy + Vz*Vz); 
}

//get magnitude for VAL3 types (a VAL3 type can be identified by checking if an INT3 can be converted to it.
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
auto GetMagnitude(const VAL3Type& V) -> decltype(std::declval<VAL3Type>().z)
{ 
	return sqrt(V.x*V.x + V.y*V.y + V.z*V.z);
}

//get magnitude for a std::vector
template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
double GetMagnitude(std::vector<VType>& V)
{
	double mag = 0;

#pragma omp parallel for reduction(+:mag)
	for (int idx = 0; idx < V.size(); idx++) {

		mag += V[idx] * V[idx];
	}

	return sqrt(mag);
}

//SFINAE sink-hole
inline double GetMagnitude(...) { return 0.0; }

//obtain distance between two Cartesian coordinates specified as a VAL2
template <typename VAL2Type, std::enable_if_t<std::is_convertible<INT2, VAL2Type>::value>* = nullptr>
auto get_distance(const VAL2Type& coord1, const VAL2Type& coord2) -> decltype(std::declval<VAL2Type>().y)
{
	return (coord2 - coord1).norm();
}

//obtain distance between two Cartesian coordinates specified as a VAL3
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
auto get_distance(const VAL3Type& coord1, const VAL3Type& coord2) -> decltype(std::declval<VAL3Type>().z)
{ 
	return (coord2 - coord1).norm();
}

//obtain distance between two Cartesian coordinates specified using std::vector
template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
double get_distance(std::vector<VType>& coord1, std::vector<VType>& coord2)
{
	if (coord1.size() != coord2.size()) return VType();

	double distance = 0;

#pragma omp parallel for reduction(+:distance)
	for (int idx = 0; idx < coord1.size(); idx++) {

		distance += (coord1[idx] - coord2[idx]) * (coord1[idx] - coord2[idx]);
	}

	return sqrt(distance);
}

//get sign of value, returning -1, 0 or 1.
template <typename Type>
int get_sign(const Type& value, std::true_type)
{
	return (Type(0) <= value) - (value <= Type(0));
}

template <typename Type> 
int get_sign(const Type& value, std::false_type)
{
	return Type(0) < value;
}

template <typename Type> 
int get_sign(const Type& value)
{
	return get_sign(value, std::is_signed<Type>());
}

///////////////////////////////////////////////////////////////////////////////

template <typename Type, std::enable_if_t<std::is_fundamental<Type>::value>* = nullptr>
void NormalizeVector(Type &Vx, Type &Vy, Type &Vz)
{ 	
	Type magnitude = GetMagnitude(Vx, Vy, Vz);
	
	if(IsNZ(magnitude)) {

		Vx /= magnitude; 
		Vy /= magnitude;
		Vz /= magnitude; 
	} 
}

//normalize a VAL3 if possible, else return zero value
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type normalize(const VAL3Type& val3) { return (val3 != VAL3Type() ? val3.normalized() : VAL3Type()); }

//normalize a VAL3 to the norm of another VAL3
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type normalize(const VAL3Type& val3, const VAL3Type& val3_2) { return (val3_2 != VAL3Type() ? val3 / val3_2.norm() : VAL3Type()); }

//normalize a VAL33 to the norm of a VAL3
template <typename VAL33Type, typename VAL3Type, std::enable_if_t<std::is_convertible<INT33, VAL33Type>::value>* = nullptr>
VAL33Type normalize(const VAL33Type& val33, const VAL3Type& val3) { return (val3 != VAL3Type() ? val33 / val3.norm() : VAL33Type()); }

///////////////////////////////////////////////////////////////////////////////

//Note : round is now a standard function since C++14 (in math.h), so no need to define it here as per previous C++ standards
//previous:
//inline int round(float fval) { return (fval < 0.0) ? (int)ceil(fval - 0.5) : (int)floor(fval + 0.5); }
//inline int round(double fval) { return (fval < 0.0) ? (int)ceil(fval - 0.5) : (int)floor(fval + 0.5); }

template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type round(const VAL3Type& val3) { return VAL3Type( round(val3.x), round(val3.y), round(val3.z) ); }

//"fixed" floor function.
//Note, using just the standard floor is not good enough : if the floating point value is very close to the upper integer value (closer than the defined floating point accuracy) then its value should be equal to it.
//The standard floor function will result in "wrong" behaviour by returning the lower integer
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type floor_epsilon(const Type& fval) { return floor(fval + EPSILON_ROUNDING); }

//as above but with ceil
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type ceil_epsilon(const Type& fval) { return ceil(fval - EPSILON_ROUNDING); }

//return "fixed" floor of each VAL3 component
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type floor(const VAL3Type& fval) { return VAL3Type( floor_epsilon(fval.x), floor_epsilon(fval.y), floor_epsilon(fval.z) ); }

//return "fixed" ceil of each DBL3 component.
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type ceil(const VAL3Type& fval) { return VAL3Type( ceil_epsilon(fval.x), ceil_epsilon(fval.y), ceil_epsilon(fval.z) ); }

//absolute value for a floating point number
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type mod(const Type& fval) { return fabs(fval); }

//absolute value for an integer
template <typename Type, std::enable_if_t<std::is_integral<Type>::value>* = nullptr>
Type mod(const Type& ival) { return ival * (2*(ival >= 0) - 1); }

//absolute value for a VAL3
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type mod(const VAL3Type& fval) { return VAL3Type(mod(fval.x), mod(fval.y), mod(fval.z)); }

//remainder after division (using floating point numbers) - fixed version of fmod from <cmath>
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type fmod_epsilon(const Type& fval, Type denom) { return fval - floor_epsilon(fval / denom) * denom; }

//fixed fmod for a VAL3
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type fmod_epsilon(const VAL3Type& fval, double denom) { return VAL3Type(fmod_epsilon(fval.x, denom), fmod_epsilon(fval.y, denom), fmod_epsilon(fval.z, denom)); }

///////////////////////////////////////////////////////////////////////////////

//convert from polar (magnitude, polar, azimuthal) coordinates, where angles are in degrees, to Cartesian coordinates (x, y, z)
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3Type Polar_to_Cartesian(const VAL3Type& polar) 
{ 
	return VAL3Type(
		polar.x*sin(polar.y*PI/180) * cos(polar.z*PI/180),
		polar.x*sin(polar.y*PI/180) * sin(polar.z*PI/180),
		polar.x*cos(polar.y*PI/180)); 
}

//convert from polar (magnitude, azimuthal, i.e. r, theta) coordinates, where angles are in degrees, to Cartesian coordinates (x, y)
template <typename VAL2Type, std::enable_if_t<std::is_convertible<INT2, VAL2Type>::value>* = nullptr>
VAL2Type Polar_to_Cartesian(const VAL2Type& polar)
{
	return VAL2Type(
		polar.x*cos(polar.y*PI / 180),
		polar.x*sin(polar.y*PI / 180));
}

//convert from Cartesian (x, y) coordinates, to polar coordinates (magnitude, azimuthal, i.e. r, theta) where angles are in degrees in interval [-180 , 180]
template <typename VAL2Type, std::enable_if_t<std::is_convertible<INT2, VAL2Type>::value>* = nullptr>
VAL2Type Cartesian_to_Polar(const VAL2Type& cartesian)
{
	if (cartesian.x > 0) {

		//right half : -90 to 90
		return VAL2Type(
			sqrt(cartesian.x * cartesian.x + cartesian.y * cartesian.y),
			atan(cartesian.y / cartesian.x) * 180.0 / PI
			);
	}
	else if (cartesian.x < 0) {

		//left half

		if (cartesian.y > 0) {

			//upper-left quadrant : 90 to 180
			return VAL2Type(
				sqrt(cartesian.x * cartesian.x + cartesian.y * cartesian.y),
				atan(-cartesian.x / cartesian.y) * 180.0 / PI + 90.0);
		}
		else {

			//lower-left quadrant : -180 to -90
			return VAL2Type(
				sqrt(cartesian.x * cartesian.x + cartesian.y * cartesian.y),
				atan(cartesian.y / cartesian.x) * 180.0 / PI - 180.0);
		}
	}
	else {

		if (cartesian.y > 0) {

			return VAL2Type(cartesian.y, 90.0);
		}
		else if (cartesian.y < 0) {

			return VAL2Type(cartesian.y, -90.0);
		}
		else return VAL2Type(0.0, 0.0);
	}
}

///////////////////////////////////////////////////////////////////////////////

//Kahan sum
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type sum_Kahan(std::vector<Type> &sum_this)
{
	//Kahan summation

	Type sum = Type();
	
	//A running compensation for lost low-order bits.
	Type c = Type();

	for(int i = 0; i < sum_this.size(); i++) {
		
		//new term to add to sum with correction from previous iteration
		Type y = sum_this[i] - c;

		//add term to temporary
		Type t = sum + y;
		
		//find correction for lost low-order bits -> to apply at next iteration
		c = (t - sum) - y;

		//now set in running sum
		sum = t;
	}
	
	return sum;
}

//Kahan sum modified by Neumaier
template <typename Type, std::enable_if_t<std::is_floating_point<Type>::value>* = nullptr>
Type sum_KahanNeumaier(std::vector<Type> &sum_this)
{
	//Kahan summation with Neumaier modification

	Type sum = Type();

	//A running compensation for lost low-order bits.
	Type c = Type();

	for (int i = 0; i < sum_this.size(); i++) {

		//new sum in temporary
		Type t = sum + sum_this[i];

		if (GetMagnitude(sum) > GetMagnitude(sum_this[i])) {

			//compensation for lost low-order bits in sum_this[i]
			c += (sum - t) + sum_this[i];
		}
		else {

			//compensation for lost low-order bits in sum
			c += (sum_this[i] - t) + sum;
		}

		//now set in running sum
		sum = t;
	}

	//apply accumulated correction before returning sum
	return (sum + c);
}

///////////////////////////////////////////////////////////////////////////////

//Use linear interpolation/extrapolation to obtain the y value at the given x, where
//p0 = (x0, y0) and p1 = (x1, y1)
//y = [(x1 - x) * y0 - (x0 - x) * y1] / (x1 - x0)
//p0 and p1 cannot have the same x coordinate - this is not checked here
template <typename VAL2Type, std::enable_if_t<std::is_convertible<INT2, VAL2Type>::value>* = nullptr>
auto interpolate(const VAL2Type& p0, const VAL2Type& p1, decltype(std::declval<VAL2Type>().x) x) -> decltype(std::declval<VAL2Type>().y)
{
	return ((p1.x - x) * p0.y - (p0.x - x) * p1.y) / (p1.x - p0.x);
}

//parametric interpolation where parameter = 0 gives start, parameter = 1 gives end.
template <typename VType>
VType parametric_interpolation(VType start, VType end, double parameter)
{
	return (start * (1 - parameter) + end * parameter);
}

//bi-linear interpolation at point p = (x,y) from points p11 = (x1, y1), p12 = (x1, y2), p21 = (x2, y1), p22 = (x2, y2).
//A function f takes values at these points as f11, f12, f21, f22. Return f(p).
//Must not have x1 == x2 or y1 == y2 : this is not checked here.
template <typename VType>
VType interpolate_bilinear(double x1, double y1, double x2, double y2, VType f11, VType f12, VType f21, VType f22, DBL2 p)
{
	if (IsNZ(x2 - x1) && IsNZ(y2 - y1)) {

		return ((x2 - p.x) * (f11 * (y2 - p.y) + f12 * (p.y - y1)) + (p.x - x1) * (f21 * (y2 - p.y) + f22 * (p.y - y1))) / ((x2 - x1) * (y2 - y1));
	}
	else {

		if (IsNZ(x2 - x1) && IsZ(y2 - y1)) {

			return (f21 * (p.x - x1) + f11 * (x2 - p.x)) / (x2 - x1);
		}
		else if (IsZ(x2 - x1) && IsNZ(y2 - y1)) {

			return (f12 * (p.y - y1) + f11 * (y2 - p.y)) / (y2 - y1);
		}
	}

	return f11;
}

///////////////////////////////////////////////////////////////////////////////

//solve line equation for line passing through point1 and point2, where the x coordinate has been provided in *pSolutionPoint - fill the other 2 and return true. If no solution then return false.
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
bool solve_line_equation_fixed_x(const VAL3Type& point1, const VAL3Type& point2, VAL3Type* pSolutionPoint)
{
	typedef decltype(std::declval<VAL3Type>().z) Type;

	VAL3Type direction = point2 - point1;

	if (IsNZ(direction.x)) {

		Type t = ((*pSolutionPoint).x - point1.x) / direction.x;
		*pSolutionPoint = direction * t + point1;

		return true;
	}
	else return false;
}

//solve line equation for line passing through point1 and point2, where the y coordinate has been provided in *pSolutionPoint - fill the other 2 and return true. If no solution then return false.
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
bool solve_line_equation_fixed_y(const VAL3Type& point1, const VAL3Type& point2, VAL3Type* pSolutionPoint)
{
	typedef decltype(std::declval<VAL3Type>().z) Type;

	VAL3Type direction = point2 - point1;

	if (IsNZ(direction.y)) {

		Type t = ((*pSolutionPoint).y - point1.y) / direction.y;
		*pSolutionPoint = direction * t + point1;

		return true;
	}
	else return false;
}

//solve line equation for line passing through point1 and point2, where the z coordinate has been provided in *pSolutionPoint - fill the other 2 and return true. If no solution then return false.
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
bool solve_line_equation_fixed_z(const VAL3Type& point1, const VAL3Type& point2, VAL3Type* pSolutionPoint)
{
	typedef decltype(std::declval<VAL3Type>().z) Type;

	VAL3Type direction = point2 - point1;

	if (IsNZ(direction.z)) {

		Type t = ((*pSolutionPoint).z - point1.z) / direction.z;
		*pSolutionPoint = direction * t + point1;

		return true;
	}
	else return false;
}

///////////////////////////////////////////////////////////////////////////////

//check if point p is on the right-hand-side of a given line, determined by the points s and e, in the oriented plane containing the three points with normal n
//the meaning of right-hand-side is: If looking towards the plane against the plane normal direction, if the points s, p, e (in this order!) are clockwise then we have LHS, else RHS
inline bool point_on_rhs_of_line(const DBL3& s, const DBL3& e, const DBL3& p, const DBL3& n)
{
	return ((e - s) ^ (p - s)) * n < 0;
}

//special case of the above where the plane is the x-y plane with normal z.
inline bool point_on_rhs_of_line(const DBL2& s, const DBL2& e, const DBL2& p)
{
	return ((DBL3(e.x, e.y, 0) - DBL3(s.x, s.y, 0)) ^ (DBL3(p.x, p.y, 0) - DBL3(s.x, s.y, 0))) * DBL3(0, 0, 1) < 0;
}

//check if point p is on the right-hand-side of plane defined by the points p1, p2, p3 (in this order!)
//the meaning of rhs then is if p is on the side of the plane containing the surface normal direction n, where n is along this cross product: (p2 - p1) ^ (p3 - p1)
//i.e. if you look towards the plane which has points p1 p2 p3 in counter-clockwise direction then you are on the rhs!
inline bool point_on_rhs_of_plane(const DBL3& p1, const DBL3& p2, const DBL3& p3, const DBL3& p)
{
	return ((p2 - p1) ^ (p3 - p1)) * (p - p1) > 0;
}

///////////////////////////////////////////////////////////////////////////////

template <typename VType, std::enable_if_t<std::is_floating_point<VType>::value>* = nullptr>
std::pair<DBL2, DBL2> linear_regression(std::vector<VType>& x_vals, std::vector<VType>& y_vals, int startIdx = 0, int lastIdx = 0)
{
	DBL2 gradient, intercept;

	double xy = 0, x = 0, y = 0, xsq = 0;

	int N;

	if (lastIdx <= 0) {

		N = (x_vals.size() < y_vals.size() ? x_vals.size() : y_vals.size());
	}
	else N = lastIdx - startIdx + 1;

	if (N < 2) return { gradient, intercept };

#pragma omp parallel for reduction(+:x, y, xy, xsq)
	for (int idx = startIdx; idx < startIdx + N; idx++) {

		x += x_vals[idx];
		y += y_vals[idx];

		xy += x_vals[idx] * y_vals[idx];
		xsq += x_vals[idx] * x_vals[idx];
	}

	//calculate gradient and intercept
	double denom = N * xsq - x * x;
	gradient.i = (N*xy - x * y) / denom;
	intercept.i = (y - gradient.i * x) / N;

	double y_dist_sq = 0;

#pragma omp parallel for reduction(+:y_dist_sq)
	for (int idx = startIdx; idx < startIdx + N; idx++) {

		y_dist_sq += pow(y_vals[idx] - (gradient.i * x_vals[idx] + intercept.i), 2);
	}

	//calculate uncertainties in gradient and intercept
	if (N > 2) {

		gradient.j = sqrt((N / (N - 2)) * y_dist_sq / denom);
		intercept.j = gradient.j * sqrt(xsq / N);
	}

	return { gradient, intercept };
}

///////////////////////////////////////////////////////////////////////////////

//return matrix multiplication of rank-3 unit antisymmetric tensor with a VAL3 - return type is a VAL33
template <typename VAL3Type, std::enable_if_t<std::is_convertible<INT3, VAL3Type>::value>* = nullptr>
VAL3<VAL3Type> epsilon3(const VAL3Type& val3)
{
	return VAL3<VAL3Type>(
		VAL3Type(0, val3.z, -val3.y),
		VAL3Type(-val3.z, 0, val3.x),
		VAL3Type(val3.y, -val3.x, 0)
		);
}

///////////////////////////////////////////////////////////////////////////////

//overloaded inverse method so can be used in templated methods
template <typename Type, std::enable_if_t<std::is_same<Type, double>::value>* = nullptr>
double inverse(const double& m)
{
	if (m) return 1 / m;
	else return 0.0;
}

//3x3 matrix inverse
template <typename Type, std::enable_if_t<std::is_same<Type, DBL33>::value>* = nullptr>
DBL33 inverse(const DBL33& m)
{
	//in DBL33 we use row-major notation: first index is the row, the second is the column

	double d_11 = m.j.j * m.k.k - m.j.k * m.k.j;
	double d_21 = m.j.k * m.k.i - m.j.i * m.k.k;
	double d_31 = m.j.i * m.k.j - m.j.j * m.k.i;

	double det = m.i.i * d_11 + m.i.j * d_21 + m.i.k * d_31;
	if (IsZ(det)) return DBL33();

	double d_12 = m.i.k * m.k.j - m.i.j * m.k.k;
	double d_22 = m.i.i * m.k.k - m.i.k * m.k.i;
	double d_32 = m.i.j * m.k.i - m.i.i * m.k.j;

	double d_13 = m.i.j * m.j.k - m.i.k * m.j.j;
	double d_23 = m.i.k * m.j.i - m.i.i * m.j.k;
	double d_33 = m.i.i * m.j.j - m.i.j * m.j.i;

	return DBL33(
		DBL3(d_11, d_12, d_13) / det,
		DBL3(d_21, d_22, d_23) / det,
		DBL3(d_31, d_32, d_33) / det);
}

//overloaded so can be used in templated methods : simple identity
template <typename Type, std::enable_if_t<std::is_same<Type, double>::value>* = nullptr>
double ident(void) { return 1.0; }

//3x3 matrix identity
template <typename Type, std::enable_if_t<std::is_same<Type, DBL33>::value>* = nullptr>
DBL33 ident(void)
{ 
	return DBL33(
		DBL3(1.0, 0.0, 0.0),
		DBL3(0.0, 1.0, 0.0),
		DBL3(0.0, 0.0, 1.0));
}

///////////////////////////////////////////////////////////////////////////////

//This function returns the solution of s = a * m^s + b * m^m^s + f
//i.e. solve for s, where m, s, f are DBL3, ^ is the cross product, a and b are constants; moreover m is a unit vector.
inline DBL3 solve_crossprod(double a, double b, const DBL3& m, const DBL3& f)
{
	double ab = a * a + b + b * b;

	return f + ((a*a + b*b) / (a*a + ab*ab)) * (a * (m ^ f) + ab * (m ^ (m ^ f)));
}

//This function returns the solution of s = a * m^s + b * m^m^s + f
//i.e. solve for s, where m, s, f are DBL3, ^ is the cross product, a and b are constants; moreover m is a unit vector perpendicular to f, so that m ^ m ^ f = -f
inline DBL3 solve_crossprod_perp(double a, double b, const DBL3& m, const DBL3& f)
{
	double ab = a * a + b + b * b;

	return (1.0 / (a*a + ab*ab)) * ((a*a + b*ab) * m + a * (a*a + b*b) * (m ^ f));
}

///////////////////////////////////////////////////////////////////////////////

//rotate object through angles in radians as psi: around y, theta: around x, phi: around z
inline DBL3 rotate_object_yxz(const DBL3& r, double psi, double theta, double phi)
{
	return DBL3(
		(cos(psi) * cos(phi) + sin(psi) * sin(theta) * sin(phi)) * r.x + (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)) * r.y + (cos(theta) * sin(psi)) * r.z,
		(cos(theta) * sin(phi)) * r.x + (cos(theta) * cos(phi)) * r.y - sin(theta) * r.z,
		(cos(psi) * sin(theta) * sin(phi) - cos(phi) * sin(psi)) * r.x + (cos(psi) * cos(phi) * sin(theta) + sin(psi) * sin(phi)) * r.y + (cos(psi) * cos(theta)) * r.z
	);
}

//rotate the vector r by setting the coordinate system x axis to given polar and azimuthal angles. The angles are in radians.
//i.e. if we have the unit vector n = [ cos(phi)sin(theta), sin(phi)sin(theta), cos(theta) ]
//then we rotate the coordinate system s.t. the new x axis is n : 1) rotate the x axis around z through azimuthal, obtaining new x' and y' axes. 2) rotate x' axis around y' axis through (pi/2 - polar), obtaining new x'', y'=y'', z'' axes
//then the rotated vector r' (returned vector) is that vector which has the same components in the new coordinate system as the vector r has in the original system.
inline DBL3 rotate_polar(const DBL3& r, double theta, double phi)
{
	return DBL3(
		cos(phi)*sin(theta) * r.x - sin(phi) * r.y - cos(phi)*cos(theta) * r.z,
		sin(phi)*sin(theta) * r.x + cos(phi) * r.y - sin(phi)*cos(theta) * r.z,
		cos(theta)		    * r.x - 0              + sin(theta)	 		 * r.z
	);
}

//same as above but the rotation is specified using the unit vector n
inline DBL3 rotate_polar(const DBL3& r, const DBL3& n)
{
	//if n = [cos(phi)sin(theta), sin(phi)sin(theta), cos(theta)]
	//where theta ranges in [0, PI], and phi ranges in [0, 2*PI] then:
	
	//nxy is sin(theta)
	double nxy = sqrt(n.x*n.x + n.y*n.y);

	//then sin(phi) = ny / nxy, and cos(phi) = nx / nxy
	//n.z = cos(theta)

	if (nxy > 0) {

		return DBL3(
			n.x * r.x - (n.y / nxy) * r.y - (n.x / nxy) * n.z * r.z,
			n.y * r.x + (n.x / nxy) * r.y - (n.y / nxy) * n.z * r.z,
			n.z * r.x - 0				  + nxy				  * r.z
		);
	}
	else {

		//special case where n = [0, 0, n.z]

		if (n.z > 0) {

			return DBL3(-r.z, r.y, r.x);
		}
		else {

			return DBL3(+r.z, r.y, -r.x);
		}
	}
}

//This is the inverse operation of rotate_polar.
//i.e. for rotate polar we have A * r = r_rot, where A is the rotation matrix for theta and phi angles.
//Then here we use the A^-1 matrix, so that if r_rot is passed as argument we get back the r vector (theta and phi same angles as for rotate_polar)
inline DBL3 invrotate_polar(const DBL3& r, double theta, double phi)
{
	return DBL3(
		cos(phi)*sin(theta) * r.x + sin(phi) * sin(theta) * r.y + cos(theta) * r.z,
		-sin(phi)			* r.x + cos(phi)			  * r.y + 0,
		-cos(phi)*cos(theta)* r.x - sin(phi)*cos(theta)   * r.y + sin(theta) * r.z
	);
}

//same as above but the rotation is specified using the unit vector n
inline DBL3 invrotate_polar(const DBL3& r, const DBL3& n)
{
	//if n = [cos(phi)sin(theta), sin(phi)sin(theta), cos(theta)]
	//where theta ranges in [0, PI], and phi ranges in [0, 2*PI] then:

	//nxy is sin(theta)
	double nxy = sqrt(n.x*n.x + n.y*n.y);

	//then sin(phi) = ny / nxy, and cos(phi) = nx / nxy
	//n.z = cos(theta)

	if (nxy > 0) {

		return DBL3(
			n.x * r.x + n.y * r.y + n.z	* r.z,
			-(n.y / nxy) * r.x + (n.x / nxy) * r.y,
			-(n.x / nxy) * n.z * r.x - (n.y / nxy) * n.z * r.y + nxy * r.z
		);
	}
	else {

		//special case where n = [0, 0, n.z]

		if (n.z > 0) {

			return DBL3(+r.z, r.y, -r.x);
		}
		else {

			return DBL3(-r.z, r.y, +r.x);
		}
	}
}

//return a vector which is rotated by theta (0 to PI) from the vector r and by phi around the vector r (0 to 2PI)
//s.t. if r is along the x axis, theta is a rotation in the xy plane, and phi is a rotation around the x axis (geometric sense rotations).
inline DBL3 relrotate_polar(const DBL3& r, double theta, double phi)
{
	//rotation around the x axis first
	DBL3 rdash = DBL3(cos(theta), sin(theta)*cos(phi), sin(phi)*sin(theta)) * r.norm();

	//now the vector to return is that vector obtained by rotating the x axis into r
	return rotate_polar(rdash, r.normalized());
}