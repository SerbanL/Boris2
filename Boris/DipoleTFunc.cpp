#include "stdafx.h"
#include "DipoleTFunc.h"

double DipoleTFunc::Nxx(DBL3 xyz, DBL3 abc)
{
	std::vector<double> sum_this(8, 0);

	sum_this[0] = atan(fx(xyz & DBL3(-1, 1, 1), abc));
	sum_this[1] = atan(fx(xyz & DBL3(1, 1, 1), abc));
	sum_this[2] = atan(fx(xyz & DBL3(-1, -1, 1), abc));
	sum_this[3] = atan(fx(xyz & DBL3(-1, 1, -1), abc));
	sum_this[4] = atan(fx(xyz & DBL3(1, -1, 1), abc));
	sum_this[5] = atan(fx(xyz & DBL3(-1, -1, -1), abc));
	sum_this[6] = atan(fx(xyz & DBL3(1, 1, -1), abc));
	sum_this[7] = atan(fx(xyz & DBL3(1, -1, -1), abc));

	return (sum_KahanNeumaier(sum_this) / (4 * PI));
}

double DipoleTFunc::Nyy(DBL3 xyz, DBL3 abc)
{
	std::vector<double> sum_this(8, 0);

	sum_this[0] = atan(fy(xyz & DBL3(-1, 1, 1), abc));
	sum_this[1] = atan(fy(xyz & DBL3(1, 1, 1), abc));
	sum_this[2] = atan(fy(xyz & DBL3(-1, -1, 1), abc));
	sum_this[3] = atan(fy(xyz & DBL3(-1, 1, -1), abc));
	sum_this[4] = atan(fy(xyz & DBL3(1, -1, 1), abc));
	sum_this[5] = atan(fy(xyz & DBL3(-1, -1, -1), abc));
	sum_this[6] = atan(fy(xyz & DBL3(1, 1, -1), abc));
	sum_this[7] = atan(fy(xyz & DBL3(1, -1, -1), abc));

	return (sum_KahanNeumaier(sum_this) / (4 * PI));
}

double DipoleTFunc::Nzz(DBL3 xyz, DBL3 abc)
{
	std::vector<double> sum_this(8, 0);

	sum_this[0] = atan(fz(xyz & DBL3(-1, 1, 1), abc));
	sum_this[1] = atan(fz(xyz & DBL3(1, 1, 1), abc));
	sum_this[2] = atan(fz(xyz & DBL3(-1, -1, 1), abc));
	sum_this[3] = atan(fz(xyz & DBL3(-1, 1, -1), abc));
	sum_this[4] = atan(fz(xyz & DBL3(1, -1, 1), abc));
	sum_this[5] = atan(fz(xyz & DBL3(-1, -1, -1), abc));
	sum_this[6] = atan(fz(xyz & DBL3(1, 1, -1), abc));
	sum_this[7] = atan(fz(xyz & DBL3(1, -1, -1), abc));

	return (sum_KahanNeumaier(sum_this) / (4 * PI));
}

double DipoleTFunc::Nxy(DBL3 xyz, DBL3 abc)
{
	std::vector<double> sum_this(8, 0);

	sum_this[0] = log(Fxy(xyz, abc));
	sum_this[1] = log(Fxy(xyz, abc & DBL3(-1, -1, 1)));
	sum_this[2] = log(Fxy(xyz, abc & DBL3(1, -1, -1)));
	sum_this[3] = log(Fxy(xyz, abc & DBL3(-1, 1, -1)));
	sum_this[4] = -log(Fxy(xyz, abc & DBL3(1, -1, 1)));
	sum_this[5] = -log(Fxy(xyz, abc & DBL3(-1, 1, 1)));
	sum_this[6] = -log(Fxy(xyz, abc & DBL3(1, 1, -1)));
	sum_this[7] = -log(Fxy(xyz, abc & DBL3(-1, -1, -1)));

	return (sum_KahanNeumaier(sum_this) / (4 * PI));
}

double DipoleTFunc::Nxz(DBL3 xyz, DBL3 abc)
{
	std::vector<double> sum_this(8, 0);

	sum_this[0] = log(Fxz(xyz, abc));
	sum_this[1] = log(Fxz(xyz, abc & DBL3(-1, -1, 1)));
	sum_this[2] = log(Fxz(xyz, abc & DBL3(1, -1, -1)));
	sum_this[3] = log(Fxz(xyz, abc & DBL3(-1, 1, -1)));
	sum_this[4] = -log(Fxz(xyz, abc & DBL3(1, -1, 1)));
	sum_this[5] = -log(Fxz(xyz, abc & DBL3(-1, 1, 1)));
	sum_this[6] = -log(Fxz(xyz, abc & DBL3(1, 1, -1)));
	sum_this[7] = -log(Fxz(xyz, abc & DBL3(-1, -1, -1)));

	return (sum_KahanNeumaier(sum_this) / (4 * PI));
}

double DipoleTFunc::Nyz(DBL3 xyz, DBL3 abc)
{
	std::vector<double> sum_this(8, 0);

	sum_this[0] = log(Fyz(xyz, abc));
	sum_this[1] = log(Fyz(xyz, abc & DBL3(-1, -1, 1)));
	sum_this[2] = log(Fyz(xyz, abc & DBL3(1, -1, -1)));
	sum_this[3] = log(Fyz(xyz, abc & DBL3(-1, 1, -1)));
	sum_this[4] = -log(Fyz(xyz, abc & DBL3(1, -1, 1)));
	sum_this[5] = -log(Fyz(xyz, abc & DBL3(-1, 1, 1)));
	sum_this[6] = -log(Fyz(xyz, abc & DBL3(1, 1, -1)));
	sum_this[7] = -log(Fyz(xyz, abc & DBL3(-1, -1, -1)));

	return (sum_KahanNeumaier(sum_this) / (4 * PI));
}

double DipoleTFunc::fx(DBL3 xyz, DBL3 abc)
{
	return (abc.y / 2 - xyz.y)*(abc.z / 2 - xyz.z) /
		((abc.x / 2 - xyz.x) * sqrt((abc.x / 2 - xyz.x)*(abc.x / 2 - xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z)));
}

double DipoleTFunc::fy(DBL3 xyz, DBL3 abc)
{
	return (abc.x / 2 - xyz.x)*(abc.z / 2 - xyz.z) /
		((abc.y / 2 - xyz.y) * sqrt((abc.x / 2 - xyz.x)*(abc.x / 2 - xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z)));
}

double DipoleTFunc::fz(DBL3 xyz, DBL3 abc)
{
	return (abc.y / 2 - xyz.y)*(abc.x / 2 - xyz.x) /
		((abc.z / 2 - xyz.z) * sqrt((abc.x / 2 - xyz.x)*(abc.x / 2 - xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z)));
}

double DipoleTFunc::Fxy(DBL3 xyz, DBL3 abc)
{
	return (abc.z / 2 - xyz.z) + sqrt((abc.x / 2 + xyz.x)*(abc.x / 2 + xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z));
}

double DipoleTFunc::Fxz(DBL3 xyz, DBL3 abc)
{
	return (abc.y / 2 - xyz.y) + sqrt((abc.x / 2 + xyz.x)*(abc.x / 2 + xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z));
}

double DipoleTFunc::Fyz(DBL3 xyz, DBL3 abc)
{
	return (abc.x / 2 + xyz.x) + sqrt((abc.x / 2 + xyz.x)*(abc.x / 2 + xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z));
}
