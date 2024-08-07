/******************************************************************************
* Copyright (C) Siarhei Uzunbajakau, 2023.
*
* This program is free software. You can use, modify, and redistribute it under
* the terms of the GNU Lesser General Public License as published by the Free
* Software Foundation, either version 3 or (at your option) any later version.
* This program is distributed without any warranty.
*
* Refer to COPYING.LESSER for more details.
******************************************************************************/

#include "exact_solution.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

using namespace dealii;

double ExactSolutionCVPII_T::value(const Point<2> & r,
	const unsigned int component) const
{
	if (r.norm() < SettingsCVPII::a1)
		return 0.5*(a2*a2-a1*a1);

	if (r.norm() > SettingsCVPII::a2)
		return 0.0;

	return -0.5*(r[0]*r[0] +r[1]*r[1] - a2*a2);
}

Tensor<1, 2> ExactSolutionCVPII_T::gradient(const Point<2> & r,
	const unsigned int component) const
{
	if (r.norm() < SettingsCVPII::a1)
		return Point<2>();

	if (r.norm() > SettingsCVPII::a2)
		return Point<2>();

	return -r;
}

ExactSolutionCVPII_Jf::ExactSolutionCVPII_Jf() : Function<2>(2) {}

void ExactSolutionCVPII_Jf::vector_value_list(const std::vector<Point<2>> & r,
	std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	auto v = values.begin();
	for (auto p: r)
	{
		if ((p.norm() > SettingsCVPII::a1) &&
				(p.norm() < SettingsCVPII::a2))
		{
			(*v)[0] =-p[1];
			(*v)[1] = p[0];
		}
		else
		{
			(*v)[0] = 0.0;
			(*v)[1] = 0.0;
		}

		v++;
	}
}

