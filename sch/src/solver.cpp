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

#include "solver.hpp"

template <>
void SolverSCH<2>::make_mesh()
{
	GridIn<2> gridin;
	gridin.attach_triangulation(triangulation);
	std::ifstream ifs("../../gmsh/data/circle_r" + std::to_string(r) + ".msh");
	gridin.read_msh(ifs);
	mark_materials();
}

template <>
void SolverSCH<3>::make_mesh()
{
	GridIn<3> gridin;
	gridin.attach_triangulation(triangulation);
	std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
	gridin.read_msh(ifs);
	mark_materials();
}

