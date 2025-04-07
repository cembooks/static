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

template<>
void
SolverFLCAXI<true>::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/cylinder_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  double cell_r;
  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    cell_r = cell->center()[0];

    if (cell_r < d_1) {
      cell->set_material_id(mid_1);
    } else if (cell_r < d_2) {
      cell->set_material_id(mid_3);
    } else if (cell_r < b) {
      cell->set_material_id(mid_2);
    }
  }
}

template<>
void
SolverFLCAXI<false>::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  Solver<2>::triangulation.reset_all_manifolds();

  double cell_r;
  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    cell_r = cell->center().norm();

    if (cell_r < d_1) {
      cell->set_material_id(mid_1);
    } else if (cell_r < d_2) {
      cell->set_material_id(mid_3);
    } else if (cell_r < b) {
      cell->set_material_id(mid_2);
    }
  }

  Solver<2>::triangulation.set_all_manifold_ids(1);
  Solver<2>::triangulation.set_manifold(1, sphere);
}
