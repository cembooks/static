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
SolverRHOAXI<true>::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/cylinder_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    if (std::abs(cell->center()[0]) < a) {
      cell->set_material_id(mid_2);
    } else {
      cell->set_material_id(mid_1);
    }
  }

  GridGenerator::hyper_cube(triangulation_slice, 0.0, b);
}

template<>
void
SolverRHOAXI<false>::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  Solver<2>::triangulation.reset_all_manifolds();

  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    if (cell->center().norm() < a) {
      cell->set_material_id(mid_2);

      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++) {
        double dif_norm = 0.0;
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; v++)
          dif_norm += std::abs(cell->face(f)->vertex(v).norm() - a);

        if (dif_norm < eps)
          cell->face(f)->set_all_manifold_ids(1);
      }
    } else {
      cell->set_material_id(mid_1);

      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++)
        if (cell->face(f)->at_boundary() &&
            (cell->face(f)->boundary_id() == bid))
          cell->face(f)->set_all_manifold_ids(1);
    }
  }

  Solver<2>::triangulation.set_manifold(1, sphere);

  GridGenerator::hyper_cube(triangulation_slice, 0.0, b - 0.01);
}
