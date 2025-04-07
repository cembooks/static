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
SolverSCHAXI<true>::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/cylinder_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    if (std::abs(cell->center()[0]) < a) {
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
        double dif_norm = 0.0;
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; ++v) {
          dif_norm += std::abs(cell->face(f)->vertex(v)[0] - a);
        }

        if (dif_norm < eps) {
          cell->face(f)->set_user_index(2);
          cell->set_user_index(1);
        }
      }
    }
  }
  GridGenerator::hyper_cube(triangulation_slice, 0.0, b);
}

template<>
void
SolverSCHAXI<false>::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  Solver<2>::triangulation.reset_all_manifolds();

  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
      double dif_norm = 0.0;
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; v++)
        dif_norm += std::abs(cell->face(f)->vertex(0).norm() -
                             cell->face(f)->vertex(v).norm());

      if ((dif_norm < eps) && (cell->center().norm() > rd)) {

        cell->face(f)->set_all_manifold_ids(1);

        double dif_norm_a = 0.0;
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; v++)
          dif_norm_a += std::abs(cell->face(f)->vertex(v).norm() - a);

        if ((dif_norm < eps) && (std::abs(cell->center().norm()) < a)) {
          cell->face(f)->set_user_index(2);
          cell->set_user_index(1);
        }
      }
    }
  }

  Solver<2>::triangulation.set_manifold(1, sphere);

  GridGenerator::hyper_cube(triangulation_slice, 0.0, b - 0.01);
}
