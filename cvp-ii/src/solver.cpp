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

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include "solver.hpp"
#include <fstream>

using namespace StaticScalarSolver;

void
SolverCVPII::make_mesh()
{
  GridIn<2> gridin;
  Triangulation<2> tria_tmp;

  gridin.attach_triangulation(tria_tmp);
  std::ifstream ifs("../../gmsh/data/circle_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  std::tuple<std::vector<Point<2>>, std::vector<CellData<2>>, SubCellData>
    mesh_description;

  mesh_description = GridTools::get_coarse_mesh_description(tria_tmp);

  GridTools::invert_all_negative_measure_cells(std::get<0>(mesh_description),
                                               std::get<1>(mesh_description));

  GridTools::consistently_order_cells(std::get<1>(mesh_description));

  Solver<2>::triangulation.create_triangulation(std::get<0>(mesh_description),
                                                std::get<1>(mesh_description),
                                                std::get<2>(mesh_description));

  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    cell->set_material_id(
      mid_1); // The cell is outside the current-carrying region.

    if ((cell->center().norm() > a1) && (cell->center().norm() < a2))
      cell->set_material_id(
        mid_2); // The cell is inside the current-carrying region.

    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++) {
      double dif_norm_a1 = 0.0;
      double dif_norm_a2 = 0.0;
      double dif_norm_b = 0.0;

      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; v++) {
        dif_norm_a1 += std::abs(cell->face(f)->vertex(v).norm() - a1);
        dif_norm_a2 += std::abs(cell->face(f)->vertex(v).norm() - a2);
        dif_norm_b += std::abs(cell->face(f)->vertex(v).norm() - b);
      }

      if ((dif_norm_a1 < eps) || (dif_norm_a2 < eps) || (dif_norm_b < eps))
        cell->face(f)->set_all_manifold_ids(1);
    }
  }

  Solver<2>::triangulation.set_manifold(1, sphere);

  GridOut gridout;
  GridOutFlags::Msh msh_flags(true, true);
  gridout.set_flags(msh_flags);

  std::ofstream ofs("../../gmsh/data/circle_r" + std::to_string(r) +
                    "_reordered.msh");
  gridout.write_msh(Solver<2>::triangulation, ofs);
}

void
SolverCVPII::fill_dirichlet_stack()
{
  Solver<2>::dirichlet_stack = { { bid_dirichlet, &dirichlet_bc } };
}

void
SolverCVPII::solve()
{
  ReductionControl control(
    Solver<2>::system_rhs.size(), 0.0, 1e-8, false, false);

  if (Settings::log_cg_convergence)
    control.enable_history_data();

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>> cg(control, memory);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(Solver<2>::system_matrix, 1.2);

  cg.solve(Solver<2>::system_matrix,
           Solver<2>::solution,
           Solver<2>::system_rhs,
           preconditioner);

  Solver<2>::constraints.distribute(Solver<2>::solution);

  if (SettingsCVPII::log_cg_convergence) {
    const std::vector<double> history_data = control.get_history_data();

    std::ofstream ofs(fname + "_cg_convergence.csv");

    unsigned int i = 1;
    for (auto item : history_data) {
      ofs << i << ", " << item << "\n";
      i++;
    }
    ofs.close();
  }
}
