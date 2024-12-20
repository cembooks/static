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

void
SolverSSOLI::make_mesh()
{
  GridIn<3> gridin;

  gridin.attach_triangulation(Solver1<3>::triangulation);
  std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  mark_materials();
}

void
SolverSSOLI::fill_dirichlet_stack()
{
  dirichlet_stack = {};
}

void
SolverSSOLI::solve()
{
  ReductionControl control(
    Solver1<3>::system_rhs.size(), 0.0, 1e-8, false, false);

  if (log_cg_convergence)
    control.enable_history_data();

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>> cg(control, memory);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(Solver1<3>::system_matrix, 1.2);

  cg.solve(Solver1<3>::system_matrix,
           Solver1<3>::solution,
           Solver1<3>::system_rhs,
           preconditioner);

  Solver1<3>::constraints.distribute(Solver1<3>::solution);

  if (SettingsSSOLI::log_cg_convergence) {
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

void
SolverSSOLI::mark_materials()
{
  Solver1<3>::triangulation.reset_all_manifolds();

  for (auto cell : Solver1<3>::triangulation.active_cell_iterators()) {
    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f) {
      double dif_norm_a = 0.0;

      for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face; v++)
        dif_norm_a += std::abs(cell->face(f)->vertex(v).norm() - a);

      if (dif_norm_a < eps)
        if (std::abs(cell->center().norm()) < a) {
          cell->face(f)->set_user_index(1);
          cell->set_user_index(1);
        }
    }
  }
}
