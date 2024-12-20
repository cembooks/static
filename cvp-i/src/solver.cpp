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

using namespace StaticVectorSolver;

void
SolverCVPI::make_mesh()
{
  GridIn<3> gridin;
  gridin.attach_triangulation(Solver1<3>::triangulation);

  std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  Solver1<3>::triangulation.reset_all_manifolds();

  for (auto cell : Solver1<3>::triangulation.active_cell_iterators()) {
    cell->set_material_id(mid_1); // The cell is outside the coil.

    if ((cell->center().norm() > a1) && (cell->center().norm() < a2))
      cell->set_material_id(mid_2); // The cell is inside the coil.

    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++) {
      double dif_norm_a1 = 0.0;
      double dif_norm_a2 = 0.0;
      double dif_norm_b = 0.0;

      for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face; v++) {
        dif_norm_a1 += std::abs(cell->face(f)->vertex(v).norm() - a1);
        dif_norm_a2 += std::abs(cell->face(f)->vertex(v).norm() - a2);
        dif_norm_b += std::abs(cell->face(f)->vertex(v).norm() - b);
      }

      if ((dif_norm_a1 < eps) || (dif_norm_a2 < eps) || (dif_norm_b < eps))
        cell->face(f)->set_all_manifold_ids(1);
    }
  }

  Solver1<3>::triangulation.set_manifold(1, sphere);
}

void
SolverCVPI::fill_dirichlet_stack()
{
  Solver1<3>::dirichlet_stack = { { bid_dirichlet, &dirichlet_bc } };
}

void
SolverCVPI::solve()
{
  ReductionControl control(
    Solver1<3>::system_rhs.size(), 1e-8, 1e-5, false, false);

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

  if (SettingsCVPI::log_cg_convergence) {
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
