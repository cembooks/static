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

#include "solver.hpp"

using namespace StaticScalarSolver;

void
SolverCVPII::make_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(Solver<2>::triangulation);

  std::ifstream ifs("../../gmsh/data/circle_r" + std::to_string(r) + ".msh");
  gridin.read_msh(ifs);

  for (auto cell : Solver<2>::triangulation.active_cell_iterators()) {
    cell->set_material_id(mid_1); // The cell is outside the
                                  // current-carrying region.

    if ((cell->center().norm() > a) && (cell->center().norm() < b))
      cell->set_material_id(mid_2); // The cell is inside the
                                    // current-carrying region.

    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++) {
      double dif_norm = 0.0;
      for (unsigned int v = 1; v < GeometryInfo<2>::vertices_per_face; v++)
        dif_norm += std::abs(cell->face(f)->vertex(0).norm() -
                             cell->face(f)->vertex(v).norm());

      if ((dif_norm < eps) && (cell->center().norm() > rd1))
        cell->face(f)->set_all_manifold_ids(1);
    }
  }

  Solver<2>::triangulation.set_manifold(1, sphere);
}

void
SolverCVPII::fill_dirichlet_stack()
{
  Solver<2>::dirichlet_stack = { { bid_dirichlet, &dirichlet_bc } };
}

void
SolverCVPII::solve()
{
  SolverControl control(Solver<2>::system_rhs.size(),
                        1e-12 * Solver<2>::system_rhs.l2_norm(),
                        false,
                        false);

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
