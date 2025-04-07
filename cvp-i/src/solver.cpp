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

    if ((cell->center().norm() > a) && (cell->center().norm() < b))
      cell->set_material_id(mid_2); // The cell is inside the coil.

    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++) {
      double dif_norm = 0.0;
      for (unsigned int v = 1; v < GeometryInfo<3>::vertices_per_face; v++)
        dif_norm += std::abs(cell->face(f)->vertex(0).norm() -
                             cell->face(f)->vertex(v).norm());

      if ((dif_norm < eps) && (cell->center().norm() > rd1))
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
  SolverControl control(1000 * Solver1<3>::system_rhs.size(),
                        1e-6 * Solver1<3>::system_rhs.l2_norm(),
                        false,
                        false);

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
