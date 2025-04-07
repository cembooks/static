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

void
SolverMMSVTII_T::make_mesh()
{
  GridIn<2> gridin;

  gridin.attach_triangulation(Solver<2, 0>::triangulation);
#if DOMAIN__ == 0
  std::ifstream ifs("../../gmsh/data/square_r" + std::to_string(r) + ".msh");
#elif DOMAIN__ == 1
  std::ifstream ifs("../../gmsh/data/circle_r" + std::to_string(r) + ".msh");
#endif
  gridin.read_msh(ifs);

#if DOMAIN__ == 1
  Solver<2, 0>::triangulation.reset_all_manifolds();

  double dif_norm;

  for (auto cell : Solver<2, 0>::triangulation.active_cell_iterators()) {
    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++) {
      dif_norm = std::abs(cell->face(f)->vertex(0).norm() -
                          cell->face(f)->vertex(1).norm());

      if ((dif_norm < eps) && (cell->center().norm() > rd1))
        cell->face(f)->set_all_manifold_ids(1);
    }
  }

  Solver<2, 0>::triangulation.set_manifold(1, sphere);
#endif
}

void
SolverMMSVTII_T::fill_dirichlet_stack()
{
  Solver<2, 0>::dirichlet_stack = { { bid_dirichlet, &dirichlet_bc } };
}

void
SolverMMSVTII_T::solve()
{
  SolverControl control(1000 * Solver<2, 0>::system_rhs.size(),
                        1e-12 * Solver<2, 0>::system_rhs.l2_norm(),
                        false,
                        false);

  if (log_cg_convergence)
    control.enable_history_data();

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>> cg(control, memory);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(Solver<2, 0>::system_matrix, 1.2);

  cg.solve(Solver<2, 0>::system_matrix,
           Solver<2, 0>::solution,
           Solver<2, 0>::system_rhs,
           preconditioner);

  Solver<2, 0>::constraints.distribute(Solver<2, 0>::solution);

  if (SettingsMMSVTII::log_cg_convergence) {
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

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void
SolverMMSVTII_A::fill_dirichlet_stack()
{
  Solver2<2, 2>::dirichlet_stack = { { bid_dirichlet, &dirichlet_bc } };
}

void
SolverMMSVTII_A::solve()
{
  SolverControl control(1000 * Solver2<2, 2>::system_rhs.size(),
                        1e-12 * Solver2<2, 2>::system_rhs.l2_norm(),
                        false,
                        false);

  if (log_cg_convergence)
    control.enable_history_data();

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>> cg(control, memory);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(Solver2<2, 2>::system_matrix, 1.2);

  cg.solve(Solver2<2, 2>::system_matrix,
           Solver2<2, 2>::solution,
           Solver2<2, 2>::system_rhs,
           preconditioner);

  Solver2<2, 2>::constraints.distribute(Solver2<2, 2>::solution);

  if (SettingsMMSVTII::log_cg_convergence) {
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
