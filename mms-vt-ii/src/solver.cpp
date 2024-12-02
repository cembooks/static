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
SolverMMSVTII_T::make_mesh()
{
  GridIn<2> gridin;
  Triangulation<2> tria_tmp;

  gridin.attach_triangulation(tria_tmp);
#if DOMAIN__ == 0
  std::ifstream ifs("../../gmsh/data/circle_r" + std::to_string(r) + ".msh");
#elif DOMAIN__ == 1
  std::ifstream ifs("../../gmsh/data/square_r" + std::to_string(r) + ".msh");
#endif
  gridin.read_msh(ifs);

  std::tuple<std::vector<Point<2>>, std::vector<CellData<2>>, SubCellData>
    mesh_description;

  mesh_description = GridTools::get_coarse_mesh_description(tria_tmp);

  GridTools::invert_all_negative_measure_cells(std::get<0>(mesh_description),
                                               std::get<1>(mesh_description));

  GridTools::consistently_order_cells(std::get<1>(mesh_description));

  Solver<2, 0>::triangulation.create_triangulation(
    std::get<0>(mesh_description),
    std::get<1>(mesh_description),
    std::get<2>(mesh_description));

  GridOut gridout;
  GridOutFlags::Msh msh_flags(true, true);
  gridout.set_flags(msh_flags);

#if DOMAIN__ == 0
  std::ofstream ofs("../../gmsh/data/circle_r" + std::to_string(r) +
                    "_reordered.msh");
#elif DOMAIN__ == 1
  std::ofstream ofs("../../gmsh/data/square_r" + std::to_string(r) +
                    "_reordered.msh");
#endif
  gridout.write_msh(Solver<2, 0>::triangulation, ofs);
}

void
SolverMMSVTII_T::fill_dirichlet_stack()
{
  Solver<2, 0>::dirichlet_stack = { { bid_dirichlet, &dirichlet_bc } };
}

void
SolverMMSVTII_T::solve()
{
  ReductionControl control(
    Solver<2, 0>::system_rhs.size(), 0.0, 1e-12, false, false);

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
  ReductionControl control(
    Solver2<2, 2>::system_rhs.size(), 0.0, 1e-12, false, false);

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
