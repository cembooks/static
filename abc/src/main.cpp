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

#include <deal.II/base/multithread_info.h>

#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"

using namespace Misc;

/**
 * \brief An extended version of the convergence table used in
 * *Asymptotic boundary condition* [(abc/)](@ref page_abc)
 * numerical experiment.
 *****************************************************************************/
class MainOutputTableABC : public MainOutputTable
{
public:
  MainOutputTableABC() = delete;
  MainOutputTableABC(int dim)
    : MainOutputTable(dim)
  {
    std::vector<std::string> order = { "d",      "b",     "m",  "p", "r",
                                       "ncells", "ndofs", "L2", "H1" };

    set_new_order(order);
  }

  virtual void format() override
  {
    MainOutputTable::format();

    set_tex_caption("d", "d");
    set_tex_caption("b", "b");
    set_tex_caption("m", "m");
  }
};

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * *Asymptotic boundary condition* [(abc/)](@ref page_abc)
 * numerical experiment.
 *****************************************************************************/
class BatchABC : public SettingsABC
{
public:
  void run()
  {
    if (nr_threads_max > 0)
      MultithreadInfo::set_thread_limit(nr_threads_max);

    Assert(DIMENSION__ > 1, ExcInternalError());
    Assert(DIMENSION__ < 4, ExcInternalError());

    const std::vector<std::string> bc_string = {
      "neumann",   // 0
      "dirichlet", // 1
      "abc"        // 2
    };

    std::string dir;

    dir = (DIMENSION__ == 2) ? "Data/ppc-" + bc_string[BC_TYPE__] + "/"
                             : "Data/shell-" + bc_string[BC_TYPE__] + "/";

    std::cout << "Program: abc\n"
              << "Dimensions: " << DIMENSION__ << "\n"
              << "Boundary condition: " << BC_TYPE__ << "\n"
              << "Writing to: " << dir << "\n";

    MainOutputTableABC table_PHI(DIMENSION__);
    for (unsigned int m = 1; m < 6; m++) {
      std::cout << "--------------------------------------------\n"
                << "--------------------------------------------\n";

      for (unsigned int p = 1; p < 4; p++) {
        table_PHI.clear();
        for (unsigned int r = 1; r < 4; r++) {
          table_PHI.add_value("d", DIMENSION__);
          table_PHI.add_value("b", BC_TYPE__);
          table_PHI.add_value("m", m);
          table_PHI.add_value("r", r);
          table_PHI.add_value("p", p);

          SolverABC<DIMENSION__> problem(
            m,
            p,
            2,
            r,
            dir + "solution_p" + std::to_string(p) + "_m" + std::to_string(m) +
              "_r" + std::to_string(r));
          problem.get_L2_norm();
          table_PHI.add_value("ndofs", problem.get_n_dofs());
          table_PHI.add_value("ncells", problem.get_n_cells());
          table_PHI.add_value("L2", problem.get_L2_norm());
          table_PHI.add_value("H1", problem.get_H1_norm());
        }
        table_PHI.save(dir + "table_PHI_m" + std::to_string(m) + "_p" +
                       std::to_string(p));
      }
    }
  }
};

int
main()
{
  try {
    BatchABC batch;
    batch.run();
  } catch (std::exception& exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  } catch (...) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
