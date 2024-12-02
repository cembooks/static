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
#include <deal.II/base/timer.h>

#include <iostream>
#include <string>

#include "misc.hpp"
#include "project_A_to_B.hpp"
#include "project_Axy_to_Bz.hpp"
#include "settings.hpp"
#include "solver.hpp"

using namespace Misc;
using namespace StaticVectorSolver;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Method of manufactured solutions, vector potential (mms-v/)](@ref
 *page_mms_v) numerical experiment.
 *****************************************************************************/
class BatchMMSV : public SettingsMMSV
{
public:
  void run()
  {
    if (nr_threads_max > 0)
      MultithreadInfo::set_thread_limit(nr_threads_max);

    Assert(DIMENSION__ < 4, ExcInternalError());
    Assert(DIMENSION__ > 1, ExcInternalError());

    std::string dir;
    std::string fname;

    if (HYPERCUBE__ == 1) {
      dir = (DIMENSION__ == 2) ? "Data/square/" : "Data/cube/";
    } else {
      dir = (DIMENSION__ == 2) ? "Data/circle/" : "Data/sphere/";
    }

    std::cout << "Program: mms-v\n"
              << "Dimensions: " << DIMENSION__ << "\n"
              << "Writing to: " << dir << "\n";

    MainOutputTable table(DIMENSION__);

    for (unsigned int p = 0; p < 3; p++) {
      table.clear();
      for (unsigned int r = 5; r < 9; r++) {
        std::cout << "Solving for A ...\n";

        fname =
          dir + "solution_A_p" + std::to_string(p) + "_r" + std::to_string(r);

        table.add_value("r", r);
        table.add_value("p", p);

        SolverMMSV<DIMENSION__> problem(p, r, fname);

        std::cout << "Projecting A in H(curl) to B in H(div) ...\n";

        ExactSolutionMMSV_B<DIMENSION__> exact_solution;

#if DIMENSION__ == 3
        fname =
          dir + "solution_B_p" + std::to_string(p) + "_r" + std::to_string(r);

        ProjectAtoB projector(p,
                              1,
                              problem.get_tria(),
                              problem.get_dof_handler(),
                              problem.get_solution(),
                              fname,
                              &exact_solution,
                              Settings::print_time_tables,
                              Settings::project_exact_solution,
                              Settings::log_cg_convergence);
#endif
#if DIMENSION__ == 2
        fname =
          dir + "solution_B_p" + std::to_string(p) + "_r" + std::to_string(r);

        ProjectAxyToBz projector(p,
                                 1,
                                 problem.get_tria(),
                                 problem.get_dof_handler(),
                                 problem.get_solution(),
                                 fname,
                                 &exact_solution,
                                 Settings::print_time_tables,
                                 Settings::project_exact_solution,
                                 Settings::log_cg_convergence);
#endif
        table.add_value("ndofs", problem.get_n_dofs());
        table.add_value("ncells", problem.get_n_cells());
        table.add_value("L2", projector.get_L2_norm());
        table.add_value("H1", 0.0);
      }
      table.save(dir + "main_table_p" + std::to_string(p));
    }
  }
};

int
main()
{
  try {
    BatchMMSV batch;
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
