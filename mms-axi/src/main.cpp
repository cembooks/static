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

#include "deal.II/base/multithread_info.h"

#include <iostream>
#include <string>

#include "misc.hpp"
#include "project_PHI_to_D.hpp"
#include "project_PHI_to_E.hpp"
#include "solver.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Axisymmetric - method of manufactured solutions (mms-axi/)](@ref
 *page_mms_axi) numerical experiment.
 *****************************************************************************/
class BatchMMSAXI : public SettingsMMSAXI
{
public:
  void run()
  {
    if (nr_threads_max > 0)
      MultithreadInfo::set_thread_limit(nr_threads_max);

    Assert(DIMENSION__ > 1, ExcInternalError());
    Assert(DIMENSION__ < 4, ExcInternalError());

    std::string dir =
      (DIMENSION__ == 2) ? "Data/cylinder2d/" : "Data/cylinder3d/";
    std::string fname;

    std::cout << "Program: mms-axi\n"
              << "Dimensions: " << DIMENSION__ << "\n"
              << "Writing to: " << dir << "\n";

    MainOutputTable table_PHI(DIMENSION__);
    MainOutputTable table_E(DIMENSION__);
    MainOutputTable table_D(DIMENSION__);

    bool axisymmetric;
    for (unsigned int p = 1; p < 4; p++) {
      table_PHI.clear();
      table_E.clear();
      table_D.clear();

#if DIMENSION__ == 2
      axisymmetric = true;
      for (unsigned int r = 27; r < 31; r++)
#endif
#if DIMENSION__ == 3
        axisymmetric = false;
      for (unsigned int r = 5; r < 9; r++)
#endif
      {
        fname =
          dir + "solution_PHI_p" + std::to_string(p) + "_r" + std::to_string(r);

        table_PHI.add_value("r", r);
        table_PHI.add_value("p", p);

        SolverMMSAXI<DIMENSION__> problem(p, r, fname);

        table_PHI.add_value("ndofs", problem.get_n_dofs());
        table_PHI.add_value("ncells", problem.get_n_cells());
        table_PHI.add_value("L2", problem.get_L2_norm());
        table_PHI.add_value("H1", problem.get_H1_norm());

        problem.clear();

        {
          fname = dir + "solution_E" + "_p" + std::to_string(p) + "_r" +
                  std::to_string(r);

          table_E.add_value("r", r);
          table_E.add_value("p", p);

          ExactSolutionMMSAXI_E<DIMENSION__> exact_solution;

          ProjectPHItoE<DIMENSION__> projector(p - 1,
                                               1,
                                               problem.get_tria(),
                                               problem.get_dof_handler(),
                                               problem.get_solution(),
                                               fname,
                                               &exact_solution,
                                               axisymmetric,
                                               Settings::print_time_tables,
                                               Settings::project_exact_solution,
                                               Settings::log_cg_convergence);

          table_E.add_value("ndofs", problem.get_n_dofs());
          table_E.add_value("ncells", problem.get_n_cells());
          table_E.add_value("L2", projector.get_L2_norm());
          table_E.add_value("H1", 0.0);
        }
        {
          fname = dir + "solution_D" + "_p" + std::to_string(p) + "_r" +
                  std::to_string(r);

          table_D.add_value("r", r);
          table_D.add_value("p", p);

          ExactSolutionMMSAXI_D<DIMENSION__> exact_solution;

          ProjectPHItoD<DIMENSION__> projector(p - 1,
                                               1,
                                               problem.get_tria(),
                                               problem.get_dof_handler(),
                                               problem.get_solution(),
                                               fname,
                                               &exact_solution,
                                               axisymmetric,
                                               Settings::print_time_tables,
                                               Settings::project_exact_solution,
                                               Settings::log_cg_convergence);
          table_D.add_value("ndofs", problem.get_n_dofs());
          table_D.add_value("ncells", problem.get_n_cells());
          table_D.add_value("L2", projector.get_L2_norm() / ep_0);
          table_D.add_value("H1", 0.0);
        }
      }

      std::cout << "Table PHI\n";
      table_PHI.save(dir + "table_PHI_p" + std::to_string(p));

      std::cout << "Table E\n";
      table_E.save(dir + "table_E_p" + std::to_string(p));

      std::cout << "Table D\n";
      table_D.save(dir + "table_D_p" + std::to_string(p));
    }
  }
};

int
main()
{
  try {
    BatchMMSAXI batch;
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
