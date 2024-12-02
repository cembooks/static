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
#include "solver.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Axisymmetric - surface charge (sch-axi/)](@ref page_sch_axi)
 * numerical experiment.
 *****************************************************************************/
class BatchSCHAXI : public SettingsSCHAXI
{
public:
  void run()
  {
    if (nr_threads_max > 0)
      MultithreadInfo::set_thread_limit(nr_threads_max);

    std::string dir =
      (CYLINDER__ == 1) ? "Data/cylinder-axi/" : "Data/sphere-axi/";

    unsigned int mapping_degree = (CYLINDER__ == 1) ? 1 : 2;

    std::cout << "Program: sch-axi\n"
              << "Mapping degree = " << mapping_degree << "\n"
              << "Condition: Cylinder = " << CYLINDER__ << "\n"
              << "Writing to: " << dir << "\n";

    MainOutputTable table_PHI(2);

    for (unsigned int p = 1; p < 4; p++) {
      table_PHI.clear();
      for (unsigned int r = 7; r < 11; r++) {
        table_PHI.add_value("r", r);
        table_PHI.add_value("p", p);

        SolverSCHAXI<static_cast<bool>(CYLINDER__)> problem(
          p,
          mapping_degree,
          r,
          dir + "solution_p" + std::to_string(p) + "_r" + std::to_string(r));
        table_PHI.add_value("ndofs", problem.get_n_dofs());
        table_PHI.add_value("ncells", problem.get_n_cells());
        table_PHI.add_value("L2", problem.get_L2_norm());
        table_PHI.add_value("H1", problem.get_H1_norm());
      }
      table_PHI.save(dir + "table_PHI_p" + std::to_string(p));
    }
  }
};

int
main()
{
  try {
    BatchSCHAXI batch;
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
