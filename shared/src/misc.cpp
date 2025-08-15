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

#include "misc.hpp"

#include <iostream>
#include <fstream>

using namespace std;
using namespace Misc;

void
MainOutputTable::format()
{
  set_precision("L2", 2);
  set_precision("H1", 2);

  set_scientific("L2", true);
  set_scientific("H1", true);

  // As in Step-51 line 1385.
  evaluate_convergence_rates(
    "L2", "ncells", ConvergenceTable::reduction_rate_log2, dimensions);
  evaluate_convergence_rates(
    "H1", "ncells", ConvergenceTable::reduction_rate_log2, dimensions);

  set_tex_caption("p", "p");
  set_tex_caption("r", "r");
  set_tex_caption("ncells", "nr. cells");
  set_tex_caption("ndofs", "nr. dofs");
  set_tex_caption("L2", "L2 norm");
  set_tex_caption("H1", "H1 norm");

  set_column_order(new_order);
}

void
MainOutputTable::save(std::string fname)
{
  format();

  std::cout << "--------------------------------------------\n";
  write_text(std::cout);
  std::cout << "\n\n";

  // Save the table in text ...
  {
    std::ofstream ofs(fname + ".txt");
    write_text(ofs);
  }
  // and tex formats.
  {
    std::ofstream ofs(fname + ".tex");
    write_tex(ofs);
  }
}
