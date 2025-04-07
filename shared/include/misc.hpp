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

#ifndef MISC8457956__
#define MISC8457956__

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/tensor.h>

namespace Misc {
using namespace dealii;

/**
 * \brief The convergence table used in multiple numerical experiments.
 *****************************************************************************/
class MainOutputTable : public ConvergenceTable
{
public:
  MainOutputTable() = delete;
  /**
   * \brief The only constructor.
   *****************************************************************************/
  MainOutputTable(int dimensions)
    : ConvergenceTable()
    , dimensions(dimensions)
  {
  }
  /**
   * \brief Saves the data in text and tex formats, and prints the data on
   * screen.
   *
   * @param[in] fname - name of the file without extension; includes full path
   * of the file.
   *****************************************************************************/
  void save(std::string fname);
  /**
   * \brief Sets a new order of columns.
   *
   * @param[in] new_order_in - The new order of columns. The default order of
   * columns is: "p", "r", "ncells", "ndofs", "L2", "H1".
   *****************************************************************************/
  void set_new_order(std::vector<std::string> new_order_in)
  {
    new_order = new_order_in;
  }

  /**
   * \brief Appends a column to the table.
   *
   * @param[in] new_column - The name of the new column to be appended.
   *****************************************************************************/
  void append_new_order(std::string new_column)
  {
    new_order.push_back(new_column);
  }

  virtual void format();

  virtual ~MainOutputTable() = default;

private:
  int dimensions;
  std::vector<std::string> new_order = {
    "p", "r", "ncells", "ndofs", "L2", "H1"
  };
};

} // namespace Misc

#endif
