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

#include "solver.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

template<>
double
DirichletSLDI<2>::value(const Point<2>& r, unsigned int component) const
{
  return -H_0 * r[0];
}

template<>
double
DirichletSLDI<3>::value(const Point<3>& r, unsigned int component) const
{
  return -H_0 * r[2];
}

#pragma GCC diagnostic pop
