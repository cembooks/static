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

#include "constants.hpp"
using namespace Constants;

template<>
unsigned int
QuadratureTableScalar<2>::sim() const
{
  return Table2DSimulation.at(indx);
}

template<>
unsigned int
QuadratureTableScalar<3>::sim() const
{
  return Table3DSimulation.at(indx);
}

template<>
unsigned int
QuadratureTableScalar<2>::enorm() const
{
  return Table2DErrorNorm.at(indx);
}

template<>
unsigned int
QuadratureTableScalar<3>::enorm() const
{
  return Table3DErrorNorm.at(indx);
}

template<>
unsigned int
QuadratureTableVector<2>::sim() const
{
  return Table2DSimulation.at(indx);
}

template<>
unsigned int
QuadratureTableVector<3>::sim() const
{
  return Table3DSimulation.at(indx);
}

template<>
unsigned int
QuadratureTableVector<2>::enorm() const
{
  return Table2DErrorNorm.at(indx);
}

template<>
unsigned int
QuadratureTableVector<3>::enorm() const
{
  return Table3DErrorNorm.at(indx);
}
