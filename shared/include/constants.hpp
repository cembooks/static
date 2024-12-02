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

#ifndef Constants_H__
#define Constants_H__

#include <deal.II/base/exceptions.h>

#include <iostream>
#include <math.h>
#include <vector>

namespace Constants {

/**
 * \brief Lists physical constants.
 *
 *******************************************************************************/
class Physics
{
public:
  Physics(){};

  /**
   * \brief The ratio between the circumference and the diameter of any circle,
   * \f$\pi\f$.
   *****************************************************************************/
  const double pi =
    3.141592653589793238462643383279502884197169399375105820974944592307816406286;

  /**
   * \brief The speed of light in free space, \f$c\f$.
   *****************************************************************************/
  const double c = 299792458.0;

  /**
   * \brief The permeability of the free space, \f$\mu_0\f$.
   ****************************************************************************/
  const double permeability_fs = 4.0 * pi * 1.0e-7;

  /**
   * \brief The permittivity of the free space, \f$\epsilon_0\f$.
   * *************************************************************************/
  const double permittivity_fs = 1.0 / (std::pow(c, 2) * permeability_fs);
};

/**
 * \brief The tables that contain the amount of quadrature points used in the
 * scalar problems.
 *******************************************************************************/
template<int dim>
class QuadratureTableScalar
{
public:
  QuadratureTableScalar(unsigned int p)
    : indx(p - 1)
  {
    Assert(indx < 4, dealii::ExcInternalError());
  }

  QuadratureTableScalar(const QuadratureTableScalar<dim>& qt)
    : indx(qt.indx)
  {
    Assert(indx < 4, dealii::ExcInternalError());
  }

  /**
   * \brief Returns the amount of quadrature points used when assembling
   * system if linear equations.
   *******************************************************************************/
  unsigned int sim() const;

  /**
   * \brief Returns the amount of quadrature points used when calculating the
   * error norms.
   *******************************************************************************/
  unsigned int enorm() const;

private:
  const unsigned int indx;
  const std::vector<unsigned int> Table2DSimulation = { 2, 3, 4, 5 };
  const std::vector<unsigned int> Table3DSimulation = { 2, 3, 4, 5 };
  const std::vector<unsigned int> Table2DErrorNorm = { 4, 5, 6, 7 };
  const std::vector<unsigned int> Table3DErrorNorm = { 4, 5, 6, 7 };
};

/**
 * \brief The tables that contain the amount of quadrature points used in
 * vector problems.
 *******************************************************************************/
template<int dim>
class QuadratureTableVector
{
public:
  QuadratureTableVector(unsigned int p)
    : indx(p)
  {
    Assert(indx < 3, dealii::ExcInternalError());
  }

  QuadratureTableVector(const QuadratureTableVector<dim>& qt)
    : indx(qt.indx)
  {
    Assert(indx < 3, dealii::ExcInternalError());
  }

  /**
   * \brief Returns the amount of quadrature points used when assembling
   * system if linear equations.
   *******************************************************************************/
  unsigned int sim() const;

  /**
   * \brief Returns the amount of quadrature points used when calculating the
   * error norms.
   *******************************************************************************/
  unsigned int enorm() const;

private:
  const unsigned int indx;
  const std::vector<unsigned int> Table2DSimulation = { 2, 3, 4, 5 };
  const std::vector<unsigned int> Table3DSimulation = { 2, 3, 4, 5 };
  const std::vector<unsigned int> Table2DErrorNorm = { 4, 5, 6, 7 };
  const std::vector<unsigned int> Table3DErrorNorm = { 4, 5, 6, 7 };
};

} // namespace Constants
#endif
