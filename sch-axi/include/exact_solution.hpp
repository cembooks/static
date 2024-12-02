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

#ifndef ExactSolutionSCHAXI_H__
#define ExactSolutionSCHAXI_H__

#include "constants.hpp"
#include "settings.hpp"
#include <deal.II/base/function.h>

using namespace dealii;

/**
 * \brief Describes exact solution, \f$\Phi\f$, of the
 * [Axisymmetric - surface charge (sch-axi/)](@ref page_sch_axi)
 * numerical experiment.
 *****************************************************************************/
template<bool is_cylinder>
class ExactSolutionSCHAXI_PHI
  : public Function<2>
  , public SettingsSCHAXI
{
public:
  ExactSolutionSCHAXI_PHI(){};

  virtual double value(const Point<2>& p,
                       const unsigned int component = 0) const override final;

  virtual Tensor<1, 2> gradient(
    const Point<2>& p,
    const unsigned int component = 0) const override final;
};

#endif
