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

#ifndef SolverSSOLIII_H__
#define SolverSSOLIII_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>

#include <string>

#include "exact_solution.hpp"
#include "settings.hpp"
#include "static_vector_solver_i.hpp"
#include "static_vector_solver_ii.hpp"

using namespace StaticVectorSolver;

/**
 * \brief Solves for the current vector potential, \f$\vec{T}\f$, in the
 * [Thick spherical coil with magnetic core (ssol-iii/)](@ref page_ssol_iii)
 * numerical experiment.
 *****************************************************************************/
class SolverSSOLIII_T
  : public SettingsSSOLIII
  , public Solver1<3, 0>
{
public:
  SolverSSOLIII_T() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   *used for mapping. Setting it to 1 will do in the most of the cases. Note,
   *that it makes sense to attach a meaningful manifold to the triangulation if
   *this parameter is greater than 1.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in ssol-iii/gmsh/build. This
   *parameter is used to compose the name of the mesh file to be uploaded from
   * ssol-iii/gmsh/data/.
   * @param[in] fname - The name of the vtk file without extension to save
   * the data.
   *****************************************************************************/
  SolverSSOLIII_T(unsigned int p,
                  unsigned int mapping_degree,
                  unsigned int r,
                  std::string fname)
    : Solver1<3, 0>(p,
                    mapping_degree,
                    3.0,
                    0.0,
                    fname,
                    nullptr,
                    SettingsSSOLIII::print_time_tables,
                    false)
    , r(r)
    , fname(fname)
  {
    Solver1<3, 0>::run();
  }

  ~SolverSSOLIII_T() = default;

private:
  const unsigned int r;
  const std::string fname;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;
};

/**
 * \brief Solves for the magnetic vector potential, \f$\vec{A}\f$, in the
 * [Thick spherical coil with magnetic core (ssol-iii/)](@ref page_ssol_iii)
 * numerical experiment.
 *****************************************************************************/
class SolverSSOLIII_A
  : public SettingsSSOLIII
  , public Solver2<3, 2>
{
public:
  SolverSSOLIII_A() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   *used for mapping. Setting it to 1 will do in the most of the cases. Note,
   *that it makes sense to attach a meaningful manifold to the triangulation if
   *this parameter is greater than 1.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in ssol-iii/gmsh/build. This
   *parameter is used to compose the name of the mesh file to be uploaded from
   * ssol-iii/gmsh/data/.
   * @param[in] triangulation_T - The triangulation created at the 0-th stage of
   * the simulation.
   * @param[in] dof_handler_T - The dof handler created at the 0-th stage of the
   * simulation.
   * @param[in] solution_T - The degrees of freedom that describe the current
   * vector potential computed at the 0-th stage of the simulation.
   * @param[in] fname - The name of the vtk file without extension to save
   * the data.
   *****************************************************************************/
  SolverSSOLIII_A(unsigned int p,
                  unsigned int mapping_degree,
                  unsigned int r,
                  const Triangulation<3>& triangulation_T,
                  const DoFHandler<3>& dof_handler_T,
                  const Vector<double>& solution_T,
                  std::string fname)
    : Solver2<3, 2>(p,
                    mapping_degree,
                    triangulation_T,
                    dof_handler_T,
                    solution_T,
                    0.001 / mu_0,
                    fname,
                    nullptr,
                    SettingsSSOLIII::print_time_tables,
                    false)
    , r(r)
    , fname(fname)
  {
    Solver2<3, 2>::run();
  }

  ~SolverSSOLIII_A() = default;

private:
  const unsigned int r;
  const std::string fname;

  // Solver2 parasites on triangulations of other solvers.
  // No make_mesh() this time around.
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;
};

#endif
