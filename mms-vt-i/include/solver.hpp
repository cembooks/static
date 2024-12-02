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

#ifndef SolverMMSVTI_H__
#define SolverMMSVTI_H__

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
 * \brief Implements the solver for current vector potential, \f$\vec{T}\f$, in
 *the [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref
 *page_mms_vt_i) numerical experiment.
 *****************************************************************************/
class SolverMMSVTI_T
  : public SettingsMMSVTI
  , public Solver1<3, 0>
{
public:
  SolverMMSVTI_T() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   *used for mapping. Setting it to 1 will do in the most of the cases. Note,
   *that it makes sense to attach a meaningful manifold to the triangulation if
   *this parameter is greater than 1.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in mms-vt-i/gmsh/build. This
   *parameter is used to compose the name of the mesh file to be uploaded from
   * mms-vt-i/gmsh/data/.
   * @param[in] fname - The name of the vtk file without extension to save
   * the data.
   *****************************************************************************/
  SolverMMSVTI_T(unsigned int p,
                 unsigned int mapping_degree,
                 unsigned int r,
                 std::string fname)
    : Solver1<3, 0>(p,
                    mapping_degree,
                    3,
                    0.0,
                    fname,
                    nullptr,
                    SettingsMMSVTI::print_time_tables,
                    false)
    , fname(fname)
    , r(r)
  {
    Solver1<3, 0>::run();
  }

  ~SolverMMSVTI_T() = default;

private:
  const std::string fname;
  const unsigned int r;

  DirichletBC_T dirichlet_bc;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;
};

/**
 * \brief Implements the solver for magnetic vector potential, \f$\vec{A}\f$, in
 *the [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref
 *page_mms_vt_i) numerical experiment.
 *****************************************************************************/
class SolverMMSVTI_A
  : public SettingsMMSVTI
  , public Solver2<3, 2>
{
public:
  SolverMMSVTI_A() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   *used for mapping. Setting it to 1 will do in the most of the cases. Note,
   *that it makes sense to attach a meaningful manifold to the triangulation if
   *this parameter is greater than 1.
   * @param[in] triangulation_T - The triangulation created at the 0-th stage of
   * the simulation.
   * @param[in] dof_handler_T - The dof handler created at the 0-th stage of the
   * simulation.
   * @param[in] solution_T - The degrees of freedom that describe the current
   * vector potential computed at the 0-th stage of the simulation.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in mms-vt-i/gmsh/build. This
   *parameter is used to compose the name of the mesh file to be uploaded from
   * mms-vt-i/gmsh/data/.
   * @param[in] fname - The name of the vtk file without extension to save
   * the data.
   *****************************************************************************/
  SolverMMSVTI_A(unsigned int p,
                 unsigned int mapping_degree,
                 const Triangulation<3>& triangulation_T,
                 const DoFHandler<3>& dof_handler_T,
                 const Vector<double>& solution_T,
                 unsigned int r,
                 std::string fname)
    : Solver2<3, 2>(p,
                    mapping_degree,
                    triangulation_T,
                    dof_handler_T,
                    solution_T,
                    0.0,
                    fname,
                    &dirichlet_bc,
                    SettingsMMSVTI::print_time_tables,
                    SettingsMMSVTI::project_exact_solution)
    , fname(fname)
    , r(r)
  {
    Solver2<3, 2>::run();
  }

  ~SolverMMSVTI_A() = default;

private:
  const std::string fname;

  const unsigned int r;

  DirichletBC_A dirichlet_bc;
  // Solver2 parasites on triangulations of other solvers.
  // No make_mesh() this time around.
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;
};

#endif
