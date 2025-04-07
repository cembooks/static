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

#ifndef SolverSSOLII_H__
#define SolverSSOLII_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <string>

#include "exact_solution.hpp"
#include "settings.hpp"
#include "static_vector_solver_i.hpp"
#include "static_vector_solver_ii.hpp"

using namespace StaticVectorSolver;

/**
 * \brief Solves for the current vector potential, \f$\vec{T}\f$, in the
 * *Thick spherical coil* [(ssol-ii/)](@ref page_ssol_ii) numerical experiment.
 *****************************************************************************/
class SolverSSOLII_T
  : public SettingsSSOLII
  , public Solver1<3, 0>
{
public:
  SolverSSOLII_T() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in <br> ssol-ii/gmsh/build. This
   * parameter is used to compose the name of the mesh file to be uploaded
   * from <br> ssol-ii/gmsh/data/.
   * @param[in] fname - The name of the vtu file without extension to save
   * the data.
   *****************************************************************************/
  SolverSSOLII_T(unsigned int p,
                 unsigned int mapping_degree,
                 unsigned int r,
                 std::string fname)
    : Solver1<3, 0>(p,
                    mapping_degree,
                    2,
                    0.0,
                    fname,
                    nullptr,
                    SettingsSSOLII::print_time_tables,
                    false,
                    true)
    , r(r)
    , fname(fname)
  {
    Solver1<3, 0>::run();
  }

  ~SolverSSOLII_T() = default;

private:
  const unsigned int r;
  const std::string fname;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;

  const SphericalManifold<3> sphere;
  const Functions::ZeroFunction<3> dirichlet_bc{ 3 };
};

/**
 * \brief Solves for the magnetic vector potential, \f$\vec{A}\f$, in the
 * *Thick spherical coil* [(ssol-ii/)](@ref page_ssol_ii) numerical experiment.
 *****************************************************************************/
class SolverSSOLII_A
  : public SettingsSSOLII
  , public Solver2<3, 2>
{
public:
  SolverSSOLII_A() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in <br> ssol-ii/gmsh/build. This
   * parameter is used to compose the name of the mesh file to be uploaded
   * from <br> ssol-ii/gmsh/data/.
   * @param[in] triangulation_T - The triangulation created at the 0-th stage of
   * the simulation.
   * @param[in] dof_handler_T - The dof handler created at the 0-th stage of the
   * simulation.
   * @param[in] solution_T - The degrees of freedom that describe the current
   * vector potential computed at the 0-th stage of the simulation.
   * @param[in] fname - The name of the vtu file without extension to save
   * the data.
   *****************************************************************************/
  SolverSSOLII_A(unsigned int p,
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
                    2,
                    0.0,
                    fname,
                    nullptr,
                    SettingsSSOLII::print_time_tables,
                    false,
                    true)
    , r(r)
    , fname(fname)
  {
    for (auto cell : Solver2<3, 2>::triangulation_T.active_cell_iterators()) {
      for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++) {
        if (cell->face(f)->at_boundary() && (cell->face(f)->boundary_id() == 1))
          cell->face(f)->set_boundary_id(2);
      }
    }

    Solver2<3, 2>::run();
  }

  ~SolverSSOLII_A() = default;

private:
  const unsigned int r;
  const std::string fname;

  // Solver2 reuses triangulations of other solvers.
  // No make_mesh() this time around.
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;
};

#endif
