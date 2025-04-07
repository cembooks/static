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

#ifndef SolverSSOLI_H__
#define SolverSSOLI_H__

#include <deal.II/base/tensor_function.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <string>

#include "exact_solution.hpp"
#include "settings.hpp"
#include "static_vector_solver_i.hpp"

using namespace StaticVectorSolver;

/**
 * \brief Solves for the magnetic vector potential, \f$\vec{A}\f$, in the
 * *Thin spherical coil* [(ssol-i/)](@ref page_ssol_i)
 * numerical experiment.
 *****************************************************************************/
class SolverSSOLI
  : public SettingsSSOLI
  , public Solver1<3>
{
public:
  SolverSSOLI() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the Nedelec finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in <br> ssol-i/gmsh/build. This
   * parameter is used to compose the name of the mesh file to be uploaded
   * from <br> ssol-i/gmsh/data/.
   * @param[in] fname - The name of the vtu file without extension to save
   * the data.
   *****************************************************************************/
  SolverSSOLI(unsigned int p,
              unsigned int mapping_degree,
              unsigned int r,
              std::string fname)
    : Solver1<3>(p,
                 mapping_degree,
                 0,
                 0.0000001 / mu_0,
                 fname,
                 nullptr,
                 SettingsSSOLI::print_time_tables,
                 false,
                 true)
    , r(r)
    , fname(fname)
  {
    Solver1<3>::run();
  }

  ~SolverSSOLI() = default;

private:
  const unsigned int r;
  const std::string fname;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;

  void mark_materials();

  const SphericalManifold<3> sphere;
};

#endif
