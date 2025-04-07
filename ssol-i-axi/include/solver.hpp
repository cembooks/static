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

#ifndef SolverSSOLIAXI_H__
#define SolverSSOLIAXI_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/numerics/fe_field_function.h>

#include "exact_solution.hpp"
#include "settings.hpp"
#include "static_scalar_solver.hpp"

#define TMR(__name) TimerOutput::Scope timer_section(timer, __name)

using namespace StaticScalarSolver;

/**
 * \brief Implements the
 * *Axisymmetric - thin spherical coil* [(ssol-i-axi)](@ref page_ssol_i_axi)
 * numerical experiment.
 *****************************************************************************/
class SolverSSOLIAXI
  : public SettingsSSOLIAXI
  , public Solver<2>
{
public:
  SolverSSOLIAXI() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the interpolating polynomials of the Lagrange
   * finite elements,
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in ssol-i-axi/gmsh/build. This
   * parameter is used to compose the name of the mesh file to be uploaded
   * from ssol-i-axi/gmsh/data/.
   * @param[in] fname - The name of the vtu file without extension to save
   * the data.
   *****************************************************************************/
  SolverSSOLIAXI(unsigned int p,
                 unsigned int mapping_degree,
                 unsigned int r,
                 std::string fname)
    : Solver<2>(p,
                mapping_degree,
                0,
                fname,
                &exact_solution,
                true,
                true,
                SettingsSSOLIAXI::print_time_tables,
                SettingsSSOLIAXI::project_exact_solution,
                true)
    , r(r)
    , fname(fname)
    , fe_slice(1)
  {
    TimerOutput::OutputFrequency tf = (SettingsSSOLIAXI::print_time_tables)
                                        ? TimerOutput::summary
                                        : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Solver run");
      Solver<2>::run();
    }
    {
      TMR("Data slice");
      data_slice(fname);
    }
  }

  ~SolverSSOLIAXI() = default;

private:
  const unsigned int r;
  const std::string fname;

  const ExactSolutionSSOLIAXI_A exact_solution;
  const dealii::Functions::ZeroFunction<2> dirichlet_function;

  // The amount of global mesh refinements that need to be done to the
  // one-dimensional mesh used for the plot of potential vs. x coordinate.
  const unsigned int nr_slice_global_refs = 10;

  // These four data members are needed for making the plot of potential
  // vs. \f$x\f$ coordinate.
  Triangulation<1, 2> triangulation_slice;
  FE_Q<1, 2> fe_slice;
  DoFHandler<1, 2> dof_handler_slice;
  Vector<double> solution_slice;

  SphericalManifold<2> sphere;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;

  void mark_materials();

  // This function makes the plot of potential vs. \f$x\f$ coordinate.
  void data_slice(std::string fname);
};

#endif
