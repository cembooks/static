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

#ifndef SolverPLS_H__
#define SolverPLS_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/numerics/fe_field_function.h>

#include "exact_solution.hpp"
#include "settings.hpp"
#include "static_scalar_solver.hpp"

#define TMR(__name) TimerOutput::Scope timer_section(timer, __name)

using namespace StaticScalarSolver;

/**
 * \brief Implements the
 * [Planes of symmetry (pls/)](@ref page_pls)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverPLS
  : public SettingsPLS
  , public Solver<dim>
{
public:
  SolverPLS() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the interpolating polynomials of the Lagrange
   * finite elements,
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in pls/gmsh/build. This parameter
   * is used to compose the name of the mesh file to be uploaded from
   * pls/gmsh/data/.
   * @param[in] fname - The name of the vtk file without extension to save
   * the data.
   *****************************************************************************/
  SolverPLS(unsigned int p, unsigned int r, std::string fname)
    : Solver<dim>(p,
                  1,
                  1,
                  fname,
                  &exact_solution,
                  false,
                  false,
                  print_time_tables,
                  project_exact_solution)
    , r(r)
    , fname(fname)
  {
    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Solver run");
      Solver<dim>::run();
    }
  }

  ~SolverPLS() = default;

private:
  const unsigned int r;
  const std::string fname;

  const ExactSolutionPLS_PHI<dim> exact_solution;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;
};

template<int dim>
void
SolverPLS<dim>::solve()
{
  ReductionControl control(
    Solver<dim>::system_rhs.size(), 0.0, 1e-8, false, false);

  if (log_cg_convergence)
    control.enable_history_data();

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>> cg(control, memory);

  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(Solver<dim>::system_matrix, 1.0);

  cg.solve(Solver<dim>::system_matrix,
           Solver<dim>::solution,
           Solver<dim>::system_rhs,
           preconditioner);

  Solver<dim>::constraints.distribute(Solver<dim>::solution);

  if (log_cg_convergence) {
    const std::vector<double> history_data = control.get_history_data();

    std::ofstream ofs(fname + "_cg_convergence.csv");

    unsigned int i = 1;
    for (auto item : history_data) {
      ofs << i << ", " << item << "\n";
      i++;
    }
    ofs.close();
  }
}
#endif
