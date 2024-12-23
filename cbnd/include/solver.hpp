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

#ifndef SolverCBND_H__
#define SolverCBND_H__

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
 * \brief Implements the solver of the Effect of curved boundaries
 * [(cbnd/)](@ref page_cbnd) numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverCBND
  : public SettingsCBND
  , public Solver<dim>
{
public:
  SolverCBND() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the interpolating polynomials of the Lagrange
   * finite elements,
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping. Setting it to 1 will do in the most of the cases. Note,
   * that it makes sense to attach a meaningful manifold to the triangulation
   * if this parameter is greater than 1.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in cbnd/gmsh/build. This parameter
   * is used to compose the name of the mesh file to be uploaded from
   * cbnd/gmsh/data/.
   * @param[in] fname - The name of the vtk file without extension to save
   * the data.
   *****************************************************************************/
  SolverCBND(unsigned int p,
             unsigned int mappipng_degree,
             unsigned int r,
             std::string fname)
    : Solver<dim>(p,
                  mappipng_degree,
                  1, // The right-hand side is free-current density.
                  fname,
                  &exact_solution,
                  false, // Is axisymmetric.
                  false, // Is vector potential.
                  print_time_tables,
                  project_exact_solution)
    , r(r)
    , fname(fname)
    , dirichlet_function_in(1.0)
    , fe_slice(1)
  {
    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Solver run");
      Solver<dim>::run();
    }
    {
      TMR("Data slice");
      data_slice(fname);
    }
  }

  unsigned int n_cells;

  ~SolverCBND() = default;

private:
  const unsigned int r;
  const std::string fname;
  const ExactSolutionCBND_PHI<dim> exact_solution;
  const dealii::Functions::ZeroFunction<dim> dirichlet_function_out;
  const dealii::Functions::ConstantFunction<dim> dirichlet_function_in;

  // The amount of global mesh refinements that need to be done to the
  // one- dimensional mesh used for the plot of potential vs. \f$x\f$
  // coordinate.
  const unsigned int nr_slice_global_refs = 10;

  // These four data members are needed for making the plot of
  // potential vs. x coordinate.

  Triangulation<1, dim> triangulation_slice;
  FE_Q<1, dim> fe_slice;
  DoFHandler<1, dim> dof_handler_slice;
  Vector<double> solution_slice;

  //  SphericalManifold<dim> sphere;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;

  // This function makes the plot of potential vs. x coordinate.
  void data_slice(std::string fname);

  void attach_manifold_refine();
};

template<int dim>
void
SolverCBND<dim>::fill_dirichlet_stack()
{
#if IS_BC_EXACT__ == 1
  Solver<dim>::dirichlet_stack = { { bid_in, &exact_solution },
                                   { bid_out, &exact_solution } };
#endif
#if IS_BC_EXACT__ == 0
  Solver<dim>::dirichlet_stack = { { bid_in, &dirichlet_function_in },
                                   { bid_out, &dirichlet_function_out } };
#endif
}

template<int dim>
void
SolverCBND<dim>::data_slice(std::string fname)
{

  if (dim == 2) {
    GridGenerator::hyper_cube(triangulation_slice, a + eps, b - eps);
  } else {
    GridGenerator::hyper_cube(triangulation_slice, a + eps, b - eps);
  }

  triangulation_slice.refine_global(nr_slice_global_refs);

  dof_handler_slice.reinit(triangulation_slice);
  dof_handler_slice.distribute_dofs(fe_slice);
  solution_slice.reinit(dof_handler_slice.n_dofs());

  Functions::FEFieldFunction<dim> potential(Solver<dim>::dof_handler,
                                            Solver<dim>::solution);

  VectorTools::interpolate(dof_handler_slice, potential, solution_slice);

  DataOut<1, dim> data_out;

  data_out.attach_dof_handler(dof_handler_slice);
  data_out.add_data_vector(solution_slice, "solution_slice");
  data_out.build_patches();

  std::ofstream out(fname + "_slice" + ".gpi");

  data_out.write_gnuplot(out);
  out.close();
}

template<int dim>
void
SolverCBND<dim>::attach_manifold_refine()
{
  //  Solver<dim>::triangulation.set_all_manifold_ids(0);
  //  Solver<dim>::triangulation.set_manifold(0,sphere);

  //  for (unsigned int i = 0; i < 3; ++i)
  //  {
  //    for (auto &cell: Solver<dim>::triangulation.active_cell_iterators())
  //     if (cell->at_boundary())
  //       cell->set_refine_flag();

  //  Solver<dim>::triangulation.execute_coarsening_and_refinement();
  //  Solver<dim>::triangulation.set_all_manifold_ids(0);
  //  Solver<dim>::triangulation.set_manifold(0,sphere);
  //  }
}

template<int dim>
void
SolverCBND<dim>::solve()
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
