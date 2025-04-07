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

#ifndef SolverFLCAXI_H__
#define SolverFLCAXI_H__

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
 * \brief Implements the solver of the
 * *Axisymmetric - floating conductor* [(flc-axi/)](@ref page_flc_axi)
 * numerical experiment.
 *****************************************************************************/
template<bool is_cylinder>
class SolverFLCAXI
  : public SettingsFLCAXI
  , public Solver<2>
{
public:
  SolverFLCAXI() = delete;

  /**
   * The constructor.
   *
   * @param[in] p - The degree of the interpolating polynomials of the Lagrange
   * finite elements,
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping.
   * @param[in] r - The parameter that encodes the degree of mesh refinement.
   * Must coincide with one of the values set in flc-axi/gmsh/build. This
   * parameter is used to compose the name of the mesh file to be uploaded
   * from flc-axi/gmsh/data/.
   * @param[in] fname - The name of the vtu file without extension to save
   * the data.
   *****************************************************************************/
  SolverFLCAXI(unsigned int p,
               unsigned int mapping_degree,
               unsigned int r,
               std::string fname)
    : Solver<2>(p,
                mapping_degree,
                0,
                fname,
                &exact_solution,
                true,
                false,
                SettingsFLCAXI::print_time_tables,
                SettingsFLCAXI::project_exact_solution,
                true)
    , r(r)
    , fname(fname)
    , dirichlet_function_in(1.0)
    , fe_slice(1)
  {
    TimerOutput::OutputFrequency tf = (SettingsFLCAXI::print_time_tables)
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

  ~SolverFLCAXI() = default;

private:
  const unsigned int r;
  const std::string fname;

  const ExactSolutionFLCAXI_PHI<is_cylinder> exact_solution;
  const dealii::Functions::ZeroFunction<2> dirichlet_function_out;
  const dealii::Functions::ConstantFunction<2> dirichlet_function_in;

  // The amount of global mesh refinements that need to be done to the
  // one- dimensional mesh used for the plot of potential vs. x coordinate.
  const unsigned int nr_slice_global_refs = 3;

  // These four data members are needed for making the plot of potential
  // vs. x coordinate.
  Triangulation<1, 2> triangulation_slice;
  FE_Q<1, 2> fe_slice;
  DoFHandler<1, 2> dof_handler_slice;
  Vector<double> solution_slice;

  SphericalManifold<2> sphere;

  virtual void make_mesh() override final;
  virtual void fill_dirichlet_stack() override final;
  virtual void solve() override final;

  // This function makes the plot of potential vs. x coordinate.
  void data_slice(std::string fname);
};

template<bool is_cylinder>
void
SolverFLCAXI<is_cylinder>::fill_dirichlet_stack()
{
  Solver<2>::dirichlet_stack = { { bid_in, &dirichlet_function_in },
                                 { bid_out, &dirichlet_function_out } };
}

template<bool is_cylinder>
void
SolverFLCAXI<is_cylinder>::data_slice(std::string fname)
{
  GridGenerator::hyper_cube(triangulation_slice, a + eps, b - eps);
  triangulation_slice.refine_global(nr_slice_global_refs);

  dof_handler_slice.reinit(triangulation_slice);
  dof_handler_slice.distribute_dofs(fe_slice);
  solution_slice.reinit(dof_handler_slice.n_dofs());

  Functions::FEFieldFunction<2> potential(Solver<2>::dof_handler,
                                          Solver<2>::solution);

  VectorTools::interpolate(dof_handler_slice, potential, solution_slice);

  //  DataOut<1,DoFHandler<1,2>> data_out;
  DataOut<1, 2> data_out;

  data_out.attach_dof_handler(dof_handler_slice);
  data_out.add_data_vector(solution_slice, "solution_slice");
  data_out.build_patches();

  std::ofstream out(fname + "_slice" + ".gpi");

  data_out.write_gnuplot(out);
  out.close();
}

template<bool is_cylinder>
void
SolverFLCAXI<is_cylinder>::solve()
{
  SolverControl control(Solver<2>::system_rhs.size(),
                        1e-8 * Solver<2>::system_rhs.l2_norm(),
                        false,
                        false);

  if (log_cg_convergence)
    control.enable_history_data();

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>> cg(control, memory);

  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(Solver<2>::system_matrix, 1.0);

  cg.solve(Solver<2>::system_matrix,
           Solver<2>::solution,
           Solver<2>::system_rhs,
           preconditioner);

  Solver<2>::constraints.distribute(Solver<2>::solution);

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
