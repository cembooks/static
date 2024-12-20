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

#ifndef StaticVectorSolverII_H__
#define StaticVectorSolverII_H__

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/numerics/vector_tools_project.h>

#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>

#include "constants.hpp"
#include "settings.hpp"
#include "static_vector_input.hpp"

#define VE scratch_data.ve
#define SE scratch_data.se

#define TMR(__name) TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

namespace StaticVectorSolver {
/**
 * \brief Solves
 * [static vector boundary value problem](@ref page_veibvp).
 *
 * Implements the following recipes:
 * - (1) Recipe for static vector solver in 3D
 * - (2) Recipe for static vector solver in 2D
 *
 * This class template is intended to be a general solver for problems in
 * magnetostatics that can be formulated in terms of the magnetic vector
 * potential, \f$\vec{A}\f$. The Bossavit's diagrams below illustrate the
 * partial differential equations that can be solved with a help of this
 * class template.
 *
 * ![](svst_svst/diagram_svst.svg)
 *
 * This class template is very similar to StaticVectorSolver::Solver1. The main
 * difference between StaticVectorSolver::Solver1 and
 * StaticVectorSolver::Solver2 can be described as the following. All inputs
 * to StaticVectorSolver::Solver1 must be given in a form of an analytical
 * expression so they can be codded in class templates derived from
 * [Function](https://www.dealii.org/current/doxygen/deal.II/classFunction.html#a8c6e33c27ac2c3c2be40af1f954b71a7)
 * class template of deal.II. The same holds for StaticVectorSolver::Solver2
 * with an exception of the input that describes the right-hand side of the
 * partial differential equation. The input on the right-hand side of the
 * partial differential equation fed into the StaticVectorSolver::Solver2 class
 * template must be a current vector potential, \f$\vec{T}\f$, expressed in a
 * form of a field function. That is, \f$\vec{T}\f$ fed into
 * StaticVectorSolver::Solver2 is a result of another deal.II simulation. This
 * allows to solve the curl-curl equation for \f$\vec{A}\f$ observing the
 * compatibility condition. The diagram below illustrates how this can be done
 * in two- and three- dimensions.
 *
 * @anchor diagram_svst2
 * ![](svst_svst/diagram_svst2.svg)
 *
 * The magnetic vector potential is modeled by the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements.
 *
 * A user of this class is supposed to do the following.
 *
 * - Derive a class from StaticVectorSolver::Solver2.
 *
 * - Override the virtual destructor.
 *
 * - Call constructor StaticVectorSolver::Solver2.
 *
 * - Override member functions
 *   @code
 *   virtual void fill_dirichlet_stack() = 0;
 *   virtual void solve() = 0;
 *   @endcode
 *
 * - Implement the member functions
 *   @code
 *   void value_list(...);
 *   @endcode
 *   of the classes StaticVectorSolver::TheCoefficient,
 *   StaticVectorSolver::Gamma, StaticVectorSolver::RobinRhs, and
 *   StaticVectorSolver::FreeSurfaceCurrent.
 *
 * - Implement member function
 *   @code
 *   double value(...) const;
 *   @endcode
 *   of the class StaticVectorSolver::Weight.
 *
 * - Explicitly call
 *   @code
 *   StaticVectorSolver::run();
 *   @endcode
 *   Alternatively, the user may call individual member functions (such as
 *   fill_dirichlet_stack(), setup(), etc.) in a proper order.
 *
 * According to the deal.II documentation of
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements, several aspects of the implementation of the Nedelec
 * elements are still experimental. At this moment only globally refined meshes
 * with consistent orientation of faces are allowed. This class utilizes
 * FE_Nedelec finite elements. Consequently, all restrictions applied to the
 * FE_Nedelec finite elements apply to this class template.
 *
 * The boundaries of the mesh must be labeled such that the
 * <code>boundary_id()</code> member function of a face object returns the
 * corresponding boundary ID. The boundary ID's must obey the following
 * convention.
 *
 *  @anchor veis_bnd_convention
 * - The Dirichlet boundary conditions are applied on the boundaries with odd
 *   boundary ID's.
 * - The Robin boundary conditions are applied on the boundaries with even
 *   boundary ID's. The boundary ID's in this case must be greater than zero.
 * - No boundary condition is applied on a boundary with zero ID. Applying no
 *   boundary condition is as good as applying the homogeneous Neumann boundary
 *   condition,
 *   \f$\big(1/\mu\big)\hat{n}\times\big(\vec{\nabla}\times\vec{A}\big)=0\f$,
 *   as implicitly implied by the first term of the
 *   [functional](@ref veibvp_functional).
 *   This boundary condition can also be imposed
 *   by assigning to a boundary an even ID greater than zero, and
 *   setting \f$ \gamma \f$ and \f$\vec{Q}\f$ to zero in the
 *   <code>value_list</code> methods of the classes
 *   StaticVectorSolver::Gamma and StaticVectorSolver::RobinRhs.
 *
 * It is assumed that the object that has been used to calculate \f$\vec{T}\f$
 * (or \f$T\f$ in 2D) is still in the computer memory such that the
 * triangulation, the degrees of freedom, and the handler of the degrees of
 * freedom are accessible while \f$\vec{A}\f$ is computed. An object of
 * the StaticVectorSolver::Solver2 class template reuses the triangulation by
 * creating an additional dof handler associated with the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. That is, \f$\vec{T}\f$ and \f$\vec{A}\f$ share the same
 * triangulation. Two separate DoFHandler objects are used for
 * \f$\vec{T}\f$ and \f$\vec{A}\f$.
 *
 * The algorithm walks synchronously through both DoFHandler objects and
 * assembles the system matrix and the right-hand side.
 * The algorithm utilizes the
 * [WorkStream](https://www.dealii.org/current/doxygen/deal.II/namespaceWorkStream.html)
 * technology of deal.II. The amount of threads used can be limited as the
 * following.
 * @code
 * #include <deal.II/base/multithread_info.h>
 * ...
 * MultithreadInfo::set_thread_limit(nr_threads_max);
 * @endcode
 *
 * @note Application examples:
 * - [mms-vt-i/](@ref page_mms_vt_i),
 * [mms-vt-ii/](@ref page_mms_vt_ii),
 * [ssol-ii/](@ref page_ssol_ii),
 * [ssol-iii/](@ref page_ssol_iii).
 *****************************************************************************/
template<int dim, int stage = 1>
class Solver2
{
public:
  Solver2() = delete;

  /**
   * \brief The only constructor.
   *
   * @param[in] p - Degree of the
   * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
   * finite elements.
   * @param[in] mapping_degree - The degree of the interpolating Lagrange
   * polynomials used for mapping. Setting it to 1 will do in the most of the
   * cases. Note, that it makes sense to attach a meaningful manifold to the
   * triangulation if this parameter is greater than 1.
   * @param[in] triangulation_T - A reference to the mesh on which the current
   * vector potential, \f$\vec{T}\f$, has been computed.
   * @param[in] dof_handler_T - A reference to the DoFHandler object that
   * describes the current vector potential, \f$\vec{T}\f$.
   * @param[in] solution_T - A reference to the degrees of freedom that describe
   * the current vector potential, \f$\vec{T}\f$.
   * @param[in] eta_squared - The gauging parameter, \f$\eta^2\f$, in the
   * [partial differential equation](@ref page_veibvp).
   * @param[in] fname - The name of the output files without extension. Names of
   * the output files will be generated by appending simulation conditions to
   * this string.
   * @param[in] exact_solution - Points to an object that describes the exact
   * solution to the problem. It is needed for calculating error norms. It is a
   * responsibility of the user to make sure that the object exists at the time
   * of the execution of run() or compute_error_norms().
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the Nedelec finite elements and saves
   * the result into the vtk file next to the solution. This may be useful for
   * debugging purposes as a comparison between the projected exact solution and
   * the solution to the boundary value problem can yield a hint on where to
   * search for bugs.
   * @param[in] write_higher_order_cells - Switches between the two modes of
   * operation of the save() function, see the description of save().
   *****************************************************************************/
  Solver2(unsigned int p,
          unsigned int mapping_degree,
          const Triangulation<dim>& triangulation_T,
          const DoFHandler<dim>& dof_handler_T,
          const Vector<double>& solution_T,
          double eta_squared = 0.0,
          std::string fname = "data",
          const Function<dim>* exact_solution = nullptr,
          bool print_time_tables = false,
          bool project_exact_solution = false,
          bool write_higher_order_cells = false)
    : triangulation_T(triangulation_T)
    , dof_handler_T(dof_handler_T)
    , solution_T(solution_T)
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // The public attribute fe.degree is the maximal polynomial degree of a
    // shape function in a single coordinate direction, not the degree of
    // the finite element. For FE_Nedelec and FE_RaviartThomas degree of
    // the finite element is: degree_of_element = fe.degree - 1.
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //        fe(dof_handler_T.get_fe().degree-1),
    , fe(p)
    , mapping_degree(mapping_degree)
    , eta_squared(eta_squared)
    , fname(fname)
    , exact_solution(exact_solution)
    , print_time_tables(print_time_tables)
    , project_exact_solution(project_exact_solution)
    , write_higher_order_cells(write_higher_order_cells)
  {
  }

  /**
   * \brief Initializes the data member
   * StaticVectorSolver::Solver2::dirichlet_stack.
   *
   * This function must be overridden by the user. It must initialize the stack
   * of the Dirichlet boundary conditions. For example:
   * @code
   *
   * using namespace dealii;
   *
   * const types::boundary_id boundary_id_1 = 1;
   * const types::boundary_id boundary_id_2 = 3;
   *
   * const DirichletFunction1<dim> dirichlet_bc_1;
   * const DirichletFunction2<dim> dirichlet_bc_2;
   *
   * template<>
   * void SolverMyProblem<dim>::fill_dirichlet_stack()
   * {
   *   Solver<dim, stage>::dirichlet_stack =
   *     {{boundary_id_1, & dirichlet_bc_1},
   *      {boundary_id_2, & dirichlet_bc_2}};
   * }
   * @endcode
   *
   * The boundary IDs must be odd numbers, see above the
   * [convention](@ref veis_bnd_convention)
   * on the boundary IDs.
   *****************************************************************************/
  virtual void fill_dirichlet_stack() = 0;

  /**
   * \brief Solves the system of linear equations.
   *****************************************************************************/
  virtual void solve() = 0;

  /**
   * \brief Initializes system matrix and the right-hand side vector.
   *
   * Initialises
   * StaticVectorSolver::Solver2::system_matrix,
   * StaticVectorSolver::Solver2::system_rhs and other arrays.
   * Applies the Dirichlet boundary conditions by constraining the system
   * matrix. Distributes degrees of freedom.
   *****************************************************************************/
  void setup();

  /**
   * \brief Assembles the system matrix and the right-hand side vector.
   *****************************************************************************/
  void assemble();

  /**
   * \brief Computes error norms.
   *****************************************************************************/
  void compute_error_norms();

  /**
   * \brief Projects exact solution.
   *
   * The mesh and the finite elements, are the same as are used for the
   * numerical solution of the boundary vale problem. The exact solution will be
   * saved in the vtk file next to the numerical solution to the boundary value
   * problem. This function works properly only if the exact solution is
   * submitted to the constructor via the input parameter
   * <code>exact_solution</code> and <code>project_exact_solution=true</code>.
   *****************************************************************************/
  void project_exact_solution_fcn();

  /**
   * \brief Saves simulation results into a vtk file.
   *
   * The following data are saved:
   * - The calculated potential under the name "VectorField".
   * - The \f$L^2\f$ error norm associated with the calculated potential under
   *   the name "L2norm". One value per mesh cell is saved.
   * - The exact solution expressed as a linear combination of the shape
   *   functions of the
   *   [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
   *   finite elements is saved under the name "VectorFieldExact". The
   *   "VectorField" and "VectorFieldExact" are modeled by exactly the same
   *   finite elements.
   *
   * The "L2norm" and "VectorFieldExact" are saved only if an exact solution
   * is submitted to the constructor. Moreover, "VectorFieldExact" is
   * calculated and saved only if <code>project_exact_solution = true</code>.
   *
   * If <code>write_higher_order_cells = false</code>, the name of the file is
   * computed by appending ".vtk" to the string contained by the parameter
   * <code>fname</code> passed to the constructor. The vtk file can be
   * inspected with a help of [Visit](https://visit.llnl.gov)
   * or [Paraview](www.paraview.org). Higher-order cells are not saved.
   * If <code>write_higher_order_cells = true</code>, the data is saved into
   * fname.vtu file preserving the higher-order cells. The file can be viewed
   * with a help of [Paraview](www.paraview.org) version 5.5.0 or higher.
   ****************************************************************************/
  void save() const;

  /**
   * \brief Releases computer memory associated with system matrix and
   * right-hand side.
   *****************************************************************************/
  void clear()
  {
    system_matrix.clear();
    system_rhs.reinit(0);
  }

  /**
   * \brief Returns a reference to a dof handler.
   *****************************************************************************/
  const DoFHandler<dim>& get_dof_handler() const { return dof_handler; }

  /**
   * \brief Returns a reference to solution.
   *****************************************************************************/
  const Vector<double>& get_solution() const { return solution; }

  /**
   * \brief Returns the number of active cells in the mesh.
   *****************************************************************************/
  unsigned int get_n_cells() const
  {
    return static_cast<unsigned int>(triangulation_T.n_active_cells());
  }

  /**
   * \brief Returns the total amount of the degrees of freedom.
   *****************************************************************************/
  unsigned int get_n_dofs() const
  {
    return static_cast<unsigned int>(dof_handler.n_dofs());
  }

  /**
   * \brief Returns \f$L^2\f$ error norm.
   *****************************************************************************/
  double get_L2_norm() const { return L2_norm; }

  /**
   * \brief Returns \f$L^{\infty}\f$ error norm.
   *****************************************************************************/
  double get_Linfty_norm() const { return L2_norm; }

  /**
   * \brief Returns degree of the interpolating Lagrange polynomials used
   * for mapping from the reference cell to the real mesh cell and back.
   *****************************************************************************/
  unsigned int get_mapping_degree() const { return mapping_degree; }

  /**
   * \brief Runs the simulation.
   *
   * Executes the following member functions in a proper order:
   * fill_dirichlet_stack(), setup(), assemble(), solve(),
   * project_exact_solution_fcn(), compute_error_norms(), save().
   *****************************************************************************/
  void run()
  {
    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Fill Dirichlet stack");
      fill_dirichlet_stack();
    }
    {
      TMR("Setup");
      setup();
    }
    {
      TMR("Assemble");
      assemble();
    }
    {
      TMR("Solve");
      solve();
    }

    if (exact_solution) {
      if (project_exact_solution) {
        TMR("Project exact solution");
        project_exact_solution_fcn();
      }

      {
        TMR("Compute error norms");
        compute_error_norms();
      }
    }

    {
      TMR("Save");
      save();
    }
  };

  virtual ~Solver2() = default;

protected:
  /**
   * \brief Reference to the mesh.
   *****************************************************************************/
  const Triangulation<dim>& triangulation_T;

  /**
   * \brief Reference to the dof handler that describes the current vector
   * potential, \f$\vec{T}\f$.
   *****************************************************************************/
  const DoFHandler<dim>& dof_handler_T;

  /**
   * \brief Reference to the degrees of freedom that describe the current vector
   * potential, \f$\vec{T}\f$.
   *****************************************************************************/
  const Vector<double>& solution_T;

  /**
   * \brief A map that contains pairs of boundary IDs and the corresponding
   * Dirichlet boundary conditions. All boundary IDs must be odd.
   *****************************************************************************/
  std::map<types::boundary_id, const Function<dim>*> dirichlet_stack;

  /**
   * \brief The finite elements.
   *****************************************************************************/
  const FE_Nedelec<dim> fe;

  /**
   * \brief The finite elements.
   *****************************************************************************/
  DoFHandler<dim> dof_handler;

  /**
   * \brief The solution vector, that is, degrees of freedom yielded by the
   * simulation.
   *****************************************************************************/
  Vector<double> solution;

  /**
   * \brief The projected exact solution vector.
   *****************************************************************************/
  Vector<double> projected_exact_solution;

  /**
   * \brief The constraints associated with the Dirichlet boundary conditions.
   ****************************************************************************/
  AffineConstraints<double> constraints;

  /**
   * \brief The sparsity pattern of the system matrix.
   *****************************************************************************/
  SparsityPattern sparsity_pattern;

  /**
   * \brief The system matrix.
   *****************************************************************************/
  SparseMatrix<double> system_matrix;

  /**
   * \brief The system right-hand side vector.
   *****************************************************************************/
  Vector<double> system_rhs;

  /**
   * \brief The \f$L^2\f$ norm.
   *****************************************************************************/
  double L2_norm;

  /**
   * \brief The \f$L^{\infty}\f$ norm.
   *****************************************************************************/
  double Linfty_norm;

private:
  const unsigned int mapping_degree;
  const double eta_squared;
  const std::string fname;
  const Function<dim>* exact_solution;
  const bool print_time_tables;
  const bool project_exact_solution;
  const bool write_higher_order_cells;

  Vector<float> L2_per_cell;
  Vector<float> Linfty_per_cell;

  // ----------------------------------------------------------------------------
  // These structures and functions are related to the Work Stream algorithm.
  // See article "WorkStream – A Design Pattern for Multicore-Enabled Finite
  // Element Computations." by BRUNO TURCKSIN, MARTIN KRONBICHLER,
  // WOLFGANG BANGERTH for more details.
  // ----------------------------------------------------------------------------

  using IteratorTuple =
    std::tuple<typename DoFHandler<dim>::active_cell_iterator,
               typename DoFHandler<dim>::active_cell_iterator>;

  using IteratorPair = SynchronousIterators<IteratorTuple>;

  struct AssemblyScratchData
  {
    AssemblyScratchData(const FiniteElement<dim>& fe,
                        const DoFHandler<dim>& dof_hand_T,
                        const Vector<double>& dofs_T,
                        unsigned int mapping_degree,
                        double eta_squared);

    AssemblyScratchData(const AssemblyScratchData& scratch_data);

    TheCoefficient<dim, stage> the_coefficient;
    Gamma<dim, stage> gamma;
    RobinRhs<dim, stage> robin_rhs;
    FreeSurfaceCurrent<dim, stage> free_surface_current;

    MappingQ<dim> mapping;
    Constants::QuadratureTableVector<dim> qt;
    FEValues<dim> fe_values;
    FEFaceValues<dim> fe_face_values;

    FEValues<dim> fe_values_T;
    FEFaceValues<dim> fe_face_values_T;

    const unsigned int dofs_per_cell;
    const unsigned int n_q_points;
    const unsigned int n_q_points_face;

    std::vector<double> the_coefficient_list;

    std::vector<Tensor<1, dim>> vector_values;
    std::vector<Tensor<1, dim>> vector_values_face;

    std::vector<double> values;
    std::vector<double> values_face;

    std::vector<double> gamma_list;
    std::vector<Tensor<1, dim>> robin_rhs_list;
    std::vector<Tensor<1, dim>> free_surface_current_list;

    const FEValuesExtractors::Vector ve;
    const FEValuesExtractors::Scalar se;

    const DoFHandler<dim>& dof_hand_T;
    const Vector<double>& dofs_T;

    bool do_robin;
    bool do_K;
    bool do_T_on_boundary;
    const double eta_squared;
  };

  struct AssemblyCopyData
  {
    FullMatrix<double> cell_matrix;
    Vector<double> cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };

  void system_matrix_local(const IteratorPair& IP,
                           AssemblyScratchData& scratch_data,
                           AssemblyCopyData& copy_data);

  void copy_local_to_global(const AssemblyCopyData& copy_data);
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
};

template<int dim, int stage>
void
Solver2<dim, stage>::setup()
{
  dof_handler.reinit(triangulation_T);
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
  for (auto item : dirichlet_stack) {
    Assert(item.first % 2 == 1, ExcInternalError());
  }
#pragma GCC diagnostic pop

  for (auto item : dirichlet_stack)
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0,            // first vector component
      *item.second, // boundary function
      item.first,   // boundary id
      constraints,  // constraints
      MappingQ<dim>(mapping_degree));

  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  if (project_exact_solution)
    projected_exact_solution.reinit(dof_handler.n_dofs());

  if (exact_solution) {
    L2_per_cell.reinit(triangulation_T.n_active_cells());
    Linfty_per_cell.reinit(triangulation_T.n_active_cells());
  }
}

template<int dim, int stage>
void
Solver2<dim, stage>::assemble()
{
  WorkStream::run(
    IteratorPair(
      IteratorTuple(dof_handler.begin_active(), dof_handler_T.begin_active())),
    IteratorPair(IteratorTuple(dof_handler.end(), dof_handler_T.end())),
    *this,
    &Solver2::system_matrix_local,
    &Solver2::copy_local_to_global,
    AssemblyScratchData(
      fe, dof_handler_T, solution_T, mapping_degree, eta_squared),
    AssemblyCopyData());
}

template<int dim, int stage>
Solver2<dim, stage>::AssemblyScratchData::AssemblyScratchData(
  const FiniteElement<dim>& fe,
  const DoFHandler<dim>& dof_hand_T,
  const Vector<double>& dofs_T,
  unsigned int mapping_degree,
  double eta_squared)
  : mapping(mapping_degree)
  , qt(fe.degree - 1)
  , fe_values(mapping,
              fe,
              QGauss<dim>(qt.sim()),
              update_gradients | update_values | update_quadrature_points |
                update_JxW_values)
  , fe_face_values(mapping,
                   fe,
                   QGauss<dim - 1>(qt.sim()),
                   update_values | update_normal_vectors |
                     update_quadrature_points | update_JxW_values)
  , fe_values_T(mapping,
                dof_hand_T.get_fe(),
                QGauss<dim>(qt.sim()),
                update_values | update_gradients)
  , fe_face_values_T(mapping,
                     dof_hand_T.get_fe(),
                     QGauss<dim - 1>(qt.sim()),
                     update_values | update_gradients)
  , dofs_per_cell(fe_values.dofs_per_cell)
  , n_q_points(fe_values.get_quadrature().size())
  , n_q_points_face(fe_face_values.get_quadrature().size())
  , the_coefficient_list(n_q_points)
  , vector_values(n_q_points)
  , vector_values_face(n_q_points_face)
  , values(n_q_points)
  , values_face(n_q_points_face)
  , gamma_list(n_q_points_face)
  , robin_rhs_list(n_q_points_face)
  , free_surface_current_list(n_q_points_face)
  , ve(0)
  , se(0)
  , dof_hand_T(dof_hand_T)
  , dofs_T(dofs_T)
  , eta_squared(eta_squared)
{
}

template<int dim, int stage>
Solver2<dim, stage>::AssemblyScratchData::AssemblyScratchData(
  const AssemblyScratchData& scratch_data)
  : mapping(scratch_data.mapping.get_degree())
  , qt(scratch_data.qt)
  , fe_values(mapping,
              scratch_data.fe_values.get_fe(),
              scratch_data.fe_values.get_quadrature(),
              update_gradients | update_values | update_quadrature_points |
                update_JxW_values)
  , fe_face_values(mapping,
                   scratch_data.fe_face_values.get_fe(),
                   scratch_data.fe_face_values.get_quadrature(),
                   update_values | update_normal_vectors |
                     update_quadrature_points | update_JxW_values)
  , fe_values_T(mapping,
                scratch_data.fe_values_T.get_fe(),
                scratch_data.fe_values_T.get_quadrature(),
                update_values | update_gradients)
  , fe_face_values_T(mapping,
                     scratch_data.fe_face_values_T.get_fe(),
                     scratch_data.fe_face_values_T.get_quadrature(),
                     update_values | update_gradients)
  , dofs_per_cell(fe_values.dofs_per_cell)
  , n_q_points(fe_values.get_quadrature().size())
  , n_q_points_face(fe_face_values.get_quadrature().size())
  , the_coefficient_list(n_q_points)
  , vector_values(n_q_points)
  , vector_values_face(n_q_points_face)
  , values(n_q_points)
  , values_face(n_q_points_face)
  , gamma_list(n_q_points_face)
  , robin_rhs_list(n_q_points_face)
  , free_surface_current_list(n_q_points_face)
  , ve(0)
  , se(0)
  , dof_hand_T(scratch_data.dof_hand_T)
  , dofs_T(scratch_data.dofs_T)
  , eta_squared(scratch_data.eta_squared)
{
}

template<int dim, int stage>
void
Solver2<dim, stage>::system_matrix_local(const IteratorPair& IP,
                                         AssemblyScratchData& scratch_data,
                                         AssemblyCopyData& copy_data)
{
  // See the following boxes:
  // (1) Recipe for static vector solver in 3D
  // (2) Recipe for static vector solver in 2D

  // The comments below refer to these recipes by number, i.e., recipe (1)
  // and recipe (2).

  copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                               scratch_data.dofs_per_cell);

  copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

  copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

  auto cell = std::get<0>(*IP);
  auto cell_T = std::get<1>(*IP);

  scratch_data.fe_values.reinit(cell);
  scratch_data.fe_values_T.reinit(cell_T);

  scratch_data.the_coefficient.value_list(
    scratch_data.fe_values.get_quadrature_points(),
    cell->material_id(),
    cell->user_index(),
    scratch_data.the_coefficient_list);

  if (dim == 2) {
    scratch_data.fe_values_T[SE].get_function_values(scratch_data.dofs_T,
                                                     scratch_data.values);
  } else if (dim == 3) {
    scratch_data.fe_values_T[VE].get_function_values(
      scratch_data.dofs_T, scratch_data.vector_values);
  }

  for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index) {
    for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < scratch_data.dofs_per_cell;
           ++j) { // Integral I_a1+I_a3 in recipes (1) and (2).
        copy_data.cell_matrix(i, j) +=
          ((1 / scratch_data.the_coefficient_list[q_index]) * // 1 / mu
             scratch_data.fe_values[VE].curl(i, q_index) *    // curl N_i
             scratch_data.fe_values[VE].curl(j, q_index)      // curl N_j
           + scratch_data.eta_squared *                       // eta^2
               scratch_data.fe_values[VE].value(i, q_index) * // N_i
               scratch_data.fe_values[VE].value(j, q_index)   // N_j
           ) *
          scratch_data.fe_values.JxW(q_index); // dV (dS in 2D)
      }

      if (dim == 2) { // Integral I_b3-1 in recipe (2).
        copy_data.cell_rhs(i) +=
          scratch_data.values.at(q_index) *                // T
          scratch_data.fe_values[VE].curl(i, q_index)[0] * // curl_s N_i
          scratch_data.fe_values.JxW(q_index);             // dS
      } else if (dim == 3) { // Integral I_b3-1 in recipe (1).
        copy_data.cell_rhs(i) +=
          (scratch_data.vector_values[q_index][0] *
             scratch_data.fe_values[VE].curl(i, q_index)[0] +
           scratch_data.vector_values[q_index][1] *
             scratch_data.fe_values[VE].curl(i, q_index)[1] +
           scratch_data.vector_values[q_index][2] *
             scratch_data.fe_values[VE].curl(i, q_index)[2] // T
           ) *                                              // curl N_i
          scratch_data.fe_values.JxW(q_index);              // dV
      }
    }
  }

  for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
    scratch_data.do_robin = ((cell->face(f)->at_boundary()) &&
                             (cell->face(f)->boundary_id() % 2 == 0) &&
                             (cell->face(f)->boundary_id() != 0));

    scratch_data.do_K =
      ((cell->user_index() > 0) && (cell->face(f)->user_index() > 0));

    scratch_data.do_T_on_boundary = (cell->face(f)->at_boundary());

    Assert(!(scratch_data.do_robin && scratch_data.do_K), ExcInternalError());

    if (scratch_data.do_robin || scratch_data.do_K ||
        scratch_data.do_T_on_boundary) {
      scratch_data.fe_face_values.reinit(cell, f);
      scratch_data.fe_face_values_T.reinit(cell_T, f);

      if (scratch_data.do_robin) {
        scratch_data.gamma.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          scratch_data.fe_face_values.get_normal_vectors(),
          cell->face(f)->boundary_id(),
          cell->material_id(),
          cell->user_index(),
          cell->face(f)->user_index(),
          scratch_data.gamma_list);

        scratch_data.robin_rhs.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          scratch_data.fe_face_values.get_normal_vectors(),
          cell->face(f)->boundary_id(),
          cell->material_id(),
          cell->user_index(),
          cell->face(f)->user_index(),
          scratch_data.robin_rhs_list);
      }

      if (scratch_data.do_K) {
        scratch_data.free_surface_current.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          scratch_data.fe_face_values.get_normal_vectors(),
          cell->material_id(),
          cell->user_index(),
          cell->face(f)->user_index(),
          scratch_data.free_surface_current_list);
      }

      if (scratch_data.do_T_on_boundary) {
        if (dim == 2) {
          scratch_data.fe_face_values_T.get_function_values(
            scratch_data.dofs_T, scratch_data.values_face);
        } else if (dim == 3) {
          scratch_data.fe_face_values_T[VE].get_function_values(
            scratch_data.dofs_T, scratch_data.vector_values_face);
        }
      }

      for (unsigned int q_index_face = 0;
           q_index_face < scratch_data.n_q_points_face;
           ++q_index_face) {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i) {
          if (scratch_data.do_robin) {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j) {
              if (dim == 2) { // Integral I_a2 in recipe (2)
                copy_data.cell_matrix(i, j) +=
                  scratch_data.gamma_list[q_index_face] * // gamma
                  (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                     scratch_data.fe_face_values[VE].value(i, q_index_face)[1] -
                   scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                     scratch_data.fe_face_values[VE].value(i,
                                                           q_index_face)[0]) *
                  (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                     scratch_data.fe_face_values[VE].value(j, q_index_face)[1] -
                   scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                     scratch_data.fe_face_values[VE].value(
                       j, q_index_face)[0]) * //    s           s
                  scratch_data.fe_face_values.JxW(
                    q_index_face);   // (n x N_i) . (n x N_j) dl
              } else if (dim == 3) { // Integral I_a2 in recipe (1).
                copy_data.cell_matrix(i, j) +=
                  scratch_data.gamma_list[q_index_face] * // gamma
                  ((scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[2] -
                    scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[1]) *
                     (scratch_data.fe_face_values.normal_vector(
                        q_index_face)[1] *
                        scratch_data.fe_face_values[VE].value(j,
                                                              q_index_face)[2] -
                      scratch_data.fe_face_values.normal_vector(
                        q_index_face)[2] *
                        scratch_data.fe_face_values[VE].value(
                          j, q_index_face)[1]) +
                   (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[2] -
                    scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[0]) *
                     (scratch_data.fe_face_values.normal_vector(
                        q_index_face)[0] *
                        scratch_data.fe_face_values[VE].value(j,
                                                              q_index_face)[2] -
                      scratch_data.fe_face_values.normal_vector(
                        q_index_face)[2] *
                        scratch_data.fe_face_values[VE].value(
                          j, q_index_face)[0]) +
                   (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[1] -
                    scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[0]) *
                     (scratch_data.fe_face_values.normal_vector(
                        q_index_face)[0] *
                        scratch_data.fe_face_values[VE].value(j,
                                                              q_index_face)[1] -
                      scratch_data.fe_face_values.normal_vector(
                        q_index_face)[1] *
                        scratch_data.fe_face_values[VE].value(
                          j, q_index_face)[0])) // (n x N_i) . (n x N_j)
                  * scratch_data.fe_face_values.JxW(q_index_face); // dS
              } else {
                Assert(false, ExcInternalError());
              }
            }

            // Integral I_b1 in recipes (1) and (2).
            copy_data.cell_rhs(i) +=
              -scratch_data.robin_rhs_list[q_index_face] *             // Q
              scratch_data.fe_face_values[VE].value(i, q_index_face) * // N_i
              scratch_data.fe_face_values.JxW(q_index_face); // dS (dl in 2D)
          } // if (scratch_data.do_robin)

          if (scratch_data.do_K) { // Integral I_b2 in recipes (1) and (2).
            copy_data.cell_rhs(i) +=
              scratch_data.free_surface_current_list[q_index_face] *   // K_f
              scratch_data.fe_face_values[VE].value(i, q_index_face) * // N_i
              scratch_data.fe_face_values.JxW(q_index_face); // dS (dl in 2D)
          }

          if (scratch_data.do_T_on_boundary) {
            if (dim == 2) { // Integral I_b3-2 in recipe (2).
              copy_data.cell_rhs(i) -=
                scratch_data.values.at(q_index_face) * // T
                (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                   scratch_data.fe_face_values[VE].value(i, q_index_face)[1] -
                 scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                   scratch_data.fe_face_values[VE].value(
                     i, q_index_face)[0]) *                    //    s
                scratch_data.fe_face_values.JxW(q_index_face); // (n x N_i) dl
            } else if (dim == 3) { // Integral I_b3-2 in recipe (1).
              copy_data.cell_rhs(i) -=
                (scratch_data.vector_values_face[q_index_face][0] *
                   (scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[2] -
                    scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[1]) -
                 scratch_data.vector_values_face[q_index_face][1] *
                   (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[2] -
                    scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[0]) +
                 scratch_data.vector_values_face[q_index_face][2] *
                   (scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
                      scratch_data.fe_face_values[VE].value(i,
                                                            q_index_face)[1] -
                    scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
                      scratch_data.fe_face_values[VE].value(
                        i, q_index_face)[0]) // T . (n x N_i )
                 ) *
                scratch_data.fe_face_values.JxW(q_index_face); // dS
            } else {
              Assert(false, ExcInternalError());
            }
          }
        }
      }
    }
  }
  cell->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim, int stage>
void
Solver2<dim, stage>::copy_local_to_global(const AssemblyCopyData& copy_data)
{
  constraints.distribute_local_to_global(copy_data.cell_matrix,
                                         copy_data.cell_rhs,
                                         copy_data.local_dof_indices,
                                         system_matrix,
                                         system_rhs);
}

template<int dim, int stage>
void
Solver2<dim, stage>::compute_error_norms()
{
  Weight<dim, stage> weight;
  const Function<dim, double>* mask = &weight;

  Constants::QuadratureTableVector<dim> qt(dof_handler.get_fe().degree - 1);
  QGauss<dim> quadrature(qt.enorm());

  VectorTools::integrate_difference(MappingQ<dim>(mapping_degree),
                                    dof_handler,
                                    solution,
                                    *exact_solution,
                                    L2_per_cell,
                                    quadrature,
                                    VectorTools::L2_norm,
                                    mask);

  L2_norm = VectorTools::compute_global_error(
    triangulation_T, L2_per_cell, VectorTools::L2_norm);

  VectorTools::integrate_difference(MappingQ<dim>(mapping_degree),
                                    dof_handler,
                                    solution,
                                    *exact_solution,
                                    Linfty_per_cell,
                                    QGauss<dim>(1),
                                    VectorTools::Linfty_norm,
                                    mask);

  Linfty_norm = Linfty_per_cell.linfty_norm();
}

template<int dim, int stage>
void
Solver2<dim, stage>::project_exact_solution_fcn()
{
  Constants::QuadratureTableVector<dim> qt(fe.degree - 1);

  AffineConstraints<double> constraints_empty;
  constraints_empty.close();

  VectorTools::project(MappingQ<dim>(mapping_degree),
                       dof_handler,
                       constraints_empty,
                       QGauss<dim>(qt.sim()),
                       *exact_solution,
                       projected_exact_solution);
}

template<int dim, int stage>
void
Solver2<dim, stage>::save() const
{
  std::vector<std::string> solution_names(dim, "VectorField");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim> data_out;

  data_out.add_data_vector(
    dof_handler, solution, solution_names, interpretation);

  if (project_exact_solution && exact_solution) {
    std::vector<std::string> solution_names_ex(dim, "VectorFieldExact");

    data_out.add_data_vector(
      dof_handler, projected_exact_solution, solution_names_ex, interpretation);
  }

  if (exact_solution) {
    data_out.add_data_vector(L2_per_cell, "L2norm");
    data_out.add_data_vector(Linfty_per_cell, "LinftyNorm");
  }

  std::ofstream ofs;

  if (write_higher_order_cells) {
    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<dim> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe.degree + 2,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);

    ofs.open(fname + ".vtu");
    data_out.write_vtu(ofs);

  } else {

    data_out.build_patches();

    ofs.open(fname + ".vtk");
    data_out.write_vtk(ofs);
  }

  ofs.close();
}

} // namespace StaticVectorSolver

#endif
