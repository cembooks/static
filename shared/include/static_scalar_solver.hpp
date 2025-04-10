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

#ifndef StaticScalarSolver_H__
#define StaticScalarSolver_H__

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

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
#include "static_scalar_input.hpp"

#define TMR(__name) TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

namespace StaticScalarSolver {

/**
 * \brief Solves
 * [static scalar boundary value problem](@ref page_seibvp).
 *
 * Implements the following recipes:
 * - (1) Recipe for static scalar solver in 3D
 * - (2) Recipe for static scalar solver in 2D (planar)
 * - (3) Recipe for static scalar solver in 2D (axisymmetric)
 * - (4) Recipe for static scalar solver in 2D (current vect. potential)
 *
 * This class template is intended to be a general solver for problems in
 * electro- and magnetostatics that can be formulated in therms of the
 * electrostatic scalar potential, \f$\Phi\f$, total magnetostatic scalar
 * potential,\f$\Psi\f$, reduced magnetostatic scalar potential, \f$\Theta\f$,
 * two-dimensional magnitude of vector potential, \f$A\f$, and scaled
 * two-dimensional magnitude of vector potential, \f$A'\f$. It can also be
 * used to solve for the current vector potential, \f$T\f$, in planar
 * two-dimensional problems. The calculated potential is saved in a vtk file,
 * see function save() for more details. The Bossavit's diagrams below
 * illustrate the partial differential equations that can be solved with a
 * help of this class template. Note, that in all cases listed below the
 * potential belongs to the \f$H(\text{grad})\f$ function space, i.e., is
 * modeled by the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements.
 *
 * ![](sss_sss/diagram_sss.svg)
 *
 * The table below lists the recommended settings for switching between
 * different types of problems. The letters in the first column of the table
 * correspond to the seven diagrams above. The `dim` parameter is the input
 * parameter of the class template. The other three parameters,
 * `type_of_pde_rhs`, `axisymmetric`, and `vector_potential`, are passed as
 * input parameters to the constructor of the class.
 *
 * Insert | dim | type_of_pde_rhs | axisymmetric | vector_potential |
 * -------|-----|-----------------|--------------|------------------|
 * A), C) |  3  |     0 or 1      |    false     |     false        |
 * B), D) |  2  |     0 or 1      |true or false |     false        |
 *   E)   |  2  |     0 or 1      |    false     |     true         |
 *   F)   |  2  |     0 or 1      |    false     |     true         |
 *   G)   |  2  |     2 or 3      |    false     |     true         |
 *
 * A user of this class is supposed to do the following.
 *
 * - Derive a class or class template from StaticScalarSolver::Solver.
 *
 * - Override the virtual destructor.
 *
 * - Call constructor StaticScalarSolver::Solver(...) and pass to it the degree
 *   of interpolating polynomials and other arguments.
 *
 * - Override the following member functions
 *   @code
 *   virtual void make_mesh() = 0;
 *   virtual void fill_dirichlet_stack() = 0;
 *   virtual void solve() = 0;
 *   @endcode
 *
 * - Instantiate class templates
 *   StaticScalarSolver::TheCoefficient, StaticScalarSolver::PdeRhs
 *   (or StaticScalarSolver::PdeRhsCvp if `type_of_pde_rhs=2`),
 *   StaticScalarSolver::Gamma, StaticScalarSolver::RobinRhs, and
 *   StaticScalarSolver::FreeSurfaceCharge, StaticScalarSolver::Weight.
 *
 * - Explicitly call
 *   @code
 *   StaticScalarSolver::Solver.run();
 *   @endcode
 *   Alternatively, the user can call individual member functions (make_mesh(),
 *   fill_dirichlet_stack(), setup(), solve(), save()) in a proper order.
 *
 * @anchor seis_bnd_convention
 * The boundaries of the mesh must be labeled such that the `boundary_id` member
 * function of a face object returns the corresponding boundary ID. The boundary
 * ID's must obey the following convention.
 * - The Dirichlet boundary conditions are applied on the boundaries with odd
 *   boundary ID's.
 * - The Robin boundary conditions are applied on the boundaries with
 *   even boundary ID's. The boundary ID's in this case must be greater than
 *   zero. The Neumann boundary conditions are considered to be special cases of
 *   the Robin boundary conditions with \f$\gamma=0\f$.
 * - No boundary condition is applied on a boundary with zero ID. Applying no
 *   boundary condition is as good as applying the homogeneous Neumann boundary
 *   condition, \f$ \hat{n}\cdot\vec{\nabla}\Phi = 0\f$, as it is implicitly
 *   implied by the first term of the
 *   [functional](@ref seibvp_functional).
 *   This boundary condition can also be imposed by assigning to a boundary an
 *   even ID greater than zero, and setting \f$ \sigma \f$ and \f$ \gamma \f$
 *   to zero in the `value_list(...)` methods of the class templates
 *   StaticScalarSolver::RobinRhs and StaticScalarSolver::Gamma.
 *
 * The constructor's argument `type_of_pde_rhs` switches the operation of the
 * class template between following four modes:
 *
 * - `type_of_pde_rhs = 1`. In this mode the right-hand side of the partial
 * differential equation is assumed to be a scalar field. Let us for the sake
 * of illustration assume that we are computing the electric scalar potential,
 * \f$\Phi\f$, and that the right-hand side is the free-charge density,
 * \f$\rho_f\f$. Then the partial differential equation reads
 * \f[
 * - \vec{\nabla} \cdot \big( \epsilon \vec{\nabla} \Phi \big)= \rho_f.
 * \f]
 * The corresponding integral in the variational formulation
 * reads
 * \f[
 * \iiint_{\Omega} \rho_f \Phi dV
 * \f]
 * in three dimensions and
 * \f[
 * \iint_{\Omega} \rho_f \Phi dS
 * \f]
 * in two dimensions. In this mode the values of \f$\rho_f\f$ at quadrature
 * points are computed by calling StaticScalarSolver::PdeRhs::value_list.
 *
 * - `type_pde_rhs = 0`. Setting `type_pde_rhs = 0` is the same as setting
 * `type_pde_rhs = 1` and \f$\rho_f = 0\f$. In this mode algorithm saves some
 * time on calling StaticScalarSolver::PdeRhs::value_list and evaluating the
 * two integrals above.
 *
 * - `type_pde_rhs = 3`. In this mode the class template computes the
 * two-dimensional current vector potential, \f$T\f$, by solving the following
 * partial differential equation:
 * \f[
 * - \vec{\nabla} \cdot \big(\vec{\nabla} T \big)=
 *   \vec{\nabla} \overset{S}{\times}\vec{J}_f.
 * \f]
 * This mode works only in two dimensions (the class template
 * StaticVectorSolver::Solver1 must be used for calculating the current vector
 * potential in three dimensions). The following two integrals represent the
 * right-hand side of the partial differential equation in the functional:
 * \f[
 * \iint_{\Omega}\vec{J}_f\cdot\bigg(\vec{\nabla}\overset{V}{\times}T\bigg)dS-
 * \underbrace{\oint_{\Gamma} \vec{J}_f \cdot \bigg(\hat{n} \overset{V}{\times}
 * T \bigg) dl}_{\text{Boundary integral}}. \f] That is to say, in this mode
 * the source on the right-hand side of the partial differential equation is
 * not a scalar field (such as \f$\rho_f\f$), but a two-dimensional vector
 * field, \f$\vec{J}_f\f$. The class template calls an object of the type
 * StaticScalarSolver::PdeRhsCvp to evaluate the values of \f$\vec{J}_f\f$ at
 * quadrature points.
 *
 * - `type_pde_rhs = 2`. There is only one difference between this mode and the
 * mode `type_pde_rhs = 3`. In this mode the boundary integral, see above,
 * is not computed. This can save simulation time if \f$\vec{J}_f = 0\f$ on the
 * boundary by definition of the problem, see [(cvp-ii)](@ref page_cvp_ii)
 * numerical experiment for an example.
 *
 * This class template utilizes the
 * [WorkStream](https://www.dealii.org/current/doxygen/deal.II/namespaceWorkStream.html)
 * technology of deal.II. The amount of threads used can be limited as the
 * following
 * @code
 * #include <deal.II/base/multithread_info.h>
 *
 * MultithreadInfo::set_thread_limit(nr_threads_max);
 * @endcode
 *
 * @note Application examples:
 * - [flc/](@ref page_flc), [rho/](@ref page_rho), [sch/](@ref page_sch) -
 *   Electric scalar potential, \f$\Phi\f$, in three-dimensional and planar
 *   two-dimensional problems.
 * - [flc-axi/](@ref page_flc_axi), [sch-axi/](@ref page_sch_axi) - Electric
 *   scalar potential, \f$\Phi\f$, in axisymmetric two-dimensional problems.
 * - [sld-i/](@ref page_sld_i) - Total scalar magnetic potential, \f$\Psi\f$,
 *   in three-dimensional and planar two-dimensional problems.
 * - [sld-ii/](@ref page_sld_ii) - Total scalar magnetic potential,
 *   \f$\Theta\f$, in three-dimensional and planar two-dimensional problems.
 * - [mwr/](@ref page_mwr) - Magnetic vector potential, \f$A\f$, in planar
 *   two-dimensional problems.
 * - [ssol-i-axi/](@ref page_ssol_i_axi),
 *   [ssol-ii-axi/](@ref page_ssol_ii_axi),
 *   [ssol-iii-axi/](@ref page_ssol_iii_axi), - Scaled magnetic vector
 *   potential, \f$A'\f$, in axisymmetric two-dimensional problems.
 * - [cvp-ii/](@ref page_cvp_ii), [mms-vt-ii](@ref page_mms_vt_ii) - Current
 *   vector potential, \f$T\f$, in planar two-dimensional problems.
 *****************************************************************************/
template<int dim, int stage = 1>
class Solver
{
public:
  Solver() = delete;
  /**
   * \brief The only constructor.
   *
   * @param[in] p - The degree of the interpolating Lagrange polynomials in
   * finite elements that model the potential.
   * @param[in] mapping_degree - The degree of the interpolating Lagrange
   * polynomials used for mapping. Setting it to 1 will do in the most of the
   * cases. Note, that it makes sense to attach a meaningful manifold to the
   * triangulation if this parameter is greater than 1.
   * @param[in] type_of_pde_rhs - Switches between four modes of operation,
   * see above.
   * @param[in] fname - The name of the output files without extension.
   * @param[in] exact_solution - Points to an object that describes the exact
   * solution to the problem. It is needed for calculating error norms. It is a
   * responsibility of the user to make sure that the object exists at the time
   * of the execution of run() or compute_error_norms().
   * @param[in] axisymmetric - If true, assumes that the problem is
   * axisymmetric. If `axisymmetric = true`, `dim` must equal 2.
   * @param[in] vector_potential - If true, assumes that the problem is
   * two-dimensional and formulated in terms of the magnitude of vector
   * potential, \f$A\f$, or in terms of the scaled magnitude of vector
   * potential, \f$A'\f$, or current vector potential, \f$T\f$. If
   * `vector_potential = true`, `dim` must equal 2.
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the Lagrange finite elements (FE_Q) and saves the
   * result into the output file next to the solution. This may be useful for
   * debugging purposes as a comparison between the projected exact solution and
   * the solution to the boundary value problem can yield a hint on where to
   * search for bugs.
   * @param[in] write_higher_order_cells - Switches between the two modes of
   * operation of the save() function, see the description of save().
   *****************************************************************************/
  Solver(unsigned int p,
         unsigned int mapping_degree,
         unsigned int type_of_pde_rhs,
         std::string fname = "data",
         const Function<dim>* exact_solution = nullptr,
         bool axisymmetric = false,
         bool vector_potential = false,
         bool print_time_tables = false,
         bool project_exact_solution = false,
         bool write_higher_order_cells = false)
    : fe(p)
    , mapping_degree(mapping_degree)
    , type_of_pde_rhs(type_of_pde_rhs)
    , fname(fname)
    , exact_solution(exact_solution)
    , axisymmetric(axisymmetric)
    , vector_potential(vector_potential)
    , print_time_tables(print_time_tables)
    , project_exact_solution(project_exact_solution)
    , write_higher_order_cells(write_higher_order_cells)
  {
    Assert(((dim == 2) || (dim == 3)), ExcInternalError());
    Assert(p < 6, ExcInternalError());
    Assert(type_of_pde_rhs < 4, ExcInternalError());

    if (axisymmetric) {
      Assert(
        dim == 2,
        ExcMessage("The setting axisymmetric=true is only allowed if dim=2."));

      Assert(
        type_of_pde_rhs < 2,
        ExcMessage(
          "The settings axisymmetric=true and type_of_pde_rhs>1 (corresponds to \
the modes for computing the current vector potential) are not compatible. \
An out-of-plane scalar current vector potential implies that the free-current \
density is an in-plane vector. Strictly speaking, it is impossible to setup a \
curl-curl equation on an axisymmetric problem domain if the free-current density \
is an in-plane vector. Charge conservation will fail."));
    }

    if (vector_potential) {
      Assert(dim == 2,
             ExcMessage(
               "The setting vector_potential=true can only be used if dim=2."));
      Assert(((type_of_pde_rhs == 2) || (type_of_pde_rhs == 3)),
             ExcMessage("The setting vector_potential=true can only be used if \
type_of_pde_rhs=2 or type_of_pde_rhs=3."));
    }

    if ((type_of_pde_rhs == 0) || (type_of_pde_rhs == 1)) {
      Assert(!vector_potential,
             ExcMessage(
               "The settings type_of_pde_rhs=0 and type_of_pde_rhs=1 can only \
be used if vector_potential=false."));
    }

    if ((type_of_pde_rhs == 2) || (type_of_pde_rhs == 3)) {
      Assert(dim == 2,
             ExcMessage(
               "The settings type_of_pde_rhs=2 and type_of_pde_rhs=3 can only \
be used if dim=2."));
      Assert(vector_potential,
             ExcMessage(
               "The settings type_of_pde_rhs=2 and type_of_pde_rhs=3 can only \
be used if vector_potential=true"));
    }
  }

  /**
   * \brief Initializes the data member
   * StaticScalarSolver::Solver::triangulation.
   *
   * This function must be overridden by the user. It must generate or load
   * the mesh, label the boundaries, and, if necessary, assign user IDs. This
   * function is an ideal place for binding manifolds to the mesh. The last is
   * a reasonable thing to do if `mapping_degree > 1`. The mesh must be stored
   * in the data member of this class,
   * StaticScalarSolver::Solver::triangulation.
   *****************************************************************************/
  virtual void make_mesh() = 0;

  /**
   * \brief Initializes the data member
   * StaticScalarSolver::Solver::dirichlet_stack.
   *
   * This function must be overridden by the user. It must initialize the stack
   * of the Dirichlet boundary conditions. For example,
   * @code
   *
   * using namespace dealii;
   *
   * const types::boundary_id boundary_id_1 = 1;
   * const types::boundary_id boundary_id_2 = 3;
   *
   * const Functions::ZeroFunction<dim> dirichlet_bc_1;
   * const Functions::ConstantFunction<dim> dirichlet_bc_2(1.0);
   *
   * template<int dim>
   * void SolverMyProblem<dim>::fill_dirichlet_stack()
   * {
   *   Solver<dim>::dirichlet_stack =
   *     {{boundary_id_1, & dirichlet_bc_1},
   *      {boundary_id_2, & dirichlet_bc_2}};
   * }
   * @endcode
   *
   * The boundary IDs must be odd numbers, see above the
   * [convention](@ref seis_bnd_convention)
   * on the boundary IDs.
   *****************************************************************************/
  virtual void fill_dirichlet_stack() = 0;

  /**
   * \brief Solves the system of linear equations.
   *****************************************************************************/
  virtual void solve() = 0;

  /**
   * \brief Initializes system matrix and the right-hand side vector, etc.
   *
   * Initialises
   * StaticScalarSolver::Solver::system_matrix,
   * StaticScalarSolver::Solver::system_rhs, and some private arrays.
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
   * The mesh and the finite elements are the same as are used for the
   * numerical solution of the boundary value problem. The exact solution will
   * be saved in the output file next to the numerical solution to the boundary
   * value problem. This function works properly only if the exact solution
   * is submitted to the constructor via the input parameter
   * `exact_solution` and `project_exact_solution = true`.
   *****************************************************************************/
  void project_exact_solution_fcn();

  /**
   * \brief Saves simulation results into a vtk or vtu file.
   *
   * The following data are saved:
   * - The calculated potential under the name "ScalarField".
   * - The \f$L^2\f$ error norm associated with the calculated potential under
   *   the name "L2norm". One value per mesh cell is saved.
   * - The \f$H^1\f$ error norm associated with the calculated potential under
   *   the name "H1seminorm". One value per mesh cell is saved.
   * - The \f$L^{\infty}\f$ error norm associated with the calculated potential
   *   under the name "LinftyNorm". One value per mesh cell is saved.
   * - The exact solution expressed as a linear combination of the shape
   *   functions of the
   *   [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
   *   `finite elements is saved under the name "ScalarFieldExact". The
   *   "Scalarfield" and "SclalarFieldExact" are modeled by exactly the same
   *   finite elements.
   *
   * The "L2norm", "H1seminorm", "LinftyNorm", and "ScalarFieldExact" are
   * saved only if an exact solution is submitted to the constructor. Moreover,
   * "ScalarFieldExact" is calculated and saved only
   * if `project_exact_solution = true`.
   *
   * If `write_higher_order_cells = false`, the name of the file is computed by
   * appending ".vtk" to the string contained by the parameter `fname` passed
   * to the constructor. The file can be inspected with a help of
   * [VisIt](https://visit.llnl.gov) or [Paraview](www.paraview.org).
   * Higher-order cells are saved as regular quadrilaterals and hexahedra.
   * If `write_higher_order_cells = true`, the name of the file is computed by
   * appending ".vtu" to the string contained by the parameter `fname`. The
   * data is saved preserving the higher-order cells. The file can be viewed
   * with a help of [Paraview](www.paraview.org) version 5.5.0 or higher.
   * [VisIt](https://visit.llnl.gov) cannot load higher-order cells.
   *****************************************************************************/
  void save() const;

  /**
   * \brief Saves the system matrix and the right-hand side into a csv file.
   *
   * All the zeros are included into the csv files. This is a very dumb and
   * inefficient way of saving sparse matrices. On the positive side - it is
   * very easy and straightforward to read the csv files. This function may be
   * useful for debugging. One can assemble the system on a coarse mesh
   * (so there are a few mesh cells and the system matrix is small) and export
   * the system matrix together with the right-hand side into another program
   * such as Matlab or GNU Octave for an analysis.
   *
   * @param[in] fname - A stem of the names of the output files. The matrix
   * will be saved into fname_matrix.csv file. The right-hand side will be save
   * into fname_rhs.csv file.
   *****************************************************************************/
  void save_matrix_and_rhs_to_csv(std::string fname) const;

  /**
   * \brief Releases computer memory associated with the system matrix and
   * right-hand side.
   *****************************************************************************/
  void clear()
  {
    system_matrix.clear();
    system_rhs.reinit(0);
  }

  /**
   * \brief Returns a reference to triangulation.
   *****************************************************************************/
  const Triangulation<dim>& get_tria() const { return triangulation; }

  /**
   * \brief Returns a reference to dof handler.
   *****************************************************************************/
  const DoFHandler<dim>& get_dof_handler() const { return dof_handler; }

  /**
   * \brief Returns a reference to the solution.
   *****************************************************************************/
  const Vector<double>& get_solution() const { return solution; }

  /**
   * \brief Returns the number of active cells in the mesh.
   *****************************************************************************/
  unsigned int get_n_cells() const
  {
    return static_cast<unsigned int>(triangulation.n_active_cells());
  }

  /**
   * \brief Returns the number of vertices.
   *****************************************************************************/
  unsigned int get_n_vertices() const
  {
    return static_cast<unsigned int>(triangulation.n_vertices());
  }

  /**
   * \brief Returns the number of used vertices.
   *****************************************************************************/
  unsigned int get_n_used_vertices() const
  {
    return static_cast<unsigned int>(triangulation.n_used_vertices());
  }

  /**
   * \brief Returns the number of lines.
   *****************************************************************************/
  unsigned int get_n_lines() const
  {
    return static_cast<unsigned int>(triangulation.n_lines());
  }

  /**
   * \brief Returns the total amount of the degrees of freedom.
   *****************************************************************************/
  unsigned int get_n_dofs() const
  {
    return static_cast<unsigned int>(dof_handler.n_dofs());
  }

  /**
   * \brief Returns the value of `type_of_pde_rhs`.
   *****************************************************************************/
  unsigned int get_type_of_pde_rhs() const { return type_of_pde_rhs; }

  /**
   * \brief Returns \f$L^2\f$ error norm.
   *****************************************************************************/
  double get_L2_norm() const { return L2_norm; }

  /**
   * \brief Returns \f$H^1\f$ error norm.
   *****************************************************************************/
  double get_H1_norm() const { return H1_norm; }

  /**
   * \brief Returns \f$L^{\infty}\f$ error norm.
   *****************************************************************************/
  double get_Linfty_norm() const { return Linfty_norm; }

  /**
   * \brief Returns degree of the interpolating Lagrange polynomials used
   * for mapping from the reference cell to the real mesh cell and back.
   *****************************************************************************/
  unsigned int get_mapping_degree() const { return mapping_degree; }

  /**
   * \brief Runs the simulation.
   *
   * Executes the following member functions in a proper order: make_mesh();
   * fill_dirichlet_stack(); setup(); assemble(); solve();
   * project_exact_solution_fcn(); compute_error_norms(); save();
   *****************************************************************************/
  void run()
  {
    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Make mesh");
      make_mesh();
    }
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

  virtual ~Solver() = default;

protected:
  /**
   * \brief A map that contains pairs of boundary IDs and the corresponding
   * Dirichlet boundary conditions.
   *
   * All boundary IDs must be odd numbers,
   * see the [convention](@ref seis_bnd_convention) above. The algorithm will
   * loop through the map and apply the boundary conditions one-by-one.
   *****************************************************************************/
  std::map<types::boundary_id, const Function<dim>*> dirichlet_stack;

  /**
   * \brief The mesh.
   *****************************************************************************/
  Triangulation<dim> triangulation;

  /**
   * \brief The finite elements.
   *****************************************************************************/
  const FE_Q<dim> fe;

  /**
   * \brief The degrees-of-freedom handler.
   *****************************************************************************/
  DoFHandler<dim> dof_handler;

  /**
   * \brief The solution vector, i.e., degrees of freedom yielded by the
   * simulation.
   *****************************************************************************/
  Vector<double> solution;

  /**
   * \brief The projected exact solution vector.
   *****************************************************************************/
  Vector<double> projected_exact_solution;

  /**
   * \brief The constraints associated with the Dirichlet boundary conditions.
   *****************************************************************************/
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
   * \brief The \f$L^2\f$ error norm.
   *****************************************************************************/
  double L2_norm;

  /**
   * \brief The \f$L^{\infty}\f$ error norm.
   *****************************************************************************/
  double Linfty_norm;

  /**
   * \brief The \f$H^1\f$ error semi-norm.
   *****************************************************************************/
  double H1_norm;

private:
  const unsigned int mapping_degree;
  const unsigned int type_of_pde_rhs;
  const std::string fname;
  const Function<dim>* exact_solution;
  const bool axisymmetric;
  const bool vector_potential;
  const bool print_time_tables;
  const bool project_exact_solution;
  const bool write_higher_order_cells;

  Vector<float> L2_per_cell;
  Vector<float> Linfty_per_cell;
  Vector<float> H1_per_cell;

  // ----------------------------------------------------------------------------
  // These structures and functions are related to the Work Stream algorithm.
  // See article "WorkStream – A Design Pattern for Multicore-Enabled Finite
  // Element Computations." by BRUNO TURCKSIN, MARTIN KRONBICHLER,
  // WOLFGANG BANGERTH for more details.
  // ----------------------------------------------------------------------------
  struct AssemblyScratchData
  {
    AssemblyScratchData(const FiniteElement<dim>& fe,
                        unsigned int type_of_pde_rhs,
                        bool axisymmetric,
                        bool vector_potential,
                        unsigned int mapping_degree);

    AssemblyScratchData(const AssemblyScratchData& scratch_data);

    TheCoefficient<dim, stage> the_coefficient;
    PdeRhs<dim, stage> pde_rhs;
    PdeRhsCvp<dim, stage> pde_rhs_cvp;
    Gamma<dim, stage> gamma;
    RobinRhs<dim, stage> robin_rhs;
    FreeSurfaceCharge<dim, stage> free_surface_charge;

    MappingQ<dim> mapping;
    Constants::QuadratureTableScalar<dim> qt;
    FEValues<dim> fe_values;
    FEFaceValues<dim> fe_face_values;

    const unsigned int dofs_per_cell;
    const unsigned int n_q_points;
    const unsigned int n_q_points_face;

    std::vector<double> the_coefficient_list;
    std::vector<double> pde_rhs_list;
    std::vector<Tensor<1, dim>> pde_rhs_cvp_list;      // Used in 2D only
    std::vector<Tensor<1, dim>> pde_rhs_cvp_list_face; // Used in 2D only
    std::vector<double> gamma_list;
    std::vector<double> robin_rhs_list;
    std::vector<double> free_surface_charge_list;

    const unsigned int type_of_pde_rhs;
    const bool axisymmetric;
    const bool vector_potential;
    double axi_mult; // Equals the distance to the axis of rotation symmetry, r,
                     // in the "Recipe for static scalar solver in 2D
                     // (axisymmetric)". Equals 1.0 in all other recipes. All
                     // integrands are multiplied by this multiplier.

    // Four auxiliary variables.
    bool do_robin;
    bool do_kappa;
    bool do_Jf_on_boundary;
    double robin_rhs_or_kappa;
  };

  struct AssemblyCopyData
  {
    FullMatrix<double> cell_matrix;
    Vector<double> cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };

  void system_matrix_local(
    const typename DoFHandler<dim>::active_cell_iterator& cell,
    AssemblyScratchData& scratch_data,
    AssemblyCopyData& copy_data);

  void copy_local_to_global(const AssemblyCopyData& copy_data);
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
};

template<int dim, int stage>
void
Solver<dim, stage>::setup()
{
  dof_handler.reinit(triangulation);
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
  for (auto item : dirichlet_stack) {
    Assert(item.first % 2 == 1, ExcInternalError());
  }
#pragma GCC diagnostic pop

  VectorTools::interpolate_boundary_values(
    MappingQ<dim>(mapping_degree), dof_handler, dirichlet_stack, constraints);

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
    L2_per_cell.reinit(triangulation.n_active_cells());
    H1_per_cell.reinit(triangulation.n_active_cells());
    Linfty_per_cell.reinit(triangulation.n_active_cells());
  }
}

template<int dim, int stage>
void
Solver<dim, stage>::assemble()
{
  WorkStream::run(
    dof_handler.begin_active(),
    dof_handler.end(),
    *this,
    &Solver::system_matrix_local,
    &Solver::copy_local_to_global,
    AssemblyScratchData(
      fe, type_of_pde_rhs, axisymmetric, vector_potential, mapping_degree),
    AssemblyCopyData());
}

template<int dim, int stage>
Solver<dim, stage>::AssemblyScratchData::AssemblyScratchData(
  const FiniteElement<dim>& fe,
  unsigned int type_of_pde_rhs,
  bool axisymmetric,
  bool vector_potential,
  unsigned int mapping_degree)
  : mapping(mapping_degree)
  , qt(fe.degree)
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
  , dofs_per_cell(fe_values.dofs_per_cell)
  , n_q_points(fe_values.get_quadrature().size())
  , n_q_points_face(fe_face_values.get_quadrature().size())
  , the_coefficient_list(n_q_points)
  , pde_rhs_list(n_q_points)
  , pde_rhs_cvp_list(n_q_points, Tensor<1, dim>())
  , pde_rhs_cvp_list_face(n_q_points_face, Tensor<1, dim>())
  , gamma_list(n_q_points_face)
  , robin_rhs_list(n_q_points_face)
  , free_surface_charge_list(n_q_points_face)
  , type_of_pde_rhs(type_of_pde_rhs)
  , axisymmetric(axisymmetric)
  , vector_potential(vector_potential)
  , axi_mult(1.0)
{
}

template<int dim, int stage>
Solver<dim, stage>::AssemblyScratchData::AssemblyScratchData(
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
  , dofs_per_cell(fe_values.dofs_per_cell)
  , n_q_points(fe_values.get_quadrature().size())
  , n_q_points_face(fe_face_values.get_quadrature().size())
  , the_coefficient_list(n_q_points)
  , pde_rhs_list(n_q_points)
  , pde_rhs_cvp_list(n_q_points, Tensor<1, dim>())
  , pde_rhs_cvp_list_face(n_q_points_face, Tensor<1, dim>())
  , gamma_list(n_q_points_face)
  , robin_rhs_list(n_q_points_face)
  , free_surface_charge_list(n_q_points_face)
  , type_of_pde_rhs(scratch_data.type_of_pde_rhs)
  , axisymmetric(scratch_data.axisymmetric)
  , vector_potential(scratch_data.vector_potential)
  , axi_mult(1.0)
{
}

template<int dim, int stage>
void
Solver<dim, stage>::system_matrix_local(
  const typename DoFHandler<dim>::active_cell_iterator& cell,
  AssemblyScratchData& scratch_data,
  AssemblyCopyData& copy_data)
{
  // See the following boxes:
  // (1) Recipe for static scalar solver in 3D
  // (2) Recipe for static scalar solver in 2D (planar)
  // (3) Recipe for static scalar solver in 2D (axisymmetric)
  // (4) Recipe for static scalar solver in 2D (current vect. potential)

  // The comments below refer to these recipes by number, i.e., recipe (1),
  // recipe (2), recipe (3), and recipe (4).

  copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                               scratch_data.dofs_per_cell);

  copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

  copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

  scratch_data.fe_values.reinit(cell);

  if ((scratch_data.type_of_pde_rhs == 2) ||
      (scratch_data.type_of_pde_rhs == 3)) {
    // The coefficient equals 1.0 if the current vector potential, T, is
    // computed.
    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      scratch_data.the_coefficient_list[q_index] = 1.0;
  } else {
    scratch_data.the_coefficient.value_list(
      scratch_data.fe_values.get_quadrature_points(),
      cell->material_id(),
      cell->user_index(),
      scratch_data.the_coefficient_list);
  }

  if (scratch_data.type_of_pde_rhs == 1)
    scratch_data.pde_rhs.value_list(
      scratch_data.fe_values.get_quadrature_points(),
      cell->material_id(),
      cell->user_index(),
      scratch_data.pde_rhs_list);

  if ((scratch_data.type_of_pde_rhs == 2) ||
      (scratch_data.type_of_pde_rhs == 3))
    scratch_data.pde_rhs_cvp.value_list(
      scratch_data.fe_values.get_quadrature_points(),
      cell->material_id(),
      cell->user_index(),
      scratch_data.pde_rhs_cvp_list);

  for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index) {
    if ((scratch_data.axisymmetric) && (!scratch_data.vector_potential)) {
      scratch_data.axi_mult =
        scratch_data.fe_values.quadrature_point(q_index)[0];
    }

    for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j) {
        // Integral I_a1 in recipes (1), (2), (3), and (4).
        copy_data.cell_matrix(i, j) +=
          scratch_data.axi_mult *
          scratch_data.the_coefficient_list[q_index] *    // epsilon
          scratch_data.fe_values.shape_grad(i, q_index) * // grad N_i
          scratch_data.fe_values.shape_grad(j, q_index) * // grad N_j
          scratch_data.fe_values.JxW(q_index);            // dV
      }

      switch (scratch_data.type_of_pde_rhs) {
        case 0:
          // Integral I_b3 in recipes (1), (2), and (3) with rho_f=0.
          copy_data.cell_rhs(i) = 0.0;
          break;
        case 1:
          // Integral I_b3 in recipes (1), (2), and (3).
          copy_data.cell_rhs(i) +=
            scratch_data.axi_mult * scratch_data.pde_rhs_list[q_index] * // rho
            scratch_data.fe_values.shape_value(i, q_index) *             // N_i
            scratch_data.fe_values.JxW(q_index);                         // dV
          break;
        case 2:
        case 3:
          // Integral I_b3-1 in recipe (4).
          copy_data.cell_rhs(i) +=
            (scratch_data.pde_rhs_cvp_list[q_index][0] *
               scratch_data.fe_values.shape_grad(i, q_index)[1] -
             scratch_data.pde_rhs_cvp_list[q_index][1] *
               scratch_data.fe_values.shape_grad(i, q_index)[0]) *
            scratch_data.fe_values.JxW(q_index); // (Jf . curl_v N_i) dS
          break;
        default:
          Assert(false,
                 ExcMessage(
                   "The parameter type_of_pde_rhs equals " +
                   std::to_string(scratch_data.type_of_pde_rhs) +
                   ". Only the following values are allowed: 0, 1, and 2."));
          break;
      }
    }
  }

  for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
    scratch_data.do_robin = (cell->face(f)->at_boundary() &&
                             (cell->face(f)->boundary_id() % 2 == 0) &&
                             (cell->face(f)->boundary_id() != 0));

    scratch_data.do_kappa =
      ((cell->user_index() > 0) && (cell->face(f)->user_index() > 0) &&
       (scratch_data.type_of_pde_rhs < 2));

    scratch_data.do_Jf_on_boundary =
      (cell->face(f)->at_boundary() && (scratch_data.type_of_pde_rhs == 3));

    Assert(
      !(scratch_data.do_robin && scratch_data.do_kappa),
      ExcMessage(
        "Robin boundary condition is applied on a boundary. The surface free-current \
charge, kappa_f, exists only on interfaces. No interface is a boundary. \
Therefore, do_robin and do_kappa are mutually exclusive."));

    Assert(
      !(scratch_data.do_kappa && scratch_data.do_Jf_on_boundary),
      ExcMessage(
        "When computing the current vector potential, type_of_pde = 2, J_f is \
integrated over the boundary. The surface free-current charge, kappa_f \
exists only on interfaces. No interface is a boundary. Therefore, \
do_Jf_on_boundary and do_kappa are mutually exclusive."));

    if (scratch_data.do_robin || scratch_data.do_kappa ||
        scratch_data.do_Jf_on_boundary) {
      scratch_data.fe_face_values.reinit(cell, f);

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

      if (scratch_data.do_kappa)
        scratch_data.free_surface_charge.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          scratch_data.fe_face_values.get_normal_vectors(),
          cell->material_id(),
          cell->user_index(),
          cell->face(f)->user_index(),
          scratch_data.free_surface_charge_list);

      if (scratch_data.do_Jf_on_boundary)
        scratch_data.pde_rhs_cvp.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          cell->material_id(),
          cell->user_index(),
          scratch_data.pde_rhs_cvp_list_face);

      for (unsigned int q_index_face = 0;
           q_index_face < scratch_data.n_q_points_face;
           ++q_index_face) {
        if ((scratch_data.axisymmetric) && (!scratch_data.vector_potential))
          scratch_data.axi_mult =
            scratch_data.fe_face_values.quadrature_point(q_index_face)[0];

        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i) {
          if (scratch_data.do_robin) {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j) {
              // Integral I_a2 in recipes (1), (2), and (3).
              copy_data.cell_matrix(i, j) +=
                scratch_data.axi_mult *
                scratch_data.gamma_list[q_index_face] * // gamma
                scratch_data.fe_face_values.shape_value(i,
                                                        q_index_face) * // N_i
                scratch_data.fe_face_values.shape_value(j,
                                                        q_index_face) * // N_j
                scratch_data.fe_face_values.JxW(q_index_face);          // dS
            }
          }

          if (scratch_data.do_robin || scratch_data.do_kappa) {
            scratch_data.robin_rhs_or_kappa = 0.0;

            // If true, the integral below is I_b1.
            if (scratch_data.do_robin)
              scratch_data.robin_rhs_or_kappa =
                scratch_data.robin_rhs_list[q_index_face];

            // If true, the integral below is I_b2.
            if (scratch_data.do_kappa)
              scratch_data.robin_rhs_or_kappa =
                scratch_data.free_surface_charge_list[q_index_face];

            // Depending on the current context:
            // integral I_b1 in recipes (1), (2), (3), and (4)
            // or
            // integral I_b2 in recipes (1), (2), and (3).
            copy_data.cell_rhs(i) +=
              scratch_data.axi_mult *
              scratch_data.robin_rhs_or_kappa * // sigma or kappa_f
              scratch_data.fe_face_values.shape_value(i, q_index_face) * // N_i
              scratch_data.fe_face_values.JxW(q_index_face);             // dS
          }

          if (scratch_data.do_Jf_on_boundary) {
            // Integral I_b3-2 in recipe (4).
            copy_data.cell_rhs(i) -=
              scratch_data.fe_face_values.shape_value(i, q_index_face) *
              (scratch_data.pde_rhs_cvp_list_face[q_index_face][0] *
                 scratch_data.fe_face_values.normal_vector(q_index_face)[1] -
               scratch_data.pde_rhs_cvp_list_face[q_index_face][1] *
                 scratch_data.fe_face_values.normal_vector(
                   q_index_face)[0] //        V
               ) *
              scratch_data.fe_face_values.JxW(
                q_index_face); // J_f.(n x N_i)dl =
                               //           S
                               // = N_i(J_f x n)dl
          }
        } // for (unsigned int i = 0; ...
      }   // for (unsigned int q_index_face = 0; ...
    }     // if (scratch_data.do_robin || scratch_data.do_kappa ||
          // scratch_data.do_Jf_on_boundary)
  }       // for (unsigned int f = 0; ...

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim, int stage>
void
Solver<dim, stage>::copy_local_to_global(const AssemblyCopyData& copy_data)
{
  constraints.distribute_local_to_global(copy_data.cell_matrix,
                                         copy_data.cell_rhs,
                                         copy_data.local_dof_indices,
                                         system_matrix,
                                         system_rhs);
}

template<int dim, int stage>
void
compute_L2_error_norm(const DoFHandler<dim>& dof_handler,
                      const Triangulation<dim>& triangulation,
                      const Vector<double>& solution,
                      const Function<dim>* exact_solution,
                      Vector<float>& L2_per_cell,
                      double& L2_norm,
                      unsigned int mapping_degree)
{
  Weight<dim, stage> weight;
  const Function<dim, double>* mask = &weight;

  Constants::QuadratureTableScalar<dim> qt(dof_handler.get_fe().degree);
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
    triangulation, L2_per_cell, VectorTools::L2_norm);
}

template<int dim, int stage>
void
compute_H1_error_norm(const DoFHandler<dim>& dof_handler,
                      const Triangulation<dim>& triangulation,
                      const Vector<double>& solution,
                      const Function<dim>* exact_solution,
                      Vector<float>& H1_per_cell,
                      double& H1_norm,
                      unsigned int mapping_degree)
{
  Weight<dim, stage> weight;
  const Function<dim, double>* mask = &weight;

  Constants::QuadratureTableScalar<dim> qt(dof_handler.get_fe().degree);
  QGauss<dim> quadrature(qt.enorm());

  VectorTools::integrate_difference(MappingQ<dim>(mapping_degree),
                                    dof_handler,
                                    solution,
                                    *exact_solution,
                                    H1_per_cell,
                                    quadrature,
                                    VectorTools::H1_seminorm,
                                    mask);

  H1_norm = VectorTools::compute_global_error(
    triangulation, H1_per_cell, VectorTools::H1_seminorm);
}

template<int dim, int stage>
void
compute_Linfty_error_norm(const DoFHandler<dim>& dof_handler,
                          const Triangulation<dim>& triangulation,
                          const Vector<double>& solution,
                          const Function<dim>* exact_solution,
                          Vector<float>& Linfty_per_cell,
                          double& Linfty_norm,
                          unsigned int mapping_degree)
{
  Weight<dim, stage> weight;
  const Function<dim, double>* mask = &weight;

  Constants::QuadratureTableScalar<dim> qt(dof_handler.get_fe().degree);
  QGauss<dim> quadrature(qt.enorm());

  VectorTools::integrate_difference(MappingQ<dim>(mapping_degree),
                                    dof_handler,
                                    solution,
                                    *exact_solution,
                                    Linfty_per_cell,
                                    QGauss<dim>(1),
                                    VectorTools::Linfty_norm,
                                    mask);

  Linfty_norm = VectorTools::compute_global_error(
    triangulation, Linfty_per_cell, VectorTools::Linfty_norm);
}

template<int dim, int stage>
void
Solver<dim, stage>::compute_error_norms()
{
  if (exact_solution) {
    Threads::Task<void> task_l2 =
      Threads::new_task(&compute_L2_error_norm<dim, stage>,
                        dof_handler,
                        triangulation,
                        solution,
                        exact_solution,
                        L2_per_cell,
                        L2_norm,
                        mapping_degree);

    Threads::Task<void> task_h1 =
      Threads::new_task(&compute_H1_error_norm<dim, stage>,
                        dof_handler,
                        triangulation,
                        solution,
                        exact_solution,
                        H1_per_cell,
                        H1_norm,
                        mapping_degree);

    Threads::Task<void> task_linfty =
      Threads::new_task(&compute_Linfty_error_norm<dim, stage>,
                        dof_handler,
                        triangulation,
                        solution,
                        exact_solution,
                        Linfty_per_cell,
                        Linfty_norm,
                        mapping_degree);
  }
}

template<int dim, int stage>
void
Solver<dim, stage>::project_exact_solution_fcn()
{
  Constants::QuadratureTableScalar<dim> qt(fe.degree);

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
Solver<dim, stage>::save() const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "ScalarField");

  if (exact_solution) {
    data_out.add_data_vector(L2_per_cell, "L2norm");
    data_out.add_data_vector(Linfty_per_cell, "LinftyNorm");
    data_out.add_data_vector(H1_per_cell, "H1seminorm");

    if (project_exact_solution)
      data_out.add_data_vector(projected_exact_solution, "ScalarFieldExact");
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

template<int dim, int stage>
void
Solver<dim, stage>::save_matrix_and_rhs_to_csv(std::string fname) const
{
  std::ofstream ofs_matrix(fname + "_matrix.csv");
  std::ofstream ofs_rhs(fname + "_rhs.csv");

  for (unsigned int i = 0; i < system_matrix.m(); ++i) {
    ofs_rhs << system_rhs(i);
    if (i < (system_matrix.m() - 1))
      ofs_rhs << "\n";

    for (unsigned int j = 0; j < system_matrix.n(); ++j) {
      ofs_matrix << std::scientific << std::setprecision(16)
                 << system_matrix.el(i, j);

      if (j < (system_matrix.m() - 1))
        ofs_matrix << ", ";
    }
    if (i < (system_matrix.m() - 1))
      ofs_matrix << "\n";
  }

  ofs_rhs.close();
  ofs_matrix.close();
}

} // namespace StaticScalarSolver

#endif
