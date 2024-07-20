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

#ifndef StaticVectorSolverI_H__
#define StaticVectorSolverI_H__

#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/numerics/vector_tools_project.h>

#include <fstream>
#include <ios>
#include <iostream>
#include <iomanip>
#include <map>

#include "constants.hpp"
#include "settings.hpp"
#include "static_vector_input.hpp"

#define VE scratch_data.ve

#define TMR(__name) \
	TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

namespace StaticVectorSolver
{
/**
 * \brief Solves
 * [static vector boundary value problem](@ref page_veibvp).
 *
 * Implements the following recipes:
 * - (1) Recipe for static vector solver in 3D
 * - (2) Recipe for static vector solver in 2D
 * - (3) Recipe for static vector solver in 3D (current vect. potential)
 *
 * This class template is intended to be a general solver for problems in
 * magnetostatics that can be formulated in terms of the magnetic vector
 * potential, \f$\vec{A}\f$. It can also be used for calculating the current
 * vector potential, \f$\vec{T}\f$, i.e., converting a closed-form analytical
 * expression for \f$\vec{J}_f\f$ into \f$\vec{T}\f$ expressed as a finite-element 
 * field-function. Such calculated \f$\vec{T}\f$ can be used as an input
 * for StaticVectorSolver::Solver2. The Bossavit's diagrams below illustrate the
 * partial differential equations that can be solved with a help of this class
 * template. In all five cases the the vector potential is modeled by the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements.
 *
 *![](svs_svs/diagram_svs.svg)
 *
 * A user of this class is supposed to do the following.
 *
 * ![](diagrams/diagram_svs.svg)
 *
 * - Derive a class from StaticVectorSolver::Solver1.
 *
 * - Override the virtual destructor.
 *
 * - Call constructor StaticVectorSolver::Solver1(...) and pass to it the degree
 *   of interpolating polynomials and other arguments.
 *
 * - Override member functions
 *   @code
 *	 virtual void make_mesh() = 0;
 *	 virtual void fill_dirichlet_stack() = 0;
 *	 virtual void solve() = 0;
 *	 @endcode
 *
 * - Implement the member functions
 *   @code
 *   void value_list(...);
 *	 @endcode
 *   of the classes StaticVectorSolver::TheCoefficient, StaticVectorSolver::PdeRhs,
 *	 StaticVectorSolver::Gamma, StaticVectorSolver::RobinRhs, and
 *	 StaticVectorSolver::FreeSurfaceCurrent.
 *
 * - Implement member function
 *	 @code
 *	 double value(...) const;
 *	 @endcode
 *	 of the class StaticVectorSolver::Weight.
 *
 * - Explicitly call
 *	 @code
 *   StaticVectorSolver::run();
 *	 @endcode
 *	 Alternatively, the user may call individual member functions (such as
 *	 make_mesh, fill_dirichlet_stack, setup, etc. ) in a proper order.
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
 * boundary_id() member function of a face object returns the corresponding
 * boundary ID.
 *
 * @anchor veis_bnd_convention
 *
 * The boundary ID's must obey the following convention.
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
 *   setting \f$ \gamma \f$ and \f$\vec{Q}\f$ to zero in the value_list(...)
 *   methods of the classes
 *   StaticVectorSolver::Gamma and StaticVectorSolver::RobinRhs.
 *
 * @note Application examples:
 * - [cvp-i/](@ref page_cvp_i),
 * [mms-v/](@ref page_mms_v),
 * [mms-vt-i/](@ref page_mms_vt_i),
 * [ssol-i/](@ref page_ssol_i),
 * [ssol-ii/](@ref page_ssol_ii),
 * [ssol-iii/](@ref page_ssol_iii).
 *****************************************************************************/
template<int dim, int stage = 1>
class Solver1
{
public:
	Solver1() = delete;

/**
 * \brief The only constructor.
 *
 * The meanings of the input parameters of this constructor, except for the
 * second parameter, are relatively straightforward. They are described below.
 * The second input parameter, type_of_pde_rhs, defines how an object
 * derived from this class template interprets the data provided by an object
 * derived from StaticVectorSolver::PdeRhs.
 *
 * @anchor svsi-type_of_pde_rhs
 * When solving for the magnetic vector potential, \f$\vec{A}\f$ the following
 * options are available:
 *
 * - type_of_pde_rhs = 0. There is no free-current density in the problem domain. 
 * The right-hand side of the partial differential equation is zero,i.e.,
 * \f[
 * \vec{\nabla}\times\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\times\vec{A}\bigg)
 * + \eta^2 \vec{A} = 0.
 * \f]
 * in a three-dimensional space and
 * \f[
 * \vec{\nabla}\overset{V}{\times}\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\overset{S}{\times}\vec{A}\bigg)+ \eta^2 \vec{A} = 0.
 * \f]
 * in a two-dimensional space. The data provided by StaticVectorSolver::PdeRhs is not used.
 *
 * - type_of_pde_rhs = 1. The data provided by StaticVectorSolver::PdeRhs is interpreted as the 
 * free-current density, \f$\vec{J}_f\f$, i.e.,  
 * \f[
 * \vec{\nabla}\times\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\times\vec{A}\bigg)
 * + \eta^2 \vec{A} = \vec{J}_f.
 * \f]
 * in a three-dimensional space and
 * \f[
 * \vec{\nabla}\overset{V}{\times}\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\overset{S}{\times}\vec{A}\bigg)+ \eta^2 \vec{A} = \vec{J}_f.
 * \f]
 * in a two-dimensional space.
 * The corresponding therm of the
 * functional is
 * \f[
 * \iiint_{\Omega} \vec{J}_f \cdot \vec{A} dV
 * \f]
 * in a three-dimensional space and
 * \f[
 * \iint_{\Omega} \vec{J}_f \cdot \vec{A} dS
 * \f]
 * in a two-dimensional space.
 *
 * - type_of_pde_rhs = 2. The data provided by StaticVectorSolver::PdeRhs is interpreted as
 * vector current potential, \f$\vec{T}\f$, i.e.,
 * \f[
 * \vec{\nabla}\times\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\times\vec{A}\bigg)
 * + \eta^2 \vec{A} = \vec{\nabla}\times\vec{T}.
 * \f]
 * in a three-dimensional space and
 * \f[
 * \vec{\nabla}\overset{V}{\times}\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\overset{S}{\times}\vec{A}\bigg)+ \eta^2 \vec{A} = \vec{\nabla}\overset{V}{\times} T.
 * \f]
 * in a two-dimensional space. Then the corresponding term
 * of the functional is
 * \f[
 * \iiint_{\Omega} \vec{T} \cdot \bigg(\vec{\nabla}\times\vec{A}\bigg) dV
 * \f]
 * in a three-dimensional space and
 * \f[
 * \iint_{\Omega} T \bigg(\vec{\nabla}\overset{S}{\times}\vec{A}\bigg) dS
 * \f]
 * in a two-dimensional space.
 *
 * - type_of_pde_rhs = 3. The data provided by StaticVectorSolver::PdeRhs is interpreted as
 * vector current potential, \f$\vec{T}\f$, i.e.,
 * \f[
 * \vec{\nabla}\times\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\times\vec{A}\bigg)
 * + \eta^2 \vec{A} = \vec{\nabla}\times\vec{T}.
 * \f]
 * in a three-dimensional space and
 * \f[
 * \vec{\nabla}\overset{V}{\times}\bigg(\dfrac{1}{\mu_0} \vec{\nabla}\overset{S}{\times}\vec{A}\bigg)+ \eta^2 \vec{A} = \vec{\nabla}\overset{V}{\times} T.
 * \f]
 * in a two-dimensional space. Then the corresponding
 * terms of the functional are
 * \f[
 * \iiint_{\Omega} \vec{T}\cdot\bigg(\vec{\nabla}\times\vec{A}\bigg) dV -
 * \iint_{\Gamma_{\Omega}}\vec{T}\cdot\bigg(\hat{n}\times\vec{A}\bigg)dS
 * \f]
 * in a three dimensional space and
 * \f[
 * \iint_{\Omega} T \bigg(\vec{\nabla}\overset{S}{\times}\vec{A}\bigg) dS -
 * \int_{\Gamma_{\Omega}} T \bigg( \hat{n} \overset{S}{\times} \vec{A}\bigg)dl
 * \f]
 * in a two-dimensional space.
 * 
 * The option type_of_pde_rhs = 2 differs from the option type_of_pde_rhs = 3
 * at one point only: the boundary integral is not calculated
 * if type_of_pde_rhs = 2. This can help to reduce simulation time a bit if the
 * current vector potential equals zero on the boundary due to a Dirichlet 
 * boundary condition. 
 *
 * The same four options are available when solving for the current 
 * vector potential, \f$\vec{T}\f$. In the case of \f$\vec{T}\f$, however,
 * the first two options, i.e., type_of_pde_rhs=0 and type_of_pde_rhs=1, 
 * do not make much sense: the curl of the free-current density must be 
 * present on the right-hand side of the partial differential equation. 
 * The program does not filter out this two options. When solving for 
 * the current vector potential, \f$\vec{T}\f$, the following two options 
 * make sense:
 *
 * - type_of_pde_rhs = 2. The data provided by StaticVectorSolver::PdeRhs 
 * is interpreted as free-current density, \f$\vec{J}_f\f$, i.e.,
 * \f[
 * \vec{\nabla}\times\bigg(\vec{\nabla}\times\vec{T}\bigg)
 * + \eta^2 \vec{T} = \vec{\nabla}\times\vec{J}_f.
 * \f]
 * in a three-dimensional space and
 * \f[
 * \vec{\nabla}\overset{V}{\times}\bigg(\vec{\nabla}\overset{S}{\times}\vec{T}\bigg)
 * + \eta^2 \vec{T} = \vec{\nabla}\overset{V}{\times} J_f.
 * \f]
 * in a two-dimensional space. Then the corresponding terms of the 
 * functional are
 * \f[
 * \iiint_{\Omega} \vec{J}_f\cdot\bigg(\vec{\nabla}\times\vec{T}\bigg) dV 
 * \f]
 * in a three dimensional space and
 * \f[
 * \iint_{\Omega} J_f \bigg(\vec{\nabla}\overset{S}{\times}\vec{T}\bigg) dS 
 * \f]
 * in a two-dimensional space.
 *
 * - type_of_pde_rhs = 3. The data provided by StaticVectorSolver::PdeRhs 
 * is interpreted as free-current density, \f$\vec{J}_f\f$.
 * \f[
 * \vec{\nabla}\times\bigg(\vec{\nabla}\times\vec{T}\bigg)
 * + \eta^2 \vec{T} = \vec{\nabla}\times\vec{J}_f.
 * \f]
 * in a three-dimensional space and
 * \f[
 * \vec{\nabla}\overset{V}{\times}\bigg(\vec{\nabla}\overset{S}{\times}\vec{T}\bigg)
 * + \eta^2 \vec{T} = \vec{\nabla}\overset{V}{\times} J_f.
 * \f]
 * in a two-dimensional space. Then the 
 * corresponding terms of the functional are
 * \f[
 * \iiint_{\Omega} \vec{J}_f\cdot\bigg(\vec{\nabla}\times\vec{T}\bigg) dV -
 * \iint_{\Gamma_{\Omega}}\vec{J}_f\cdot\bigg(\hat{n}\times\vec{T}\bigg)dS
 * \f]
 * in a three dimensional space and
 * \f[
 * \iint_{\Omega} J_f \bigg(\vec{\nabla}\overset{S}{\times}\vec{T}\bigg) dS -
 * \int_{\Gamma_{\Omega}} J_f \bigg( \hat{n} \overset{S}{\times} \vec{T}\bigg)dl
 * \f]
 * in a two-dimensional space.
 *
 * Here again the option type_of_pde_rhs = 3 differs from the option 
 * type_of_pde_rhs = 2 by the boundary integral in the functional.
 *
 * Note, the code that implements the solver for \f$\vec{A}\f$ is identical
 * to the code that implements the solver for \f$\vec{T}\f$ as the list of
 * arguments of the constructor of this class template does not contain an
 * argument that toggles the code between two modes: "solving for A"
 * mode and "solving for T" mode. The user toggles the code between
 * these two modes by feeding the right data through StaticVectorSolver::PdeRhs
 * and StaticVectorSolver::TheCoefficient and ignoring options 
 * type_of_pde_rhs = 0 and type_of_pde_rhs = 1 when solving for \f$\vec{T}\f$.
 * If the magnetic vector potential, \f$\vec{A}\f$, is being computed, 
 * StaticVectorSolver::PdeRhs::value_list must return the values of \f$\vec{T}\f$ 
 * (type_of_pde_rhs = 2, type_of_pde_rhs = 3) or values of \f$\vec{J}_f\f$
 * (type_of_pde_rhs = 1). If the current vector potential, \f$\vec{T}\f$, is being 
 * computed, StaticVectorSolver::PdeRhs::value_list must return the values of 
 * \f$\vec{J}_f\f$ (type_of_pde_rhs = 2, type_of_pde_rhs = 3). 
 * The StaticVectorSolver::TheCoefficient must return 1.0 when solving for 
 * \f$\vec{T}\f$ and \f$\mu\f$ when solving for \f$\vec{A}\f$.
 *
 * @param[in] p - Degree of the FE_Nedelec finite elements.
 * @param[in] mapping_degree - The degree of the interpolating Lagrange
 * polynomials used for mapping. Setting it to 1 will do in the most of the
 * cases. Note, that it makes sense to attach a meaningful manifold to the
 * triangulation if this parameter is greater than 1.
 * @param[in] type_of_pde_rhs - Defines how this class interprets the data
 * provided by StaticVectorSolver::PdeRhs, see above.
 * @param[in] eta_squared - The gauging parameter, \f$\eta^2\f$, in the
 * [partial differential equation](@ref page_veibvp).
 * @param[in] fname - The name of the output files without extension. Names of
 * the output files will be generated by appending simulation conditions to
 * this string.
 * @param[in] exact_solution - Points to an object that describes the exact
 * solution to the problem. It is needed for calculating error norms. It is a
 * responsibility of the user to make sure that the object exists at the time of
 * the execution of run() or compute_error_norms().
 * @param[in] print_time_tables - If true, prints time tables on the screen.
 * @param[in] project_exact_solution - If true, projects the exact solution
 * onto the space spanned by the Nedelec finite elements and saves 
 * the result into the vtk file next to the solution. This may be useful for 
 * debugging purposes as a comparison between the projected exact solution and 
 * the solution to the boundary value problem can yield a hint on where to 
 * search for bugs.
 *****************************************************************************/
	Solver1(
			unsigned int p,
			unsigned int mapping_degree,
			unsigned int type_of_pde_rhs = 3,
			double eta_squared = 0.0,
			std::string fname = "data",
			const Function<dim> * exact_solution = nullptr,
			bool print_time_tables = false,
			bool project_exact_solution = false):
				fe(p),
				mapping_degree(mapping_degree),
				type_of_pde_rhs(type_of_pde_rhs),
				eta_squared(eta_squared),
				fname(fname),
				exact_solution(exact_solution),
				print_time_tables(print_time_tables),
				project_exact_solution(project_exact_solution)
	{
		Assert( p < 3, ExcInternalError());
	}

/**
 * \brief Initializes the data member StaticVectorSolver::Solver1::triangulation.
 *
 * This function must be overridden by the user. It must generate the mesh,
 * label the boundaries, and, if necessary, assign user IDs. The mesh must be
 * stored in the data member of this class
 * StaticVectorSolver::Solver1::triangulation.
 *****************************************************************************/
	virtual void make_mesh() = 0;

/**
 *	\brief Initializes the data member StaticVectorSolver::Solver1::dirichlet_stack.
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
 * 	Solver<dim, stage>::dirichlet_stack =
 * 		{{boundary_id_1, & dirichlet_bc_1},
 * 		 {boundary_id_2, & dirichlet_bc_2}};
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
 * StaticVectorSolver::Solver1::system_matrix,
 * StaticVectorSolver::Solver1::system_rhs and other arrays.
 * Applies the Dirichlet boundary conditions by constraining the system
 * matrix. Distributes degrees of freedom.
 *****************************************************************************/
	void setup();

/**
 * \brief Assembles the system matrix and the right-hand side vector.
 *****************************************************************************/
	void assemble();

/**
 * \brief Computes \f$L^2\f$ error norms.
 *****************************************************************************/
	void compute_error_norms();

/**
 * \brief Projects exact solution.
 *
 * The mesh and the finite elements, are the same as are used for the numerical
 * solution of the boundary vale problem. The exact solution will be saved in
 * the vtk file next to the numerical solution to the boundary value problem.
 * This function works properly only if the exact solution is submitted to the
 * constructor via the input parameter exact_solution and
 * project_exact_solution=true.
 *****************************************************************************/
	void project_exact_solution_fcn();

/**
 * \brief Saves the result of the simulation into a vtk file.
 *
 * The calculated vector potential is saved into a vtk file under a name
 * "VectorField"
 ****************************************************************************/
	void save() const;

/**
 * \brief Releases computer memory associated with system matrix and
 * right-hand side.
 *****************************************************************************/
	void clear() {system_matrix.clear(); system_rhs.reinit(0);}

/**
 * \brief Returns a reference to triangulation.
 *****************************************************************************/
	const Triangulation<dim> & get_tria() const {return triangulation;}

/**
 * \brief Returns a reference to dof handler.
 *****************************************************************************/
	const DoFHandler<dim> & get_dof_handler() const {return dof_handler;}

/**
 * \brief Returns a reference to solution.
 *****************************************************************************/
	const Vector<double> & get_solution() const {return solution;}

/**
 * \brief Returns the number of active cells in the mesh.
 *****************************************************************************/
	unsigned int get_n_cells() const
		{return static_cast<unsigned int>(triangulation.n_active_cells());}

/**
 * \brief Returns the total amount of the degrees of freedom.
 *****************************************************************************/
	unsigned int get_n_dofs() const
		{return static_cast<unsigned int>(dof_handler.n_dofs());}

/**
 * \brief Returns the value of **type_of_pde_rhs**.
 *****************************************************************************/
	unsigned int get_rhs_type() const {return type_of_pde_rhs;}

/**
 * \brief Returns \f$L^2\f$ error norm.
 *****************************************************************************/
	double get_L2_norm() const {return L2_norm;}

/**
 * \brief Returns \f$L^{infty}\f$ error norm.
 *****************************************************************************/
	double get_Linfty_norm() const {return Linfty_norm;}

/**
 * \brief Returns degree of the interpolating Lagrange polynomials used
 * for mapping from the reference cell to the real mesh cell and back.
 *****************************************************************************/
	unsigned int get_mapping_degree() const {return mapping_degree;}

/**
 * \brief Runs the simulation.
 *
 * Executes the following member functions in a proper order: make_mesh();
 * fill_dirichlet_stack(); setup(); assemble(); solve(); save();
 *****************************************************************************/
	void run()
	{
		TimerOutput::OutputFrequency tf =
			(print_time_tables) ? TimerOutput::summary : TimerOutput::never;

		TimerOutput timer(
			std::cout,
			tf,
			TimerOutput::cpu_and_wall_times_grouped);

		{TMR("Make mesh"); make_mesh();}
		{TMR("Fill Dirichlet stack"); fill_dirichlet_stack();}
		{TMR("Setup"); setup();}
		{TMR("Assemble"); assemble();}
		{TMR("Solve"); solve();}

		if (exact_solution)
		{
			if (project_exact_solution)
				{TMR("Project exact solution"); project_exact_solution_fcn();}

			{TMR("Compute error norms"); compute_error_norms();}
		}

		{TMR("Save"); save();}
	};

	virtual ~Solver1() = default;

protected:

	/**
	 * \brief A map that contains pairs of boundary IDs and the corresponding
	 * Dirichlet boundary conditions. All boundary IDs must be odd.
	 ****************************************************************************/
	std::map<types::boundary_id, const Function<dim> *> dirichlet_stack;

	/**
	 * \brief The mesh.
	 ****************************************************************************/
	Triangulation<dim> triangulation;

	/**
	 * \brief The finite elements.
	 ****************************************************************************/
	const FE_Nedelec<dim> fe;

	/**
	 * \brief The dof handler.
	 ****************************************************************************/
	DoFHandler<dim> dof_handler;

	/**
	 * \brief The solution vector, that is, degrees of freedom yielded by the
	 * simulation.
	 ****************************************************************************/
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
	 ****************************************************************************/
	SparsityPattern sparsity_pattern;

	/**
	 * \brief The system matrix.
	 ****************************************************************************/
	SparseMatrix<double> system_matrix;

	/**
	 * \brief The system right-hand side vector.
	 ****************************************************************************/
	Vector<double> system_rhs;

/**
 * \brief The \f$L^2\f$ norm.
 *****************************************************************************/
	double L2_norm;

/**
 * \brief The \f$L^{infty}\f$ norm.
 *****************************************************************************/
	double Linfty_norm;

private:
	const unsigned int mapping_degree;
	const unsigned int type_of_pde_rhs;
	const double eta_squared;
	const std::string fname;
	const Function<dim> * exact_solution;
  const bool print_time_tables;
	const bool project_exact_solution;

	Vector<float> L2_per_cell;
	Vector<float> Linfty_per_cell;

// ----------------------------------------------------------------------------
// These structures and functions are related to the Work Stream algorithm.
// See article "WorkStream – A Design Pattern for Multicore-Enabled Finite
// Element Computations." by BRUNO TURCKSIN, MARTIN KRONBICHLER,
// WOLFGANG BANGERTH for more details.
// ----------------------------------------------------------------------------
	struct AssemblyScratchData
	{
		AssemblyScratchData(const FiniteElement<dim> &fe,
			const unsigned int type_of_pde_rhs,
			unsigned int mapping_degree);

		AssemblyScratchData(const AssemblyScratchData & scratch_data);

	  TheCoefficient<dim, stage> the_coefficient;
		PdeRhs<dim, stage> pde_rhs;
		Gamma<dim, stage> gamma;
		RobinRhs<dim, stage> robin_rhs;
		FreeSurfaceCurrent<dim, stage> free_surface_current;

		MappingQ<dim> mapping;
		Constants::QuadratureTableVector<dim> qt;
		FEValues<dim> fe_values;
		FEFaceValues<dim> fe_face_values;

		const unsigned int dofs_per_cell;
		const unsigned int n_q_points;
		const unsigned int n_q_points_face;

		std::vector<double> the_coefficient_list;
		std::vector<Tensor<1,dim>> pde_rhs_list;
		std::vector<Tensor<1,dim>> pde_rhs_list_face;
		std::vector<double> gamma_list;
		std::vector<Tensor<1,dim>> robin_rhs_list;
		std::vector<Tensor<1,dim>> free_surface_current_list;

		const unsigned int type_of_pde_rhs;
		const FEValuesExtractors::Vector ve;
		bool do_robin;
		bool do_K;
		bool do_T_on_boundary;
	};

	struct AssemblyCopyData
	{
		FullMatrix<double>	cell_matrix;
		Vector<double>	cell_rhs;
		std::vector<types::global_dof_index>	local_dof_indices;
	};

	void system_matrix_local(
		const typename DoFHandler<dim>::active_cell_iterator &cell,
		AssemblyScratchData & scratch_data,
		AssemblyCopyData & copy_data);

	void copy_local_to_global(const AssemblyCopyData &copy_data);
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
};

template<int dim, int stage>
void Solver1<dim, stage>::setup()
{
	dof_handler.reinit(triangulation);
	dof_handler.distribute_dofs(fe);

	constraints.clear();
	DoFTools::make_hanging_node_constraints(dof_handler, constraints);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
	for (auto item : dirichlet_stack)
	{
		Assert(item.first %2 == 1, ExcInternalError());
	}
#pragma GCC diagnostic pop

	for (auto item : dirichlet_stack)
		VectorTools::project_boundary_values_curl_conforming_l2(
			dof_handler,
			0, // first vector component
			*item.second, // boundary function
			item.first, // boundary id
			constraints, // constraints
			MappingQ<dim> (mapping_degree)
			);

	constraints.close();

	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(
		dof_handler,
		dsp,
		constraints,
		false);

	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());

	if (project_exact_solution)
		projected_exact_solution.reinit(dof_handler.n_dofs());

	if (exact_solution)
	{
		L2_per_cell.reinit(triangulation.n_active_cells());
		Linfty_per_cell.reinit(triangulation.n_active_cells());
	}
}

template<int dim, int stage>
void Solver1<dim, stage>::assemble()
{
	WorkStream::run(
		dof_handler.begin_active(),
		dof_handler.end(),
		*this,
		&Solver1::system_matrix_local,
		&Solver1::copy_local_to_global,
		AssemblyScratchData(fe, type_of_pde_rhs, mapping_degree),
		AssemblyCopyData());
}

template<int dim, int stage>
Solver1<dim, stage>::AssemblyScratchData::AssemblyScratchData(
		const FiniteElement<dim> &fe, const unsigned int type_of_pde_rhs,
		unsigned int mapping_degree):
		mapping(mapping_degree),
		qt(fe.degree-1),
		fe_values(mapping,
			fe,
			QGauss<dim>(qt.sim()),
			update_gradients |
			update_values |
			update_quadrature_points |
			update_JxW_values),
		fe_face_values(mapping,
			fe,
			QGauss<dim - 1>(qt.sim()),
			update_values |
			update_normal_vectors |
			update_quadrature_points |
			update_JxW_values),
		dofs_per_cell(fe_values.dofs_per_cell),
		n_q_points(fe_values.get_quadrature().size()),
		n_q_points_face(fe_face_values.get_quadrature().size()),
		the_coefficient_list(n_q_points),
		pde_rhs_list(n_q_points),
		pde_rhs_list_face(n_q_points_face),
		gamma_list(n_q_points_face),
		robin_rhs_list(n_q_points_face),
		free_surface_current_list(n_q_points_face),
		type_of_pde_rhs(type_of_pde_rhs),
		ve(0)
{
}

template<int dim, int stage>
Solver1<dim, stage>::AssemblyScratchData::AssemblyScratchData(
		const AssemblyScratchData & scratch_data):
		mapping(scratch_data.mapping.get_degree()),
	  qt(scratch_data.qt),
		fe_values(mapping,
			scratch_data.fe_values.get_fe(),
			scratch_data.fe_values.get_quadrature(),
			update_gradients |
			update_values |
			update_quadrature_points |
			update_JxW_values),
		fe_face_values(mapping,
			scratch_data.fe_face_values.get_fe(),
			scratch_data.fe_face_values.get_quadrature(),
			update_values |
			update_normal_vectors |
			update_quadrature_points |
			update_JxW_values),
		dofs_per_cell(fe_values.dofs_per_cell),
		n_q_points(fe_values.get_quadrature().size()),
		n_q_points_face(fe_face_values.get_quadrature().size()),
		the_coefficient_list(n_q_points),
		pde_rhs_list(n_q_points),
		pde_rhs_list_face(n_q_points_face),
		gamma_list(n_q_points_face),
		robin_rhs_list(n_q_points_face),
		free_surface_current_list(n_q_points_face),
		type_of_pde_rhs(scratch_data.type_of_pde_rhs),
		ve(0)
{
}

template<int dim, int stage>
void Solver1<dim, stage>::system_matrix_local(
	const typename DoFHandler<dim>::active_cell_iterator & cell,
	AssemblyScratchData & scratch_data, AssemblyCopyData & copy_data)
{
// See the following boxes:
// (1) Recipe for static vector solver in 3D
// (2) Recipe for static vector solver in 2D
// (3) Recipe for static vector solver in 3D (current vect. potential)

// The comments below refer to these recipes by number, i.e., recipe (1),
// recipe (2), and recipe (3).

	copy_data.cell_matrix.reinit(
		scratch_data.dofs_per_cell, scratch_data.dofs_per_cell);

	copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

	copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

	scratch_data.fe_values.reinit(cell);

	scratch_data.the_coefficient.value_list(
			scratch_data.fe_values.get_quadrature_points(),
			cell->material_id(),
			cell->user_index(),
			scratch_data.the_coefficient_list);

	scratch_data.pde_rhs.value_list(
			scratch_data.fe_values.get_quadrature_points(),
			cell->material_id(),
			cell->user_index(),
			scratch_data.pde_rhs_list);

// The curl of a Nedelec shape function is Tensor<1,3> in three dimensions and
// Tensor<1,1> in two dimensions. It is somewhat difficult to fit this into the
// "dim" class template paradigm of deal.II. Consequently, we have to use the
// if (dim==2) and if (dim==3) filters or instantiate this function template
// explicitly for two and three dimensions. We choose the former.
//
// If the current vector potential, T, is used in recipes (1) and (2), it must
// be implemented by the same class that implements the right-hand side of the
// curl-curl equation. Neglect the fact that the right-hand side is the curl
// of the current vector potential, not the current vector potential itself.
// So, "pde_rhs" is a misnomer if the current vector potential is used on the
// right-hand side of the PDE. The current vector potential must be Tensor<1,3>
// in three dimensions. In two dimensions it must fill the first component of
// Tensor<1,2>. The second component is ignored:
//
// Tensor<1,2> a;
// a[0] = T;
// a[1] = whatever;
//
// The same holds for the free-current density, J_f, if the right-hand side of
// the PDE is the curl of J_f, i.e., recipe (3).

	for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
	{
		for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
		{
			for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
			{
				copy_data.cell_matrix(i, j) +=
					(// Integral I_a1+I_a3 in recipes (1), (2), and (3).
						(1 / scratch_data.the_coefficient_list[q_index]) *  // 1 / mu
						scratch_data.fe_values[VE].curl(i, q_index) *       // curl N_i
						scratch_data.fe_values[VE].curl(j, q_index)         // curl N_j
						+
						eta_squared *                                       // eta^2
						scratch_data.fe_values[VE].value(i, q_index) *      // N_i
						scratch_data.fe_values[VE].value(j, q_index)        // N_j
						) * scratch_data.fe_values.JxW(q_index);            // dV (dS in 2D)
			}

			switch(scratch_data.type_of_pde_rhs)
			{
				case 0:
					// Integral I_b3 in recipes (1) and (2) with J_f=0.
					copy_data.cell_rhs(i) = 0.0;
					break;
				case 1:
					// Integral I_b3 in recipes (1) and (2).
					copy_data.cell_rhs(i) +=
						scratch_data.pde_rhs_list[q_index] *                  // J_f
						scratch_data.fe_values[VE].value(i, q_index) *        // N_i
						scratch_data.fe_values.JxW(q_index);                  // dV (dS in 2D)
					break;
				case 2: case 3:
					if (dim == 2)
					{ // Integral I_b3-1 in recipe (2)
						copy_data.cell_rhs(i) +=
							scratch_data.pde_rhs_list[q_index][0] *             // T
							scratch_data.fe_values[VE].curl(i, q_index)[0] *    // curl_s N_i
							scratch_data.fe_values.JxW(q_index);                // dS
					}
					else if (dim == 3)
					{ // Integral I_b3-1 in recipe (1) and (3).
						copy_data.cell_rhs(i) +=
						(
							scratch_data.pde_rhs_list[q_index][0] *
							scratch_data.fe_values[VE].curl(i, q_index)[0]
							+
							scratch_data.pde_rhs_list[q_index][1] *
							scratch_data.fe_values[VE].curl(i, q_index)[1]
							+
							scratch_data.pde_rhs_list[q_index][2] *
							scratch_data.fe_values[VE].curl(i, q_index)[2]
							)*
							scratch_data.fe_values.JxW(q_index);  // If recipe (1): T.(curl N_i)dV.
							                                      // If recipe (3): J_f.(curl N_i)dV.
						// So, the scratch_data.pde_rhs_list must contain the values of T, or the
						// values of J_f, depending on the context.
					}
					else
					{
						Assert(false, ExcInternalError());
					}
					break;
				default:
					Assert(false, ExcInternalError());
					break;
			}
		}
	}

	for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
	{
		scratch_data.do_robin =
			((cell->face(f)->at_boundary()) &&
			(cell->face(f)->boundary_id() % 2 == 0) &&
			(cell->face(f)->boundary_id() != 0));

		scratch_data.do_K =
			((cell->user_index() > 0) &&
			(cell->face(f)->user_index() > 0));

		scratch_data.do_T_on_boundary =
			((cell->face(f)->at_boundary()) &&
			(scratch_data.type_of_pde_rhs == 3));

		Assert( !(scratch_data.do_robin && scratch_data.do_K), ExcInternalError() );

		if (scratch_data.do_robin || scratch_data.do_K || scratch_data.do_T_on_boundary)
		{

			scratch_data.fe_face_values.reinit(cell, f);

			if (scratch_data.do_robin)
			{
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

			if (scratch_data.do_K)
			{
				scratch_data.free_surface_current.value_list(
					scratch_data.fe_face_values.get_quadrature_points(),
					scratch_data.fe_face_values.get_normal_vectors(),
					cell->material_id(),
					cell->user_index(),
					cell->face(f)->user_index(),
					scratch_data.free_surface_current_list);
			}

			if (scratch_data.do_T_on_boundary)
			{
				scratch_data.pde_rhs.value_list(
					scratch_data.fe_face_values.get_quadrature_points(),
					cell->material_id(),
					cell->user_index(),
					scratch_data.pde_rhs_list_face);
			}

			for (unsigned int q_index_face = 0;
				q_index_face < scratch_data.n_q_points_face; ++q_index_face)
			{
				for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
				{
					if (scratch_data.do_robin)
					{
						for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
						{
							if (dim == 2)
							{
								// Integral I_a2 in recipe (2).
								copy_data.cell_matrix(i, j) +=
								scratch_data.gamma_list[q_index_face] *                    // gamma
								(
									scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[1] -
									scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[0]  //   S
									)*                                                         // n x N_i
								(
									scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
									scratch_data.fe_face_values[VE].value(j, q_index_face)[1] -
									scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
									scratch_data.fe_face_values[VE].value(j, q_index_face)[0]   //   S
									) *                                                         // n x N_j
									scratch_data.fe_face_values.JxW(q_index_face);              // dl
							}
							else if (dim == 3)
							{
								// Integral I_a2 in recipe (1).
								copy_data.cell_matrix(i, j) +=
								scratch_data.gamma_list[q_index_face] *                     //gamma
								(
									(
										scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
										scratch_data.fe_face_values[VE].value(i, q_index_face)[2] -
										scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
										scratch_data.fe_face_values[VE].value(i, q_index_face)[1]
									)*
									(
										scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
										scratch_data.fe_face_values[VE].value(j, q_index_face)[2] -
										scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
										scratch_data.fe_face_values[VE].value(j, q_index_face)[1]
									) +
									(
										scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
										scratch_data.fe_face_values[VE].value(i, q_index_face)[2] -
										scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
										scratch_data.fe_face_values[VE].value(i, q_index_face)[0]
									) *
									(
										scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
										scratch_data.fe_face_values[VE].value(j, q_index_face)[2] -
										scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
										scratch_data.fe_face_values[VE].value(j, q_index_face)[0]
									) +
									(
										scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
										scratch_data.fe_face_values[VE].value(i, q_index_face)[1] -
										scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
										scratch_data.fe_face_values[VE].value(i, q_index_face)[0]
									) *
									(
										scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
										scratch_data.fe_face_values[VE].value(j, q_index_face)[1] -
										scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
										scratch_data.fe_face_values[VE].value(j, q_index_face)[0]
									)
								)
								* scratch_data.fe_face_values.JxW(q_index_face); // (n x N_i) . (n x N_j) dS
							}else
							{
								Assert(false, ExcInternalError());
							}
						}

						// Integral I_b1 in recipes (1), (2), and (3).
						copy_data.cell_rhs(i) +=
						-scratch_data.robin_rhs_list[q_index_face]*                   // Q
						 scratch_data.fe_face_values[VE].value(i, q_index_face) *     // N_i
						 scratch_data.fe_face_values.JxW(q_index_face);               // dS (dl in 2D)
					} // if (scratch_data.do_robin)

					if ( scratch_data.do_K )
					{
						// Integral I_b2 in recipes (1) and (2).
						copy_data.cell_rhs(i) +=
						 scratch_data.free_surface_current_list[q_index_face]*        // K_f
						 scratch_data.fe_face_values[VE].value(i, q_index_face)*      // N_i
						 scratch_data.fe_face_values.JxW(q_index_face);               // dS (dl in 2D)
					}
/**/
					if (scratch_data.do_T_on_boundary)
					{
						if (dim == 2)
						{
							// Integral I_b3-2 in recipe (2).
							copy_data.cell_rhs(i) -=
								scratch_data.pde_rhs_list_face[q_index_face][0] *           // T
							(
								scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
								scratch_data.fe_face_values[VE].value(i, q_index_face)[1] -
								scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
								scratch_data.fe_face_values[VE].value(i, q_index_face)[0]   //   S
								) *                                                         // n x N_i
								scratch_data.fe_face_values.JxW(q_index_face);              // dl
						}
						else if (dim == 3)
						{
							// Integral I_b3-2 in recipes (1) and (3).
							copy_data.cell_rhs(i) -=
							(
									scratch_data.pde_rhs_list_face[q_index_face][0] *
								(
									scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[2] -
									scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[1]
								)
								- scratch_data.pde_rhs_list_face[q_index_face][1] *
								(
									scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[2] -
									scratch_data.fe_face_values.normal_vector(q_index_face)[2] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[0]
								)
								+ scratch_data.pde_rhs_list_face[q_index_face][2] *
								(
									scratch_data.fe_face_values.normal_vector(q_index_face)[0] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[1] -
									scratch_data.fe_face_values.normal_vector(q_index_face)[1] *
									scratch_data.fe_face_values[VE].value(i, q_index_face)[0]
								)
							) * scratch_data.fe_face_values.JxW(q_index_face); // If recipe (1): T.(n x N_i)dS.
					}                                                    // If recipe (3): J_f.(n x N_i)dS.
						// So, the scratch_data.pde_rhs_list_face must contain the values of
						// T or the values of J_f, depending on the context.
						else
						{
							Assert(false, ExcInternalError());
						}
					}
/**/
				}
			}
		}
	}
	cell->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim, int stage>
void Solver1<dim, stage>::copy_local_to_global(const AssemblyCopyData &copy_data)
{
	constraints.distribute_local_to_global(
		copy_data.cell_matrix,
		copy_data.cell_rhs,
		copy_data.local_dof_indices,
		system_matrix,
		system_rhs);
}

template<int dim, int stage>
void Solver1<dim, stage>::compute_error_norms()
{
	Weight<dim,stage> weight;
	const Function<dim, double> * mask = & weight;

	Constants::QuadratureTableVector<dim> qt(dof_handler.get_fe().degree - 1);
	QGauss<dim> quadrature(qt.enorm());

	VectorTools::integrate_difference(
		MappingQ<dim> (mapping_degree),
		dof_handler,
		solution,
		* exact_solution,
		L2_per_cell,
		quadrature,
		VectorTools::L2_norm,
		mask);

	L2_norm = VectorTools::compute_global_error(
		triangulation,
		L2_per_cell,
		VectorTools::L2_norm);

	VectorTools::integrate_difference(
		MappingQ<dim> (mapping_degree),
		dof_handler,
		solution,
		* exact_solution,
		Linfty_per_cell,
		QGauss<dim>(1),
		VectorTools::Linfty_norm,
		mask);

	Linfty_norm = Linfty_per_cell.linfty_norm();
}

template<int dim, int stage>
void Solver1<dim, stage>::project_exact_solution_fcn()
{
	Constants::QuadratureTableVector<dim> qt(fe.degree - 1);

	AffineConstraints<double> constraints_empty;
	constraints_empty.close();

	VectorTools::project(
		MappingQ<dim>(mapping_degree),
		dof_handler,
		constraints_empty,
		QGauss<dim> (qt.sim()),
		* exact_solution,
		projected_exact_solution);
}

template<int dim, int stage>
void Solver1<dim, stage>::save() const
{
	std::vector<std::string> solution_names(dim, "VectorField");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  interpretation(dim,
                 DataComponentInterpretation::component_is_part_of_vector);

	DataOut<dim> data_out;

	data_out.add_data_vector(dof_handler,
                           solution,
                           solution_names,
                           interpretation);

	if (project_exact_solution && exact_solution)
	{
		std::vector<std::string> solution_names_ex(dim, "VectorFieldExact");

		data_out.add_data_vector(dof_handler,
		projected_exact_solution,
		solution_names_ex,
		interpretation);
	}

	if (exact_solution)
	{
		data_out.add_data_vector(L2_per_cell, "L2norm");
		data_out.add_data_vector(Linfty_per_cell, "LinftyNorm");
	}

	data_out.build_patches();

	std::ofstream out(fname + ".vtk");

	data_out.write_vtk(out);
	out.close();
}

} // namespace StaticVectorSolver

#endif

