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

#ifndef ProjectPhiToE_H__
#define ProjectPhiToE_H__

#include "project_Hgrad_to_Hcurl.hpp"

namespace StaticScalarSolver
{
/**
 * \brief Computes the electric field \f$\vec{E}\f$ as a negative
 * gradient of the scalar electric potential, \f$\Phi\f$, i.e.,
 * \f$\vec{E} = - \vec{\nabla} \Phi\f$.
 *
 * The problem is assumed to be two- or three- dimensional. This class template
 * is meant to be used in conjunction with StaticScalarSolver::Solver.
 *
 * The scalar electric potential, \f$\Phi\f$, is assumed to be an output of
 * an object of the class StaticScalarSolver::Solver. That is, \f$\Phi\f$ is
 * expressed as a linear combination of the shape functions of the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements. In other words, \f$\Phi\f$ is in the \f$H(\text{grad})\f$
 * function space. This class template computes electric field \f$\vec{E}\f$
 * as a negative gradient of the scalar electric potential, \f$\Phi\f$, such
 * that \f$\vec{E}\f$ is expressed as a linear combination of the shape functions
 * of the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. That is, \f$\vec{E}\f$ is in the \f$H(\text{curl})\f$
 * function space. In other words, this class template takes a scalar function
 * \f$\Phi\f$ from the \f$H(\text{grad})\f$ function space, calculates its gradient,
 * multiplies the gradient by \f$-1\f$, and projects the result into
 * \f$H(\text{curl})\f$ function space, which is a proper space for the electric
 * field \f$\vec{E}\f$. The operation of projection is illustrated by
 * means of the Bossavit's diagrams below.
 *
 * @anchor txt_prj_PHI_to_E_dia
 *![](sss_prj/diagram_prj_phi_to_e.svg)
 *
 * It is assumed that the object of the class StaticScalarSolver::Solver that
 * yielded \f$\Phi\f$ is still in the computer memory such that the
 * triangulation, the degrees of freedom, and the handler of the degrees of
 * freedom are accessible while \f$\vec{E}\f$ is computed. An object of
 * the StaticScalarSolver::ProjectPHItoE class template reuses the triangulation by
 * creating an additional dof handler associated with the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. That is, \f$\vec{E}\f$ and \f$\Phi\f$ share the same triangulation
 * but are modeled by different finite elements:
 * \f$\vec{E}\f$ - by
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * and \f$\Phi\f$ - by
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
 * In this disposition there are two DoFHandler objects associated with one
 * Triangulation object. The algorithm walks synchronously through both
 * DoFHandler objects and assembles the system matrix and the right-hand side.
 * The algorithm utilizes the
 * [WorkStream](https://www.dealii.org/current/doxygen/deal.II/namespaceWorkStream.html)
 * technology of deal.II. The amount of threads used can be limited as the
 * following
 * @code
 * #include <deal.II/base/multithread_info.h>
 * ...
 * MultithreadInfo::set_thread_limit(nr_threads_max);
 * @endcode
 *
 * The output data are saved into a vtk file. The following data are saved:
 * - The calculated electrostatic field, \f$\vec{E} = - \vec{\nabla} \Phi\f$,
 *   under the name "VectorField".
 * - The \f$L^2\f$ error norm associated with the calculated electrostatic field
 *   under the name "L2norm". One value per mesh cell is saved.
 * - The exact solution projected onto \f$H(\text{curl})\f$ function space is
 *   saved under the name "VectorFieldExact". The "VectorField"
 *   and "VectorFieldExact" are modeled by exactly the same finite elements,
 *   i.e., FE_Nedelec .
 *
 * The "L2norm", and "VectorFieldExact" are saved only if an exact solution is
 * passes to the constructor.  Moreover, "VectorFieldExact" is calculated and
 * saved only if the input parameter project_exact_solution passed to the
 * constructor equals true.
 *
 * The name of the output file is computed by appending ".vtk" to the string
 * contained by the input parameter fname of the constructor. The vtk file
 * can be viewed with a help of
 * [Visit](https://visit.llnl.gov)
 * software package of the Lawrence Livermore National Laboratory.
 *
 * The user is supposed to derive an object from this class template. All usual
 * computations, i.e., setup, assembling the linear system, etc., happen
 * automatically.
 *
 * @note Application examples:
 * - [mms/](@ref page_mms), [int/](@ref page_int) - Projection nr. 1 and
 *   projection nr. 3 (planar).
 * - [mms-axi/](@ref page_mms_axi), [int-axi/](@ref page_int_axi) - Projection
 *   nr. 3 (axisymmetric).
 *****************************************************************************/
template <int dim, int stage = 1>
class ProjectPHItoE : public ProjectHgradToHcurl<dim, stage>
{
	public:
		ProjectPHItoE() = delete;
/**
 * \brief The only constructor.
 *
 * Constructs the object and runs the calculations. That is, there is no need
 * to call other functions such as run().
 * @param[in] p - The degree of the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements.
 * @param[in] mapping_degree - The degree of the interpolating polynomials used
 * for mapping. Setting it to 1 will do in the most of the cases. Note, that it
 * makes sense to attach a meaningful manifold to the triangulation if this
 * parameter is greater than 1.
 * @param[in] triangulation_PHI - Reference to the triangulation inside the
 * object of the class StaticScalarSolver::Solver that has yielded the
 * scalar electric potential, \f$\Phi\f$.
 * @param[in] dof_handler_PHI - Reference to the dof handler inside the class
 * StaticScalarSolver::Solver that has yielded the scalar electric
 * potential, \f$\Phi\f$.
 * @param[in] solution_PHI - Vector filled with the degrees of freedom that
 * together with the shape functions of the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements model
 * the total scalar electric potential, \f$\Phi\f$.
 * @param[in] fname - The name of the output files without extension.
 * @param[in] exact_solution - Points to an object that describes the exact
 * solution to the problem. It is needed for calculating error norms. The
 * exact solution object must exist in the computer memory at the time of
 * calling this constructor.
 * @param[in] axisymmetric - If true, assumes that the problem is
 * axisymmetric. In this case, of course, the problem must be two- dimensional.
 * @param[in] print_time_tables - If true, prints time tables on the screen.
 * @param[in] project_exact_solution - If true, projects the exact solution
 * onto the space spanned by the Nedelec finite elements and saves the
 * result into the vtk file next to \f$\vec{E}\f$.
 * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
 * gradient solver into a file. The name of the file is generated by appending
 * "_cg_convergence.csv" to fname.
 *****************************************************************************/
		ProjectPHItoE(
			unsigned int p,
			unsigned int mapping_degree,
			const Triangulation<dim> & triangulation_PHI,
			const DoFHandler<dim> & dof_handler_PHI,
			const Vector<double> & solutioin_PHI,
			std::string fname = "E",
			const Function<dim> * exact_solution = nullptr,
			bool axisymmetric = false,
			bool print_time_tables = false,
			bool project_exact_solution = false,
			bool log_cg_convergence = false) :
				ProjectHgradToHcurl<dim,stage>(
				p,
				mapping_degree,
				triangulation_PHI,
				dof_handler_PHI,
				solutioin_PHI,
				fname,
				exact_solution,
				axisymmetric,
				false,
				print_time_tables,
				project_exact_solution,
				log_cg_convergence){}
};

} // namespace StaticScalarSolver

#endif

