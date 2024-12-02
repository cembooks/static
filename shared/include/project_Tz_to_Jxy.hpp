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

#ifndef ProjectTzToJxy_H__
#define ProjectTzToJxy_H__

#include "project_Hgrad_to_Hdiv.hpp"

namespace StaticScalarSolver {
/**
 * \brief Computes the two-dimensional free-current density \f$\vec{J}_f\f$
 * as a vector curl of the current vector potential, \f$T\f$, i.e.,
 * \f$\vec{J}_f = \vec{\nabla} \overset{V}{\times} T\f$.
 *
 * The problem is assumed to be planar two-dimensional. The free-current density
 * \f$\vec{J}_f\f$ is an in-plane vector. The current vector potential, \f$T\f$,
 * is an out-of-plane vector. This class template is meant to be used in
 *conjunction with StaticScalarSolver::Solver.
 *
 * The current vector potential, \f$T\f$, is assumed to be an output of
 * an object of the class StaticScalarSolver::Solver. That is, \f$T\f$ is
 * expressed as a linear combination of the shape functions of the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements. In other words, \f$T\f$ is in the \f$H(\text{grad})\f$
 * function space. This class template computes free-current density
 * \f$\vec{J}_f\f$ such that \f$\vec{J}_f\f$ is expressed as a linear
 *combination of the shape functions of the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{J}_f\f$ is in the \f$H(\text{div})\f$
 * function space. In other words, this class template takes a scalar function
 * \f$T\f$ (an out-of-plane vector is, essentially, a scalar) from the
 * \f$H(\text{grad})\f$ function space, calculates its vector curl and projects
 * the result into \f$H(\text{div})\f$ function space, which is a proper space
 * for the free-current density \f$\vec{J}_f\f$. The operation of projection is
 * illustrated by means of the
 * [Bossavit's diagrams](@ref sec_bossavit_dia2)
 * below.
 *
 *![](sss_prj/diagram_prj_Tz_to_Jxy.svg)
 *
 * It is assumed that the object of the class StaticScalarSolver::Solver that
 * yielded \f$T\f$ is still in the computer memory such that the
 * triangulation, the degrees of freedom, and the handler of the degrees of
 * freedom are accessible while \f$\vec{J}_f\f$ is computed. An object of
 * the StaticScalarSolver::ProjectTzToJxy class template reuses the
 *triangulation by creating an additional dof handler associated with the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{J}_f\f$ and \f$T\f$ share the same
 *triangulation but are modeled by different finite elements: \f$\vec{J}_f\f$ -
 *by
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * and \f$T\f$ - by
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
 *
 * MultithreadInfo::set_thread_limit(nr_threads_max);
 * @endcode
 *
 * The output data are saved into a vtk file. The following data are saved:
 * - The calculated free-current density, \f$\vec{J}_f = \vec{\nabla}
 *\overset{V}{\times} T\f$, under the name "VectorField".
 * - The \f$L^2\f$ error norm associated with the calculated displacement
 *   under the name "L2norm". One value per mesh cell is saved.
 * - The exact solution projected onto \f$H(\text{div})\f$ function space is
 *   saved under the name "VectorFieldExact". The "VectorField"
 *   and "VectorFieldExact" are modeled by exactly the same finite elements,
 *   i.e., FE_RaviartThomas.
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
 * @note Application examples:
 *
 * - [cvp-ii](@ref page_cvp_ii),
 * [mms-vt-ii](@ref page_mms_vt_ii).
 *****************************************************************************/
template<int stage = 1>
class ProjectTzToJxy : public ProjectHgradToHdiv<2, stage>
{
public:
  ProjectTzToJxy() = delete;
  /**
   * \brief The only constructor.
   *
   * Constructs the object and runs the calculations. That is, there is no need
   * to call other functions such as run().
   *
   * @param[in] p - The degree of the
   * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
   * finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   *used for mapping. Setting it to 1 will do in the most of the cases. Note,
   *that it makes sense to attach a meaningful manifold to the triangulation if
   *this parameter is greater than 1.
   * @param[in] triangulation_T - Reference to the triangulation inside an
   * object of the class StaticScalarSolver::Solver that has yielded the
   * current vector potential, \f$T\f$.
   * @param[in] dof_handler_T - Reference to the dof handler inside the class
   * StaticScalarSolver::Solver that has yielded the current vector
   * potential, \f$T\f$.
   * @param[in] solution_T - Vector filled with the degrees of freedom that
   * together with the shape functions of the
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
   * finite elements model the current vector potential, \f$T\f$.
   * @param[in] fname - The name of the output files without extension.
   * @param[in] exact_solution - Points to an object that describes the exact
   * solution to the problem. It is needed for calculating error norms. The
   * exact solution object must exist in the computer memory at the time of
   * calling this constructor.
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the Raviart-Thomas finite elements and saves the
   * result into the vtk file next to \f$\vec{J}_f\f$.
   * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
   * gradient solver into a file. The name of the file is generated by appending
   * "_cg_convergence.csv" to fname.
   *****************************************************************************/
  ProjectTzToJxy(unsigned int p,
                 unsigned int mapping_degree,
                 const Triangulation<2>& triangulation_T,
                 const DoFHandler<2>& dof_handler_T,
                 const Vector<double>& solutioin_T,
                 std::string fname = "Jxy",
                 const Function<2>* exact_solution = nullptr,
                 bool print_time_tables = false,
                 bool project_exact_solution = false,
                 bool log_cg_convergence = false)
    : ProjectHgradToHdiv<2, stage>(p,
                                   mapping_degree,
                                   triangulation_T,
                                   dof_handler_T,
                                   solutioin_T,
                                   fname,
                                   exact_solution,
                                   false, // axisymmetric
                                   true,  // vector potential
                                   print_time_tables,
                                   project_exact_solution,
                                   log_cg_convergence)
  {
  }
};

} // namespace StaticScalarSolver

#endif
