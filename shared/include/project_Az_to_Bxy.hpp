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

#ifndef ProjectAzToBxy_H__
#define ProjectAzToBxy_H__

#include "project_Hgrad_to_Hdiv.hpp"

namespace StaticScalarSolver {
/**
 * \brief Computes the two-dimensional magnetic field \f$\vec{B}\f$
 * as a vector curl of the magnetic vector potential, \f$A\f$, i.e.,
 * \f$\vec{B} = \vec{\nabla} \overset{V}{\times} A\f$.
 *
 * The problem is assumed to be planar two-dimensional. The magnetic field
 * \f$\vec{B}\f$ is an in-plane vector. The magnetic vector potential, \f$A\f$,
 * is an out-of-plane vector.
 *
 * The magnetic vector potential, \f$A\f$, is assumed to be an output of
 * an object of the class StaticScalrSolver::Solver. That is, \f$A\f$ is
 * expressed as a linear combination of the shape functions of the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements. In other words, \f$A\f$ is in the \f$H(\text{grad})\f$
 * function space. This class template computes magnetic field \f$\vec{B}\f$
 * such that \f$\vec{B}\f$ is expressed as a linear combination of the shape
 * functions of the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{B}\f$ is in the \f$H(\text{div})\f$
 * function space. In other words, this class template takes a scalar function
 * \f$A\f$ (an out-of-plane vector is, essentially, a scalar) from the
 * \f$H(\text{grad})\f$ function space, calculates its vector curl and projects
 * the result into \f$H(\text{div})\f$ function space, which is a proper space
 * for the magnetic field \f$\vec{B}\f$. The operation of projection is
 * illustrated by means of the Bossavit's diagram below.
 *
 *![](sss_prj/diagram_prj_Az_to_Bxy.svg)
 *
 * It is assumed that the object of the class StaticScalarSolver::Solver that
 * yielded \f$A\f$ is still in the computer memory such that the
 * triangulation, the degrees of freedom, and the handler of the degrees of
 * freedom are accessible while \f$\vec{B}\f$ is computed. An object of
 * the StaticScalarSolver::ProjectAzToBxy class template reuses the
 *triangulation by creating an additional dof handler associated with the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{B}\f$ and \f$A\f$ share the same
 *triangulation but are modeled by different finite elements: \f$\vec{B}\f$ - by
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * and \f$A\f$ - by
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
 * - The calculated magnetic field, \f$\vec{B} = \vec{\nabla}
 *\overset{V}{\times}A\f$, under the name "VectorField".
 * - The \f$L^2\f$ error norm associated with the calculated magnetic field
 *   under the name "L2norm". One value per mesh cell is saved.
 * - The exact solution projected onto \f$H(\text{div})\f$ function space is
 *   saved under the name "VectorFieldExact". The "VectorField"
 *   and "VectorFieldExact" are modeled by exactly the same finite elements,
 *   i.e., FE_RaviartThomas.
 *
 * The "L2norm", and "VectorFieldExact" are saved only if an exact solution is
 * passed to the constructor.  Moreover, "VectorFieldExact" is calculated and
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
 * - [mwr/](@ref page_mwr)
 *****************************************************************************/
template<int stage = 1>
class ProjectAzToBxy : public ProjectHgradToHdiv<2, stage>
{
public:
  ProjectAzToBxy() = delete;
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
   * @param[in] triangulation_A - Reference to the triangulation inside an
   * object of the class StaticScalarSolver::Solver that has yielded the
   * magnetic vector potential, \f$A\f$.
   * @param[in] dof_handler_A - Reference to the dof handler inside the class
   * StaticScalarSolver::Solver that has yielded the magnetic vector
   * potential, \f$A\f$.
   * @param[in] solution_A - Vector filled with the degrees of freedom that
   * together with the shape functions of the
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
   * finite elements model the magnetic vector potential, \f$A\f$.
   * @param[in] fname - The name of the output files without extension.
   * @param[in] exact_solution - Points to an object that describes the exact
   * solution to the problem. It is needed for calculating error norms. The
   * exact solution object must exist in the computer memory at the time of
   * calling this constructor.
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the Raviart-Thomas finite elements and saves the
   * result into the vtk file next to \f$\vec{B}\f$.
   * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
   * gradient solver into a file. The name of the file is generated by appending
   * "_cg_convergence.csv" to fname.
   *****************************************************************************/
  ProjectAzToBxy(unsigned int p,
                 unsigned int mapping_degree,
                 const Triangulation<2>& triangulation_A,
                 const DoFHandler<2>& dof_handler_A,
                 const Vector<double>& solutioin_A,
                 std::string fname = "Bxy",
                 const Function<2>* exact_solution = nullptr,
                 bool print_time_tables = false,
                 bool project_exact_solution = false,
                 bool log_cg_convergence = false)
    : ProjectHgradToHdiv<2, stage>(p,
                                   mapping_degree,
                                   triangulation_A,
                                   dof_handler_A,
                                   solutioin_A,
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
