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

#ifndef ProjectPsiToH_H__
#define ProjectPsiToH_H__

#include "project_Hgrad_to_Hcurl.hpp"

namespace StaticScalarSolver {
/**
 * \brief Computes the auxiliary vector field \f$\vec{H}\f$ as a negative
 * gradient of the total scalar magnetic potential, \f$\Psi\f$, i.e.,
 * \f$\vec{H} = - \vec{\nabla} \Psi\f$.
 *
 * The problem is assumed to be two- or three- dimensional. This class template
 * is meant to be used in conjunction with StaticScalarSolver::Solver.
 *
 * The total scalar magnetic potental, \f$\Psi\f$, is assumed to be an output of
 * an object of the class StaticScalarSolver::Solver. That is, \f$\Psi\f$ is
 * expressed as a linear combination of the shape functions of the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements. In other words, \f$\Psi\f$ is in the \f$H(\text{grad})\f$
 * function space. This class template computes auxiliary vector field
 * \f$\vec{H}\f$ as a negative gradient of the total scalar magnetic potential,
 * \f$\Psi\f$, such that \f$\vec{H}\f$ is expressed as a linear combination
 * of the shape functions of the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. That is, \f$\vec{H}\f$ is in the \f$H(\text{curl})\f$
 * function space. In other words, this class template takes a scalar function
 * \f$\Psi\f$ from the \f$H(\text{grad})\f$ function space, calculates its
 * gradient, multiplies the gradient by \f$-1\f$, and projects the result into
 * \f$H(\text{curl})\f$ function space, which is a proper space for the
 * auxiliary vector field \f$\vec{H}\f$. The operation of projection is
 * illustrated by means of the Bossavit's diagrams below.
 *
 * ![](sss_prj/diagram_prj_psi_to_h.svg)
 *
 * It is assumed that the object of the class StaticScalarSolver::Solver that
 * yielded \f$\Psi\f$ is still in the computer memory such that the
 * triangulation, the degrees of freedom, and the handler of the degrees of
 * freedom are accessible while \f$\vec{H}\f$ is computed. An object of
 * the StaticScalarSolver::ProjectPSItoH class template reuses the triangulation
 * by creating an additional dof handler associated with the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. That is, \f$\vec{H}\f$ and \f$\Psi\f$ share the same
 * triangulation but are modeled by different
 * finite elements: \f$\vec{H}\f$ - by
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * and \f$\Psi\f$ - by
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
 * - The calculated field, \f$\vec{H} = - \vec{\nabla} \Psi\f$,
 *   under the name "VectorField".
 * - The \f$L^2\f$ and \f$L^{\infty}\f$ error norms associated with the
 *   calculated \f$\vec{H}\f$ field under the names "L2norm" and "LinftyNorm".
 *   One value per mesh cell is saved.
 * - The exact solution projected onto \f$H(\text{curl})\f$ function space is
 *   saved under the name "VectorFieldExact". The "VectorField"
 *   and "VectorFieldExact" are modeled by exactly the same finite elements,
 *   i.e., FE_Nedelec .
 *
 * The "L2norm", "LinftyNorm", and "VectorFieldExact" are saved only if an
 * exact solution is passed to the constructor.  Moreover, "VectorFieldExact"
 * is calculated and saved only if the input parameter `project_exact_solution`
 * passed to the constructor equals true.
 *
 * The user is supposed to derive an object from this class template. All usual
 * computations, i.e., setup, assembling the linear system, etc., happen
 * automatically.
 *
 * @note Application examples:
 * - [sld-i/](@ref page_sld_i) - Projection nr. 2 and
 *   projection nr. 4 (planar).
 *****************************************************************************/
template<int dim, int stage = 1>
class ProjectPSItoH : public ProjectHgradToHcurl<dim, stage>
{
public:
  ProjectPSItoH() = delete;
  /**
   * \brief The only constructor.
   *
   * Constructs the object and runs the calculations. That is, there is no need
   * to call other functions such as run().
   * @param[in] p - The degree of the
   * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
   * finite elements.
   * @param[in] mapping_degree - The degree of the interpolating polynomials
   * used for mapping. Setting it to 1 will do in the most of the cases. Note,
   * that it makes sense to attach a meaningful manifold to the triangulation if
   * this parameter is greater than 1.
   * @param[in] triangulation_PSI - Reference to the triangulation inside an
   * object of the class StaticScalarSolver::Solver that has yielded the
   * total scalar magnetic potential, \f$\Psi\f$.
   * @param[in] dof_handler_PSI - Reference to the dof handler inside the class
   * StaticScalarSolver::Solver that has yielded the total scalar magnetic
   * potential, \f$\Psi\f$.
   * @param[in] solution_PSI - Vector filled with the degrees of freedom that
   * together with the shape functions of the
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
   * finite elements model
   * the total scalar magnetic potential, \f$\Psi\f$.
   * @param[in] fname - The name of the output files without extension.
   * @param[in] exact_solution - Points to an object that describes the exact
   * solution to the problem. It is needed for calculating error norms. The
   * exact solution object must exist in the computer memory at the time of
   * calling this constructor.
   * @param[in] axisymmetric - If true, assumes that the problem is
   * axisymmetric. In this case, of course, the problem must be two-
   * dimensional.
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the Nedelec finite elements and saves the
   * result into the output file next to \f$\vec{H}\f$.
   * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
   * gradient solver into a file. The name of the file is generated by appending
   * "_cg_convergence.csv" to fname.
   * @param[in] write_higher_order_cells - If false, the data is saved in the
   * file fname.vtk. Higher-order cells are not saved. If true, the data is
   * saved into fname.vtu file preserving the higher-order cells. In this case
   * the file can be viewed with a help of
   * [Paraview](www.paraview.org) version 5.5.0 or higher.
   *****************************************************************************/
  ProjectPSItoH(unsigned int p,
                unsigned int mapping_degree,
                const Triangulation<dim>& triangulation_PSI,
                const DoFHandler<dim>& dof_handler_PSI,
                const Vector<double>& solution_PSI,
                std::string fname = "H",
                const Function<dim>* exact_solution = nullptr,
                bool axisymmetric = false,
                bool print_time_tables = false,
                bool project_exact_solution = false,
                bool log_cg_convergence = false,
                bool write_higher_order_cells = false)
    : ProjectHgradToHcurl<dim, stage>(p,
                                      mapping_degree,
                                      triangulation_PSI,
                                      dof_handler_PSI,
                                      solution_PSI,
                                      fname,
                                      exact_solution,
                                      axisymmetric,
                                      false,
                                      print_time_tables,
                                      project_exact_solution,
                                      log_cg_convergence,
                                      write_higher_order_cells)
  {
  }
};

} // namespace StaticScalarSolver

#endif
