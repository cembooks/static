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

#ifndef ProjectAtoB_H__
#define ProjectAtoB_H__

#include "project_Hcurl_to_Hdiv.hpp"

namespace StaticVectorSolver {
/**
 * \brief Computes the magnetic field \f$\vec{B}\f$ as a curl of the magnetic
 * vector potential, \f$\vec{A}\f$,
 * i.e., \f$\vec{B} = \vec{\nabla} \times \vec{A}\f$.
 *
 * The problem is assumed to be three-dimensional. This class template
 * is meant to be used in conjunction with StaticVectorSolver::Solver1 or
 * StaticVectorSolver::Solver2.
 *
 * The magnetic vector potential, \f$\vec{A}\f$, is assumed to be an output of
 * an object of the class
 * StaticVectorSolver::Solver1 or StaticVectorSolver::Solver2 .  That is,
 * \f$\vec{A}\f$ is expressed as a linear combination of the shape functions
 * of the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. In other words, \f$\vec{A}\f$ is in the \f$H(\text{curl})\f$
 * function space. This class template computes magnetic field \f$\vec{B}\f$
 * as a curl of the magnetic vector potential, \f$\vec{A}\f$, such
 * that \f$\vec{B}\f$ is expressed as a linear combination of the shape
 * functions of the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{B}\f$ is in the \f$H(\text{div})\f$
 * function space. In other words, this class template takes a vector field
 * \f$\vec{A}\f$ from the \f$H(\text{curl})\f$ function space, calculates its
 * curl, and projects the result into \f$H(\text{div})\f$ function space, which
 * is a proper space for the magnetic field \f$\vec{B}\f$. The operation of
 * projection is illustrated by  means of the Bossavit's diagram below.
 *
 * ![](sss_prj/diagram_prj_A_to_B.svg)
 *
 * It is assumed that the object of the class StaticVectorSolver::Solver1 or
 * StaticVectorSolver::Solver2 that yielded \f$\vec{A}\f$ is still in the
 * computer memory such that the triangulation, the degrees of freedom, and the
 * handler of the degrees of freedom are accessible while \f$\vec{B}\f$ is
 * computed. An object of the StaticVectorSolver::ProjectAtoB class template
 * reuses the triangulation by creating an additional dof handler associated
 * with the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{A}\f$ and \f$\vec{B}\f$ share the same
 * triangulation but are modeled by different
 * finite elements: \f$\vec{A}\f$ - by
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * and
 * \f$\vec{B}\f$ - by
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html).
 * In this disposition there are two DoFHandler objects associated with one
 * Triangulation object. The algorithm walks synchronously through both
 * DoFHandler objects and assembles the system matrix and the right-hand side.
 *
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
 * - The calculated magnetic field, \f$\vec{B} = \vec{\nabla} \times \vec{A}\f$,
 *   under the name "VectorField".
 * - The \f$L^2\f$ error norm associated with the calculated magnetic field
 *   under the name "L2norm". One value per mesh cell is saved.
 * - The exact solution projected onto \f$H(\text{div})\f$ function space is
 *   saved under the name "VectorFieldExact". The "VectorField"
 *   and "VectorFieldExact" are modeled by exactly the same finite elements,
 *   i.e., FE_RaviartThomas.
 *
 * The "L2norm", and "VectorFieldExact" are saved only if an exact solution is
 * passes to the constructor. Moreover, "VectorFieldExact" is calculated and
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
 * - [mms-v/](@ref page_mms_v),
 * [mms-vt-i/](@ref page_mms_vt_i),
 * [ssol-i/](@ref page_ssol_i),
 * [ssol-ii/](@ref page_ssol_ii),
 * [ssol-iii/](@ref page_ssol_iii).
 *****************************************************************************/
template<int stage = 1>
class ProjectAtoB : public ProjectHcurlToHdiv<stage>
{
public:
  ProjectAtoB() = delete;

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
   * used for mapping. Setting it to 1 will do in the most of the cases. Note,
   * that it makes sense to attach a meaningful manifold to the triangulation if
   * this parameter is greater than 1.
   * @param[in] triangulation_A - Reference to the triangulation inside an
   * object of the class StaticVectorSolver::Solver1 or
   * StaticVectorSolver::Solver2 that has yielded the magnetic vector potential,
   * \f$\vec{A}\f$.
   * @param[in] dof_handler_A - Reference to the dof handler inside the class
   * StaticVectorSolver::Solver1 or StaticVectorSolver::Solver2that has yielded
   * the magnetic vector potential, \f$\vec{A}\f$.
   * @param[in] solution_A - Vector filled with the degrees of freedom that
   * together with the shape functions of the
   * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
   * finite elements model the magnetic vector potential, \f$\vec{A}\f$.
   * @param[in] fname - The name of the output files without extension.
   * @param[in] exact_solution - Points to an object that describes the exact
   * solution to the problem. It is needed for calculating error norms. The
   * exact solution object must exist in the computer memory at the time of
   * calling this constructor.
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the FE_Nedelec finite elements and saves the
   * result into the vtk file next to \f$\vec{B}\f$.
   * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
   * gradient solver into a file. The name of the file is generated by appending
   * "_cg_convergence.csv" to fname.
   * @param[in] write_higher_order_cells - If false, the data is saved in the
   * file fname.vtk. Higher-order cells are not saved. If true, the data is
   * saved into fname.vtu file preserving the higher-order cells. In this case
   * the file can be viewed with a help of
   * [Paraview](www.paraview.org) version 5.5.0 or higher.
   *****************************************************************************/
  ProjectAtoB(unsigned int p,
              unsigned int mapping_degree,
              const Triangulation<3>& triangulation_A,
              const DoFHandler<3>& dof_handler_A,
              const Vector<double>& solutioin_A,
              std::string fname = "B",
              const Function<3>* exact_solution = nullptr,
              bool print_time_tables = false,
              bool project_exact_solution = false,
              bool log_cg_convergence = false,
              bool write_higher_order_cells = false)
    : ProjectHcurlToHdiv<stage>(p,
                                mapping_degree,
                                triangulation_A,
                                dof_handler_A,
                                solutioin_A,
                                fname,
                                exact_solution,
                                print_time_tables,
                                project_exact_solution,
                                log_cg_convergence,
                                write_higher_order_cells)
  {
  }
};

} // namespace StaticScalarSolver

#endif
