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

#ifndef ProjectTtoJ_H__
#define ProjectTtoJ_H__

#include "project_Hcurl_to_Hdiv.hpp"

using namespace StaticVectorSolver;

namespace StaticVectorSolver {
/**
 * \brief Computes the volume free-current density \f$\vec{J}_f\f$ as a curl
 * of the current vector potential, \f$\vec{T}\f$,
 * i.e., \f$\vec{J}_f = \vec{\nabla} \times \vec{T}\f$.
 *
 * The problem is assumed to be three-dimensional. This class template
 * is meant to be used in conjunction with StaticVectorSolver::Solver1.
 *
 * The current vector potential, \f$\vec{T}\f$, is assumed to be an output of
 * an object of the class
 * StaticVectorSolver::Solver1. That is, \f$\vec{T}\f$ is expressed as a
 * linear combination of the shape functions of the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. In other words, \f$\vec{T}\f$ is in the \f$H(\text{curl})\f$
 * function space. This class template computes volume free-current density
 * \f$\vec{J}_f\f$ as a curl of the current vector potential, \f$\vec{T}\f$,
 * such that \f$\vec{J}_f\f$ is expressed as a linear combination of the shape
 * functions of the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{J}_f\f$ is in the \f$H(\text{div})\f$
 * function space. In other words, this class template takes a vector field
 * \f$\vec{T}\f$ from the \f$H(\text{curl})\f$ function space, calculates its
 * curl, and projects the result into \f$H(\text{div})\f$ function space, which
 * is a proper space for the volume free-current density \f$\vec{J}_f\f$. The
 * operation of projection is illustrated by  means of the Bossavit's diagram
 * below.
 *
 * ![](sss_prj/diagram_prj_T_to_J.svg)
 *
 * It is assumed that the object of the class StaticVectorSolver::Solver1
 * that yielded \f$\vec{T}\f$ is still in the computer
 * memory such that the triangulation, the degrees of freedom, and the handler
 * of the degrees of freedom are accessible while \f$\vec{J}_f\f$ is computed.
 * An object of the StaticVectorSolver::ProjectTtoJ class template reuses
 * the triangulation by creating an additional dof handler associated with the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements. That is, \f$\vec{T}\f$ and \f$\vec{J}_f\f$ share the same
 * triangulation but are modeled by different finite elements:
 * \f$\vec{T}\f$ - by
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * and
 * \f$\vec{J}_f\f$ - by
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
 * ...
 * MultithreadInfo::set_thread_limit(nr_threads_max);
 * @endcode
 *
 * The output data are saved into a vtk file. The following data are saved:
 * - The calculated volume free-current density,
 *   \f$\vec{J}_f = \vec{\nabla} \times \vec{T}\f$, under the name
 *   "VectorField".
 * - The \f$L^2\f$ error norm associated with the calculated volume free-current
 *   density under the name "L2norm". One value per mesh cell is saved.
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
 * - [cvp-i/](@ref page_cvp_i),
 * [mms-vt-i/](@ref page_mms_vt_i),
 * [ssol-ii/](@ref page_ssol_ii),
 * [ssol-iii/](@ref page_ssol_iii).
 *****************************************************************************/
template<int stage = 1>
class ProjectTtoJ : public ProjectHcurlToHdiv<stage>
{
public:
  ProjectTtoJ() = delete;

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
   * @param[in] triangulation_T - Reference to the triangulation inside an
   * object of the class StaticVectorSolver::Solver1 that has yielded the
   * magnetic vector potential, \f$T\f$.
   * @param[in] dof_handler_T - Reference to the dof handler inside the class
   * StaticVectorSolver::Solver1 that has yielded the magnetic vector
   * potential, \f$T\f$.
   * @param[in] solution_T - Vector filled with the degrees of freedom that
   * together with the shape functions of the
   * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
   * finite elements model the magnetic vector potential, \f$T\f$.
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
   * @param[in] write_higher_order_cells - If false, the data is saved in the
   * file fname.vtk. Higher-order cells are not saved. If true, the data is
   * saved into fname.vtu file preserving the higher-order cells. In this case
   * the file can be viewed with a help of
   * [Paraview](www.paraview.org) version 5.5.0 or higher.
   *****************************************************************************/
  ProjectTtoJ(unsigned int p,
              unsigned int mapping_degree,
              const Triangulation<3>& triangulation_T,
              const DoFHandler<3>& dof_handler_T,
              const Vector<double>& solutioin_T,
              std::string fname = "J",
              const Function<3>* exact_solution = nullptr,
              bool print_time_tables = false,
              bool project_exact_solution = false,
              bool log_cg_convergence = false,
              bool write_higher_order_cells = false)
    : ProjectHcurlToHdiv<stage>(p,
                                mapping_degree,
                                triangulation_T,
                                dof_handler_T,
                                solutioin_T,
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
