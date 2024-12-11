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

#ifndef ProjectAxyToBz_H__
#define ProjectAxyToBz_H__

#include "project_Hcurl_to_L2.hpp"

namespace StaticVectorSolver {
/**
 * \brief Computes an out-of-plane magnetic field, \f$B\f$, as a scalar curl
 * of an in-plane magnetic vector potential, \f$\vec{A}\f$,
 * i.e., \f$B = \vec{\nabla} \overset{S}{\times} \vec{A}\f$.
 *
 * The problem is assumed to be planar two-dimensional. In this configuration
 * the magnetic vector potential, \f$\vec{A}\f$, is an in-plane vector and the
 * magnetic field, \f$B\f$, is an out-of-plane vector. That is, in the initial
 * three-dimensional problem form which the two-dimensional problem is derived
 * the magnetic field has only \f$z\f$ component and the magnetic vector
 * potential has only \f$x\f$ and \f$y\f$ components.
 *
 * This class template
 * is meant to be used in conjunction with StaticVectorSolver::Solver1 or
 * StaticVectorSolver::Solver2.
 *
 * The magnetic vector potential, \f$\vec{A}\f$, is assumed to be an output of
 * an object of the class
 * StaticVectorSolver::Solver1 or StaticVectorSolver::Solver2.  That is,
 * \f$\vec{A}\f$ is expressed as a linear combination of the shape functions
 * of the
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * finite elements. In other words, \f$\vec{A}\f$ is in the \f$H(\text{curl})\f$
 * function space. This class template computes magnetic field \f$B\f$
 * as a scalar curl of the magnetic vector potential, \f$\vec{A}\f$, such
 * that \f$B\f$ is expressed as a linear combination of the shape functions of
 * the
 * [FE_DGQ](https://www.dealii.org/current/doxygen/deal.II/classFE__DGQ.html)
 * finite elements. That is, \f$B\f$ is in the \f$L^2\f$
 * function space. In other words, this class template takes a vector field
 * \f$\vec{A}\f$ from the \f$H(\text{curl})\f$ function space, calculates its
 * scalar curl, and projects the result into \f$L^2\f$ function space, which is
 * a proper space for the out-of-plane magnetic field \f$B\f$. The operation of
 * projection is illustrated by  means of the Bossavit's diagram below.
 *
 * ![](svs_prj/diagram_prj_Axy_to_Bz.svg)
 *
 * It is assumed that the object of the class StaticVectorSolver::Solver1 or
 * StaticVectorSolver::Solver2 that yielded \f$\vec{A}\f$ is still in the
 * computer memory such that the triangulation, the degrees of freedom, and the
 * handler of the degrees of freedom are accessible while \f$B\f$ is computed.
 * An object of the StaticVectorSolver::ProjectAxyToBz class template reuses
 * the triangulation by creating an additional dof handler associated with the
 * [FE_DGQ](https://www.dealii.org/current/doxygen/deal.II/classFE__DGQ.html)
 * finite elements. That is, \f$\vec{A}\f$ and \f$B\f$ share the same
 * triangulation but are modeled by different
 * finite elements: \f$\vec{A}\f$ - by
 * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
 * and
 * \f$B\f$ - by
 * [FE_DGQ](https://www.dealii.org/current/doxygen/deal.II/classFE__DGQ.html)
 * In this disposition there are two DoFHandler objects associated with one
 * Triangulation object. The algorithm walks synchronously through both
 * DoFHandler objects and assembles the system matrix and the right-hand side.
 *
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
 * The output data are saved into a vtk file. The following data are saved:
 * - The calculated magnetic field, \f$B = \vec{\nabla}\overset{S}{\times}
 *   \vec{A}\f$, under the name "ScalarField".
 * - The \f$L^2\f$ and \f$L^{\infty}\f$ error norms associated with the
 *   calculated electrostatic field under the names "L2norm" and "LinftyNorm".
 *   One value per mesh cell is saved.
 * - The exact solution projected onto \f$L^2\f$ function space is
 *   saved under the name "ScalarFieldExact". The "ScalarField"
 *   and "ScalarFieldExact" are modeled by exactly the same finite elements,
 *   i.e., FE_DGQ.
 *
 * The "L2norm", "LinftyNorm", and "VectorFieldExact" are saved only if an
 * exact solution is passed to the constructor.  Moreover, "VectorFieldExact"
 * is calculated and saved only if the input parameter
 * <code>project_exact_solution</code> passed to the constructor equals true.
 *
 * The user is supposed to derive an object from this class template. All usual
 * computations, i.e., setup, assembling the linear system, etc., happen
 * automatically.
 *
 * @note Application examples:
 *
 * - [mms-v/](@ref page_mms_v),
 * [mms-vt-ii/](@ref page_mms_vt_ii).
 *
 *****************************************************************************/
template<int stage = 1>
class ProjectAxyToBz : public ProjectHcurlToL2<stage>
{
public:
  ProjectAxyToBz() = delete;

  /**
   * \brief The only constructor.
   *
   * Constructs the object and runs the calculations. That is, there is no need
   * to call other functions such as run().
   *
   * @param[in] p - The degree of the
   * [FE_DGQ](https://www.dealii.org/current/doxygen/deal.II/classFE__DGQ.html)
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
   * StaticVectorSolver::Solver1 or StaticVectorSolver::Solver2 that has yielded
   * the magnetic vector potential, \f$\vec{A}\f$.
   * @param[in] solution_A - Vector filled with the degrees of freedom that
   * together with the shape functions of the
   * [FE_Nedelec](https://www.dealii.org/current/doxygen/deal.II/classFE__Nedelec.html)
   * finite elements model the magnetic vector potential, \f$\vec{A}\f$.
   * @param[in] fname - The name of the output files without extension.
   * @param[in] exact_solution - Points to an object that describes the exact
   * \f$B\f$. It is needed for calculating error norms. The
   * exact solution object must exist in the computer memory at the time of
   * calling this constructor.
   * @param[in] print_time_tables - If true, prints time tables on the screen.
   * @param[in] project_exact_solution - If true, projects the exact solution
   * onto the space spanned by the FE_DGQ finite elements and saves the
   * result into the vtk file next to \f$B\f$.
   * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
   * gradient solver into a file. The name of the file is generated by appending
   * "_cg_convergence.csv" to fname
   * @param[in] write_higher_order_cells - If false, the data is saved in the
   * file fname.vtk. Higher-order cells are not saved. If true, the data is
   * saved into fname.vtu file preserving the higher-order cells. In this case
   * the file can be viewed with a help of
   * [Paraview](www.paraview.org) version 5.5.0 or higher.
   *****************************************************************************/
  ProjectAxyToBz(unsigned int p,
                 unsigned int mapping_degree,
                 const Triangulation<2>& triangulation_A,
                 const DoFHandler<2>& dof_handler_A,
                 const Vector<double>& solutioin_A,
                 std::string fname = "Bz",
                 const Function<2>* exact_solution = nullptr,
                 bool print_time_tables = false,
                 bool project_exact_solution = false,
                 bool log_cg_convergence = false,
                 bool write_higher_order_cells = false)
    : ProjectHcurlToL2<stage>(p,
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
