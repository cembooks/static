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

#ifndef StaticVectorInput_H__
#define StaticVectorInput_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

namespace StaticVectorSolver {

/**
 * \brief Implements the permeability, \f$\mu\f$, in the partial differential
 * equation (i) of the [vector boundary value problem](@ref veibvp_bvp).
 *
 *
 * If the magnetic vector potential \f$\vec{A}\f$ is being computed, the
 * coefficient must equal the permeability, \f$\mu\f$. If the current vector
 * potential \f$\vec{T}\f$ is being computed, the coefficient must equal 1.
 *
 * This class is declared in shared/include/static_vector_input.hpp but
 * must be implemented in xyz/src/static_vector_input.cpp, where xyz is
 * the directory of the current numerical experiment. That is, the declaration
 * is shared between all numerical experiments while implementation is specific
 * for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The `dim` template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the `stage` template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class.
 *****************************************************************************/
template<int dim, int stage = 1>
class TheCoefficient : private Settings
{
public:
  /**
   * \brief Computes values of permeability, \f$\mu\f$, at quadrature
   * points.
   *
   * This function is called by the StaticVectorSolver::Solver1 and
   * StaticVectorSolver::Solver2 classes during the  assembly of the system
   * matrix. It is called ones per each cell. This function must fill the vector
   * `values` with the values of \f$\mu\f$. The `values[i]` is interpreted as
   * the value of \f$\mu\f$ at the quadrature point `r[i]`. The code snippet
   * below provides an example in which permeability equals
   *
   * \f$\mu = \mu_0 \big(x^2+y^2+1\big)\f$,
   *
   * where \f$\mu_0\f$ is the permeability of free space.
   * @code
   * #pragma GCC diagnostic push
   * #pragma GCC diagnostic ignored "-Wunused-parameter"
   * template<>
   * void TheCoefficient<3>::value_list(
   * const std::vector<Point<3>> &r,
   * types::material_id mid,
   * unsigned int cuid,
   * std::vector<double> & values) const
   * {
   *   Assert(r.size() == values.size(),
   *   ExcDimensionMismatch(r.size(), values.size()));
   *
   *   auto v = values.begin();
   *   for (auto p: r)
   *   {
   *     *v = MU0 * ( pow(p[0],2) + pow(p[1],2) + 1 );
   *     v++;
   *   }
   * }
   * #pragma GCC diagnostic pop
   * @endcode
   *
   * @param[in] r - A vector that contains the quadrature points of the cell
   * being processed.
   * @param[in] mid - The material ID.
   * @param[in] cuid - The cell user ID.
   * @param[out] values - The output data.
   ***************************************************************************/
  void value_list(const std::vector<Point<dim>>& r,
                  types::material_id mid,
                  unsigned int cuid,
                  std::vector<double>& values) const;
};

/**
 * \brief Implements the source vector field on the right-hand side of the
 * partial differential equation (i) of the
 * [vector boundary value problem](@ref veibvp_bvp).
 *
 * How the source vector field modeled by this function is
 * [interpreted](@ref svsi-type_of_pde_rhs)
 * by the StaticVectorSolver::Solver1 class template depends on the second
 * parameter passed to the constructor of StaticVectorSolver::Solver1, i.e.,
 * the parameter `type_of_pde_rhs`. It can be interpreted as current density,
 * \f$\vec{J}_f\f$, or as the current vector potential, \f$\vec{T}\f$. See the
 * documentation of StaticVectorSolver::Solver1 for more details. The
 * StaticVectorSolver::Solver2 class template does not use
 * StaticVectorSolver::PdeRhs as the right-hand-side of the PDE is fed to the
 * StaticVectorSolver::Solver2 via the input parameters of the constructor
 * `dof_handler` and `solution`.
 *
 * This class is declared in shared/include/static_vector_input.hpp but must be
 * implemented in xyz/src/static_vector_input.cpp, where xyz is the directory of
 * the current numerical experiment. That is, the declaration is shared between
 * all numerical experiments while implementation is specific for each
 * individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The `dim` template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the `stage` template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class.
 ******************************************************************************/
template<int dim, int stage = 1>
class PdeRhs : private Settings
{
public:
  /**
   * \brief Computes the vector field on the right-hand side of the partial
   * differential equation
   *
   * This function is called by the StaticVectorSolver::Solver1 class during
   * the assembly of the right-hand side. It is called ones per each cell. This
   * function must fill `values` vector with the source vectors sampled at
   * quadrature points. The `values[i]` is interpreted as the source vector at
   * quadrature point `r[i]`. The code snippet below provides an example in
   * which the source field equals
   *
   * \f[
   * \vec{J}_f = (1 / \mu_0)(x^2+y^2) \hat{k},
   * \f]
   *
   * where \f$\mu_0\f$ is the permeability of free space.
   *
   * @code
   * #pragma GCC diagnostic push
   * #pragma GCC diagnostic ignored "-Wunused-parameter"
   * template<>
   * void PdeRhs<3>::value_list(
   *   const std::vector<Point<3>> &r,
   *   types::material_id mid,
   *   unsigned int cuid,
   *   std::vector<Tensor<1,3>> & values) const
   * {
   *   Assert(r.size() == values.size(),
   *     ExcDimensionMismatch(r.size(), values.size()));
   *
   *   auto v = values.begin();
   *   for (auto p: r)
   *   {
   *     (*v)[0] = 0.0;
   *     (*v)[1] = 0.0;
   *     (*v)[2] = (1 / MU0) * ( pow(p[0],2) + pow(p[1],2) );
   *     v++;
   *   }
   * }
   * #pragma GCC diagnostic pop
   * @endcode
   *
   * @param[in] r - A vector that contains the quadrature points of the cell
   * being processed.
   * @param[in] mid - The material ID.
   * @param[in] cuid - The cell user ID.
   * @param[out] values - The output data.
   *****************************************************************************/
  void value_list(const std::vector<Point<dim>>& r,
                  types::material_id mid,
                  unsigned int cuid,
                  std::vector<Tensor<1, dim>>& values) const;
};

/**
 * \brief Implements the coefficient \f$\gamma\f$ on the left-hand side of the
 * Robin boundary condition (iii) of the
 * [static vector boundary value problem](@ref veibvp_bvp).
 *
 * This class is declared in shared/include/static_vector_input.hpp but must be
 * implemented in xyz/src/static_vector_input.cpp, where xyz is the directory of
 * the current numerical experiment. That is, the declaration is shared between
 * all numerical experiments while implementation is specific for each
 * individual numerical experiment. See
 * [the structure of the code](@ref page_code)
 * for more details.
 *
 * The `dim` template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the `stage` template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class.
 *****************************************************************************/
template<int dim, int stage = 1>
class Gamma : private Settings
{
public:
  /**
   * \brief Computes values of the parameter \f$\gamma\f$ of the Robin boundary
   * condition at quadrature points.
   *
   * This function is called by the StaticVectorSolver::Solver1 and
   * StaticVectorSolver::Solver2 classes during the assembly of the system
   *matrix. It is called ones per each boundary face if the face ID is an even
   *integer greater than zero, see [boundary ID convention](@ref
   *veis_bnd_convention). This function must fill the vector `values` with the
   *values of \f$\gamma\f$. The `values[i]` is interpreted as \f$\gamma\f$ at
   *quadrature point `r[i]`. The code snippet below provides an example in which
   * \f$\gamma\f$ equals
   *
   * \f[
   * \gamma = (1 / \mu) ( \sqrt{x^2+y^2}+2 )
   * \f],
   *
   * where \f$\mu\f$ is the permeability.
   * @code
   * #pragma GCC diagnostic push
   * #pragma GCC diagnostic ignored "-Wunused-parameter"
   * template<>
   * void Gamma<3>::value_list(
   *   const std::vector<Point<3>>& r,
   *   const std::vector<Tensor<1, 3>>& n,
   *   types::boundary_id bid,
   *   types::material_id mid,
   *   unsigned int cuid,
   *   unsigned int fuid,
   *   std::vector<double>& values) const
   * {
   *   Assert(r.size() == values.size(),
   *     ExcDimensionMismatch(r.size(), values.size()));
   *
   *   auto v = values.begin();
   *   for (auto p: r)
   *   {
   *     double mu = MU0 * (pow(p[0],2) + pow(p[1],2) + 1);
   *     *v = (1 / mu) * (sqrt(pow(p[0],2) + pow(p[1],2)) + 2);
   *     v++;
   *   }
   * }
   * #pragma GCC diagnostic pop
   * @endcode
   *
   * @param[in] r - A vector that contains the quadrature points of the cell
   * being processed.
   * @param[in] n - A vector that contains the normal vectors at quadrature
   * points.
   * @param[in] bid - The boundary ID.
   * @param[in] mid - The material ID.
   * @param[in] cuid - The cell user ID.
   * @param[in] fuid - The face user ID.
   * @param[out] values - The output data.
   ***************************************************************************/
  void value_list(const std::vector<Point<dim>>& r,
                  const std::vector<Tensor<1, dim>>& n,
                  types::boundary_id bid,
                  types::material_id mid,
                  unsigned int cuid,
                  unsigned int fuid,
                  std::vector<double>& values) const;
};

/**
 * \brief Implements the vector field, \f$\vec{Q}\f$, on the right-hand side
 * of the Robin boundary condition (iii) of the
 * [static vector boundary value problem](@ref veibvp_bvp).
 *
 * This class is declared in shared/include/static_vector_input.hpp but must be
 * implemented in xyz/src/static_vector_input.cpp, where xyz is the directory of
 * the current numerical experiment. That is, the declaration is shared between
 * all numerical experiments while implementation is specific for each
 * individual numerical experiment. See
 * [the structure of the code](@ref page_code)
 * for more details.
 *
 * The `dim` template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the `stage` template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class.
 *****************************************************************************/
template<int dim, int stage = 1>
class RobinRhs : private Settings
{
public:
  /**
   * \brief Computes values of the vector field \f$\vec{Q}\f$ on the right-hand
   * side of the Robin boundary condition at quadrature points.
   *
   * This function is called by the StaticVectorSolver::Solver1 and
   * StaticVectorSolver::Solver2 classes during the assembly of the system
   * right-hand side. It is called ones per each boundary face if the face
   * ID is an even integer greater than zero, see
   * [boundary ID convention](@ref veis_bnd_convention).
   * This function must fill the vector `values` with the values of
   * \f$\vec{Q}\f$. The `values[i]` is interpreted as \f$\vec{Q}\f$ at
   * quadrature point `r[i]`. The code snippet below provides an example in
   * which \f$\vec{Q}\f$ equals zero.
   *
   * @code
   * #pragma GCC diagnostic push
   * #pragma GCC diagnostic ignored "-Wunused-parameter"
   * template<>
   * void RobinRhs<3>::value_list(
   * const std::vector<Point<3>>& r,
   * const std::vector<Tensor<1, 3>>& n,
   * types::boundary_id bid,
   * types::material_id mid,
   * unsigned int cuid,
   * unsigned int fuid,
   * std::vector<Tensor<1,3>>& values) const
   * {
   *   Assert(r.size() == values.size(),
   *     ExcDimensionMismatch(r.size(), values.size()));
   *
   *  for (unsigned int i = 0 ; i < r.size(); i++)
   *  {
   *    values[i][0] = 0.0;
   *    values[i][1] = 0.0;
   *    values[i][2] = 0.0;
   *  }
   * }
   * #pragma GCC diagnostic pop
   * @endcode
   *
   * @param[in] r - A vector that contains the quadrature points of the cell
   * being processed.
   * @param[in] n - A vector that contains the normal vectors at quadrature
   * points.
   * @param[in] bid - The boundary ID.
   * @param[in] mid - The material ID.
   * @param[in] cuid - The cell user ID.
   * @param[in] fuid - The face user ID.
   * @param[out] values - The output data.
   ***************************************************************************/
  void value_list(const std::vector<Point<dim>>& r,
                  const std::vector<Tensor<1, dim>>& n,
                  types::boundary_id bid,
                  types::material_id mid,
                  unsigned int cuid,
                  unsigned int fuid,
                  std::vector<Tensor<1, dim>>& values) const;
};

/**
 * \brief Implements the surface free-current density, \f$\vec{K}_f\f$,
 * on the right-hand side of the continuity condition (v) of the
 * [static vector boundary value problem](@ref veibvp_bvp).
 *
 * This class is declared in shared/include/static_vector_input.hpp but must
 * be implemented in xyz/src/static_vector_input.cpp, where xyz is the
 * directory of the current numerical experiment. That is, the declaration
 * is shared between all numerical experiments while implementation is
 * specific for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The `dim` template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the `stage` template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)

 * @endcode
 * member function of this class.
 *****************************************************************************/
template<int dim, int stage = 1>
class FreeSurfaceCurrent : private Settings
{
public:
  /**
   * \brief Computes values of the surface free-current density,
   * \f$\vec{K}_f\f$, on the right-hand side of the continuity condition at
   * quadrature points.
   *
   * This function is called by the StaticVectorSolver::Solver1 and
   * StaticVectorSolver::Solver2 classes during the assembly of the system
   * right-hand side. It is called ones per each face if both user IDs, the
   * cell user ID and the face user ID, are greater than zero. This function
   * must fill the vector `values` with the  values of \f$\vec{K}_f\f$. The
   * `values[i]` is interpreted as \f$\vec{K}_f\f$ at quadrature point `r[i]`.
   * The code snippet below provides an example in which \f$\vec{K}_f\f$
   * equals zero.
   *
   * @code
   * #pragma GCC diagnostic push
   * #pragma GCC diagnostic ignored "-Wunused-parameter"
   * template<>
   * void FreeSurfaceCurrent<3>::value_list(
   * const std::vector<Point<3>>& r,
   * const std::vector<Tensor<1, 3>>& n,
   * types::material_id mid,
   * unsigned int cuid,
   * unsigned int fuid,
   * std::vector<Tensor<1,3>>& values) const
   * {
   *   Assert(r.size() == values.size(),
   *     ExcDimensionMismatch(r.size(), values.size()));
   *
   *  for (unsigned int i = 0 ; i < r.size(); i++)
   *  {
   *    values[i][0] = 0.0;
   *    values[i][1] = 0.0;
   *    values[i][2] = 0.0;
   *  }
   * }
   * #pragma GCC diagnostic pop
   * @endcode
   *
   * @param[in] r - A vector that contains the quadrature points of the cell
   * being processed.
   * @param[in] n - A vector that contains the normal vectors at quadrature
   * points.
   * @param[in] mid - The material ID.
   * @param[in] cuid - The cell user ID.
   * @param[in] fuid - The face user ID.
   * @param[out] values - The output data.
   ***************************************************************************/
  void value_list(const std::vector<Point<dim>>& r,
                  const std::vector<Tensor<1, dim>>& n,
                  types::material_id mid,
                  unsigned int cuid,
                  unsigned int fuid,
                  std::vector<Tensor<1, dim>>& values) const;
};

/**
 * \brief Implements the weight function for calculating the \f$L^2\f$ and
 * \f$L^{\infty}\f$ error norms.
 *
 * This class is declared in shared/include/static_vector_input.hpp but must
 * be implemented in xyz/src/static_vector_input.cpp, where xyz is the
 * directory of the current numerical experiment. That is, the declaration is
 * shared between all numerical experiments while implementation is specific
 * for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The `dim` template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the `stage` template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *****************************************************************************/
template<int dim, int stage = 1>
class Weight
  : public Function<dim>
  , private Settings
{
public:
  /**
   * \brief Returns the value of weight at point `r`. All error norms,
   * \f$L^2\f$, \f$L^{\infty}\f$, and \f$H^1\f$, at point `r` will be multiplied
   * by this value.
   *
   * @param [in] r - The quadrature point at which the error norm is currently
   * evaluated.
   * @param [in] component - The component of the function. Must equal zero in
   *the case of a scalar function.
   *****************************************************************************/
  virtual double value(const Point<dim>& r,
                       const unsigned int component = 0) const override final;
};

} // namespace StaticVectorSolver

#endif
