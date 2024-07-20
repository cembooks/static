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

#ifndef StaticScalarInput_H__
#define StaticScalarInput_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include "constants.hpp"
#include "settings.hpp"

namespace StaticScalarSolver
{

/**
 * \brief Implements the coefficient (\f$\epsilon\f$, \f$\mu\f$, etc.) in the
 * div-grad partial differential equation (i) of the
 * [scalar boundary value problem](@ref seibvp_bvp).
 *
 * This class implements the coefficient of the
 * [scalar boundary value problem](@ref seibvp_bvp)
 * that reflects the property of the materials. In the case of a problem
 * formulated in terms of the electrostatic potential, \f$\Phi\f$, this
 * coefficient equals permittivity, \f$\epsilon\f$. In the case of a problem
 * formulated in terms of the total magnetic potential, \f$\Psi\f$, or reduced
 * magnetic potential, \f$\Theta\f$, this coefficient equals permeability,
 * \f$\mu\f$. In the case of a planar two-dimensional problem formulated in
 * terms of the vector potential, \f$A\f$, this coefficient equals
 * \f$1 / \mu\f$. In the case of an axisymmetric two-dimensional problem
 * formulated in terms of the scaled vector potential, \f$A'\f$, this
 * coefficient equals \f$1 / \mu'\f$. If the current vector potential,
 * \f$T\f$, is being calculated, this coefficient equals 1 at all points of
 * the problem domain. Simply put, this is the coefficient in between the
 * divergence and the gradient on the left-hand side of the div-grad partial
 * differential equation.
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but
 * must be implemented in xyz/src/static_scalar_input.cpp, where xyz is
 * the directory of the current numerical experiment. That is, the declaration
 * is shared between all numerical experiments while implementation is specific
 * for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class template.
 *****************************************************************************/
template<int dim, int stage = 1>
class TheCoefficient : private Settings
{
public:

/**
 * \brief Computes the values of the coefficient at quadrature points.
 *
 * This function is called by objects derived from the StaticScalarSolver::Solver
 * template during the assembly of the system matrix. It is called ones per each
 * cell. This function must fill the vector values[...] with the values of the
 * coefficient. The values[i] is interpreted as the value of the coefficient
 * at quadrature point r[i]. The code snippet below provides an example in
 * which the permittivity equals the permittivity of the free space at each
 * quadrature point.
 *
 * @code
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wunused-parameter"
 *
 * template<>
 * void TheCoefficient<2>::value_list(
 * 	const std::vector<Point<2>> &r,
 * 	types::material_id mid,
 * 	unsigned int cuid,
 * 	std::vector<double> & values) const
 * 	{
 * 		Assert(r.size() == values.size(),
 * 		ExcDimensionMismatch(r.size(), values.size()));
 *
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = ep_0;
 * 	}
 *
 * #pragma GCC diagnostic pop
 * @endcode
 *
 * @param[in] r - A vector that contains the quadrature points of the cell
 * being processed.
 * @param[in] mid - The material ID.
 * @param[in] cuid - The cell user ID.
 * @param[out] values - The output data.
 ***************************************************************************/
	void value_list(
	const std::vector<Point<dim>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const;
};

/**
 * \brief Implements the right-hand side
 * (\f$\rho_f\f$, \f$\rho_t\f$, \f$\rho_r\f$, or \f$J_f\f$) in the div-grad
 * partial differential equation (i) of the
 * [scalar boundary value problem](@ref seibvp_bvp).

 * In the case of a problem formulated in terms of the electrostatic scalar
 * potential, \f$\Phi\f$, the right-hand side is the volume free-charge density,
 * \f$\rho_f\f$. In the case of a problem formulated in terms of the total
 * magnetic scalar potential, \f$\Psi\f$ the right-hand side is the volume
 * magnetic-charge density, \f$\rho_t\f$. (Note, that normally in the literature
 * \f$\rho_t / \mu_0\f$ is designated as volume magnetic-charge density.)
 * In the case of a problem formulated in terms of the reduced magnetic scalar
 * potential, \f$\Theta\f$, the right-hand side is the reduced volume magnetic-
 * charge density, \f$\rho_r\f$. In the cases of a planar two-dimensional
 * problem formulated in terms of the magnitude of the vector potential, \f$A\f$,
 * and an axisymmetric two-dimensional problem formulated in terms of the scaled
 * magnitude of the vector potential, \f$A'\f$, the right-hand side is the
 * magnitude of the free-current density, \f$J_f\f$. In the case of the planar
 * two-dimensional problem formulated in terms of the magnitude of the current
 * vector potential, \f$T\f$, the StaticScalarSolver::PdeRhsCvp is used instead
 * of this class template.
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but
 * must be implemented in xyz/src/static_scalar_input.cpp, where xyz is
 * the directory of the current numerical experiment. That is, the declaration
 * is shared between all numerical experiments while implementation is specific
 * for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class template.
 *****************************************************************************/
template<int dim, int stage = 1>
class PdeRhs : private Settings
{
public:
/**
 * \brief Computes the right-hand side of the div-grad partial differential
 * equation at quadrature points.
 *
 * This function is called by objects derived from the StaticScalarSolver::Solver
 * template during the assembly of the system right-hand side. It is called ones per each
 * cell. This function must fill the vector values[...] with the values of the
 * right-hand side of the partial differential equation. The values[i] is
 * interpreted as the right-hand side of the partial differential equation at
 * quadrature point r[i]. The right- hand side of the partial differential
 * equation may vary from point-to-point, or may have the same value at all
 * points. The following code snippet provides an example.
 * @code
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wunused-parameter"
 *
 * template<>
 * void PdeRhs<2>::value_list(
 * 	const std::vector<Point<2>> &r,
 * 	types::material_id mid,
 * 	unsigned int cuid,
 * 	std::vector<double> & values) const
 * {
 * 	Assert(r.size() == values.size(),
 * 		ExcDimensionMismatch(r.size(), values.size()));
 *
 * 	auto v = values.begin();
 * 	for (auto p: r)
 * 	{
 * 		*v = ep_0 * (
 * 			2*p[0]*pow(p[1],2)*sin(k*p[0]) +
 * 			2*p[1]*pow(p[0],2)*sin(k*p[1]) +
 * 			k*(pow(p[0],2)*pow(p[1],2) + 1) *
 * 			(cos(k*p[0]) + cos(k*p[1]))
 * 			);
 * 		v++;
 * 	}
 * }
 *
 * #pragma GCC diagnostic pop
 * @endcode
 *
 * @param[in] r - A vector that contains the quadrature points of the cell
 * being processed.
 * @param[in] mid - The material ID.
 * @param[in] cuid - The cell user ID.
 * @param[out] values - The output data.
 ***************************************************************************/
	void value_list(
	const std::vector<Point<dim>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const;
};

/**
 * \brief Implements the two-dimensional free-current density \f$\vec{J}_f\f$
 * on the right-hand side of the partial differential equation (i) of the
 * [scalar boundary value problem](@ref seibvp_bvp).
 *
 * If a planar two-dimensional problem is formulated in terms of the magnitude
 * of the current vector potential, \f$T\f$, the div-grad equation (i) in
 * [scalar boundary value problem](@ref seibvp_bvp) is replaced with
 * \f[
 * -\vec{\nabla}\cdot \bigg( \vec{\nabla} T \bigg) =
 * \vec{\nabla}\overset{S}{\times} \vec{J}_f.
 * \f]
 * The purpose of this class template is to implement \f$\vec{J}_f\f$ in the
 * last equation. This class template must be ignored if the problem formulated in
 * terms of \f$\Phi\f$, \f$\Psi\f$, \f$\Theta\f$, \f$A\f$, and \f$A'\f$.
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but
 * must be implemented in xyz/src/static_scalar_input.cpp, where xyz is
 * the directory of the current numerical experiment. That is, the declaration
 * is shared between all numerical experiments while implementation is specific
 * for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class template.
 *
 * @note This class template can only be used in two-dimensional problems, so,
 * strictly peaking, the template parameter dim is redundant. This template
 * parameter, however, makes the code more uniform.
 *****************************************************************************/
template<int dim, int stage = 1>
class PdeRhsCvp : private Settings
{
public:
/**
 * \brief Computes the two-dimensional free-current density \f$\vec{J}_f\f$
 * on the right-hand side of the partial differential equation at quadrature
 * points.
 *
 * This function is called by objects derived from the StaticScalarSolver::Solver
 * template during the assembly of the system right-hand side. It is called ones
 * per each cell and each face located on the boundary. This function must fill
 * the vector values with the values of \f$\vec{J}_f\f$. The following code
 * snippet provides an example.
 * @code
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wunused-parameter"
 *
 * template<>
 * void StaticScalarSolver::PdeRhsCvp<2>::value_list(
 *	const std::vector<Point<2>> &r,
 *	types::material_id mid,
 *	unsigned int cuid,
 *	std::vector<Tensor<1, 2>> & values) const
 * {
 *	Assert(r.size() == values.size(),
 *		ExcDimensionMismatch(r.size(), values.size()));
 *
 *		auto v = values.begin();
 *		for (auto p: r)
 *		{
 *			if (mid == mid_2)
 *			{
 *				(*v)[0] =-p[1];
 *				(*v)[1] = p[0];
 *			}
 *			else
 *			{
 *				(*v)[0] = 0.0;
 *				(*v)[1] = 0.0;
 *			}
 *
 *			v++;
 *		}
 *  }
 *
 * #pragma GCC diagnostic pop
 * @endcode
 *
 * @param[in] r - A vector that contains the quadrature points of the cell
 * being processed.
 * @param[in] mid - The material ID.
 * @param[in] cuid - The cell user ID.
 * @param[out] values - The output data.
 ***************************************************************************/
	void value_list(
	const std::vector<Point<dim>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1, dim>> & values) const;
};

/**
 * \brief Implements the coefficient \f$\gamma\f$ in the left-hand side of the
 * Robin boundary condition (iii) of the
 * [static scalar boundary value problem](@ref seibvp_bvp).
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but
 * must be implemented in xyz/src/static_scalar_input.cpp, where xyz
 * is the directory of the current numerical experiment. That is, the
 * declaration is shared between all numerical experiments while implementation
 * is specific for each individual numerical experiment. See
 * [the structure of the code](@ref page_code)
 * for more details.
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class template.
 *****************************************************************************/
template<int dim, int stage = 1>
class Gamma : private Settings
{
public:
/**
 * \brief Computes the coefficient \f$\gamma\f$ at quadrature points.
 *
 * This function is called by objects derived from the StaticScalarSolver::Solver
 * template during the assembly of the system matrix.
 * It is called once for each face located on the boundary if the cell ID is
 * even and is greater than zero. This function must fill the vector values[...]
 * with the values of \f$\gamma\f$. The values[i] is interpreted as
 * \f$\gamma\f$ at quadrature point r[i]. The \f$\gamma\f$ may vary from
 * point-to-point, or may have the same value at all points. The following code
 * snippet provides an example.
 * @code
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wunused-parameter"
 * template<>
 * void Gamma<2>::value_list(
 * 	const std::vector<Point<2>> &r,
 * 	const std::vector<Tensor<1, 2>> & n,
 * 	types::boundary_id bid,
 * 	types::material_id mid,
 * 	unsigned int cuid,
 * 	unsigned int fuid,
 * 	std::vector<double> & values) const
 * {
 * 	Assert(r.size() == values.size(),
 * 		ExcDimensionMismatch(r.size(), values.size()));
 *
 * 	if (bid == 2)
 * 	{
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = ep_0 / r[i].norm();
 * 	}
 * 	else
 * 	{
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = 0.0;
 * 	}
 * }
 *
 * #pragma GCC diagnostic pop
 * @endcode
 *
 * @param[in] r - A vector that contains the quadrature points of the cell
 * being processed.
 * @param[in] n - A vector that contains normal vectors at the quadrature
 * points of the cell being processed.
 * @param[in] bid - The boundary ID.
 * @param[in] mid - The material ID.
 * @param[in] cuid - The cell user ID.
 * @param[in] fuid - The face user ID.
 * @param[out] values - The output data.
 ***************************************************************************/
	void value_list(
	const std::vector<Point<dim>> &r,
	const std::vector<Tensor<1, dim>> & n,
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const;
};

/**
 * \brief Implements the right-hand side (\f$\sigma\f$ or \f$Q\f$) of the Robin
 * boundary condition (iii) in the
 * [static scalar boundary value problem](@ref seibvp_bvp).
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but
 * must be implemented in xyz/src/static_scalar_input.cpp, where xyz
 * is the directory of the current numerical experiment. That is, the
 * declaration is shared between all numerical experiments while implementation
 * is specific for each individual numerical experiment. See
 * [the structure of the code](@ref page_code)
 * for more details.
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class template.
 *****************************************************************************/
template<int dim, int stage = 1>
class RobinRhs : private Settings
{
public:
/**
 * \brief Computes the right-hand side of the Robin boundary condition
 * (\f$\sigma\f$ or \f$Q\f$).
 *
 * This function is called by objects derived from the StaticScalarSolver::Solver
 * template during the assembly of the system right-hand side. It is
 * called once for each face located on the boundary if the face ID is even
 * and is greater than zero. This function must fill the vector values[...]
 * with the values of the right-hand side of the Robin boundary condition.
 * The values[i] is interpreted as the right-hand side at quadrature
 * point r[i]. The value of the right-hand side may vary from
 * point-to-point, or may have the same value at all points. The following
 * code snippet provides an example.
 * @code
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wunused-parameter"
 *
 * template<>
 * void RobinRhs<2>::value_list(
 * 	const std::vector<Point<2>> &r,
 * 	const std::vector<Tensor<1, 2>> & n,
 * 	types::boundary_id bid,
 * 	types::material_id mid,
 * 	unsigned int cuid,
 * 	unsigned int fuid,
 * 	std::vector<double> & values) const
 * {
 * 	Assert(r.size() == values.size(),
 * 		ExcDimensionMismatch(r.size(), values.size()));
 *
 * 	if (bid == 2)
 * 	{
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = sigma;
 * 	}
 * 	else
 * 	{
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = 0.0;
 * 	}
 * }
 *
 * #pragma GCC diagnostic pop
 * @endcode
 *
 * @param[in] r - A vector that contains the quadrature points of the cell
 * being processed.
 * @param[in] n - A vector that contains normal vectors at the quadrature
 * points of the cell being processed.
 * @param[in] bid - The boundary ID.
 * @param[in] mid - The material ID.
 * @param[in] cuid - The cell user ID.
 * @param[in] fuid - The call user ID.
 * @param[out] values - The output data.
 ***************************************************************************/
	void value_list(
	const std::vector<Point<dim>> &r,
	const std::vector<Tensor<1, dim>> & n,
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const;
};

/**
 * \brief Implements the right-hand side (\f$\kappa_f\f$, \f$\kappa_t\f$,
 * \f$\kappa_r\f$, or \f$K_f\f$), of the continuity condition (v) in the
 * [static scalar boundary value problem](@ref seibvp_bvp).
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but
 * must be implemented in xyz/src/static_scalar_input.cpp, where xyz
 * is the directory of the current numerical experiment. That is, the
 * declaration is shared between all numerical experiments while implementation
 * is specific for each individual numerical experiment. See
 * [the structure of the code](@ref page_code)
 * for more details.
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to implement
 * @code
 * void value_list(...)
 * @endcode
 * member function of this class template.
 *****************************************************************************/
template<int dim, int stage = 1>
class FreeSurfaceCharge : private Settings
{
public:
/**
 * \brief Computes the surface free-charge density
 * (\f$\kappa_f\f$ or \f$K_f\f$) on the right-hand side in the continuity
 * condition.
 *
 * This function is called by the StaticScalarSolver::Solver class during the
 * assembly of the right-hand side of the system of linear equations. It is
 * called once for each face if both user IDs, the cell user ID and the face
 * user ID, are greater than zero. This function must fill the vector
 * values[...] with the values of the surface free-charge density. The
 * values[i] is interpreted as surface free-charge density at quadrature
 * point r[i]. The value of the surface free-charge density may vary from
 * point-to-point, or may have the same value at all points. The following
 * code snippet provides an example.
 * @code
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wunused-parameter"
 *
 *template<>
 *void FreeSurfaceCharge<2>::value_list(
 * 	const std::vector<Point<2>> &r,
 * 	const std::vector<Tensor<1,2>> & n,
 * 	types::material_id mid,
 * 	unsigned int cuid,
 * 	unsigned int fuid,
 * 	std::vector<double> & values) const
 * {
 * 	Assert(r.size() == values.size(),
 * 		ExcDimensionMismatch(r.size(), values.size()));
 *
 * 	if ( (cuid == 1) && (fuid == 2) )
 * 	{
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = 1.0;
 * 	}
 * 	else
 * 	{
 * 		for (unsigned int i = 0; i < values.size(); i++)
 * 			values[i] = 0.0;
 * 	}
 * }
 * #pragma GCC diagnostic pop
 * @endcode
 *
 * @param[in] r - A vector that contains the quadrature points of the cell being
 * processed.
 * @param[in] n - The vector normal to the face.
 * @param[in] mid - The material ID.
 * @param[in] cuid - The cell user ID.
 * @param[in] fuid - The call user ID.
 * @param[out] values - The output data.
 ***************************************************************************/
	void value_list(
	const std::vector<Point<dim>> &r,
	const std::vector<Tensor<1, dim>> & n,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const;
};

/**
 * \brief Implements the weight function for calculating the \f$L^2\f$
 * and \f$H^1\f$ error norms.
 *
 * Objects of this class template are passed to
 * [VectorTools::integrate_difference](https://www.dealii.org/current/doxygen/deal.II/namespaceVectorTools.html).
 *
 * This class template is used for selecting regions of the problem domain in
 * which the error norms are calculated.
 *
 * This class template is declared in shared/include/static_scalar_input.hpp but must
 * be implemented in xyz/src/static_scalar_input.cpp, where xyz is the
 * directory of the current numerical experiment. That is, the declaration is
 * shared between all numerical experiments while implementation is specific
 * for each individual numerical experiment. See
 * [the structure of the code](@ref page_code) for more details.
 *
 * The user is supposed to implement
 * @code
 * double value(...)
 * @endcode
 * member function of this class template.
 *****************************************************************************/
template<int dim, int stage = 1>
class Weight : public Function<dim>,  private Settings
{
public:
/**
 * \brief Return the value of weight at point r. Both error norms, \f$L^2\f$
 * and \f$H^1\f$, at point \f$r\f$ will be multiplied by this value.
 *
 * @param r - The quadrature point at which the error norm is currently
 * evaluated.
 *****************************************************************************/
	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;
};

} // namespace StaticScalarSolver

#endif

