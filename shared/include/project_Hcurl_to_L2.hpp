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

#ifndef ProjectHcurlToL2_H__
#define ProjectHcurlToL2_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <string>
#include <fstream>
#include <ios>
#include <iostream>
#include <iomanip>

#include "constants.hpp"

#define SE scratch_data.se

#define TMR(__name) \
	TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

namespace StaticVectorSolver
{
/**
 * \brief Projects from \f$H(\text{curl})\f$ to \f$L^2\f$.
 *
 * This class template is not supposed to be used directly. Instead
 * the wrap-around class template StaticVectorSolver::ProjectAxyToBz
 * must be used. The name of the wrap-around class template is assumed to be
 * more familiar to readers in electromagnetics.
 *
 * Implements the following recipes:
 * - (1) Recipe for projection from H(curl) to \f$L^2\f$ nr. 14.
 *
 * This class template is supposed to be used in pair with
 * StaticVectorSolver::Solver1 or StaticVectorSolver::Solver2. In all planar
 * two-dimensional problems formulated in terms of the magnetic vector potential,
 * \f$\vec{A}\f$, the numerically calculated potential needs to be converted into
 * magnetic field, \f$B\f$, as
 * \f[
 * B = \vec{\nabla} \overset{S}{\times} \vec{A}.
 * \f]
 * The magnetic vector potential belongs to the H(curl) function space. The
 * magnetic field belongs to the \f$L^2\f$ function space. Therefore, one needs to
 * compute the equation above such that the input, \f$\vec{A}\f$, is in H(curl)
 * and the output, \f$B\f$, is in \f$L^2\f$. Such computation can be
 * envisioned as a some kind of projection from one function space, H(curl),
 * into another, \f$L^2\f$.
 *
 * This class template does this projection. The Bossavit's diagram
 * below illustrates it.
 *
 * ![](svs_prj/diagram_prj_Axy_to_Bz.svg)
 *
 * The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to derive an object from this class template. All usual
 * computations, i.e., setup, assembling the linear system, etc., happen
 * automatically.
 *****************************************************************************/
template<int stage = 1>
class ProjectHcurlToL2
{
public:
	ProjectHcurlToL2() = delete;

/**
 * \brief The only constructor.
 *
 * @param[in] p - The degree of the discontinuous Lagrange finite elements,
 * FE_DGQ, that span the \f$L^2\f$ function space.
 * @param[in] mapping_degree - The degree of the interpolating polynomials used
 * for mapping. Setting it to 1 will do in the most of the cases. Note, that it
 * makes sense to attach a meaningful manifold to the triangulation if this
 * parameter is greater than 1.
 * @param[in] triangulation_Hcurl - Reference to the triangulation inside the
 * object of the class StaticVectorSolver::Solver1 or
 * StaticVectorSolver::Solver2 that has yielded the magnetic vector potential
 * in the  H(curl) function space.
 * @param[in] dof_handler_Hcurl - Reference to the dof handler inside the class
 * StaticVectorSolver::Solver1 or StaticVectorSolver::Solver2 that has yielded
 * the magnetic vector potential in the H(curl) function space.
 * @param[in] solution_Hcurl - Vector filled with the degrees of freedom that
 * together with the shape functions of the Nedelec finite elements model the
 * magnetic vector potential in the H(curl) function space.
 * @param[in] fname - The name of the output files without extension.
 * @param[in] exact_solution - Points to an object that describes the exact
 * solution to the problem. It is needed for calculating error norms.
 * @param[in] print_time_tables - If true, prints time tables on the screen.
 * @param[in] project_exact_solution - If true, projects the exact solution
 * onto the space spanned by the Raviart-Thomas finite elements and saves the
 * result into the vtk file next to the solution. This may be useful for
 * debugging purposes.
 * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
 * gradient solver into a file. The name of the file is generated by appending
 * "_cg_convergence.csv" to fname.
 ******************************************************************************/
	ProjectHcurlToL2(
		unsigned int p,
		unsigned int mapping_degree,
		const Triangulation<2> & triangulation_Hcurl,
		const DoFHandler<2> & dof_handler_Hcurl,
		const Vector<double> & solutioin_Hcurl,
		std::string fname = "L2",
		const Function<2> * exact_solution = nullptr,
		bool print_time_tables = false,
		bool project_exact_solution = false,
		bool log_cg_convergence = false
		);

/**
 * \brief Returns \f$L^2\f$ error norm.
 *****************************************************************************/
	double get_L2_norm() {return L2_norm;};

/**
 * \brief Returns \f$L^{\infty}\f$ error norm.
 *****************************************************************************/
	double get_Linfty_norm(){return Linfty_norm;}

/**
 * \brief Returns the number of active cells in the mesh.
 *****************************************************************************/
	unsigned int get_n_cells() const
		{return static_cast<unsigned int>(triangulation_Hcurl.n_active_cells());}

/**
 * \brief Returns the total amount of the degrees of freedom.
 *****************************************************************************/
	unsigned int get_n_dofs() const
		{return static_cast<unsigned int>(dof_handler_L2.n_dofs());}

/**
 * \brief Releases computer memory associated with system matrix and
 * right-hand side.
 *****************************************************************************/
	void clear() {system_matrix.clear(); system_rhs.reinit(0);}

/**
 * \brief Returns a reference to triangulation.
 *****************************************************************************/
	const Triangulation<2> & get_tria() const {return triangulation_Hcurl;}

/**
 * \brief Returns a reference to dof handler associated with the Raviart-Thomas
 * finite elements.
 *****************************************************************************/
	const DoFHandler<2> & get_dof_handler() const {return dof_handler_L2;}

/**
 * \brief Returns a reference to solution, i.e., the result of the projection.
 ******************************************************************************/
	const Vector<double> & get_solution() const {return solution_L2;}

/**
 * \brief Saves the system matrix and the right-hand side into a csv file.
 *
 * All the zeros included into the csv files. This is a very dumb and
 * inefficient way of saving sparse matrices. On the positive side - it is very
 * easy and straightforward to read the csv files. This function may be useful
 * for debugging. One can assemble the system on a coarse mesh (so there are
 * a few mesh cells and the system matrix is small) and export the system matrix
 * together with the right-hand side into another program such as Matlab or
 * GNU Octave for an analysis.
 *
 * @param[in] fname - A stem of the names of the output files. The matrix will
 * be saved into fname_matrix.csv file. The right-hand side will be save into
 * fname_rhs.csv file.
 *****************************************************************************/
	void save_matrix_and_rhs_to_csv(std::string fname) const;

private:

	void setup();
	void assemble();
	void solve();
	void save() const;
	void compute_error_norms();
	void project_exact_solution_fcn();

	const std::string fname;

	const DoFHandler<2> & dof_handler_Hcurl;
	const Vector<double> & solution_Hcurl;

	const Triangulation<2> & triangulation_Hcurl;
	const FE_DGQ<2> fe_L2;
	DoFHandler<2> dof_handler_L2;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double> solution_L2;
	Vector<double> system_rhs;

	Vector<double> projected_exact_solution;

	AffineConstraints<double> constraints;

	const Function<2> * exact_solution;

	const unsigned int mapping_degree;
	const bool project_exact_solution;
	const bool log_cg_convergence;

	Vector<double> L2_per_cell;
	double L2_norm;

	Vector<double> Linfty_per_cell;
	double Linfty_norm;

// ----------------------------------------------------------------------------
// These structures and functions are related to the Work Stream algorithm.
// See article "WorkStream – A Design Pattern for Multicore-Enabled Finite
// Element Computations." by BRUNO TURCKSIN, MARTIN KRONBICHLER,
// WOLFGANG BANGERTH for more details.
// ----------------------------------------------------------------------------

	using IteratorTuple =
		std::tuple<typename DoFHandler<2>::active_cell_iterator,
		typename DoFHandler<2>::active_cell_iterator>;

	using IteratorPair = SynchronousIterators<IteratorTuple>;

	struct AssemblyScratchData
	{
		AssemblyScratchData(const FiniteElement<2> &fe,
			const DoFHandler<2> & dof_handr_Hcurl,
			const Vector<double> & dofs_Hcurl,
			unsigned int mapping_degree);

		AssemblyScratchData(const AssemblyScratchData & scratch_data);

		MappingQ<2> mapping;
		Constants::QuadratureTableVector<2> qt;
		FEValues<2> fe_values_L2;
		FEValues<2> fe_values_Hcurl;

		const unsigned int dofs_per_cell;
		const unsigned int n_q_points;

		std::vector<std::vector<Tensor<1,2>>> vector_gradients;
		Tensor<1,1> curl_vec_in_Hcurl;

		const FEValuesExtractors::Scalar se;

		const DoFHandler<2> & dof_hand_Hcurl;
		const Vector<double> & dofs_Hcurl;
	};

	struct AssemblyCopyData
	{
		FullMatrix<double>	cell_matrix;
		Vector<double>	cell_rhs;
		std::vector<types::global_dof_index> local_dof_indices;
	};

	void system_matrix_local(
		const IteratorPair & IP,
		AssemblyScratchData & scratch_data,
		AssemblyCopyData & copy_data);

	void copy_local_to_global(const AssemblyCopyData &copy_data);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
};

template<int stage>
ProjectHcurlToL2<stage>::ProjectHcurlToL2(
	unsigned int p,
	unsigned int mapping_degree,
	const Triangulation<2> & triangulation_Hcurl,
	const DoFHandler<2> & dof_handler_Hcurl,
	const Vector<double> & solution_Hcurl,
	std::string fname,
	const Function<2> * exact_solution,
	bool print_time_tables,
	bool project_exact_solution,
	bool log_cg_convergence
	):
		fname(fname),
		dof_handler_Hcurl(dof_handler_Hcurl),
		solution_Hcurl(solution_Hcurl),
		triangulation_Hcurl(triangulation_Hcurl),
		fe_L2(p),
		exact_solution(exact_solution),
		mapping_degree(mapping_degree),
		project_exact_solution(project_exact_solution),
		log_cg_convergence(log_cg_convergence)
{
	TimerOutput::OutputFrequency tf =
		(print_time_tables) ? TimerOutput::summary : TimerOutput::never;

	TimerOutput timer(
		std::cout,
		tf,
		TimerOutput::cpu_and_wall_times_grouped);

	{TMR("Setup"); setup();}
	{TMR("Assemble"); assemble();}
	{TMR("Solve"); solve();}

	if (exact_solution)
	{
		{TMR("Compute error norms"); compute_error_norms();}

		if (project_exact_solution)
		{
			{TMR("Project exact solution"); project_exact_solution_fcn();}
		}
	}

	{TMR("Save"); save();}
}

template<int stage>
void ProjectHcurlToL2<stage>::setup()
{
	constraints.close();

	dof_handler_L2.reinit(triangulation_Hcurl);
	dof_handler_L2.distribute_dofs(fe_L2);
	constraints.close();

	DynamicSparsityPattern dsp(dof_handler_L2.n_dofs(), dof_handler_L2.n_dofs());
	DoFTools::make_sparsity_pattern(
		dof_handler_L2,
		dsp,
		constraints,
		false);

	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	solution_L2.reinit(dof_handler_L2.n_dofs());
	system_rhs.reinit(dof_handler_L2.n_dofs());

	if (project_exact_solution)
		projected_exact_solution.reinit(dof_handler_L2.n_dofs());

	if (exact_solution)
	{
		L2_per_cell.reinit(triangulation_Hcurl.n_active_cells());
		Linfty_per_cell.reinit(triangulation_Hcurl.n_active_cells());
	}
}

template<int stage>
void ProjectHcurlToL2<stage>::solve()
{
	ReductionControl control(system_rhs.size(), 0.0, 1e-8, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

	PreconditionSSOR<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(system_matrix, 1.2);

	cg.solve(system_matrix, solution_L2, system_rhs, preconditioner);

	if (log_cg_convergence)
	{
		const std::vector<double> history_data = control.get_history_data();

		std::ofstream ofs(fname + "_cg_convergence.csv");

		unsigned int i = 1;
		for (auto item : history_data)
		{
			ofs << i << ", " << item  << "\n";
			i++;
		}
		ofs.close();
	}
}

template<int stage>
void ProjectHcurlToL2<stage>::save() const
{
	std::vector<std::string> solution_names(1, "ScalarField");
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
		data_component_interpretation(
			1, DataComponentInterpretation::component_is_scalar);

	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler_L2);
	data_out.add_data_vector(
			solution_L2,
			solution_names,
			DataOut<2>::type_dof_data,
			data_component_interpretation);

	if (project_exact_solution && exact_solution)
	{
		solution_names.at(0) = "ScalarFieldExact";

		data_out.add_data_vector(
			projected_exact_solution,
			solution_names,
			DataOut<2>::type_dof_data,
			data_component_interpretation);
	}


	if (exact_solution)
	{
		solution_names.at(0) = "L2norm";

		data_out.add_data_vector(
			L2_per_cell,
			solution_names,
			DataOut<2>::type_cell_data,
			data_component_interpretation);

		solution_names.at(0) = "LinftyNorm";

		data_out.add_data_vector(
			Linfty_per_cell,
			solution_names,
			DataOut<2>::type_cell_data,
			data_component_interpretation);
	}


//	if (exact_solution)
//		data_out.add_data_vector(L2_per_cell, "L2norm");

	data_out.build_patches();

	std::ofstream out(fname + ".vtk");

	data_out.write_vtk(out);
	out.close();
}

template<int stage>
void ProjectHcurlToL2<stage>::compute_error_norms()
{
	Weight<2,stage> weight;
	const Function<2, double> * mask = & weight;

	Constants::QuadratureTableScalar<2> qt(dof_handler_L2.get_fe().degree + 1);

	VectorTools::integrate_difference(
		MappingQ<2> (mapping_degree),
		dof_handler_L2,
		solution_L2,
		* exact_solution,
		L2_per_cell,
		QGauss<2>(qt.enorm()),
		VectorTools::L2_norm,
		mask);

	L2_norm = VectorTools::compute_global_error(
		triangulation_Hcurl,
		L2_per_cell,
		VectorTools::L2_norm);

	VectorTools::integrate_difference(
		MappingQ<2> (mapping_degree),
		dof_handler_L2,
		solution_L2,
		* exact_solution,
		Linfty_per_cell,
		QGauss<2>(1),
		VectorTools::Linfty_norm,
		mask // & B_mask
		);

	Linfty_norm = Linfty_per_cell.linfty_norm();
}

template<int stage>
void ProjectHcurlToL2<stage>::project_exact_solution_fcn()
{
	Constants::QuadratureTableScalar<2> qt(fe_L2.degree + 1);

	AffineConstraints<double> constraints_empty;
	constraints_empty.close();

	VectorTools::project(
		MappingQ<2> (mapping_degree),
		dof_handler_L2,
		constraints_empty,
		QGauss<2> (qt.sim()),
		* exact_solution,
		projected_exact_solution);
}

template<int stage>
ProjectHcurlToL2<stage>::AssemblyScratchData::AssemblyScratchData(
	const FiniteElement<2> &fe,
	const DoFHandler<2> & dof_hand_Hcurl,
	const Vector<double> & dofs_Hcurl,
	unsigned int mapping_degree):
		mapping(mapping_degree),
		qt(fe.degree),
		fe_values_L2(mapping, fe, QGauss<2>(qt.sim()),
			update_values |
			update_quadrature_points |
			update_JxW_values),
		fe_values_Hcurl(mapping, dof_hand_Hcurl.get_fe(), QGauss<2>(qt.sim()),
			update_gradients),
		dofs_per_cell(fe_values_L2.dofs_per_cell),
		n_q_points(fe_values_L2.get_quadrature().size()),
		vector_gradients(n_q_points,std::vector<Tensor<1,2>>(2)),
		se(0),
		dof_hand_Hcurl(dof_hand_Hcurl),
		dofs_Hcurl(dofs_Hcurl)
{
}

template<int stage>
ProjectHcurlToL2<stage>::AssemblyScratchData::AssemblyScratchData(
		const AssemblyScratchData & scratch_data):
		mapping(scratch_data.mapping.get_degree()),
		qt(scratch_data.qt),
		fe_values_L2(mapping, scratch_data.fe_values_L2.get_fe(),
			scratch_data.fe_values_L2.get_quadrature(),
			update_values |
			update_quadrature_points |
			update_JxW_values),
		fe_values_Hcurl(mapping, scratch_data.fe_values_Hcurl.get_fe(),
			scratch_data.fe_values_Hcurl.get_quadrature(),
			update_gradients),
		dofs_per_cell(fe_values_L2.dofs_per_cell),
		n_q_points(fe_values_L2.get_quadrature().size()),
		vector_gradients(n_q_points,std::vector<Tensor<1,2>>(2)),
		se(0),
		dof_hand_Hcurl(scratch_data.dof_hand_Hcurl),
		dofs_Hcurl(scratch_data.dofs_Hcurl)
{
}

template<int stage>
void ProjectHcurlToL2<stage>::system_matrix_local(
	const IteratorPair & IP,
	AssemblyScratchData & scratch_data, AssemblyCopyData & copy_data)
{
// See the color box
// Recipe for projections from H(curl) to L^2 nr. 14

	copy_data.cell_matrix.reinit(
			scratch_data.dofs_per_cell, scratch_data.dofs_per_cell);

	copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

	copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

	scratch_data.fe_values_L2.reinit(std::get<0>(*IP));
	scratch_data.fe_values_Hcurl.reinit(std::get<1>(*IP));

	scratch_data.fe_values_Hcurl.get_function_gradients(
		scratch_data.dofs_Hcurl,
		scratch_data.vector_gradients);

	for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
	{
		for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
		{
			for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
			{
				copy_data.cell_matrix(i, j) +=
					scratch_data.fe_values_L2[SE].value(i,q_index) *
					scratch_data.fe_values_L2[SE].value(j,q_index) *
					scratch_data.fe_values_L2.JxW(q_index);
			}
			scratch_data.curl_vec_in_Hcurl[0] =
				scratch_data.vector_gradients.at(q_index).at(1)[0] -
				scratch_data.vector_gradients.at(q_index).at(0)[1];

			copy_data.cell_rhs(i) +=
				scratch_data.curl_vec_in_Hcurl[0] *
				scratch_data.fe_values_L2[SE].value(i,q_index) *
				scratch_data.fe_values_L2.JxW(q_index);
		}
	}

	std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
}

template<int stage>
void ProjectHcurlToL2<stage>::copy_local_to_global(const AssemblyCopyData &copy_data)
{
	constraints.distribute_local_to_global(
		copy_data.cell_matrix,
		copy_data.cell_rhs,
		copy_data.local_dof_indices,
		system_matrix,
		system_rhs);
}

template<int stage>
void ProjectHcurlToL2<stage>::assemble()
{
	WorkStream::run(
		IteratorPair(IteratorTuple(dof_handler_L2.begin_active(),
			dof_handler_Hcurl.begin_active())),
		IteratorPair(IteratorTuple(dof_handler_L2.end(),
			dof_handler_Hcurl.end())),
		*this,
		&ProjectHcurlToL2<stage>::system_matrix_local,
		&ProjectHcurlToL2<stage>::copy_local_to_global,
		AssemblyScratchData(fe_L2,dof_handler_Hcurl,solution_Hcurl, mapping_degree),
		AssemblyCopyData());
}

} // namespace StaticVectorSolver

#endif

