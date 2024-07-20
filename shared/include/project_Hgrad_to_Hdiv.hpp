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

#ifndef ProjectHgradToHdiv_H__
#define ProjectHgradToHdiv_H__

#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_raviart_thomas.h>

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

#define VE scratch_data.ve

#define TMR(__name) \
	TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

namespace StaticScalarSolver
{
/**
 * \brief Projects from \f$H(\text{grad})\f$ to \f$H(\text{div})\f$
 *
 * This class template is not supposed to be used directly. Instead, one of the
 * wrap-around class templates, i.e., StaticScalarSolver::ProjectPHItoD,
 * StaticScalarSolver::ProjectPSItoB, etc., must be used. The names of the
 * wrap-around class templates are assumed to be more familiar to readers in
 * electromagnetics.
 *
 * Implements the following recipes:
 * - (1) Recipes for projections from H(grad) to H(div) nr. 6 and 7
 * - (2) Recipes for projections from H(grad) to H(div) nr. 8 and 9 (planar)
 * - (3) Recipes for projections from H(grad) to H(div) nr. 8 and 9 (axisym.)
 * - (4) Recipes for projections from H(grad) to H(div) nr. 10 and 11 (planar)
 * - (5) Recipes for projections from H(grad) to H(div) nr. 10 (axisym.)
 *
 * This class template is supposed to be used in pair with
 * StaticScalarSolver::Solver. In some problems the numerically calculated
 * scalar potential needs to be converted into a vector field. For example,
 * an electric scalar potential, \f$\Phi_h\f$, calculated with a help of
 * StaticScalarSolver::Solver may need to be converted into displacement
 * as
 * \f[
 * \vec{D}_h = - \epsilon \vec{\nabla} \Phi_h.
 * \f]
 * The electric scalar potential belongs to the H(grad) function space. The
 * displacement belongs to the H(div) function space. Therefore, one needs to
 * compute the equation above such that the input, \f$\Phi_h\f$, is in H(grad)
 * and the output, \f$\vec{D}_h\f$, is in H(div). Such computation can be
 * envisioned as a some kind of projection from one function space into another.
 * This class template does such projections. The Bossavit's diagrams
 * below illustrate the projections that can be made with a help of this class
 * template.
 *
 * ![](sss_prj/diagram_prj_Hgrad_to_Hdiv.svg)
 *
 * The table below lists the recommended settings for switching between
 * different types of projections. The six letters in the first column of the
 * table correspond to the six diagrams above. The dim parameter is the
 * input parameter of the class template. The other two parameters,
 * axisymmetric and vector_potential, are passed as input parameters
 * to the constructor of the class.
 *
 * Insert         | dim | axisymmetric | vector_potential |
 *----------------|-----|--------------|------------------|
 * A), B), C), D) |  3  |    false     |     false        |
 * E), F)         |  2  |  true/false  |     true         |
 *
 * The dim template parameter is, as per usual, the amount of spatial
 * dimensions. The purpose of the stage template parameter is discussed in
 * [here](@ref txt_stage_parameter).
 *
 * The user is supposed to derive an object from this class template. All usual
 * computations, i.e., setup, assembling the linear system, etc., happen
 * automatically.
 *****************************************************************************/
template<int dim, int stage = 1>
class ProjectHgradToHdiv
{
public:
	ProjectHgradToHdiv() = delete;

/**
 * \brief The only constructor.
 *
 * @param[in] p - The degree of the
 * [FE_RaviartThomas](https://www.dealii.org/current/doxygen/deal.II/classFE__RaviartThomas.html)
 * finite elements.
 * @param[in] mapping_degree - The degree of the interpolating polynomials used
 * for mapping. Setting it to 1 will do in the most of the cases. Note, that it
 * makes sense to attach a meaningful manifold to the triangulation if this
 * parameter is greater than 1.
 * @param[in] triangulation_Hgrad - Reference to the triangulation inside the
 * object of the class StaticScalarSolver::Solver that has yielded the
 * scalar potential in the H(grad) function space.
 * @param[in] dof_handler_Hgrad - Reference to the dof handler inside the class
 * StaticScalarSolver::Solver that has yielded the scalar potential in the
 * H(grad) function space.
 * @param[in] solution_Hgrad - Vector filled with the degrees of freedom that
 * together with the shape functions of the
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html)
 * finite elements model the scalar potential in the H(grad) function space.
 * @param[in] fname - The name of the output files without extension.
 * @param[in] exact_solution - Points to an object that describes the exact
 * solution to the problem. It is needed for calculating error norms.
 * @param[in] axisymmetric - If true, assumes that the problem is axisymmetric.
 * If axisymmetric=true, dim must equal 2.
 * @param[in] vector_potential - If true, assumes that the problem is
 * two-dimensional and formulated in terms of the magnitude of vector potential,
 * \f$A\f$, or in terms of the scaled magnitude of vector potential, \f$A'\f$,
 * or current vector potential, \f$T\f$. If vector_potential=true, dim
 * must equal 2.
 * @param[in] print_time_tables - If true, prints time tables on the screen.
 * @param[in] project_exact_solution - If true, projects the exact solution
 * onto the space spanned by the Raviart-Thomas finite elements and saves the
 * result into the vtk file next to the solution. This may be useful for
 * debugging purposes.
 * @param[in] log_cg_convergence - If true, logs convergence of the conjugate
 * gradient solver into a file. The name of the file is generated by appending
 * "_cg_convergence.csv" to fname.
 ******************************************************************************/
	ProjectHgradToHdiv(
		unsigned int p,
		unsigned int mapping_degree,
		const Triangulation<dim> & triangulation_Hgrad,
		const DoFHandler<dim> & dof_handler_Hgrad,
		const Vector<double> & solutioin_Hgrad,
		std::string fname = "Hdiv",
		const Function<dim> * exact_solution = nullptr,
		bool axisymmetric = false,
		bool vector_potential = false,
		bool print_time_tables = false,
		bool project_exact_solution = false,
		bool log_cg_convergence = false);

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
		{return static_cast<unsigned int>(triangulation_Hgrad.n_active_cells());}

/**
 * \brief Returns the total amount of the degrees of freedom.
 *****************************************************************************/
	unsigned int get_n_dofs() const
		{return static_cast<unsigned int>(dof_handler_Hdiv.n_dofs());}

/**
 * \brief Releases computer memory associated with system matrix and
 * right-hand side.
 *****************************************************************************/
	void clear() {system_matrix.clear(); system_rhs.reinit(0);}

/**
 * \brief Returns a reference to triangulation.
 *****************************************************************************/
	const Triangulation<dim> & get_tria() const {return triangulation_Hgrad;}

/**
 * \brief Returns a reference to dof handler associated with the Raviart-Thomas
 * finite elements.
 *****************************************************************************/
	const DoFHandler<dim> & get_dof_handler() const {return dof_handler_Hdiv;}

/**
 * \brief Returns a reference to solution, i.e., the result of the projection.
 ******************************************************************************/
	const Vector<double> & get_solution() const {return solution_Hdiv;}

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

	const DoFHandler<dim> & dof_handler_Hgrad;
	const Vector<double> & solution_Hgrad;

	const Triangulation<dim> & triangulation_Hgrad;
	const FE_RaviartThomas<dim> fe_Hdiv;
	DoFHandler<dim> dof_handler_Hdiv;

	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double> solution_Hdiv;
	Vector<double> system_rhs;

	Vector<double> projected_exact_solution;

	AffineConstraints<double> constraints;

	const Function<dim> * exact_solution;

	const unsigned int mapping_degree;
	const bool axisymmetric;
	const bool vector_potential;
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
		std::tuple<typename DoFHandler<dim>::active_cell_iterator,
		           typename DoFHandler<dim>::active_cell_iterator>;

	using IteratorPair = SynchronousIterators<IteratorTuple>;

	struct AssemblyScratchData
	{
		AssemblyScratchData(const FiniteElement<dim> &fe,
			const DoFHandler<dim> & dof_handr_Hgrad,
			const Vector<double> & dofs_Hgrad,
			bool axisymmetric,
			bool vector_potential,
			unsigned int mapping_degree);

		AssemblyScratchData(const AssemblyScratchData & scratch_data);

		MappingQ<dim> mapping;
		Constants::QuadratureTableVector<dim> qt;
		FEValues<dim> fe_values_Hdiv;
		FEValues<dim> fe_values_Hgrad;

		const unsigned int dofs_per_cell;
		const unsigned int n_q_points;

		TheCoefficient<dim,stage> the_coefficient;
		std::vector<double> the_coefficient_list;

		std::vector<Tensor<1,dim>> vector_gradients;

		// Two-dimensional vector curl of an out-of-plane (oop) vector.
		std::vector<Tensor<1,dim>> nabla_xV_oopvector;

		const FEValuesExtractors::Vector ve;

		const DoFHandler<dim> & dof_hand_Hgrad;
		const Vector<double> & dofs_Hgrad;

		const bool axisymmetric;
		const bool vector_potential;

		double axi_mult; // Equals the distance to the axis of rotation symmetry, r,
		                 // in the recipes for axisymmetric projections. Equals 1.0
		                 // in all other recipes. All integrands are multiplied by
		                 // this multiplier.
	};

	struct AssemblyCopyData
	{
		FullMatrix<double> cell_matrix;
		Vector<double> cell_rhs;
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

template<int dim, int stage>
ProjectHgradToHdiv<dim,stage>::ProjectHgradToHdiv(
	unsigned int p,
	unsigned int mapping_degree,
	const Triangulation<dim> & triangulation_Hgrad,
	const DoFHandler<dim> & dof_handler_Hgrad,
	const Vector<double> & solution_Hgrad,
	std::string fname,
	const Function<dim> * exact_solution,
	bool axisymmetric,
	bool vector_potential,
	bool print_time_tables,
	bool project_exact_solution,
	bool log_cg_convergence):
		fname(fname),
		dof_handler_Hgrad(dof_handler_Hgrad),
		solution_Hgrad(solution_Hgrad),
		triangulation_Hgrad(triangulation_Hgrad),
 // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 // The public attribute fe.degree is the maximal polynomial degree of a
 // shape function in a single coordinate direction, not the degree of
 // the finite element. For FE_Nedelec and FE_RaviartThomas degree of
 // the finite element is: degree_of_element = fe.degree - 1.
 // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 //		fe_Hdiv(dof_handler_Hgrad.get_fe().degree - 1 ),
		fe_Hdiv(p),
		exact_solution(exact_solution),
		mapping_degree(mapping_degree),
		axisymmetric(axisymmetric),
		vector_potential(vector_potential),
		project_exact_solution(project_exact_solution),
		log_cg_convergence(log_cg_convergence)
{
	if (axisymmetric)
	{
			Assert( dim == 2, ExcMessage(
"The setting axisymmetric=true is only allowed if dim=2."));
	}

	if (vector_potential)
	{
				Assert( dim == 2, ExcMessage(
"The setting vector_potential=true can only be used if dim=2."));
	}

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

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::setup()
{
	constraints.close();

	dof_handler_Hdiv.reinit(triangulation_Hgrad);
	dof_handler_Hdiv.distribute_dofs(fe_Hdiv);

	DynamicSparsityPattern dsp(dof_handler_Hdiv.n_dofs(), dof_handler_Hdiv.n_dofs());
	DoFTools::make_sparsity_pattern(
		dof_handler_Hdiv,
		dsp,
		constraints,
		false);

	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	solution_Hdiv.reinit(dof_handler_Hdiv.n_dofs());
	system_rhs.reinit(dof_handler_Hdiv.n_dofs());

	if (project_exact_solution)
		projected_exact_solution.reinit(dof_handler_Hdiv.n_dofs());

	if (exact_solution)
	{
		L2_per_cell.reinit(triangulation_Hgrad.n_active_cells());
		Linfty_per_cell.reinit(triangulation_Hgrad.n_active_cells());
	}
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::solve()
{
	ReductionControl control(system_rhs.size(), 0.0, 1e-12, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

	PreconditionJacobi<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(system_matrix, 1.0);

	cg.solve(system_matrix, solution_Hdiv, system_rhs, preconditioner);

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

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::save() const
{
	std::vector<std::string> solution_names(dim, "VectorField");
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	interpretation(dim,
		DataComponentInterpretation::component_is_part_of_vector);

	DataOut<dim> data_out;

	data_out.add_data_vector(dof_handler_Hdiv,
		solution_Hdiv,
		solution_names,
		interpretation);

	if (project_exact_solution)
	{
		std::vector<std::string> solution_names_ex(dim, "VectorFieldExact");

		data_out.add_data_vector(dof_handler_Hdiv,
		projected_exact_solution,
		solution_names_ex,
		interpretation);
	}

	data_out.add_data_vector(L2_per_cell, "L2norm");
	data_out.add_data_vector(Linfty_per_cell, "LinftyNorm");

	data_out.build_patches();

	std::ofstream out(fname + ".vtk");

	data_out.write_vtk(out);
	out.close();
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::compute_error_norms()
{
	Weight<dim,stage> weight;
	const Function<dim, double> * mask = & weight;

	Constants::QuadratureTableVector<dim> qt(dof_handler_Hdiv.get_fe().degree - 1);
	QGauss<dim> quadrature(qt.enorm());

	VectorTools::integrate_difference(
		MappingQ<dim> (mapping_degree),
		dof_handler_Hdiv,
		solution_Hdiv,
		* exact_solution,
		L2_per_cell,
		quadrature,
		VectorTools::L2_norm,
		mask);

	L2_norm = VectorTools::compute_global_error(
		triangulation_Hgrad,
		L2_per_cell,
		VectorTools::L2_norm);

	VectorTools::integrate_difference(
		MappingQ<dim> (mapping_degree),
		dof_handler_Hdiv,
		solution_Hdiv,
		* exact_solution,
		Linfty_per_cell,
		QGauss<dim>(1),
		VectorTools::Linfty_norm,
		mask // & B_mask
		);

	Linfty_norm = Linfty_per_cell.linfty_norm();
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::project_exact_solution_fcn()
{
	Constants::QuadratureTableVector<dim> qt(fe_Hdiv.degree - 1);

	AffineConstraints<double> constraints_empty;
	constraints_empty.close();

	VectorTools::project(
		MappingQ<dim> (mapping_degree),
		dof_handler_Hdiv,
		constraints_empty,
		QGauss<dim> (qt.sim()),
		* exact_solution,
		projected_exact_solution);
}

template<int dim, int stage>
ProjectHgradToHdiv<dim,stage>::AssemblyScratchData::AssemblyScratchData(
	const FiniteElement<dim> &fe,
	const DoFHandler<dim> & dof_hand_Hgrad,
	const Vector<double> & dofs_Hgrad,
	bool axisymmetric,
	bool vector_potential,
	unsigned int mapping_degree):
		mapping(mapping_degree),
		qt(fe.degree - 1),
		fe_values_Hdiv(mapping, fe, QGauss<dim>(qt.sim()),
			update_values |
			update_quadrature_points |
			update_JxW_values),
		fe_values_Hgrad(mapping, dof_hand_Hgrad.get_fe(), QGauss<dim>(qt.sim()),
			update_gradients),
		dofs_per_cell(fe_values_Hdiv.dofs_per_cell),
		n_q_points(fe_values_Hdiv.get_quadrature().size()),
		the_coefficient_list(n_q_points),
		vector_gradients(n_q_points,Tensor<1,dim>()),
		nabla_xV_oopvector(n_q_points,Tensor<1,dim>()),
		ve(0),
		dof_hand_Hgrad(dof_hand_Hgrad),
		dofs_Hgrad(dofs_Hgrad),
		axisymmetric(axisymmetric),
		vector_potential(vector_potential),
		axi_mult(1.0)
{
}

template<int dim, int stage>
ProjectHgradToHdiv<dim,stage>::AssemblyScratchData::AssemblyScratchData(
		const AssemblyScratchData & scratch_data):
		mapping(scratch_data.mapping.get_degree()),
		qt(scratch_data.qt),
			fe_values_Hdiv(mapping, scratch_data.fe_values_Hdiv.get_fe(),
				scratch_data.fe_values_Hdiv.get_quadrature(),
				update_values |
				update_quadrature_points |
				update_JxW_values),
			fe_values_Hgrad(mapping,scratch_data.fe_values_Hgrad.get_fe(),
				scratch_data.fe_values_Hgrad.get_quadrature(),
				update_gradients),
			dofs_per_cell(fe_values_Hdiv.dofs_per_cell),
			n_q_points(fe_values_Hdiv.get_quadrature().size()),
			the_coefficient_list(n_q_points),
			vector_gradients(n_q_points,Tensor<1,dim>()),
			nabla_xV_oopvector(n_q_points,Tensor<1,dim>()),
			ve(0),
			dof_hand_Hgrad(scratch_data.dof_hand_Hgrad),
			dofs_Hgrad(scratch_data.dofs_Hgrad),
			axisymmetric(scratch_data.axisymmetric),
			vector_potential(scratch_data.vector_potential),
			axi_mult(1.0)
{
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::system_matrix_local(
	const IteratorPair & IP,
	AssemblyScratchData & scratch_data, AssemblyCopyData & copy_data)
{
// See the color boxes
// (1) Recipe for projections from H(grad) to H(div) nr. 6 and 7
// (2) Recipe for projections from H(grad) to H(div) nr. 8 and 9 (planar)
// (3) Recipe for projections from H(grad) to H(div) nr. 8 and 9 (axisym.)
// (4) Recipe for projections from H(grad) to H(div) nr. 10 and 11 (planar)
// (5) Recipe for projections from H(grad) to H(div) nr. 10 (axisym.)
//
// The comments below refer to these recipes by number, i.e., recipe (1),
// recipe (2), etc.

	copy_data.cell_matrix.reinit(
		scratch_data.dofs_per_cell, scratch_data.dofs_per_cell);

	copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

	copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

	scratch_data.fe_values_Hdiv.reinit(std::get<0>(*IP));
	scratch_data.fe_values_Hgrad.reinit(std::get<1>(*IP));

	scratch_data.fe_values_Hgrad.get_function_gradients(
		scratch_data.dofs_Hgrad,
		scratch_data.vector_gradients);

	if (! scratch_data.vector_potential)
		scratch_data.the_coefficient.value_list(
			scratch_data.fe_values_Hdiv.get_quadrature_points(),
			std::get<0>(*IP)->material_id(),
			std::get<0>(*IP)->user_index(),
			scratch_data.the_coefficient_list);

	for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
	{

		scratch_data.axi_mult = 1.0;

		if ((scratch_data.axisymmetric)&&(!scratch_data.vector_potential))
			scratch_data.axi_mult = scratch_data.fe_values_Hdiv.quadrature_point(q_index)[0];

		for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
		{
			for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
			{
				// Mass matrix is the same in all recipes.
				copy_data.cell_matrix(i, j) +=
					scratch_data.axi_mult *                                      // 1.0 or r
					scratch_data.fe_values_Hdiv[VE].value(i,q_index) *           // curl N_i
					scratch_data.fe_values_Hdiv[VE].value(j,q_index) *           // curl N_j
					scratch_data.fe_values_Hdiv.JxW(q_index);                    // dV (dS in 2D)
			}

			if (scratch_data.vector_potential) // A and A'.
			{
				// The option vector_potential = true exists only in 2D
				scratch_data.nabla_xV_oopvector.at(q_index)[0] =
					scratch_data.vector_gradients.at(q_index)[1];

				scratch_data.nabla_xV_oopvector.at(q_index)[1] = // nabla_xV_oopvector is the
				 -scratch_data.vector_gradients.at(q_index)[0];  // two-dimensional vector
				                                                 // field in H(div),
				                                                 // 'curl_v A' for short.
				double tmp;
				tmp =
						scratch_data.nabla_xV_oopvector.at(q_index) *               // curl_v A
						scratch_data.fe_values_Hdiv[VE].value(i,q_index) *          // N_i
						scratch_data.fe_values_Hdiv.JxW(q_index);                   // dS

				if (axisymmetric) // A'.
				{// Integral b_i in recipe (5).
					copy_data.cell_rhs(i) -= tmp;

// The recipe (5) yields B'=rB, where r is the distance to the axis of rotation symmetry.
// If you do not like to have the scaled magnetic field, B', and would rather have the magnetic
// field itself, B, replace the line above with the following line. Note, that a quadrature point
// can, in general, be at the origin. If so, the line below will yield an error.

//					copy_data.cell_rhs(i) -= tmp/scratch_data.fe_values_Hdiv.quadrature_point(q_index)[0];
				}
				else // A.
				{// Integral b_i in recipe (4).
					copy_data.cell_rhs(i) += tmp;
				}

			}else // Phi, Psi, Theta.
			{
			// Integral b_i in recipes (1), (2), and (3).
				copy_data.cell_rhs(i) -=
					scratch_data.axi_mult *                                     // 1.0 or r
					scratch_data.the_coefficient_list[q_index] *                // epsilon
					scratch_data.vector_gradients.at(q_index) *                 // grad PHI
					scratch_data.fe_values_Hdiv[VE].value(i,q_index) *          // N_i
					scratch_data.fe_values_Hdiv.JxW(q_index);                   // dV (dS in 2D)
			}
		}
	}

	std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::copy_local_to_global(const AssemblyCopyData &copy_data)
{
	constraints.distribute_local_to_global(
		copy_data.cell_matrix,
		copy_data.cell_rhs,
		copy_data.local_dof_indices,
		system_matrix,
		system_rhs);
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::assemble()
{
	WorkStream::run(
		IteratorPair(IteratorTuple(dof_handler_Hdiv.begin_active(),
     dof_handler_Hgrad.begin_active())),
		IteratorPair(IteratorTuple(dof_handler_Hdiv.end(),
     dof_handler_Hgrad.end())),
		*this,
		&ProjectHgradToHdiv::system_matrix_local,
		&ProjectHgradToHdiv::copy_local_to_global,
		AssemblyScratchData(
			fe_Hdiv,
			dof_handler_Hgrad,
			solution_Hgrad,
			axisymmetric,
			vector_potential,
			mapping_degree),
		AssemblyCopyData());
}

template<int dim, int stage>
void ProjectHgradToHdiv<dim,stage>::save_matrix_and_rhs_to_csv(std::string fname) const
{
	std::ofstream ofs_matrix(fname + "_matrix.csv");
	std::ofstream ofs_rhs(fname + "_rhs.csv");

	for (unsigned int i = 0; i < system_matrix.m(); ++i )
	{
		ofs_rhs << system_rhs(i);
		if ( i < (system_matrix.m() - 1) )
				ofs_rhs << "\n";

		for (unsigned int j = 0; j < system_matrix.n(); ++j )
		{
			ofs_matrix
				<< std::scientific
				<< std::setprecision(16)
				<< system_matrix.el(i,j);

			if ( j <  (system_matrix.m() - 1) )
				ofs_matrix << ", ";
		}
		if ( i < (system_matrix.m() - 1) )
			ofs_matrix << "\n";
	}

	ofs_rhs.close();
	ofs_matrix.close();
}

} // namespace StaticScalarSolver

#endif

