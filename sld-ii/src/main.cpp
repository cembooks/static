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

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/multithread_info.h>

#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"
#include "project_THETA_to_H.hpp"
#include "project_THETA_to_B.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Magnetostatic shield - 2 (sld-ii/)](@ref page_sld_ii)
 * numerical experiment.
 *
 * The purpose of this class is to make the main function to be similar to the
 * main function of the deal.II
 * [Step-6](https://dealii.org/developer/doxygen/deal.II/step_6.html)
 * tutorial.
 *****************************************************************************/
class BatchSLDII : public SettingsSLDII
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		Assert( DIMENSION__ > 1 , ExcInternalError())
		Assert( DIMENSION__ < 4 , ExcInternalError())

		std::string fname;
		std::string dir = (DIMENSION__ == 2) ? "Data/square/" : "Data/cube/";

		std::cout << "----------------------------------------\n"
			<< "Program: sld-ii\n"
			<< "Version: 1.0\n"
			<< "Dimension = " << DIMENSION__ << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_THETA(DIMENSION__);
		MainOutputTable table_B(DIMENSION__);
		MainOutputTable table_H(DIMENSION__);

		for (unsigned int p = 1; p < 4; p++)
		{
			table_THETA.clear();
			table_B.clear();
			table_H.clear();

#if DIMENSION__ == 2
			for (unsigned int r = 10; r < 14; r++)
#endif
#if DIMENSION__ == 3
			for (unsigned int r = 5; r < 9; r++)
#endif
			{
				table_THETA.add_value("r", r);
				table_THETA.add_value("p", p);

				fname = dir + "solution_THETA_p" + std::to_string(p)
					+ "_r" + std::to_string(r);

				SolverSLDII<DIMENSION__> problem(p, 2, r, fname);
				table_THETA.add_value("ndofs", problem.get_n_dofs());
				table_THETA.add_value("ncells", problem.get_n_cells());
				table_THETA.add_value("L2", problem.get_L2_norm());
				table_THETA.add_value("H1", problem.get_H1_norm());

				problem.clear();

				{
					table_H.add_value("r", r);
					table_H.add_value("p", p);

					fname = dir + "solution_H_p" + std::to_string(p)
						+ "_r" + std::to_string(r);

					if (SettingsSLDII::print_time_tables)
						std::cout << "Time table H \n";

					ExactSolutionSLDII_H<DIMENSION__> exact_solution;

					ProjectTHETAtoH<DIMENSION__> projector(
						p-1,
						problem.get_mapping_degree(),
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						false,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_H.add_value("ndofs", projector.get_n_dofs());
					table_H.add_value("ncells", projector.get_n_cells());
					table_H.add_value("L2", projector.get_L2_norm());
					table_H.add_value("H1", 0.0);
				}
				{
					table_B.add_value("r", r);
					table_B.add_value("p", p);

					fname = dir + "solution_B_p" + std::to_string(p)
						+ "_r" + std::to_string(r);

					if (SettingsSLDII::print_time_tables)
						std::cout << "Time table B \n";

					ExactSolutionSLDII_B<DIMENSION__> exact_solution;

					ProjectTHETAtoB<DIMENSION__> projector(
						p-1,
						problem.get_mapping_degree(),
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						false,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_B.add_value("ndofs", projector.get_n_dofs());
					table_B.add_value("ncells", projector.get_n_cells());
					table_B.add_value("L2", projector.get_L2_norm() / mu_0);
					table_B.add_value("H1", 0.0);
				}
			}
			// Saving convergence tables
			std::cout << "Table THETA\n";
			table_THETA.save(dir + "table_THETA_p" + std::to_string(p));

			std::cout << "Table H\n";
			table_H.save(dir + "table_H_p" + std::to_string(p));

			std::cout << "Table B\n";
			table_B.save(dir + "table_B_p" + std::to_string(p));
		}
	}


};

int main()
{
	try
	{
		BatchSLDII batch;
		batch.run();
		}
	catch (std::exception &exc)
	{
		std::cerr
			<< std::endl
			<< std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
			<< exc.what() << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}
	catch (...)
	{
		std::cerr
			<< std::endl
			<< std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr
			<< "Unknown exception!" << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}
	return 0;
}

