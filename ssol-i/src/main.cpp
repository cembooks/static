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

#include <deal.II/base/timer.h>
#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"
#include "project_A_to_B.hpp"

#include "deal.II/base/multithread_info.h"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Thin spherical coil (ssol-i/)](@ref page_ssol_i)
 * numerical experiment.
 *
 * The purpose of this class is to make the main function to be similar to the
 * main function of the deal.II
 * [Step-6](https://dealii.org/developer/doxygen/deal.II/step_6.html)
 * tutorial.
 *****************************************************************************/
class BatchSSOLI : public SettingsSSOLI
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		std::string dir = "Data/sphere/";
		std::string fname;

		std::cout
			<< "Program: ssol-i\n"
			<< "Dimensions: " << "3\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_B(3);

		for (unsigned int p = 0; p < 1; p++)
		{
			table_B.clear();

			for (unsigned int r = 14; r < 18; r++)
			{	
				fname = dir + "solution_A_p" + std::to_string(p) +
					"_r" + std::to_string(r);

				if (SettingsSSOLI::print_time_tables)
					std::cout << "Time table A \n";

				SolverSSOLI problem(p, 1, r, fname);

				problem.clear();
				{
					fname = dir + "solution_B_p" + std::to_string(p) +
						"_r" + std::to_string(r);

					table_B.add_value("r", r);
					table_B.add_value("p", p);

					if (SettingsSSOLI::print_time_tables)
						std::cout << "Time table B \n";

					ExactSolutionSSOLI_B exact_solution;

					ProjectAtoB projector(
						p,
						1,
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_B.add_value("ndofs", projector.get_n_dofs());
					table_B.add_value("ncells", projector.get_n_cells());
					table_B.add_value("L2", projector.get_L2_norm());
					table_B.add_value("H1", 0.0);
				}
			}
			// Saving convergence tables

			std::cout << "Table B\n";
			table_B.save(dir + "table_B_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchSSOLI batch;
		batch.run();
	}
	catch (std::exception &exc)
	{
		std::cerr
			<< std::endl
			<< std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr
			<< "Exception on processing: " << std::endl
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

