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
#include <deal.II/base/multithread_info.h>

#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"
#include "project_T_to_J.hpp"

#include "settings.hpp"

using namespace Misc;
using namespace StaticVectorSolver;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Current vector potential (cvp-i/)](@ref page_cvp_i)
 * numerical experiment.
 *****************************************************************************/
class BatchCVPI : public SettingsCVPI
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		std::string fname;
		std::string dir = "Data/sphere-linear/";

		std::cout
			<< "Program: cvp-i\n"
			<< "Dimensions: " << "3\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_J(3);

		for (unsigned int p = 0; p < 3; p++)
		{
			table_J.clear();

			for (unsigned int r = 12; r < 16; r++)
			{
				std::cout << "Solving for T. \n";

				fname = dir + "solution_T_p" + std::to_string(p)
					+ "_r" + std::to_string(r);

				if (SettingsCVPI::print_time_tables)
					std::cout << "Time table\n";

				SolverCVPI problem(p, 2, r, fname);
				problem.clear();

				std::cout << "Converting T to J. \n";

				fname = dir + "solution_J_p" + std::to_string(p)
					+ "_r" + std::to_string(r);

				table_J.add_value("r", r);
				table_J.add_value("p", p);

				ExactSolutionCVPI_Jf exact_solution;

				if (SettingsCVPI::print_time_tables)
					std::cout << "Time table\n";

				ProjectTtoJ projector(
					p,
					2,
					problem.get_tria(),
					problem.get_dof_handler(),
					problem.get_solution(),
					fname,
					& exact_solution,
					Settings::print_time_tables,
					Settings::project_exact_solution,
					Settings::log_cg_convergence);
		
					table_J.add_value("ndofs", projector.get_n_dofs());
					table_J.add_value("ncells", projector.get_n_cells());
					table_J.add_value("L2", projector.get_L2_norm()); 
					table_J.add_value("H1", 0.0);
			}
			// Saving convergence table

			std::cout << "Table J\n";
			table_J.save(dir + "table_J_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchCVPI batch;
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

