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
#include "project_Tz_to_Jxy.hpp"

#include "settings.hpp"

using namespace Misc;
using namespace StaticScalarSolver;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Current vector potential (cvp-ii/)](@ref page_cvp_ii)
 * numerical experiment.
 *****************************************************************************/
class BatchCVPII : public SettingsCVPII
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		std::string fname;
		std::string dir = "Data/circle-linear/";

		std::cout
			<< "Program: cvp-ii\n"
			<< "Dimensions: " << "2\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_T(2);
		MainOutputTable table_J(2);

		for (unsigned int p = 1; p < 4; p++)
		{
			table_T.clear();
			table_J.clear();

			for (unsigned int r = 12; r < 16; r++)
			{
				std::cout << "Solving for T ...\n";

				fname = dir + "solution_T_p" + std::to_string(p)
					+ "_r" + std::to_string(r);

				table_T.add_value("r", r);
				table_T.add_value("p", p);

				if (SettingsCVPII::print_time_tables)
					std::cout << "Time table\n";

				SolverCVPII problem(p, 2, r, fname);
				problem.clear();

				table_T.add_value("ndofs", problem.get_n_dofs());
				table_T.add_value("ncells", problem.get_n_cells());
				table_T.add_value("L2", problem.get_L2_norm());
				table_T.add_value("H1", problem.get_H1_norm());

				std::cout << "Converting T to J. \n";

				fname = dir + "solution_J_p" + std::to_string(p)
					+ "_r" + std::to_string(r);

				table_J.add_value("r", r);
				table_J.add_value("p", p);

				ExactSolutionCVPII_Jf exact_solution;

				if (SettingsCVPII::print_time_tables)
					std::cout << "Time table\n";

				ProjectTzToJxy projector(
				p-1,
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
			std::cout << "Table T\n";
			table_T.save(dir + "table_T_p" + std::to_string(p));

			std::cout << "Table J\n";
			table_J.save(dir + "table_J_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchCVPII batch;
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

